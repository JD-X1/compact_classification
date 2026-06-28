#!/usr/bin/env python3

import argparse
import math
import os
from typing import Dict, List, Tuple

import pandas as pd


RANKS = [
    "superkingdom",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
    "strain",
]


def clean_value(value) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    if text.lower() in {"", "nan", "none"}:
        return ""
    return text


def split_lineage(text: str) -> List[str]:
    return [clean_value(x) for x in clean_value(text).split(";") if clean_value(x)]


def join_lineage(parts: List[str]) -> str:
    return ";".join([x for x in parts if x])


def pick_query_col(df: pd.DataFrame) -> str:
    names = ["query", "query_name", "query_names", "name", "placement_name", "sample"]
    cols = {c.lower(): c for c in df.columns}
    for name in names:
        if name in cols:
            return cols[name]
    raise ValueError(
        "No per-query key column found in the gappa assignment table "
        f"(columns={list(df.columns)}). This usually means gappa was run without "
        "--per-query-results and only the sample-level profile.tsv was produced. "
        "Pass per_query.tsv (run `gappa examine assign --per-query-results`)."
    )


def pick_score_col(df: pd.DataFrame) -> str:
    names = ["like_weight_ratio", "likelihood_weight_ratio", "lwr", "probability", "score", "support"]
    cols = {c.lower(): c for c in df.columns}
    for name in names:
        if name in cols:
            return cols[name]
    return ""


def clamp01(value) -> float:
    """Clamp a per-query like-weight ratio to [0, 1].

    The LWR from gappa's per_query.tsv is, by definition, the marginal ML
    likelihood of a placement normalized over the query's candidate edges
    (Matsen et al. 2010, BMC Bioinformatics 11:538), so it already lies in
    [0, 1]. We clamp defensively rather than rescale: a value >1 would mean a
    summed/sample-level column (profile.tsv's LWR) was passed by mistake, and
    silently dividing such a value by 100 (the previous behaviour) crushed
    confident placements. Clamping keeps a real LWR untouched and caps an
    out-of-range value at 1 instead of mangling it.
    """
    try:
        x = float(value)
    except Exception:
        return 0.0
    if x < 0:
        return 0.0
    if x > 1:
        return 1.0
    return x


def lineage_from_row(row: pd.Series) -> Tuple[str, str]:
    values = []
    for rank in RANKS:
        if rank in row.index:
            val = clean_value(row[rank])
            if val:
                values.append((rank, val))
    if values:
        return join_lineage([x[1] for x in values]), values[-1][0]
    for col in row.index:
        val = clean_value(row[col])
        if ";" in val:
            parts = split_lineage(val)
            if parts:
                rank = RANKS[min(len(parts), len(RANKS)) - 1]
                return join_lineage(parts), rank
    return "", "unclassified"


def read_profile(path: str) -> Dict[str, Dict]:
    """Per-query placement: taxopath, deepest rank, and best-edge LWR.

    Reads gappa's per_query.tsv (keyed on the pquery `name`), keeping for each
    query the single highest-LWR taxopath. The LWR is the native, calibrated
    placement-confidence quantity and is reported as-is (clamped to [0, 1]); no
    bespoke composite score is computed.
    """
    if not os.path.exists(path):
        return {}
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        return {}
    qcol = pick_query_col(df)
    scol = pick_score_col(df)
    if scol:
        idx = df.groupby(qcol)[scol].idxmax()
        best = df.loc[idx]
    else:
        best = df.groupby(qcol, as_index=False).first()

    out = {}
    for _, row in best.iterrows():
        q = clean_value(row[qcol])
        lineage, rank = lineage_from_row(row)
        lwr = clamp01(row[scol]) if scol else 0.0
        out[q] = {"lineage": lineage, "rank": rank, "lwr": lwr}
    return out


def read_marker(marker_stats: str, kept_genes: str) -> Dict:
    out = {
        "total": 0,
        "kept": 0,
        "kept_fraction": 0.0,
        "mean_coverage": 0.0,
        "mean_gap_fraction": 1.0,
    }
    if not marker_stats or not os.path.exists(marker_stats):
        return out
    df = pd.read_csv(marker_stats, sep="\t")
    if df.empty or "gene" not in df.columns:
        return out

    keep = set()
    if kept_genes and os.path.exists(kept_genes):
        with open(kept_genes, "r", encoding="utf-8") as fh:
            keep = {clean_value(line) for line in fh if clean_value(line)}
    else:
        keep = set(df["gene"].astype(str))

    out["total"] = int(len(df))
    out["kept"] = int(sum(1 for x in df["gene"].astype(str) if x in keep))
    out["kept_fraction"] = (out["kept"] / out["total"]) if out["total"] else 0.0

    if "coverage" in df.columns:
        sub = df[df["gene"].astype(str).isin(keep)] if keep else df
        if not sub.empty:
            out["mean_coverage"] = float(sub["coverage"].astype(float).mean())
    if "gap_fraction" in df.columns:
        sub = df[df["gene"].astype(str).isin(keep)] if keep else df
        if not sub.empty:
            out["mean_gap_fraction"] = float(sub["gap_fraction"].astype(float).mean())
    return out


def taxon_at_depth(taxopath: str, depth: int) -> str:
    parts = split_lineage(taxopath)
    return parts[depth] if 0 <= depth < len(parts) else ""


def weighted_support_at_depth(
    observations: List[Tuple[str, float]], depth: int
) -> Tuple[Dict[str, float], int, float]:
    """Confidence-weighted distribution over taxa at one taxopath depth.

    Returns (support, n_markers, n_eff): support maps taxon -> weight fraction
    (sums to 1); n_markers is the count of markers carrying a taxon at this
    depth; n_eff is Kish's effective sample size (Sum w)^2 / Sum w^2, used for
    the small-sample entropy correction.
    """
    support: Dict[str, float] = {}
    total = 0.0
    sum_w2 = 0.0
    n = 0
    for taxopath, weight in observations:
        taxon = taxon_at_depth(taxopath, depth)
        if not taxon:
            continue
        w = float(weight) if weight and float(weight) > 0 else 0.0
        if w <= 0:
            continue
        support[taxon] = support.get(taxon, 0.0) + w
        total += w
        sum_w2 += w * w
        n += 1
    if total <= 0:
        return {}, 0, 0.0
    support = {t: s / total for t, s in support.items()}
    n_eff = (total * total) / sum_w2 if sum_w2 > 0 else 0.0
    return support, n, n_eff


def shannon_effective(
    support: Dict[str, float], n_eff: float
) -> Tuple[float, float, int, float]:
    """Shannon entropy and the effective number of taxa exp(H) (Hill order 1).

    H is the weighted Shannon entropy of the support distribution; H_mm applies
    the Miller-Madow upward bias correction (S-1)/(2 n_eff); the effective
    number is exp(H_mm), capped at the observed richness S (the maximum a
    diversity of order 1 can take). Returns (H, H_mm, richness, effective).
    """
    ps = [p for p in support.values() if p > 0]
    richness = len(ps)
    if richness == 0:
        return 0.0, 0.0, 0, 0.0
    H = -sum(p * math.log(p) for p in ps)
    H_mm = H + (richness - 1) / (2.0 * n_eff) if n_eff and n_eff > 0 else H
    effective = min(math.exp(H_mm), float(richness))
    return H, H_mm, richness, effective


def heterogeneity_rows(
    gappa: Dict[str, Dict], mag_name: str, min_markers: int
) -> List[Dict]:
    """Per-MAG effective number of taxa at each taxopath depth (Upgrade 1).

    Aggregates the per-marker placements (one gappa pquery each) into a
    confidence-weighted lineage distribution per depth and reports its effective
    number of taxa exp(H). ~1 => single lineage (clean); ~2 => binary chimera.
    status is 'insufficient_markers' when too few markers placed for the entropy
    estimate to be trustworthy; the numbers are still emitted for transparency.
    """
    observations = [
        (g["lineage"], g.get("lwr", 0.0))
        for g in gappa.values()
        if clean_value(g.get("lineage"))
    ]
    n_placed = len(observations)
    status = "ok" if n_placed >= min_markers else "insufficient_markers"
    max_depth = max((len(split_lineage(tp)) for tp, _ in observations), default=0)

    rows = []
    for depth in range(max_depth):
        support, n_at, n_eff = weighted_support_at_depth(observations, depth)
        if not support:
            continue
        H, H_mm, richness, effective = shannon_effective(support, n_eff)
        rows.append(
            {
                "mag": mag_name,
                "level": depth,
                "n_markers_placed": n_placed,
                "n_markers_at_level": n_at,
                "n_eff_markers": round(n_eff, 4),
                "richness": richness,
                "shannon_entropy": round(H, 6),
                "shannon_entropy_mm": round(H_mm, 6),
                "effective_num_taxa": round(effective, 4),
                "status": status,
            }
        )
    return rows


HETEROGENEITY_COLUMNS = [
    "mag",
    "level",
    "n_markers_placed",
    "n_markers_at_level",
    "n_eff_markers",
    "richness",
    "shannon_entropy",
    "shannon_entropy_mm",
    "effective_num_taxa",
    "status",
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--profile", required=True, help="gappa per_query.tsv")
    ap.add_argument("--out-report", required=True)
    ap.add_argument("--marker-stats", default="")
    ap.add_argument("--kept-genes", default="")
    ap.add_argument("--min-markers", type=int, default=10)
    ap.add_argument("--out-flag", default="")
    ap.add_argument("--out-heterogeneity", default="")
    ap.add_argument("--mag-name", default="")
    args = ap.parse_args()

    gappa = read_profile(args.profile)
    marker = read_marker(args.marker_stats, args.kept_genes)

    rows_rep = []
    for key in sorted(gappa.keys()):
        g = gappa[key]
        rows_rep.append(
            {
                "query_id": key,
                "assigned_rank": g["rank"],
                "gappa_lineage": g["lineage"],
                "lwr": round(clamp01(g["lwr"]), 6),
                "marker_kept": marker["kept"],
                "marker_total": marker["total"],
                "marker_kept_fraction": round(marker["kept_fraction"], 6),
            }
        )

    pd.DataFrame(
        rows_rep,
        columns=[
            "query_id",
            "assigned_rank",
            "gappa_lineage",
            "lwr",
            "marker_kept",
            "marker_total",
            "marker_kept_fraction",
        ],
    ).to_csv(args.out_report, sep="\t", index=False)

    if args.out_heterogeneity:
        het = heterogeneity_rows(gappa, args.mag_name, args.min_markers)
        pd.DataFrame(het, columns=HETEROGENEITY_COLUMNS).to_csv(
            args.out_heterogeneity, sep="\t", index=False
        )

    # Marker filtering must have run (total > 0) to judge the count; an empty
    # marker table means filtering was disabled and the warning does not apply.
    # A sentinel file (rather than text inside the report) keeps the report a
    # clean TSV; the flag's presence is the signal, its content the message.
    below_threshold = marker["total"] > 0 and marker["kept"] < args.min_markers
    if args.out_flag:
        if below_threshold:
            with open(args.out_flag, "w", encoding="utf-8") as fh:
                fh.write(
                    "WARNING: TOTAL MARKERS DETECTED IS BELOW THRESHOLD FOR ACCURATE CLASSIFICATION\n"
                )
        elif os.path.exists(args.out_flag):
            os.remove(args.out_flag)


if __name__ == "__main__":
    main()
