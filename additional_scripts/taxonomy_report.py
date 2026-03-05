#!/usr/bin/env python3

import argparse
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


def lineage_lca(lineages: List[str]) -> str:
    values = [split_lineage(x) for x in lineages if clean_value(x)]
    if not values:
        return ""
    out = []
    m = min(len(x) for x in values)
    for i in range(m):
        col = {x[i] for x in values}
        if len(col) != 1:
            break
        out.append(values[0][i])
    return join_lineage(out)


def lineages_agree(a: str, b: str) -> Tuple[bool, str]:
    left = split_lineage(a)
    right = split_lineage(b)
    if not left or not right:
        return False, ""
    m = min(len(left), len(right))
    if left[:m] == right[:m]:
        return True, join_lineage(left[:m])
    return False, lineage_lca([a, b])


def read_tax_tree(path: str) -> Dict[str, str]:
    out = {}
    if not path or not os.path.exists(path):
        return out
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            out[clean_value(parts[0])] = clean_value(parts[1])
    return out


def pick_query_col(df: pd.DataFrame) -> str:
    names = ["query", "query_name", "query_names", "name", "placement_name", "sample"]
    cols = {c.lower(): c for c in df.columns}
    for name in names:
        if name in cols:
            return cols[name]
    return df.columns[0]


def pick_score_col(df: pd.DataFrame) -> str:
    names = ["like_weight_ratio", "likelihood_weight_ratio", "lwr", "probability", "score", "support"]
    cols = {c.lower(): c for c in df.columns}
    for name in names:
        if name in cols:
            return cols[name]
    return ""


def normalize_score(value) -> float:
    try:
        x = float(value)
    except Exception:
        return 0.0
    if x < 0:
        return 0.0
    if x > 1 and x <= 100:
        return x / 100.0
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
        gscore = normalize_score(row[scol]) if scol else 0.0
        out[q] = {"lineage": lineage, "rank": rank, "gappa_score": gscore}
    return out


def read_pairwise(path: str, tax_map: Dict[str, str], top_k: int) -> Dict[str, Dict]:
    if not os.path.exists(path):
        return {}
    df = pd.read_csv(path, sep="\t")
    if df.empty or "query_names" not in df.columns or "distance" not in df.columns:
        return {}
    out = {}
    for q, sub in df.groupby("query_names"):
        curr = sub.sort_values("distance", ascending=True).reset_index(drop=True)
        top = curr.head(max(1, top_k))
        d1 = float(top.iloc[0]["distance"])
        d2 = float(top.iloc[1]["distance"]) if len(top) > 1 else d1
        margin = (d2 - d1) / (d2 + 1e-9) if d2 >= 0 else 0.0
        margin = max(0.0, min(1.0, margin))
        leaves = [clean_value(x) for x in top["leaf_name"].tolist() if clean_value(x)]
        lineages = [tax_map[x] for x in leaves if x in tax_map]
        out[clean_value(q)] = {
            "top1": d1,
            "margin": margin,
            "nearest_leaves": leaves,
            "similarity_lineage": lineage_lca(lineages) if lineages else "",
        }
    return out


def read_marker(marker_stats: str, kept_genes: str) -> Dict:
    out = {
        "total": 0,
        "kept": 0,
        "kept_fraction": 0.0,
        "mean_coverage": 0.0,
        "mean_gap_fraction": 1.0,
        "marker_score": 0.0,
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

    m1 = out["kept_fraction"]
    m2 = max(0.0, min(1.0, out["mean_coverage"]))
    m3 = max(0.0, min(1.0, 1.0 - out["mean_gap_fraction"]))
    out["marker_score"] = (m1 + m2 + m3) / 3.0
    return out


def tier(score: float, high: float, medium: float) -> str:
    if score >= high:
        return "high"
    if score >= medium:
        return "medium"
    return "low"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--profile", required=True)
    ap.add_argument("--pairwise", required=True)
    ap.add_argument("--tax-tree", required=True)
    ap.add_argument("--out-confidence", required=True)
    ap.add_argument("--out-consensus", required=True)
    ap.add_argument("--out-report", required=True)
    ap.add_argument("--marker-stats", default="")
    ap.add_argument("--kept-genes", default="")
    ap.add_argument("--top-k", type=int, default=5)
    ap.add_argument("--high-cutoff", type=float, default=0.70)
    ap.add_argument("--medium-cutoff", type=float, default=0.45)
    args = ap.parse_args()

    tax_map = read_tax_tree(args.tax_tree)
    gappa = read_profile(args.profile)
    pair = read_pairwise(args.pairwise, tax_map, args.top_k)
    marker = read_marker(args.marker_stats, args.kept_genes)

    rows_conf = []
    rows_cons = []
    rows_rep = []
    keys = sorted(set(gappa.keys()) | set(pair.keys()))

    for key in keys:
        g = gappa.get(key, {"lineage": "", "rank": "unclassified", "gappa_score": 0.0})
        p = pair.get(key, {"top1": 0.0, "margin": 0.0, "nearest_leaves": [], "similarity_lineage": ""})

        gscore = max(0.0, min(1.0, float(g["gappa_score"])))
        dscore = (1.0 / (1.0 + max(0.0, float(p["top1"]))) + float(p["margin"])) / 2.0
        mscore = max(0.0, min(1.0, float(marker["marker_score"])))
        score = 0.40 * gscore + 0.35 * dscore + 0.25 * mscore
        confidence = tier(score, args.high_cutoff, args.medium_cutoff)

        agree, shared = lineages_agree(g["lineage"], p["similarity_lineage"])
        agreement = "agree" if agree else ("insufficient" if not shared else "disagree")
        fallback_rank = g["rank"] if g["lineage"] else "superkingdom"

        rows_conf.append(
            {
                "query_id": key,
                "confidence_score": round(score, 6),
                "confidence_tier": confidence,
                "gappa_score": round(gscore, 6),
                "distance_score": round(dscore, 6),
                "distance_margin": round(float(p["margin"]), 6),
                "marker_score": round(mscore, 6),
            }
        )
        rows_cons.append(
            {
                "query_id": key,
                "gappa_lineage": g["lineage"],
                "similarity_lineage": p["similarity_lineage"],
                "agreement": agreement,
                "shared_prefix": shared,
            }
        )
        rows_rep.append(
            {
                "query_id": key,
                "assigned_rank": g["rank"],
                "fallback_rank": fallback_rank,
                "gappa_lineage": g["lineage"],
                "similarity_lineage": p["similarity_lineage"],
                "agreement": agreement,
                "confidence_score": round(score, 6),
                "confidence_tier": confidence,
                "marker_kept": marker["kept"],
                "marker_total": marker["total"],
                "marker_kept_fraction": round(marker["kept_fraction"], 6),
                "distance_top1": round(float(p["top1"]), 6),
                "distance_margin": round(float(p["margin"]), 6),
                "nearest_leaves_top3": ";".join(p["nearest_leaves"][:3]),
            }
        )

    pd.DataFrame(rows_conf).to_csv(args.out_confidence, sep="\t", index=False)
    pd.DataFrame(rows_cons).to_csv(args.out_consensus, sep="\t", index=False)
    pd.DataFrame(rows_rep).to_csv(args.out_report, sep="\t", index=False)


if __name__ == "__main__":
    main()
