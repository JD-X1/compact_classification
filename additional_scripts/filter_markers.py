#!/usr/bin/env python3

import argparse
import os
import sys
from typing import List, Tuple

from Bio import SeqIO


def base_name(path: str) -> str:
    name = os.path.basename(path)
    for suffix in (".trimal", ".aln.partial.fas", ".aln"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return os.path.splitext(name)[0]


def is_target_id(value: str, taxon: str) -> bool:
    text = str(value).split()[0]
    return (
        text == taxon
        or text.startswith(f"{taxon}_")
        or text.startswith(f"{taxon}..")
        or text.split("|", 1)[0] == taxon
    )


def find_query_record(records, taxon: str):
    for rec in records:
        for field in (rec.id, rec.name, rec.description):
            if field and is_target_id(field, taxon):
                return rec
    return None


def metrics(seq: str) -> Tuple[int, int, int, int]:
    total = len(seq)
    if total == 0:
        return 0, 0, 0, 0
    upper = seq.upper()
    gap_like = sum(1 for c in upper if c in {"-", "?", ".", "~", "X"})
    present = total - gap_like
    informative = sum(1 for c in upper if c not in {"-", "?", ".", "~", "X"})
    return total, gap_like, present, informative


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--taxon", required=True)
    ap.add_argument("--inputs", nargs="+", required=True)
    ap.add_argument("--out-keep", required=True)
    ap.add_argument("--out-stats", required=True)
    ap.add_argument("--min-coverage", type=float, default=0.3)
    ap.add_argument("--max-gap-fraction", type=float, default=0.7)
    ap.add_argument("--min-informative-sites", type=int, default=20)
    ap.add_argument("--min-markers", type=int, default=10)
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out_keep) or ".", exist_ok=True)
    os.makedirs(os.path.dirname(args.out_stats) or ".", exist_ok=True)

    kept: List[str] = []
    rows: List[Tuple[str, int, int, int, float, float, int]] = []

    for path in sorted(args.inputs):
        gene = base_name(path)
        records = list(SeqIO.parse(path, "fasta"))
        query = find_query_record(records, args.taxon)
        if query is None:
            rows.append((gene, 0, 0, 0, 0.0, 1.0, 0))
            continue

        total, gap_like, present, informative = metrics(str(query.seq))
        if total == 0:
            coverage = 0.0
            gap_fraction = 1.0
        else:
            coverage = present / total
            gap_fraction = gap_like / total

        keep = (
            coverage >= args.min_coverage
            and gap_fraction <= args.max_gap_fraction
            and informative >= args.min_informative_sites
        )
        if keep:
            kept.append(gene)
        rows.append((gene, total, present, gap_like, coverage, gap_fraction, informative))

    with open(args.out_keep, "w", encoding="utf-8") as out_keep:
        for gene in kept:
            out_keep.write(f"{gene}\n")

    with open(args.out_stats, "w", encoding="utf-8") as out_stats:
        out_stats.write(
            "gene\taln_len\tpresent_sites\tgap_like_sites\tcoverage\tgap_fraction\tinformative_sites\n"
        )
        for row in rows:
            out_stats.write(
                f"{row[0]}\t{row[1]}\t{row[2]}\t{row[3]}\t{row[4]:.6f}\t{row[5]:.6f}\t{row[6]}\n"
            )

    if len(kept) < args.min_markers:
        sys.stderr.write(
            f"WARNING: only {len(kept)} markers passed filter for {args.taxon}; "
            f"below the recommended minimum of {args.min_markers}. Proceeding with "
            f"concatenation, but this placement should be treated as low-confidence.\n"
        )
    if not kept:
        sys.stderr.write(
            f"WARNING: no markers passed filter for {args.taxon}; "
            f"wrote empty keep list to route this MAG to UNCLASSIFIABLE.\n"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
