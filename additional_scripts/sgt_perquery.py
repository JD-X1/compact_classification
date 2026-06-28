#!/usr/bin/env python
"""Reshape per-gene SGT best-hit profiles into a per-query table for taxonomy_report.

In SGT mode each marker gene is one query and its best-hit taxopath is that marker's
vote. taxonomy_report.read_profile wants a table with a query/name column, an LWR
column, and a taxopath; the per-gene profile.tsv has LWR + taxopath but no name. This
tags each gene's row with the gene id and concatenates them, so taxonomy_report can
compute the between-marker dispersion D1. Unplaceable genes (header-only sentinels from
sgt_place_guard) carry no data row and so cast no vote, matching the off-ramp semantics.
"""
import argparse
import os

import sgt_place_guard


def gene_from_profile_path(path):
    return os.path.basename(os.path.dirname(path))


def assemble_per_query(profile_paths):
    """Return [(gene, lwr, taxopath), ...] over placeable genes, skipping unplaceable
    sentinels. Columns are resolved by header name, not position."""
    rows = []
    for path in profile_paths:
        if sgt_place_guard.is_unplaceable_profile(path):
            continue
        gene = gene_from_profile_path(path)
        with open(path) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            try:
                lwr_i = header.index("LWR")
                tax_i = header.index("taxopath")
            except ValueError:
                continue
            for line in fh:
                if not line.strip():
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) <= max(lwr_i, tax_i):
                    continue
                rows.append((gene, f[lwr_i], f[tax_i]))
    return rows


def write_per_query_tsv(profile_paths, out_path):
    rows = assemble_per_query(profile_paths)
    with open(out_path, "w") as fh:
        fh.write("name\tLWR\ttaxopath\n")
        for gene, lwr, taxopath in rows:
            fh.write(f"{gene}\t{lwr}\t{taxopath}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--profiles", nargs="*", default=[],
                    help="every per-gene profile.tsv for the MAG")
    ap.add_argument("--out", required=True, help="per-query table for taxonomy_report")
    args = ap.parse_args()
    write_per_query_tsv(args.profiles, args.out)


if __name__ == "__main__":
    main()
