#!usr/bin/env python3

import os
import re
import argparse
import gzip
import pandas as pd
from copy import deepcopy
from pathlib import Path
from Bio import SeqIO, SearchIO
from typing import Iterable, List, Optional, Tuple

"""
Usage:
python ./additional_scripts/splitter.py -i {input} -o {output}"

where -i is an fasta file containing
a the sequences to be split and -o is the
out file name for the split sequences.

we are going to have to check the mag names
from input_metadata.tsv and use the file
name to identify the sequences we want to
extract. As a rule these are appended to
the end of the input sequence, but we want
to rename the sequences to the mag name
where necessary.

"""

FA_EXTS = (".fa", ".fasta", ".fas", ".faa", ".aln")

def infer_context(input_fa: Path, outdir_opt: Optional[str]) -> Tuple[str, str, Path]:
    name = input_fa.name
    
    for ext in sorted(FA_EXTS, key=len, reverse=True):
        if name.endswith(ext):
            name=name[:-len(ext)]
            break
    gene = name
    mag = None
    for suf in ("_mafft_out", "_working_dataset", "_fish_out"):
        if input_fa.parent.name.endswith(suf):
            mag = input_fa.parent.name[:-len(suf)]
            break 
    if not mag:
        raise SystemExit(f"Could not infer MAG ID to {input_fa}")

    outdir = None
    for parent in [input_fa.parent] + list(input_fa.parents):
        if (parent / f"{mag}_fish_out").exists() or (parent / f"{mag}_working dataset").exists():
            outdir = parent
            break
        if outdir is None: 
            outdir = input_fa.parent.parent.resolve()
    return mag, gene, outdir

def find_hmmout(outdir: Path, mag: str, gene: str) -> Optional[Path]:
    base = outdir / f"{mag}_fish_out" / "tmp" / mag
    p = base / f"{gene}.hmmout"
    if p.is_file():
        return p
    for cand in sorted(base.glob(f"{gene}*")):
        if cand.is_file():
            return cand
    return None

def pick_best_hit(hmm_path: Path, input_ids: List[str]) -> Optional[str]:
    for fmt in ("hmmer3-text", "hmmer3-domtab", "hmmer3-tab"):
        try:
            q = SearchIO.read(str(hmm_path), fmt)
            break
        except Exception:
            q = None
    if q is None:
        return None
    best_id, best_ev = None, None
    for hit in q:
        if not any(hit.id in i or i in hit.id for i in input_ids):
            continue
        evs = [hsp.evalue for hsp in hit.hsps if getattr(hsp, "evalue", None) is not None]
        if not evs:
            continue
        ev = min(evs)
        if best_ev is None or ev < best_ev:
            best_ev, best_id = ev, hit.id
    if best_id is None and len(q) > 0:
        return q[0].id
    return best_id

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="input fasta file")
    parser.add_argument("-d", "--directory", help="directory containing mag files")
    parser.add_argument("-o", "--output", help="outputfasta file name and path")
    parser.add_argument("-r", "--reference_output", help="reference output fasta file name and path")
    args = parser.parse_args()

    input_fa = Path(args.input).resolve()
    out_fa = Path(args.output).resolve()
    ref_fa = Path(args.reference_output).resolve() if args.reference_output else None

    mag, gene, outdir = infer_context(input_fa, args.directory)
    print(f"Identified MAG: {mag}, gene: {gene}, outdir: {outdir}")
    opener = gzip.open if str(input_fa).endswith(".gz") else open

    with opener(input_fa, "rt") as fh:
        records = list(SeqIO.parse(fh, "fasta"))
    if not records:
        raise SystemExit(f"No sequences found in {input_fa}")
    input_ids = [r.id for r in records]

    hmm_path = find_hmmout(outdir, mag, gene)

    if hmm_path is None:
        raise SystemExit(f"HMMER output not found for: {outdir}/{mag}_fish_out/tmp/{mag}/{gene}*.hmmout")
    best = pick_best_hit(hmm_path, input_ids) or records[0].id

    out_fa.parent.mkdir(parents=True, exist_ok=True)
    ref_out = Path(args.reference_output).resolve() if args.reference_output else None
    if ref_out:
        ref_out.parent.mkdir(parents=True, exist_ok=True)

    query_index = None
    for idx, rec in enumerate(records):
        if best in rec.id or rec.id in best:
            query_index = idx
            break
        if query_index is None:
            query_index = 0

        query_record = deepcopy(records[query_index])
        query_record.id = mag
        query_record.description = mag

        with open(out_fa, "w") as out:
            SeqIO.write(query_record, out, "fasta")

        if ref_out:
            reference_records = [records[idx] for idx in range(len(records)) if idx != query_index]
            with open(ref_out, "w") as ref_handle:
                SeqIO.write(reference_records, ref_handle, "fasta")

if __name__ == "__main__":
    main()