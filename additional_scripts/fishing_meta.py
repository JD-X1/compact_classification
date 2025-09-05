#!usr/bin/env python3
import os
import sys
import argparse
from pathlib import Path
from typing import Tuple, Optional, List

PROT_EXTS = (".faa", ".fa", ".fasta", ".faa.gz", ".fa.gz", ".fasta.gz")

def sanatize_id(s: str) -> str:
    for bad, rep in [("_",""),("@",""),("..",""),(" ",""),("*","")]:
        s = s.replace(bad, rep)
    return s

def find_proteome_in_dir(d: Path) -> Optional[Path]:
    cand = d / "translated_protein.fasta"
    if cand.is_file():
        return cand
    cand = d / "eukaryota_odb12" / "translated_protein.fasta"
    if cand.is_file():
        return cand
    
    for ext in PROT_EXTS:
        for p in d.glob(f"*{ext}"):
            if p.is_file():
                return p

    for p in d.rglob("translated_protein.fasta"):
        if p.is_file():
            return p
    return None

def resolve_location_and_file(p: Path) -> Tuple[Path, str]:
    if p.is_file():
        return p.parent, p.name
    if p.is_dir():
        found = find_proteome_in_dir(p)
        if found:
            return (found.parent, found.name)
        raise FileNotFoundError(f"No proteome file found in directory {p}")
    raise FileNotFoundError(f"Path does not exist: {p}")

def infer_unique_id(location_dir: Path) -> str:
    base = location_dir.name
    if base == "eukaryota_odb12":
        mag = location_dir.parent.name
    else:
        mag = base
    return sanatize_id(mag)

def main():
    parser = argparse.ArgumentParser(
        description="Create input_metadata.tsv lines for PhyloFisher from proteome output"
    )
    parser.add_argument(
        "-p", "--path",
        nargs="+",
        required=True,
        help="One or more paths to proteome files or directories containing them."
    )
    parser.add_argument(
        "-o", "--out",
        required=True,
        help="Output path for the input_metadata.tsv file"
    )
    args = parser.parse_args()

    out_path = Path(args.out)
    if out_path.parent.as_posix() != "":
        out_path.parent.mkdir(parents=True, exist_ok=True)
    
    lines: List[str] = []
    for raw in args.path:
        p = Path(raw).resolve()
        loc_dir, file_name = resolve_location_and_file(p)
        uid = infer_unique_id(loc_dir)
        if not uid:
            uid = sanatize_id(Path(file_name).stem)

        higher_taxonomy = "Unknownazoa"
        lower_taxonomy = "Unknownsea"
        blast_seed = "none"
        long_name = "Unknownus nameus"
        data_type = "MAG"
        source = "SequencerProbably"

        out_line = "\t".join([
            str(loc_dir),
            file_name,
            uid,
            higher_taxonomy,
            lower_taxonomy,
            blast_seed,
            long_name,
            data_type,
            source
        ])
        lines.append(out_line)

    with open(out_path, "w") as fh:
        for line in lines:
            print(line)
            fh.write(line + "\n")

if __name__ == "__main__":
    main()
