#!usr/bin/env python3

import os
import argparse
import pandas as pd
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

FA_EXT = (".fa", ".fasta", ".fas", ".faa", ".aln")

def open_text(path: Path):
    if str(os.path).endswith(".gz"):
        return gzip.open(os.path, "rt")
    return open(path, "r")

def strip_known_ext(name: str) -> str:
    for ext in sorted(FA_EXT, key=len, reverse=True):
        if name.endswith(ext):
            return name[: -len(ext)]
    return name

def infer_mag_from_header(fasta_path: Path) -> Optional[str]:
    for parent in [p.parent] + list(p.parents):
        nm = parent.name
        for suf in ("_working_dataset", "_fish_out", "_mafft_out"):
            if nm.endswith(suf):
                return nm[: -len(suf)]
    return None



# get the input and output file names
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="input fasta file")
parser.add_argument("-d", "--directory", help="directory containing mag files")
parser.add_argument("-o", "--output", help="outputfasta file name and path")
args = parser.parse_args()

# check if directory resources/q_frags/
# exists, if not create it
# if not os.path.exists("resources/q_frags/"):
#    os.makedirs("resources/q_frags/")

# read in the input fasta file
input_fasta = args.input
print(input_fasta)
output_path = "/"+ str(input_fasta).split("/")[1] + "/"
print("Output Path: " + output_path)
gene = input_fasta.split("/")[-1].split(".")[0]
print("Gene: " + gene)
output_fasta = args.output
print("Output Fasta: " + output_fasta)
species = str(args.input).split("/")[-2].replace("_working_dataset", "")
hits_path = output_path + species + "_fish_out/tmp/" + species + "/" + gene + ".hmmout"
print(hits_path)
print("HMM Path: " + hits_path)

# read in the Unique ID column 
# from input_metadata.tsv
# store as a list

def get_valid_keys(input_fasta):
    fasta_dict = SeqIO.to_dict(SeqIO.parse(input_fasta, "fasta"))
    valid_keys = []
    for key, value in fasta_dict.items():
        if species in key:
            valid_keys.append(key.split("_HMM_")[0])
    return valid_keys

def eval_filter(hits_path, input_fasta):
    hits = SearchIO.read(hits_path, "hmmer3-text")
    print(hits.id)
    valid_keys = get_valid_keys(input_fasta)
    best_hit_id = valid_keys[0]
    best_hit_iter = 0
    for i in range(1, len(hits)):
        if hits[i].evalue < hits[best_hit_iter].evalue:
            if any(hits[i].id in vkey for vkey in valid_keys):
                best_hit_iter = i
                best_hit_id = hits[i].id
    print("Best Hit Iterator: " + "\t" + str(best_hit_iter))
    print("Best Hit ID: " + "\t" + best_hit_id)
    return best_hit_id

def splitter(input_fasta, output_fasta, species, best_hit):
    fasta_dict = SeqIO.to_dict(SeqIO.parse(input_fasta, "fasta"))
    valid_keys = get_valid_keys(input_fasta)
    print("Valid Keys: ")
    print(valid_keys)
    with open(output_fasta, "w") as out:
        for key, value in fasta_dict.items():
            if best_hit in key:
                if any(best_hit in vkey for vkey in valid_keys):
                    print(key + " best hit is clear !!")
                    value.id = species
                    value.description = species
                    SeqIO.write(value, out, "fasta")



best_hit = eval_filter(hits_path, input_fasta)
print("Best Hit: " + best_hit)

splitter(input_fasta, output_fasta, species, best_hit)
