#!usr/bin/python3

import os
import argparse
import pandas as pd
from Bio import SeqIO
from Bio import SearchIO

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
gene = input_fasta.split("/")[-1].split(".")[0]
print("Gene: " + gene)
output_fasta = args.output
species = str(args.input).split("/")[1].replace("_working_dataset", "")
hits_path = "resources/" + species + "_fish_out/tmp/" + species + "/" + gene + ".hmmout"
print("HMM Path: " + hits_path)

# read in the Unique ID column 
# from input_metadata.tsv
# store as a list
mag_names = []

mag_f = os.listdir(str(args.directory))

def eval_filter(hits_path):
    hits = SearchIO.read(hits_path, "hmmer3-text")
    best_hit_id = hits[0].id
    best_hit_iter = 0
    for i in range(1, len(hits)):
        if hits[i].evalue < hits[best_hit_iter].evalue:
            best_hit_iter = i
            best_hit_id = hits[i].id
    print("Best Hit Iterator: " + "\t" + str(best_hit_iter))
    print("Best Hit ID: " + "\t" + best_hit_id)
    return best_hit_id

def splitter(input_fasta, output_fasta, mag_names, species, best_hit):
    print(input_fasta)
    fasta_dict = SeqIO.to_dict(SeqIO.parse(input_fasta, "fasta"))
    with open(output_fasta, "w") as out:
        for key, value in fasta_dict.items():
            if best_hit in key:
                print(key + " best hit is clear !!")
                value.id = species
                value.description = species
                SeqIO.write(value, out, "fasta")


best_hit = eval_filter(hits_path)
print("Best Hit: " + best_hit)

splitter(input_fasta, output_fasta, mag_names, species, best_hit)
