#!usr/bin/python3

import os
import argparse
import pandas as pd
from Bio import SeqIO
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
parser.add_argument("-o", "--output", help="outputfasta file name and path")
args = parser.parse_args()

# check if directory resources/q_frags/
# exists, if not create it
if not os.path.exists("resources/q_frags/"):
    os.makedirs("resources/q_frags/")

# read in the input fasta file
input_fasta = args.input
output_fasta = args.output

# read in the Unique ID column 
# from input_metadata.tsv
# store as a list
metadata = pd.read_csv("resources/input_metadata.tsv", sep="\t")
mag_names = metadata["Unique ID"].tolist()

"""
turn the below into a function
# read in the fasta file as a dictionary

fasta_dict = SeqIO.to_dict(SeqIO.parse(input_fasta, "fasta"))

# loop through our mag names
# and from the sequence dictionary
# extract all sequences that contain
# the Unique ID in their name
# write these all out to the output file
with open(output_fasta, "w") as out:
    for mag in mag_names:
        for key, value in fasta_dict.items():
            if mag in key:
                value.id = mag
                value.description = mag
                SeqIO.write(value, out, "fasta")
"""

def splitter(input_fasta, output_fasta, mag_names):
    print(input_fasta)
    fasta_dict = SeqIO.to_dict(SeqIO.parse(input_fasta, "fasta"))
    with open(output_fasta, "w") as out:
        for mag in mag_names:
            for key, value in fasta_dict.items():
                if mag in key:
                    value.id = mag
                    value.description = mag
                    SeqIO.write(value, out, "fasta")

splitter(input_fasta, output_fasta, mag_names)

