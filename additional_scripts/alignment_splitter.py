#!usr/bin/python

import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("-a", "--alignment", help="input fasta file")
parser.add_argument("-t", "--taxon", help="name of target taxon in input fasta file")


def split_mag_from_aln(input_fasta, taxon_name):
    records = list(SeqIO.parse(input_fasta, "fasta"))
    mag_records = []
    nonmag_records = []
    for record in records:
        if taxon_name in record.id:
            mag_records.append(record)
        else:
            nonmag_records.append(record)
    SeqIO.write(mag_records, taxon_name + "_q.aln", "fasta")
    SeqIO.write(nonmag_records, taxon_name + "_ref.aln", "fasta")

if __name__ == "__main__":
    # arg parsing
    args = parser.parse_args()
    split_mag_from_aln(args.alignment, args.taxon)
        