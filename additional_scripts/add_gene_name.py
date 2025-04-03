#!usr/bin/python

import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("-a", "--alignment", help="input fasta file")
parser.add_argument("-t", "--taxon", help="name of target taxon in input fasta file")
parser.add_argument("-g", "--gene", help="name of gene in input fasta file")
parser.add_argument("-o", "--output", help="output fasta file name and path")
args = parser.parse_args()


def add_gene_name_to_fasta(input_fasta, gene_name, output_fasta):
    records = list(SeqIO.parse(input_fasta, "fasta"))
    for record in records:
        if args.taxon in record.id:
            record.id = record.id + "_" + gene_name
            record.description =  record.description + "_" + gene_name
    SeqIO.write(records, output_fasta, "fasta")


add_gene_name_to_fasta(args.alignment, args.gene, args.output)