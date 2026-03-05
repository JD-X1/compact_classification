#!usr/bin/env python

import argparse
from Bio import SeqIO


def is_target_id(value, taxon_name):
    text = str(value).split()[0]
    return (
        text == taxon_name
        or text.startswith(f"{taxon_name}_")
        or text.startswith(f"{taxon_name}..")
        or text.split("|", 1)[0] == taxon_name
    )


def is_target_record(record, taxon_name):
    fields = [record.id, record.name, record.description]
    for field in fields:
        if field and is_target_id(field, taxon_name):
            return True
    return False


def split_mag_from_aln(input_fasta, taxon_name, output_dir, gene_name=None):
    records = list(SeqIO.parse(input_fasta, "fasta"))
    mag_records = []
    nonmag_records = []
    for record in records:
        if is_target_record(record, taxon_name):
            mag_records.append(record)
        else:
            nonmag_records.append(record)
    if gene_name == None:
        SeqIO.write(mag_records, output_dir + taxon_name + "_q.aln", "fasta")
        SeqIO.write(nonmag_records, output_dir + taxon_name + "_ref.aln", "fasta")
    else:
        SeqIO.write(mag_records, output_dir + gene_name + "_q.aln", "fasta")
        SeqIO.write(nonmag_records, output_dir + gene_name + "_ref.aln", "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--alignment", required=True, help="input fasta file")
    parser.add_argument("-t", "--taxon", required=True, help="name of target taxon in input fasta file")
    parser.add_argument("-g", "--gene", help="name of gene being targeted", default=None)
    parser.add_argument("-o", "--output_dir", help="output directory", default="/output/")
    args = parser.parse_args()
    split_mag_from_aln(args.alignment, args.taxon, args.output_dir, args.gene)
