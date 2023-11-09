from Bio import SeqIO

from collections import defaultdict
import os
import glob

def get_taxa_from_fasta(fasta_file):
    """Return a set of taxa names from a given FASTA file."""
    with open(fasta_file, "r") as handle:
        taxa = {record.id for record in SeqIO.parse(handle, "fasta")}
    return taxa

def filter_and_write_fasta(input_file, output_file, taxa_to_keep):
    """Filter the input FASTA based on taxa_to_keep and write to output_file."""
    with open(input_file, "r") as handle, open(output_file, "w") as out_handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id in taxa_to_keep:
                SeqIO.write(record, out_handle, "fasta")

def main():
    # List of FASTA files
    fasta_files = glob.glob("resources/ref_alns/*.final")

    # Get taxa from each FASTA and find the intersection
    all_taxa_sets = [get_taxa_from_fasta(f) for f in fasta_files]
    taxa_to_keep = set.intersection(*all_taxa_sets)
    # reduce taxa to keep to 30
    taxa_to_keep = set(list(taxa_to_keep)[:30])
    # Filter each FASTA and write to new file
    for f in fasta_files:
        output_file = os.path.join("resources/ref_alns_reduce", os.path.basename(f))
        print(output_file)
        filter_and_write_fasta(f, output_file, taxa_to_keep)
    
if __name__ == "__main__":
    main()
