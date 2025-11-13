#!usr/bin/env python

import argparse
from Bio import SeqIO
import dendropy as dp
import argparse

def filter_tree_by_alignment(tree_file, alignment_file, output_tree_file, purge_tip=None):
    # read in with SeqIO and extract list of taxa
    alignment = SeqIO.to_dict(SeqIO.parse(alignment_file, "fasta"))
    aln_taxa = set(alignment.keys())
    # read in the tree
    tree = dp.Tree.get_from_path(tree_file, "newick")
    tree_taxa = tree.taxon_namespace.labels()
    to_prune = []
    for taxon in tree_taxa:
        if taxon not in aln_taxa:
            to_prune.append(taxon)
    if purge_tip:
        to_prune.append(purge_tip)
    # prune the tree
    tree.prune_taxa_with_labels(to_prune)
    # write the pruned tree to the output file
    with open(output_tree_file, "w") as out_file:
        out_file.write(tree.as_string("newick"))

def main():
    parser = argparse.ArgumentParser(description="Filter a tree to only include taxa present in an alignment.")
    parser.add_argument("-a", "--alignment", required=True, help="Path to the alignment file (in FASTA format).")
    parser.add_argument("-t", "--tree", required=True, help="Path to the tree file (in Newick format).")
    parser.add_argument("-p", "--purge", required=False, default=None, help="Tip name to purge.")
    parser.add_argument("-o", "--output", required=True, help="Path to the output tree file.")
    
    args = parser.parse_args()

    filter_tree_by_alignment(args.tree, args.alignment, args.output, args.purge)

if __name__ == "__main__":
    main()