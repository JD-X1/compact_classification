#! usr/bin/env python3

import os
import Bio

# get list of alignments in directory mafft that end in .final
alignments = os.listdir('resources/q_alns/')
alignments = [x for x in alignments if x.endswith('.aln')]

# get list of taxa in each single gene alignment and create a new file that takes lines from
# tax_tree.txt that match the taxa in the alignment

for alignment in alignments:
    taxa = []
    protein = alignment.split('.')[0]
    with open('resources/ref_alns/' + protein + '.fas.aln', 'r') as f:
        for line in f:
            if line.startswith('>'):
                taxa.append(line.strip().replace('>', ''))
    for taxon in taxa:
        with open('resources/tax_tree.txt', 'r') as f:
            for line in f:
                if taxon in line:
                    with open('tax_tree_' + protein + '.txt', 'a') as g:
                        g.write(line)
            