#!usr/bin/env python3

import os
from Bio import Entrez
import pandas as pd
Entrez.email = 'jduque2@ucmerced.edu'
# check if file tax_IDs.txt exists
# tax_IDs.txt has two columns
# column 1: unique ID
# column 2: taxonomic ID
# if file doesn't exist then throw an error
# else read and store both columns in a dictionary
# dictionary key: unique ID
# dictionary value: taxonomic ID

if os.path.isfile('tax_IDs.txt'):
    tax_ids = pd.read_csv('tax_IDs.txt', sep='\t', header=None)
    tax_ids.columns = ['unique_ID', 'tax_ID']
    tax_IDs = tax_ids['tax_ID'].tolist()
    unique_IDs = tax_ids['unique_ID'].tolist()
    id_dict = dict(zip(unique_IDs, tax_IDs))
else:
    raise FileNotFoundError('tax_ids.txt not found in current directory')

# loop through the taxonomic IDs and use them to construct the lineage tree file
# lineage_tree.txt this file has two columns
# output taxonomy tree file tab-delimited with two columns
#
# Column 1: Unique ID for each species
# Column 2: Taxonomy tree for each species, string separated by semicolons

# create a list to store the lineage trees
lineage_trees = []

# check if file tax_tree.txt exists
# if file exists then read it in as a dataframe
# get the unique IDs from the first column and store them in a list
if os.path.isfile('tax_tree.txt'):
    tax_tree = pd.read_csv('tax_tree.txt', sep='\t', header=None)
    tax_tree_IDs = tax_tree[0].tolist()
    # now we can remove any taxonomic IDs from tax_tree.txt by removing them from the list
    # tax_IDs
    for tax_tree_ID in tax_tree_IDs:
        if tax_tree_ID in unique_IDs:
            unique_IDs.remove(tax_tree_ID)

    print("number of completed species: ", len(tax_tree_IDs))
    print("number of remaining species: ", len(unique_IDs))

for i, unique_id in enumerate(unique_IDs):
    handle = Entrez.efetch(db='taxonomy', id=id_dict[unique_id], retmode='xml')
    records = Entrez.read(handle)
    handle.close()
    lineage = records[0]['Lineage']
    lineage = lineage.split(';')[1:]
    lineage = [x.strip() for x in lineage]
    # join lineage with semicolons
    lineage = ';'.join(lineage) 
    line = unique_IDs[i] + '\t' + lineage
    if records[0]['ScientificName'] != None:
        line = line + ';' + records[0]['ScientificName']
    with open('tax_tree.txt', 'a') as f:
        f.write(line + '\n')


