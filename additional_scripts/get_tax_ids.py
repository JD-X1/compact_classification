#!usr/bin/env python3

import os
from Bio import Entrez
import pandas as pd

# this script is designed to construct a taxonomy tree for the species in the
# the file ../PhyloFishScratch/metadata.csv
# 
# metadata.tsv has two columns we're interested in:
# Unique ID: one per species
# Lower Taxonomy: lowest taxonomic label of interest (e.g. fungi, etc.)
# 
# output taxonomy tree file tab-delimited with two columns
# Column 1: Unique ID for each species
# Column 2: Taxonomy tree for each species, string separated by semicolons
#
#
# We will use the NCBI taxonomy database to construct the tree using the Lower Taxonomy
# as a search term. We will then use the NCBI taxonomy database to construct the tree
# from lineage informayion

# read in metadata.csv which has a header
metadata = pd.read_csv('../PhyloFishScratch/metadata.tsv', sep='\t', header=0)

# create a list of unique IDs
unique_IDs = metadata['Unique ID'].tolist()
print("number of unique IDs: ", len(unique_IDs))
# create a list of long names
long_names = metadata['Long Name'].tolist()
print("number of long names: ", len(long_names))

# create a list of lower taxonomy labels
lower_tax = metadata['Lower Taxonomy'].tolist()
print("number of lower taxa: ", len(lower_tax))
# loop through the lower taxonomy labels and use them as search terms in the NCBI taxonomy database
# to get the taxonomic IDs for each species
Entrez.email = 'jduque2@ucmerced.edu'

# create a list to store the taxonomic IDs
tax_IDs = []
lower_taxa_error = []
# print total number of species

for i in range(len(lower_tax)):
    # search the NCBI taxonomy database for the lower taxonomy label
    # first we will try to search for the long name label
    handle = Entrez.esearch(db='taxonomy', term=long_names[i])
    record = Entrez.read(handle)
    handle.close()
    # check if the search returned a record
    # if it doesn't then we will search for the lower taxonomy label
    if len(record['IdList']) > 0:
        tax_IDs.append(record['IdList'][0])
    elif len(record['IdList']) == 0:
        handle = Entrez.esearch(db='taxonomy', term=lower_tax[i])
        record = Entrez.read(handle)
        handle.close()
        if len(record['IdList']) > 0:
            tax_IDs.append(record['IdList'][0])
        else:
            tax_IDs.append('WARNING')
            lower_taxa_error.append(lower_tax[i])
    line = unique_IDs[i] + '\t' + tax_IDs[i] + '\n'
    with open('tax_IDs.txt', 'a') as f:
        f.write(line)

print('total tax Ids ' + len(tax_IDs))

# write out unique IDs that didn't have a taxonomic ID
with open('lower_taxa_error.txt', 'w') as f:
    for taxa in lower_taxa_error:
        f.write(taxa + '\n')

