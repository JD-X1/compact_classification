#!usr/bin/python3

import os
import pandas as pd

# arguements needed
# -i input directory
#   - directory containing unaligned single gene fasta
#     files
# -t taxanomic information
#    - two column tab delimited file we only need the first
#    - first column is taxon unique ID
# -o output directory

# acquire single gene unaligned fasta files
proteins = os.listdir('resources/')

# check if output directory exists