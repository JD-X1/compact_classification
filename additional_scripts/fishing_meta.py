#!usr/bin/python3

import os
import sys
import pandas as pd
import subprocess

# Script to create the input_metadata.tsv file to
# run fisher.py

# Input:
# 1. Path to the directory containing the compleasm
#    output directories
# Output:
# 1. tab-delimited file containing the following
#    columns:
#    a. Location of the proteome file
#    b. File Name
#    c. Unique ID
#        Unique IDs cannot contain underscores “_”,
#        at symbols “@”, double dots “..”, white
#        space, or asterisk “*” 
#    d. Higher Taxonomy - always "Amoebozoa"
#    e. Lower Taxonomy - always "Discosea"
#    f. Blast Seed - always "none"
#    g. Long Name - always "Janus doeus"
#    h. Data Type - always "MAG"
#    i. Source - always "SequencerProbably"

# Usage:
# python3 fishing_meta.py <path_to_compleasm_output_dir>

# get the path to the compleasm output directory
path_to_compleasm_output_dir = sys.argv[1]
comp_cores = sys.argv[2]
working_set = sys.argv[3]

# get the list of directories in the compleasm output directory
dir_list = os.listdir(path_to_compleasm_output_dir)

# build empty dataframe with the appropriate columns
df = pd.DataFrame(columns = ['Location',
                             'File Name',
                             'Unique ID',
                             'Higher Taxonomy',
                             'Lower Taxonomy',
                             'Blast Seed',
                             'Long Name',
                             'Data Type',
                             'Source'])

# iterate through the directories

for dir in dir_list:

    # get the path to the proteome file
    path_to_proteome_file = os.path.join(path_to_compleasm_output_dir,
                                         dir,
                                         'eukaryota_odb10/')
    
    # get the file name
    file_name = "translated_protein.fasta"

    # get the unique ID remove the illegal
    # characters:
        # underscores “_”
        # at symbols “@”
        # double dots “..”
        # white spaces
        # asterisk “*”
    unique_id = dir.replace("_", "-")
    unique_id = unique_id.replace("@", "-")
    unique_id = unique_id.replace("..", "-")
    unique_id = unique_id.replace(" ", "-")
    unique_id = unique_id.replace("*", "-")

    # get the higher taxonomy
    higher_taxonomy = "Amoebozoa"

    # get the lower taxonomy
    lower_taxonomy = "Discosea"

    # get the blast seed
    blast_seed = "none"

    # get the long name
    long_name = "Janus doeus"

    # get the data type
    data_type = "MAG"

    # get the source
    source = "SequencerProbably"

    # add the row to the dataframe
    df.loc[len(df.index)] = [path_to_proteome_file,
                             file_name,
                             unique_id,
                             higher_taxonomy,
                             lower_taxonomy,
                             blast_seed,
                             long_name,
                             data_type,
                             source]
    
# write the dataframe to a tab-delimited file
df.to_csv("resources/input_metadata.tsv", sep = "\t", index = False)


# check for for directory:
if not os.path.exists("resources/PhyloFishScratch"):
    os.system("cp -r resources/PhyloFisherDatabase_v1.0/database resources/PhyloFishScratch")

# run config.py command
configuration_line = "config.py -d resources/PhyloFishScratch -i resources/input_metadata.tsv"

subprocess.run(configuration_line)

# run fisher.py command with the number of cores
#fisher_line =  "fisher.py --output resources/fisher_out --threads " + comp_cores

#os.system(fisher_line)

# run the informant command

#informant_line = "informant.py -i ./resources/fisher_out"
#os.system(informant_line)

# run the working database command

#working_database_line = "working_dataset_constructor.py -i resources/fisher_out/ -o " + working_set

#os.system(working_database_line)

