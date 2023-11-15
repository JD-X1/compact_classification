#!usr/bin/bash -l

# only one argument is required:
# -t <number of threads>


# get the number of threads
while getopts t: option
do
    case "${option}" in
        t) TCores=${OPTARG};;
    esac
done



cd resources/
config.py -d PhyloFishScratch -i input_metadata.tsv
fisher.py --threads $TCores -o fish_out
informant.py -i fish_out --orthologs-only
working_dataset_constructor.py -i fish_out -o working_dataset
cd ..
