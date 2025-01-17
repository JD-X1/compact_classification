#!/bin/bash

# arguement -m is location of target mag files
# arguement -d is ref_alignment filepath PF vs EP

while getopts m: flag
do
    case "${flag}" in
        m) mag_dir=${OPTARG};;
        d) db_type=${OPTARG};;
    esac
done

echo "Checking for PhyloFisherDatabase_v1.0"
if [ ! -d resources/PhyloFisherDatabase_v1.0 ]; then
    echo "PhyloFisherDatabase_v1.0 not found. Please download and extract the database to resources/PhyloFisherDatabase_v1.0"
    echo "Download instructions in the resources directory run the following commands:\n\n"
    echo "wget https://ndownloader.figshare.com/files/29093409"
    echo "tar -xzvf 29093409"
    exit 1
fi

echo "DB Found"
echo "Generating working directory"

if [ ! -d resources/PhyloFishScratch ]; then
    cp -r resources/PhyloFisherDatabase_v1.0/database resources/PhyloFishScratch
fi


### Check for Snakemake installation

###

snakemake -s rules/taxa_class.smk \
    --cores $SLURM_NTASKS_PER_NODE \
    --config mag_dir=$mag_dir mode=$db_type\
    --use-conda -p --keep-going \
    --rerun-incomplete \
    --conda-frontend mamba 


snakemake -s rules/taxa_class.smk \
    --cores $SLURM_NTASKS_PER_NODE \
    --config mag_dir=$mag_dir \
    --use-conda -p --keep-going \
    --rerun-incomplete \
    --conda-frontend mamba 
