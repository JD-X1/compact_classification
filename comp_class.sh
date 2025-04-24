#!/bin/bash

# arguement -m is location of target mag files
# arguement -d is ref_alignment filepath PF vs EP
# arguement -t is number of threads that defaults to 1
threads=1
db_type=PF
mag_dir=resources/test/
placement_method=SGT
while getopts m:d:t:p: flag
do
    case "${flag}" in
        m) mag_dir=${OPTARG};;
        d) db_type=${OPTARG};;
        t) threads=${OPTARG};;
        p) placement_method=${OPTARG};;
    esac
done

if ${placement_method} == "SGT"; then
    $placement_method=""
elif ${placement_method} == "CONCAT"; then
    $placement_method="_concat"
else
    echo "Invalid placement method. Please use SGT or CONCAT"
    exit 1
fi


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

if [ -z "$(ls $mag_dir/*.fna 2>/dev/null)" ]; then
    echo "No .fna files found in $mag_dir"
    exit 1
fi


snakemake -s rules/taxa_class_${db_type}${placement_method}.smk \
    --cores ${threads} \
    --config mag_dir=$mag_dir mode=$db_type\
    --use-conda -p --keep-going \
    --rerun-incomplete
