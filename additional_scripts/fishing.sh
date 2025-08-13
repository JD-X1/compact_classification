#!usr/bin/env bash

# only two arguments are required:
# -t <number of threads>
# -i <input metadata file>

# get the number of threads
while getopts t:i:o: option
do
    case "${option}" in
        t) TCores=${OPTARG};;
        i) Input=${OPTARG};;
        o) Outdir=${OPTARG};;
    esac
done

# get the file basename
mag=$(basename ${Input} _input_metadata.tsv)
working_dataset=${mag}_working_dataset
log_dir=${Outdir}logs/FishingLogs
phyloscratch_dir=${Outdir}${mag}_PhyloFishScratch
echo ${phyloscratch_dir}
echo 'Fishing for ${mag}'
echo 'log outdir: ${log_dir}'
echo 'phyloscratch dir: ${phyloscratch_dir}'

# check if log directory exists
# if not, create it
echo "Gathering Bait"
if [ ! -d ${log_dir} ]
then
    mkdir -p ${log_dir}
fi

if [ ! -d ${log_dir} ]
then
    mkdir -p ${log_dir}
fi


if [ ! -d ${phyloscratch_dir} ]
then
    echo "Creating the PhyloFishScratch database"
    cp -r /compact_classification/resources/PhyloFisherDatabase_v1.0/database ${phyloscratch_dir}
    echo "PhyloFishScratch database created"
fi

echo "Casting Lines"
cd ${Outdir}
# extract basename of input file
Input=$(basename $Input)
echo $Input

# get file basename

config.py -d ${phyloscratch_dir} -i $Input
echo "Configuration of PhyloFisher Modules Complete"
echo "Waiting for the Fish to Bite"
fisher.py --threads $TCores -o ${mag}_fish_out --keep_tmp 1> logs/FishingLogs/${mag}_fisher.log 2> logs/FishingLogs/${mag}_fisher.log
echo "Fish Caught"
informant.py -i ${mag}_fish_out --orthologs_only 1> logs/FishingLogs/${mag}_informant.log 2> logs/FishingLogs/${mag}_informant.log
echo "Informant Complete"
echo "Choosing the best fish"
working_dataset_constructor.py -i ${mag}_fish_out -o ${mag}_working_dataset 1> logs/FishingLogs/${mag}_working_dataset_constructor.log 2> logs/FishingLogs/${mag}_working_dataset_constructor.log
echo "Fish on the grill"
cd ..