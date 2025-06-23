#!usr/bin/bash --login

# only two arguments are required:
# -t <number of threads>
# -i <input metadata file>

# get the number of threads
while getopts t:i: option
do
    case "${option}" in
        t) TCores=${OPTARG};;
        i) Input=${OPTARG};;
    esac
done

conda activate fisher


# get the file basename
mag=$(basename ${Input} _input_metadata.tsv)


# check if log directory exists
# if not, create it
echo "Gathering Bait"
if [ ! -d logs ]
then
    mkdir logs
fi

if [ ! -d logs/FishingLogs ]
then
    mkdir logs/FishingLogs
fi

# check if the PhyloFishScratch database exists
# if not, create it
if [ ! -d resources/${mag}_PhyloFishScratch ]
then
    echo "Creating the PhyloFishScratch database"
    cp -r resources/PhyloFisherDatabase_v1.0/database resources/${mag}_PhyloFishScratch
    echo "PhyloFishScratch database created"
fi

echo "Casting Lines"
cd resources
# extract basename of input file
Input=$(basename $Input)


# get file basename

config.py -d ${mag}_PhyloFishScratch -i $Input
echo "Configuration of PhyloFisher Modules Complete"
echo "Waiting for the Fish to Bite"
fisher.py --threads $TCores -o ${mag}_fish_out --keep_tmp 1> ../logs/FishingLogs/${mag}_fisher.log 2> ../logs/FishingLogs/${mag}_fisher.log
echo "Fish Caught"
informant.py -i ${mag}_fish_out --orthologs_only 1> ../logs/FishingLogs/${mag}_informant.log 2> ../logs/FishingLogs/${mag}_informant.log
echo "Informant Complete"
echo "Choosing the best fish"
working_dataset_constructor.py -i ${mag}_fish_out -o ${mag}_working_dataset 1> ../logs/FishingLogs/${mag}_working_dataset_constructor.log 2> ../logs/FishingLogs/${mag}_working_dataset_constructor.log
echo "Fish on the grill"
cd ..