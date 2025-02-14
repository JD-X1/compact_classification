#!usr/bin/bash -l

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

# get the file basename
file=$(basename ${Input} _input_metadata.tsv)


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
if [ ! -d resources/${file}_PhyloFishScratch ]
then
    echo "Creating the PhyloFishScratch database"
    cp -r resources/PhyloFisherDatabase_v1.0/database resources/${file}_PhyloFishScratch
    echo "PhyloFishScratch database created"
fi

echo "Casting Lines"
cd resources
# extract basename of input file
Input=$(basename $Input)


# get file basename

config.py -d ${file}_PhyloFishScratch -i $Input
echo "Configuration of PhyloFisher Modules Complete"
echo "Waiting for the Fish to Bite"
fisher.py --threads $TCores -o ${file}_fish_out 1> ../logs/FishingLogs/${file}_fisher.log 2> ../logs/FishingLogs/${file}_fisher.log
echo "Fish Caught"
informant.py -i ${file}_fish_out --orthologs_only 1> ../logs/FishingLogs/${file}_informant.log 2> ../logs/FishingLogs/${file}_informant.log
echo "Informant Complete"
echo "Choosing the best fish"
working_dataset_constructor.py -i ${file}_fish_out -o ${file}_working_dataset 1> ../logs/FishingLogs/${file}_working_dataset_constructor.log 2> ../logs/FishingLogs/${file}_working_dataset_constructor.log
echo "Fish on the grill"
cd ..
