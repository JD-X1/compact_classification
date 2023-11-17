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
if [ ! -d resources/PhyloFishScratch ]
then
    echo "Creating the PhyloFishScratch database"
    cp -r resources/PhyloFisherDatabase_v1.0/database resources/PhyloFishScratch
    echo "PhyloFishScratch database created"
fi

echo "Casting Lines"
cd resources
# extract basename of input file
Input=$(basename $Input)

config.py -d PhyloFishScratch -i $Input
echo "Configuration of PhyloFisher Modules Complete"
echo "Waiting for the Fish to Bite"
fisher.py --threads $TCores -o fish_out 1> ../logs/FishingLogs/fisher.log 2> ../logs/FishingLogs/fisher.log
echo "Fish Caught"
informant.py -i fish_out --orthologs_only 1> ../logs/FishingLogs/informant.log 2> ../logs/FishingLogs/informant.logs
echo "Informant Complete"
echo "Choosing the best fish"
working_dataset_constructor.py -i fish_out -o working_dataset 1> ../logs/FishingLogs/working_dataset_constructor.log 2> ../logs/FishingLogs/working_dataset_constructor.log
echo "Fish on the grill"
cd ..
