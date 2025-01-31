
###################################################
# Uncomment the following if you need to install  #
# mamba and snakemake                             #
###################################################

# curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
# bash Miniforge3-$(uname)-$(uname -m).sh
# conda init
# mamba create -c conda-forge -c bioconda -n snakemake snakemake


###################################################
#                   Test Script                   #
###################################################

# ADD core management flag
sh comp_class.sh -m resources/test/ -d PF -t 22

# check for any ouput
# throw loud error if no
# output is found
if [ ! -d ./Lenisialimosa_summary.csv]; then
    echo "Pipeline Failed Check Logs For Errors"
    exit 1
fi
