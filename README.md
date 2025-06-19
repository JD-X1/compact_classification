# compact_class

**compact_class** is a snakemake pipeline designed to classify recently assembled putatively eukaryotic metagenome assembled genomes (MAGs). The pipeline carries out a phylogenetic placement using a pre-generated database of 240 single gene alignments and phylogenies. The pipeline automates the identification of potential orthologs, alignment of query orthologs to reference alignments, and the phylogenetic placement itself. compact_class is intended to be used on linux HPCC environments.

## Contents
 - [Setup](#setup)
 - [Example data](#example-data)
 - [Input](#input)
 - [Output](#output)
 - [Running your own data](#running-your-own-data)
 - [Getting help](#getting-help)
 - [Citations](#citations)

## Setup

This workflow was implemented in snakemake and depends on conda to manage dependencies [dependencies] (#citations).

To install conda/mamba we reccomend following the official install guide. See details: https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html.

Once installed you will have to create an environment containing the latest version of snakemake:
```
conda create -c conda-forge -c bioconda -n snakemake snakemake compleasm
```

Clone compact_class into your working directory:
```
# get github repo
git clone git@github.com:JD-X1/compact_classification.git

# change dir
cd compact_class

# Remember to activate the conda environment before attempting to run compact_class
conda activate snakemake
```

You will also need to download both the eukaryota BUSCO database as well ass the PhyloFisher database:
```
# move into the resources directory
cd resources

# Download PhylofisherDB 
wget https://ndownloader.figshare.com/files/29093409

# Download the BUSCO database
python compleasm.py download lineages eukaryota

# Unpack Compressed Directory
tar -xzvf 29093409

# Delete Compressed version of Database
rm 29093409

# Return to bas compact_class directory
cd ..
```

## Building your own Docker Image (Under)

1. Make sure you have [Docker installed](https://www.docker.com/products/docker-desktop) according to your operating system.

2. To pull the Docker build of Extensiphy, run this command.

```bash
docker pull ????
```

3. We'll build your Extensiphy Docker container using this command.
* `-i` makes the container interactive.
* `-t` specifies the image to use as a template.
* `--name` specifies the container name.
```bash
docker run --name ep_container -i -t mctavishlab/extensiphy bash
```
Your command line prompt should change to indicate that you are now working
inside your Extensiphy container.

You can exit the docker container by typing `exit`.

To restart it and return to interactive analyses, run:

```bash
docker container restart ep_container
docker exec -it ep_container bash
```

## Example data

To test the pipeline, you can run the included example.

```
conda activate snakemake

snakemake -s rules/taxa_class_EP_concat.smk \
    --cores [available threads] --use-conda \
    --config  mag_dir=resources/test mode=EP\
    -p --keep-going --rerun-incomplete

```

## Input

compact_class is designed to directly take in MAGs in fasta format and will analyze all files in a directory containing the ".fna" extension.

Avoid the following symbols in fasta sequence/file names as certain tools in this pipeline do not interact well with them.

```
'_'
'@'
'..'
'whitespace'
'*'
```


## Output

All output is saved in the `resources/` directory with each individual MAG in your target direct generating an.

| Directory / File Path                            | Description                                                  |
|--------------------------------------------------|--------------------------------------------------------------|
| resources/{mag}_input_metadata.tsv              | Metadata file containing input details for each MAG         |
| resources/{mag}_q_frags/{gene}.fas              | Query fragment FASTA files for each single gene             |
| resources/{mag}_mafft_out/{gene}.aln            | MAFFT alignment output for each gene                        |
| resources/{mag}_mafft_out/{gene}.aln.partial.fas| Filtered/partial version of each gene alignment             |
| resources/{mag}_mafft_out/{gene}.trimal         | Trimmed gene alignment produced using Trimal                |
| resources/{mag}_epa_out/{mag}_epa_out.jplace    | Phylogenetic placement results in `.jplace` format          |
| resources/{mag}_epa_out/profile.tsv             | Taxonomic profile summary derived from placement results    |
| resources/busco_out/{mag}/summary.txt           | BUSCO summary report assessing MAG genome completeness      |


## Running your own data

To run your own data place all query mags into single directory and ensure each file uses the ".fna" extension. Assuming the necessary databases are available you can run the following from within the compact_classification directory:

```
conda activate snakemake

snakemake -s rules/taxa_class_EP_concat.smk \
    --cores [available threads] --use-conda \
    --config  mag_dir=[MY_DIRECTORY] mode=EP \
    -p --keep-going --rerun-incomplete
```

## Getting help

Open an issue or email [jduque2@ucmerced.edu].

## Citations

If using this pipeline, please cite:
 - compleasm: https://github.com/huangnengCSU/compleasm
 - minimprot: https://github.com/lh3/miniprot
 - HMMER: https://github.com/EddyRivasLab/hmmer
 - sepp: https://github.com/smirarab/sepp
 - PhyloFisher: https://github.com/TheBrownLab/PhyloFisher
 - blastp: http://blast.ncbi.nlm.nih.gov/Blast.cgi
 - cd-hit: https://github.com/weizhongli/cdhit
 - diamond: https://github.com/bbuchfink/diamond
 - MAFFT: https://github.com/GSLBiotech/mafft
 - divvier: https://github.com/simonwhelan/Divvier
 - trimal: https://github.com/inab/trimal
 - geneSticher.py: https://github.com/ballesterus/Utensils
 - EPA-NG: https://github.com/pierrebarbera/epa-ng
 - GAPPA: https://github.com/lczech/gappa

 Formatting for README.md largely taken from:
 https://github.com/o-william-white/skim2mito/
