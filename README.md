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


## Running compact_class with Singularity

1. Make sure you have [Singularity installed](https://apptainer.org/docs/user/latest/quick_start.html) according to your operating system.

2. To retrieve the singularity build of compact_class, download the `.sif` image:

```bash
 singularity pull library://jduque2/mag_classification/compact_class:latest
```

3. To set up input/output directories you'll want to have two directories in your working directory one containing your input mags and the second an empty directory for logs and output files to be written to. In the rest of this guide the input directory is 'input' and the output directory is 'output'

4. You can run compact_class container interactively using this command.

Singularity Flags:
* `shell` starts and interactive shell inside the container
* `--bind` maps local directories into the container, changesmade inside the container will be reflected outside. Format [local/directory]:[container/directory]

Container Required Arguements:
* `cores` threads to use (default = 2)
* `input` path to input folder in container (defualt: `/input/`)
* `mode` whether to perform phylogenetic placement of hits on single gene trees (sgt) or on species tree (concat) (default: concat)
* `output` path to output folder in container (default: `/output/`)
* `database` which database to use between PhyloFisher (PF) and EukProt (EP) (default: PF)

```bash
singularity shell \
    --bind $(pwd)/input:/input \
    --bind $(pwd)/output:/output \
    compact_class_latest.sif [cores] [input] [mode] [output] [database]
```

Example:

```bash
singularity shell \
    --bind $(pwd)/input/:/input \
    --bind $(pwd)/output/:/output \
    compact_class_latest.sif 4 input/ concat output/ PF
```

Your command line prompt should change to indicate that you are now working inside your compact_class container.

You can exit the Singularity container by typig exit.

5. The container can also be run non-interactively:

```bash
singularity run \
    --bind $(pwd)/input:/input \
    --bind $(pwd)/input:/input \
    compact_class.sif [cores] [input] [mode] [output] [database]
```

## Using Snakemake Directly

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

All output is saved in your output directory with each individual MAG in your target direct generating an.

| Directory / File Path                            | Description                                                  |
|--------------------------------------------------|--------------------------------------------------------------|
| output/{mag}_input_metadata.tsv              | Metadata file containing input details for each MAG         |
| output/{mag}_q_frags/{gene}.fas              | Query fragment FASTA files for each single gene             |
| output/{mag}_mafft_out/{gene}.aln            | MAFFT alignment output for each gene                        |
| output/{mag}_mafft_out/{gene}.aln.partial.fas| Filtered/partial version of each gene alignment             |
| output/{mag}_mafft_out/{gene}.trimal         | Trimmed gene alignment produced using Trimal                |
| output/{mag}_epa_out/{mag}_epa_out.jplace    | Phylogenetic placement results in `.jplace` format          |
| output/{mag}_epa_out/profile.tsv             | Taxonomic profile summary derived from placement results    |
| output/busco_out/{mag}/summary.txt           | BUSCO summary report assessing MAG genome completeness      |


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
 - https://github.com/o-william-white/skim2mito/
 - https://github.com/McTavishLab/extensiphy
