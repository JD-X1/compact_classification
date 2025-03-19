# compact_class

compact_class is a pipeline designed to classify newly generated eukaryotic Metagenome Assembled Genomes (MAGs) using phylogenetic placements of highly conserved single genes using RAxML-EPA. It is currently designed to do this using the PhyloFisher [need a link] database. (WARNING: this pipeline has not be tested outside of HPCC environments, it is not reccomended to run this pipeline locally.)

## Installation and Dependencies
1. To use compact_class, first clone the repo into your working directory. 

```
git clone git@github.com:JD-X1/compact_classification.git
```

2. Then, ensure that you have a conda or mamba installation (for help installing mamba see [link here]). To run compact_class you'll need to do so in a conda/mamba environment with snakemake package installed. To create such an environment:

```
conda create -c conda-forge -c bioconda -n snakemake snakemake

```

3. Lastly, if you intend to use PhyloFisherDB as your database (only option currently available) be sure to place the directory in the resources directory. You can download the database into the correct location by running the following lines of code from inside the compact_classification directory.
```
cd resources
wget https://ndownloader.figshare.com/files/29093409
tar -xzvf 29093409
rm 29093409
cd ..
```

# How to Run on Your Data

Assuming you are using the PhyloFisher database, you should first place all your query MAGs into a single directory. You can then run the following:

```
sh comp_class.sh -m [directory containing MAGs] -d df -t [threads count]
```

### Required Flags

```
- (-m) directory containing metagenome assembled genomes in fasta format
- (-d) indicate the database type used (PF/EP), (DEFAULT: PF)
- (-t) threads (DEFAULT: 1)
```

### Output Files:

compact_class outputs a number of files per input MAG:

```
- [MAG_NAME]_summary.csv
```
CSV file summarizing the number of placements that support membership to a particular clade starting from Eukarya and moving downward into narrower and narrower clades.

In the resources directory you will find an number of directories containing intermediary files you may find useful:

```
- resources/[MAG_NAME]_epa_out/[GENE]/
```
Contains the output from RAxML EPA including the placements in 'jplace' format. As well as the output for gappa's 'assign' command which reports highest likelihood placement in the event that multiple placements were reported. 


```
- resources/[MAG_NAME]_mafft_out/[GENE].aln
```
Contains alignments of query homologs aligned to the corresponding single-gene reference alignments.


# Testing Your Installation

Be aware no clean up options have been implemented quite yet. As a result this will produce quite a few intermediary directories and steps in the resources directory.

```
conda activate snakemake

sh test_pf.sh
```

# Ensure MAGS end in extension '.fna'

Current behavior is implemented to allow me to control which files in a folder I want to be put through the pipeline without having to delete/move files. 

```
[MAGname].fna # different fasta file extensions will not be detected 
```

Also avoid the following symbols in the name as certain tools in this pipeline do not interact well with them.

```

'_'
'@'
'..'
'whitespace'
'*'

```
