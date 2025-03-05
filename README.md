# compact_class

## Instructions

compact_class is designed to classify newly generated Metagenome Assembled Genomes (MAGs) using multiple single gene phylogenetic placement. 

### Installation
```
git clone git@github.com:JD-X1/compact_classification.git
```




### Step 1: Download PhyloFisherDB

```
cd compact_classification/resources
wget https://ndownloader.figshare.com/files/29093409
tar -xzvf 29093409
rm 29093409
cd ..
```

### Dependencies

The pipleine is designed to manage the packages for the pipeline through conda or mamba. It has only two software dependencies conda and through conda an installation of snakemake.

Dependencies:
```
conda
    - snakemake
```

If you already have conda installed simply make a new environment w/ snakemake installed. To create a conda  environment w/ snakemake:
```
conda create -n snakemake snakemake
```

### Test Your Installation

Be aware no clean up options have been implemented quite yet. As a result this will produce quite a few intermediary directories and steps in the resources directory.
```
conda activate snakemake

sh test_pf.sh
```

### Ensure MAGS end in extension '.fna'

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