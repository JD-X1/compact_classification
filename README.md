# compact_class

## Instructions

compact_class is designed to classify newly generated Metagenome Assembled Genomes (MAGs) using multiple single gene phylogenetic placement. 

**Be Aware of MAG naming conventions** 

### Step 0: Download PhyloFisherDB

```
cd resources
wget https://ndownloader.figshare.com/files/29093409
tar -xzvf 29093409
rm 29093409
cd ..
```

### Step 1: Ensure MAGS end in extension '.fna'

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

### Step 2: Install Mamba + Snakemake

The pipleine is designed to manage the packages for the pipeline through Mamba. I recommend the miniforge installation. 

For unix-like platforms:
```
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
conda init
mamba create -c conda-forge -c bioconda -n snakemake snakemake

```

### Step 3: Run pipeline Test

Be aware no clean up options have been implemented quite yet. As a result this will produce quite a few intermediary directories and steps in the resources directory.
```
mamba activate raxml-ng

sh test_pf.sh
```
