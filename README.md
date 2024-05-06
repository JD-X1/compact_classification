# compact_class

## Instructions

compact_class is designed to classify newly generated MAGs using multiple single gene phylogenetic placement. 

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
[mag_name].fna # different fasta file extensions will not be detected 
```

Also avoid the following symbols in the name as certain tools in this pipeline do not interact well with them.

```

'_'
'@'
'..'
'whitespace'
'*'

```

### Step 2: Run pipeline

Be aware no clean up options have been implemented quite yet. As a result this will produce quite a few intermediary directories and steps in the resources directory.

```
sh comp_class.sh -m [directory containing MAGs]
```
