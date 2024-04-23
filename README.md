# compact_class

## Instuctions

compact_class is designed to classify newly generated MAGs useing multiple single gene phylogenetic placement. 

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

### Step 2: Run pipeline

```
sh comp_class.sh -m [directory containing MAGs]
```
