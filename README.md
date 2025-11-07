# Compact Class

**compact_class** is a snakemake pipeline designed to classify recently assembled putatively eukaryotic metagenome assembled genomes (MAGs). The pipeline carries out a phylogenetic placement using a pre-generated database of 240 single gene alignments and phylogenies. The pipeline automates the identification of potential orthologs, alignment of query orthologs to reference alignments, and the phylogenetic placement itself. compact_class is intended to be used on linux HPCC environments.

## Contents
 - [Running compact_class with Singularity](#running-compact_class-with-singularity)
 - [Manual Setup](#setup-for-manual-execution)
 - [Using Snakemake Directly](#using-snakemake-directly)
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

3. Prepare local directories for your MAGs and the workflow outputs (for
   example, `input/` and `output/`). The output directory is used as the
   container home directory so that cached downloads (BUSCO, Compleasm,
   reference scratch space) persist between runs.

4. Launch an interactive shell if you want to inspect the container
   environment before executing the workflow:

   ```bash
   singularity shell --cleanenv \
       --bind $(pwd)/input:/input \
       --bind $(pwd)/output:/output \
       compact_class_latest.sif
   ```

   Useful flags:
   * `--bind` mounts host directories into the container. Any changes
     written inside `/input` or `/output` appear in the host directories.
   * `--cleanenv` is a required arguement, ensures host environmental variable
      don't override environmental variables set within the container

   Exit the interactive shell with `exit` when you are ready to leave the
   container.

5. Run the workflow non-interactively using the container entrypoint. The
   entrypoint accepts both short and long options:

   ```bash
   singularity run \
       --bind $(pwd)/input:/input \
       --bind $(pwd)/output:/output \
       compact_class_latest.sif \
       --cores 8 \
       --mag-dir /input \
       --mode concat \
       --outdir /output \
       --db PF
   ```

   Additional flags mirror the native Snakemake configuration:
   * `--mode` selects between concatenated placements (`concat`) and
     single-gene tree placements (`sgt`).
   * `--db` chooses the reference database (`PF` for PhyloFisher or `EP`
     for the extended EukProt database).
   * `--augustus`, `--trim`, and `--purge <UID or long name>` toggle the
     optional workflow branches. See [Running your own data](#running-your-own-data)
     for details on these switches.

## Manual Setup

1. Install a recent version of [conda or mamba](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).
   Using mamba is recommended because the workflow relies on multiple
   environments and larger bioinformatics packages.

2. Clone this repository and create the Snakemake control environment:

   ```bash
   git clone https://github.com/JD-X1/compact_classification.git
   cd compact_classification
   mamba env create -f envs/snakemake.yaml
   conda activate snakemake
   ```

   The remaining environments required by each rule are created
   automatically by Snakemake from the YAML files in `envs/`.

3. Download the reference databases into `resources/`:

   ```bash
   cd resources
   # PhyloFisher reference data (required for PF runs)
   wget https://ndownloader.figshare.com/files/29093409 -O PhyloFisherDatabase_v1.0.tgz
   tar -xzf PhyloFisherDatabase_v1.0.tgz

   # Optional: extended EukProt reference set used by the EP modes
   # (PF_ExtendedEukProtDB_v0.1 should unpack into resources/)
   ```

   The pipeline also expects the Compleasm/BUSCO downloads
   (`resources/mb_downloads` and `resources/busco_downloads`). These
   directories are created automatically the first time the workflow runs
   the `compleasm` and BUSCO download steps, but they can be staged in
   advance to avoid redundant downloads on shared systems.

   > **Note:** The extended EukProt database is not bundled with this
   > repository. Please refer to the project documentation or contact the
   > maintainers for access instructions if you plan to run the EP modes.

   Return to the repository root once the downloads are complete:

   ```bash
   cd ..
   ```

## Using Snakemake Directly

This workflow was implemented in Snakemake and depends on conda/mamba to
provide rule-specific environments. After completing the [Setup](#setup)
steps, you can run the pipeline directly or use the provided helper
script.

To invoke Snakemake directly:

```
snakemake -s rules/taxa_class_PF_concat.smk \
    --config mag_dir=resources/test outdir=resources/test_output database=PF \
    --cores 8 --use-conda --conda-frontend mamba \
    -p --keep-going --rerun-incomplete
```

Key configuration parameters:

| Config key   | Description |
|--------------|-------------|
| `mag_dir`    | Directory containing the MAG fasta files. |
| `outdir`     | Output directory. All intermediate and final files are written here. |
| `mode`       | `concat` (default) for concatenated supermatrix placement or `sgt` for single-gene placements. |
| `database`   | `PF` (PhyloFisher) or `EP` (extended EukProt). |
| `augustus`   | Set to `True` to run BUSCO with Augustus gene prediction instead of Compleasm-only mode. |
| `trim`       | Set to `True` to enable additional alignment trimming. |
| `purge`      | Optional unique ID or long taxon name to remove from the reference database prior to placement. |


## Input


compact_class accepts MAG assemblies in FASTA format. Files ending in
`.fna`, `.fa`, `.fasta` (and their gzipped equivalents) are detected
automatically. All FASTA sequences should be nucleotide assemblies; if a
proteome is provided accidentally the workflow will fail during the BUSCO
stage.


Avoid the following symbols in fasta sequence/file names as certain tools in this pipeline do not interact well with them.

```
'_'
'@'
'..'
'whitespace'
'*'
```


## Output

All output is written beneath the configured `outdir` (or `/output` when
running inside the container). For each MAG (`{mag}`) the workflow
generates the following key files and directories:

| Directory / File Path                                      | Description |
|------------------------------------------------------------|-------------|
| `{outdir}/busco_out/{mag}/summary.txt`                     | Compleasm/BUSCO summary for the assembly. |
| `{outdir}/{mag}_input_metadata.tsv`                        | Metadata describing the query assembly and matched reference taxa. |
| `{outdir}/{mag}_PhyloFishScratch/`                         | Copy of the selected reference database subset used during placement. |
| `{outdir}/{mag}_working_dataset/`                          | Gene families identified by `fishing.sh` for downstream processing. |
| `{outdir}/{mag}_q_frags/{gene}.fas`                        | Query sequences for each single-copy marker gene. |
| `{outdir}/{mag}_mafft_out/{gene}.aln(.partial.fas/.trimal)`| Raw, filtered, and trimmed alignments for each marker gene. |
| `{outdir}/{mag}_SuperMatrix.fas`                           | Concatenated alignment used for species-tree placement (concat mode). |
| `{outdir}/{mag}_q.aln` / `{outdir}/{mag}_ref.aln`          | Split query/reference alignments used by EPA-NG. |
| `{outdir}/{mag}_ref.tre`                                   | Pruned reference tree tailored to the current MAG. |
| `{outdir}/{mag}_epa_out/{mag}_epa_out.jplace`              | EPA-NG placement results in JPLACE format. |
| `{outdir}/{mag}_epa_out/profile.tsv`                       | Taxonomic assignments derived from GAPPA. |
| `{outdir}/{mag}_summary.csv` (sgt modes)                   | Aggregated per-gene placement summary (generated by SGT workflows). |
| `{outdir}/{mag}_epa_out/profile_summary.tsv` (sgt modes)   | Consolidated summary of GAPPA per-gene profiles. |

Intermediate directories can be safely removed after extracting the
final reports if storage is a concern.


## Running your own data

Place all MAG assemblies in a dedicated directory and ensure the files
have consistent extensions (see [Input](#input)). Example Snakemake
invocation that runs the PhyloFisher concatenated workflow with
additional trimming:

```
snakemake -s rules/taxa_class_PF_concat.smk \
    --cores 16 --use-conda --conda-frontend mamba \
    --config mag_dir=/path/to/mags outdir=/path/to/output \
             database=PF mode=concat trim=True \
    -p --keep-going --rerun-incomplete
```

## Getting help

Open an issue or email [jduque2@ucmerced.edu].

# Citations

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
