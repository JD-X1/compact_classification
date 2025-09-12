#!usr/bin/env python

import os
import datetime

GENOME_EXTS   = [".fna", ".fa", ".fasta", ".fna.gz", ".fa.gz", ".fasta.gz"]
PROTEOME_EXTS = [".faa", ".faa.gz", ".aa.fa", ".aa.fasta"]

print("[{:%Y-%m-%d %H:%M:%S}]: Checking for resources directory...".format(datetime.datetime.now()))
if os.path.exists("/compact_classification/resources/"):
    RESOURCES_DIR = "/compact_classification/resources/"
    print("[{:%Y-%m-%d %H:%M:%S}]: Found resources directory at /compact_classification/resources/".format(datetime.datetime.now()))
elif os.path.exists("resources/"):
    RESOURCES_DIR = "resources/"
    print("[{:%Y-%m-%d %H:%M:%S}]: Found resources directory at resources/".format(datetime.datetime.now()))
else:
    raise ValueError("[{:%Y-%m-%d %H:%M:%S}]: If running from source code and not singularity, please ensure that the 'resources' directory is present in the current working directory and be sure it contains the necessary databases.".format(datetime.datetime.now()))

print("[{:%Y-%m-%d %H:%M:%S}]: Checking for additional scripts directory...".format(datetime.datetime.now()))
if os.path.exists("/compact_classification/additional_scripts/"):
    ADDITIONAL_SCRIPTS_DIR = "/compact_classification/additional_scripts/"
    print("[{:%Y-%m-%d %H:%M:%S}]: Found additional scripts directory at /compact_classification/additional_scripts/".format(datetime.datetime.now()))
elif os.path.exists("additional_scripts/"):
    ADDITIONAL_SCRIPTS_DIR = "additional_scripts/"
    print("[{:%Y-%m-%d %H:%M:%S}]: Found additional scripts directory at additional_scripts/".format(datetime.datetime.now()))
else:
    raise ValueError("[{:%Y-%m-%d %H:%M:%S}]: If running from source code and not singularity, please ensure that the 'additional_scripts' directory is present in the current working directory and be sure it contains the necessary scripts.".format(datetime.datetime.now()))


def get_genes_from_goneFishing(mag):
    checkpoint_output = config["outdir"] + "{mag}_working_dataset"

    # Init an empty list to store genes with corresponding trees
    valid_genes = []

    # Check if the checkpoint output directory exists to handle cases where it might not
    if os.path.exists(checkpoint_output):
        # List files in the output directory
        files = os.listdir(checkpoint_output)

        # Check each file to see if a corresponding tree file exists
        for f in files:
            if f.endswith('.fas'):
                gene_name = os.path.splitext(f)[0]
                tree_file = os.path.join("/compact_classification/resources/ref_trees", gene_name, gene_name + ".raxml.support")
                if os.path.exists(tree_file):
                    valid_genes.append(gene_name)

    return valid_genes

def sanitize_gene_name(name):
    return name.replace("/", "_").replace(" ", "_").replace("\\", "_")

def get_superMatrix_targets_for_mag(mag):
    ckpt_out = checkpoints.goneFishing.get(mag=mag).output[0]
    gene_files = glob_wildcards(os.path.join(ckpt_out, "{gene}.fas")).gene
    print([sanitize_gene_name(gene) for gene in gene_files])
    return [sanitize_gene_name(gene) for gene in gene_files]

def _mag_path_candidates(base):
    return [base + ext for ext in GENOME_EXTS + PROTEOME_EXTS]

def find_mag_file(wildcards):
    base = os.path.join(config["mag_dir"], wildcards.mag)
    for cand in _mag_path_candidates(base):
        if os.path.exists(cand):
            return cand
    raise ValueError(
    f"[{datetime.datetime.now():%Y-%m-%d %H:%M:%S}]: "
    f"No input file found for {wildcards.mag} in {config['mag_dir']} "
    f"with any of the expected extensions for genomic sequence: {GENOME_EXTS} "
    f"or proteomic sequence: {PROTEOME_EXTS}"
    )

def is_proteome(path):
    return any(path.endswith(ext) for ext in PROTEOME_EXTS)

# checking if config["outdir"] ends with a slash if it doesn't add one
if not config["outdir"].endswith("/"):
    config["outdir"] += "/"
if not config["mag_dir"].endswith("/"):
    config["mag_dir"] += "/"
outdir = config["outdir"]
mag_f = os.listdir(config["mag_dir"])

# Check if config["augustus"] is set, if not set it to default
print("[{:%Y-%m-%d %H:%M:%S}]: Command invoked with the following options:".format(datetime.datetime.now()))
print("[{:%Y-%m-%d %H:%M:%S}]: Output directory: {}".format(datetime.datetime.now(), config["outdir"]))
print("[{:%Y-%m-%d %H:%M:%S}]: MAG directory: {}".format(datetime.datetime.now(), config["mag_dir"]))

augustus = False
if "augustus" in config:
    augustus = True
    print("[{:%Y-%m-%d %H:%M:%S}]: Will use Augustus for BUSCO runs.".format(datetime.datetime.now()))
else:
    print("[{:%Y-%m-%d %H:%M:%S}]: Will use Compleasm for BUSCO runs.".format(datetime.datetime.now()))
trim_alignments = False
if "trim" in config:
    trim_alignments = True
    print("[{:%Y-%m-%d %H:%M:%S}]: Trimming alignments with trimAl + divvier.".format(datetime.datetime.now()))

proteome_input = False
if "--proteome" in config:
    proteome_input = True
    print("[{:%Y-%m-%d %H:%M:%S}]: Using proteome input instead of BUSCO Output.".format(datetime.datetime.now()))

#### CHANGE THE FILE EXTENSION
#### Swap out this check for file extension
#### with a test to check that the MAG files
#### are actually nucleic acid fastas.


mag_f = [f for f in mag_f if any(f.endswith(ext) for ext in (GENOME_EXTS + PROTEOME_EXTS))]

if mag_f == []:
    raise ValueError("#################\nNo MAGs found in the specified directory.\n#################\n")

# get mag names
genes = []
mags = []
for f in mag_f:
    mag = ""
    if len(f.split(".")) >= 2:
        # get mag name
        mag = '.'.join(f.split(".")[:-1])
    # append to list
    else:
        raise ValueError(
            "#################\nMake sure MAG file names follow the folowwing format:\n\t[unique id].[file extension] \n#################\n"
        )
    if "_" in mag or " " in mag:
        mag = mag.replace("_", "").replace(" ", "")
        os.rename(
            os.path.join(config["mag_dir"], f),
            os.path.join(config["mag_dir"], mag + ".fna")
        )
        print(f"Renaming MAG file to: {mag}.fna")
    print(f"Processing MAG: " + mag)
    mags.append(mag)

rule all:
    input:
        expand(config["outdir"] + "{mag}_working_dataset", mag=mags),
        expand(
            config["outdir"] + "{mag}_q_frags/{gene}.fas",
                mag=mags,
                gene=lambda wildcards: [
                    gene
                    for mag in mags
                    for gene in get_superMatrix_targets_for_mag(mag)
                ]
            ),
        expand(
            config["outdir"] + "{mag}_mafft_out/{gene}.aln",
                mag=mags,
                gene=lambda wildcards: [
                    gene
                    for mag in mags
                    for gene in get_superMatrix_targets_for_mag(mag)
                ]
            ),
        expand(config["outdir"] + "{mag}_q.aln", mag=mags),
        expand(config["outdir"] + "{mag}_ref.aln", mag=mags),
        expand(config["outdir"] + "{mag}_SuperMatrix.fas", mag=mags),
        expand(config["outdir"] + "{mag}_epa_out/{mag}_epa_out.jplace", mag=mags),
        expand(config["outdir"] + "{mag}_epa_out/profile.tsv", mag=mags)
        
    
rule run_busco:
    input:
        find_mag_file
    output:
        branch(
            augustus,
            [config["outdir"] + "busco_out/{mag}/short_summary.specific.eukaryota_odb12.{mag}.txt", config["outdir"] + "busco_out/{mag}/eukaryota_odb12/translated_protein.fasta"],
            [config["outdir"] + "busco_out/{mag}/summary.txt", config["outdir"] + "busco_out/{mag}/eukaryota_odb12/translated_protein.fasta"]
       )
    conda:
        branch(augustus,
        "busco", # executes rule w/ buscos w/ augustus
        "compleasm" # executes rule w/ compleasm's miniprot
        )
    threads: workflow.cores
    params:
        busco_mode = lambda wildcards, input: "proteins" if is_proteome(input[0]) else "genome",
        resources_dir = RESOURCES_DIR
    log:
        config["outdir"] + "logs/busco/{mag}.log"
    shell:
        branch(augustus,
        """
        echo "Running BUSCO for {wildcards.mag} (mode: {params.busco_mode})"
        busco -i {input} -m {params.busco_mode} \
            -l {params.resources_dir}/busco_downloads/lineages/eukaryota_odb12 \
            -c {threads} \
            -f --augustus \
            -o {config[outdir]}busco_out/{wildcards.mag} 1> {log} 2> {log}
        mkdir -p {config[outdir]}busco_out/{wildcards.mag}/eukaryota_odb12/
        if [[ "{params.busco_mode}" == "proteins" ]]; then
            cat {config[outdir]}busco_out/{wildcards.mag}/run_eukaryota_odb12/augustus_output/*faa* >> {config[outdir]}busco_out/{wildcards.mag}/eukaryota_odb12/translated_protein.fasta
        elif [[ "{params.busco_mode}" == "genome" ]]; then
            cat {config[outdir]}busco_out/{wildcards.mag}/run_eukaryota_odb12/augustus_output/*faa* >> {config[outdir]}busco_out/{wildcards.mag}/eukaryota_odb12/translated_protein.fasta
        fi
        """,
        
        """
        ### compleasm route requires genomic input
        if [[ "{params.busco_mode}" = "proteins" ]]; then
            echo "Error: Proteome input provided but compleasm mode selected. Please provide genomic input for compleasm."
            exit 2
        fi
        echo "Running compleasm for {wildcards.mag}"
        compleasm run -a {input} -t {threads} \
            -l eukaryota \
            -L {params.resources_dir}/mb_downloads/ \
            -o {config[outdir]}busco_out/{wildcards.mag} \
            1> {log} 2> {log}
        """
        )

rule fishing_meta:
    input:
        branch(
            proteome_input,
            config["mag_dir"],
            config["outdir"] + "busco_out/{mag}/eukaryota_odb12/translated_protein.fasta"
            )
        
    output:
        config["outdir"] + "{mag}_input_metadata.tsv"
    conda:
        "pline_max"
    params:
        ADD_SCRIPTS=ADDITIONAL_SCRIPTS_DIR
    threads: 1
    priority: 0
    shell:
        "python {params.ADD_SCRIPTS}/fishing_meta.py -p {input} -o {output}"


checkpoint goneFishing:
    input:
        config["outdir"] + "{mag}_input_metadata.tsv"
    output:
        directory(config["outdir"] + "{mag}_working_dataset")
    conda:
        "fisher"
    params:
        ADD_SCRIPTS=ADDITIONAL_SCRIPTS_DIR,
        resources_dir=RESOURCES_DIR
    threads: workflow.cores
    priority: 0
    shell:
        "bash {params.ADD_SCRIPTS}/fishing.sh -t {threads} -i {input} -r {params.resources_dir} -o {config[outdir]}"


rule splitter:
    input:
        dir=config["outdir"] + "{mag}_working_dataset/",
        tar=config["outdir"] + "{mag}_working_dataset/{gene}.fas",
        mag_dir=config["mag_dir"]
    output:
       config["outdir"] + "{mag}_q_frags/{gene}.fas"
    conda:
        "pline_max"
    params:
        ADD_SCRIPTS=ADDITIONAL_SCRIPTS_DIR
    threads: 1
    priority: 0
    log:
        config["outdir"] + "logs/splitter/{mag}_{gene}.log"
    shell:
        "python {params.ADD_SCRIPTS}splitter.py -i {input.tar} -d {input.mag_dir} -o {output} 1> {log} 2> {log}"


rule mafft:
    input:
        query=config["outdir"] + "{mag}_q_frags/{gene}.fas",
        reference=config["outdir"] + "{mag}_PhyloFishScratch/database/alignments/{gene}.fas.aln"
    output:
        config["outdir"] + "{mag}_mafft_out/{gene}.aln"
    conda:
        "pline_max"
    threads: 22
    priority: 0
    log:
        config["outdir"] + "logs/mafft/{mag}/{mag}_{gene}_mafft.log"
    shell:
        "mafft --auto --addfragments {input.query} --keeplength --thread {threads} {input.reference} > {output} 2> {log}"

#superMatrix_targets = get_superMatrix_targets(mag)
#print(len(superMatrix_targets))

rule divvier:
    input:
        config["outdir"] + "{mag}_mafft_out/{gene}.aln"
    output:
        config["outdir"] + "{mag}_mafft_out/{gene}.aln.partial.fas",
        config["outdir"] + "{mag}_mafft_out/{gene}.aln.PP"
    log:
        config["outdir"] + "logs/divvier/{mag}/{mag}_{gene}_divvier.log"
    conda:
        "div"
    threads: 1
    shell:
        """
        divvier -mincol 4 -partial -divvygap {input} > {log} 2> {log}
        """

rule trimal:
    input:
        config["outdir"] + "{mag}_mafft_out/{gene}.aln.partial.fas"
    output:
        config["outdir"] + "{mag}_mafft_out/{gene}.trimal"
    log:
        config["outdir"] + "logs/trimal/{mag}/{mag}_{gene}_trimal.log"
    conda:
        "trimal"
    threads: 1
    shell:
        """
        trimal -in {input} -gt 0.8 -out {output} > {log} 2> {log}
        """

rule concat:
    input:
        branch(
            trim_alignments,
            lambda wildcards: [
                f"{config['outdir']}{wildcards.mag}_mafft_out/{gene}.trimal"
                for gene in get_superMatrix_targets_for_mag(wildcards.mag)
            ],
            lambda wildcards: [
                f"{config['outdir']}{wildcards.mag}_mafft_out/{gene}.aln"
                for gene in get_superMatrix_targets_for_mag(wildcards.mag)
            ]
        )
    output:
        config["outdir"] + "{mag}_SuperMatrix.fas"
    conda:
        "pythas_two"
    params:
        ADD_SCRIPTS=ADDITIONAL_SCRIPTS_DIR
    threads: 1
    priority: 0
    params:
        out_dir=config["outdir"]
    log:
        config["outdir"] + "logs/concat/{mag}.log"
    shell:
        """
        FIXED_ALNS=()
        mkdir -p {params.out_dir}{wildcards.mag}_relabeled
        for i in $(realpath {params.out_dir}{wildcards.mag}_mafft_out/*trimal)
        do
        prot=$(basename ${{i}} .trimal)
        python {params.ADD_SCRIPTS}add_gene_name.py -a ${{i}} -g ${{prot}} -t {wildcards.mag} -o {params.out_dir}{wildcards.mag}_relabeled/${{prot}}.fas
        FIXED_ALNS+=("{params.out_dir}{wildcards.mag}_relabeled/${{prot}}.fas")
        done
        python2 {params.ADD_SCRIPTS}geneStitcher.py -in ${{FIXED_ALNS[@]}} 
        mv SuperMatrix.fas {output}
        """

rule alignment_splitter:
    input:
        config["outdir"] + "{mag}_SuperMatrix.fas"
    output:
        query=config["outdir"] + "{mag}_q.aln",
        ref=config["outdir"] + "{mag}_ref.aln"
    conda:
        "pline_max"
    threads: 1
    priority: 0
    params:
        out_dir=config["outdir"],
        ADD_SCRIPTS=ADDITIONAL_SCRIPTS_DIR
    log:
        config["outdir"] + "logs/alignment_splitter/{mag}.log"
    shell:
        "python {params.ADD_SCRIPTS}alignment_splitter.py -a {input} -t {wildcards.mag} -o {params.out_dir} > {log} 2> {log}"

rule sub_tree:
    input:
        aln=config["outdir"] + "{mag}_ref.aln",
        #tree="/compact_classification/resources/ref_concat_PF_alt3.tre" # uncomment for singularity
        tree="resources/ref_concat_PF_alt3.tre"  # uncomment for direct snakemake execution
    output:
        config["outdir"] + "{mag}_ref.tre"
    conda:
        "dendropy"
    params:
        ADD_SCRIPTS=ADDITIONAL_SCRIPTS_DIR
    threads: 1
    priority: 0
    log:
        config["outdir"] + "logs/sub_tree/{mag}.log"
    shell:
        """
        python {params.ADD_SCRIPTS}sub_tree.py -a {input.aln} -t {input.tree} -o {output} > {log} 2> {log}
        """

rule raxml_epa:
    input:
        q_aln= config["outdir"] + "{mag}_q.aln",
        ref_aln= config["outdir"] + "{mag}_ref.aln",
        ref_tree= config["outdir"] + "{mag}_ref.tre"
    output:
        config["outdir"] + "{mag}_epa_out/{mag}_epa_out.jplace"
    conda:
        "pline_max"
    threads: 22
    priority: 0
    params:
        out_dir=config["outdir"]
    log:
        config["outdir"] + "logs/raxml_epa/{mag}_epa.log"
    shell:
        """
        # Ensure the output directory exists
        mkdir -p {params.out_dir}{wildcards.mag}_epa_out/
        mkdir -p {params.out_dir}logs/raxml_epa/

        # Execute RAxML with the correct path to the alignment and tree files,
        # specifying only a run name for the output (not a path).
        ulimit -n 65536
        ulimit -s unlimited
        epa-ng --ref-msa {input.ref_aln} \
         --tree {input.ref_tree} \
         --query {input.q_aln} \
         --model LG -T {threads} >{log}
        # The actual output file is named according to the RAxML naming convention,
        # incorporating the run name. Ensure this matches your output specification.
        mv epa_result.jplace {output}
        mv epa_info.log {log}
        """

rule gappa:
    input:
        config["outdir"] + "{mag}_epa_out/{mag}_epa_out.jplace"
    output:
        config["outdir"] + "{mag}_epa_out/profile.tsv"
    conda:
        "gappa"
    threads: 1
    params:
        out_dir=config["outdir"]
    priority: 0
    log:
        config["outdir"] + "logs/gappa/{mag}.log"
    shell:
        """
        gappa examine assign \
            --jplace-path {input} \
            --taxon-file /compact_classification/resources/tax_tree.txt \
            --out-dir {params.out_dir}{wildcards.mag}_epa_out \
            --allow-file-overwriting --best-hit --verbose > {log}
        """
