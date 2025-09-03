#!usr/bin/env python

import os

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

# checking if config["outdir"] ends with a slash if it doesn't add one
if not config["outdir"].endswith("/"):
    config["outdir"] += "/"
if not config["mag_dir"].endswith("/"):
    config["mag_dir"] += "/"
outdir = config["outdir"]
mag_f = os.listdir(config["mag_dir"])

# Check if config["augustus"] is set, if not set it to default
augustus = False
if "augustus" in config:
    augustus = True
    print("Using Augustus for BUSCO runs.")

trim_alignments = False
if "trim" in config:
    trim_alignments = True
    print("Trimming alignments with trimAl + divvier.")

proteome_input = False
if "--proteome" in config:
    proteome_input = True
    print("Using proteome input instead of BUSCO Output.")

#### CHANGE THE FILE EXTENSION
#### Swap out this check for file extension
#### with a test to check that the MAG files
#### are actually nucleic acid fastas.
mag_f = [f for f in mag_f if f.endswith(".fna")]
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
        config["mag_dir"] + "{mag}.fna"
    output:
        branch(
            augustus,
            [config["outdir"] + "busco_out/{mag}/short_summary.specific.eukaryota_odb12.{mag}.txt", config["outdir"] + "busco_out/{mag}/eukaryota_odb12/translated_protein.fasta"],
            [config["outdir"] + "busco_out/{mag}/summary.txt", config["outdir"] + "busco_out/{mag}/eukaryota_odb12/translated_protein.fasta"]
       )
    conda:
        branch(augustus,
        "busco", # executes rule w/ buscos Waugustus
        "compleasm" # executes rule w/ compleasm's miniprot
        )
    threads: workflow.cores
    log:
        config["outdir"] + "logs/busco/{mag}.log"
    shell:
        branch(augustus,
        """
        echo "Running BUSCO w/ augustus for {wildcards.mag}..."
        echo "Running BUSCO with available threads: {threads}"
        busco -i {input} -m genome -l /compact_classification/resources/busco_downloads/lineages/eukaryota_odb12 -c {threads} -f --augustus -o {config[outdir]}busco_out/{wildcards.mag} 1> {log} 2> {log}
        mkdir -p {config[outdir]}busco_out/{wildcards.mag}/eukaryota_odb12/
        cat {config[outdir]}busco_out/{wildcards.mag}/run_eukaryota_odb12/augustus_output/*faa* >> {config[outdir]}busco_out/{wildcards.mag}/eukaryota_odb12/translated_protein.fasta
        """,
        """
        echo "Running compleasm for {wildcards.mag}..."
        echo "Running compleasm with available threads: {threads}"
        compleasm run -a {input} -t {threads} \
            -l eukaryota \
            -L /compact_classification/resources/mb_downloads/ \
            -o {config[outdir]}busco_out/{wildcards.mag} \
            1> {log} 2> {log}
        """)

rule fishing_meta:
    input:
        branch(
            proteome_input,
            "/input/{proteome}",
            config["outdir"] + "busco_out/{mag}/eukaryota_odb12/translated_protein.fasta"
            )
        
    output:
        config["outdir"] + "{mag}_input_metadata.tsv"
    conda:
        "pline_max"
    threads: 1
    priority: 0
    shell:
        "python /compact_classification/additional_scripts/fishing_meta.py -p {input} -o {output}"


checkpoint goneFishing:
    input:
        config["outdir"] + "{mag}_input_metadata.tsv"
    output:
        directory(config["outdir"] + "{mag}_working_dataset")
    conda:
        "fisher"
    threads: 22
    priority: 0
    shell:
        "sh /compact_classification/additional_scripts/fishing.sh -t {threads} -i {input} -o {config[outdir]}"


rule splitter:
    input:
        dir=config["outdir"] + "{mag}_working_dataset/",
        tar=config["outdir"] + "{mag}_working_dataset/{gene}.fas",
        mag_dir=config["mag_dir"]
    output:
       config["outdir"] + "{mag}_q_frags/{gene}.fas"
    conda:
        "pline_max"
    threads: 1
    priority: 0
    log:
        config["outdir"] + "logs/splitter/{mag}_{gene}.log"
    shell:
        "python /compact_classification/additional_scripts/splitter.py -i {input.tar} -d {input.mag_dir} -o {output} 1> {log} 2> {log}"


rule mafft:
    input:
        query=config["outdir"] + "{mag}_q_frags/{gene}.fas",
        reference=config["outdir"] + "{mag}_PhyloFishScratch/alignments/{gene}.fas.aln"
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
        python /compact_classification/additional_scripts/add_gene_name.py -a ${{i}} -g ${{prot}} -t {wildcards.mag} -o {params.out_dir}{wildcards.mag}_relabeled/${{prot}}.fas
        FIXED_ALNS+=("{params.out_dir}{wildcards.mag}_relabeled/${{prot}}.fas")
        done
        python2 /compact_classification/additional_scripts/geneStitcher.py -in ${{FIXED_ALNS[@]}} 
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
        out_dir=config["outdir"]
    log:
        config["outdir"] + "logs/alignment_splitter/{mag}.log"
    shell:
        "python /compact_classification/additional_scripts/alignment_splitter.py -a {input} -t {wildcards.mag} -o {params.out_dir} > {log} 2> {log}"

rule sub_tree:
    input:
        aln=config["outdir"] + "{mag}_ref.aln",
        tree="/compact_classification/resources/ref_concat_PF_alt3.tre" # uncomment for singularity
        #tree="resources/ref_concat_PF_alt3.tre"  # uncomment for direct snakemake execution
    output:
        config["outdir"] + "{mag}_ref.tre"
    conda:
        "dendropy"
    threads: 1
    priority: 0
    log:
        config["outdir"] + "logs/sub_tree/{mag}.log"
    shell:
        """
        python /compact_classification/additional_scripts/sub_tree.py -a {input.aln} -t {input.tree} -o {output} > {log} 2> {log}
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
