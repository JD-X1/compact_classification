#!usr/bin/env python

import os
import pandas as pd


# resources mags/

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
                tree_file = os.path.join("resources/ref_trees", gene_name, gene_name + ".raxml.support")
                if os.path.exists(tree_file):
                    valid_genes.append(gene_name)

    return valid_genes

def get_superMatrix_targets(wildcards):
    ckpt_out = checkpoints.goneFishing.get(mag=wildcards.mag).output[0]
    genes = glob_wildcards(os.path.join(ckpt_out, "{gene}.fas")).gene
    return genes

def sanitize_gene_name(name):
    return name.replace("/", "_").replace(" ", "_").replace("\\", "_")

def get_superMatrix_targets_for_mag(mag):
    ckpt_out = checkpoints.goneFishing.get(mag=mag).output[0]
    gene_files = glob_wildcards(os.path.join(ckpt_out, "{gene}.fas")).gene
    print([sanitize_gene_name(gene) for gene in gene_files])
    return [sanitize_gene_name(gene) for gene in gene_files]


def all_input():
    # Initialize an empty list to collect all inputs
    all_inputs = []
    
    # Iterate over each MAG to generate inputs
    for mag in mags:
        # Dynamically get genes for the current MAG after goneFishing
        # Check if the checkpoint directory for this MAG exists to infer completion
        checkpoint_dir = config["outdir"] + "{mag}_working_dataset"
        mag_inputs = [
            config["outdir"] + "busco_out/{mag}/summary.txt",
            config["outdir"] + "busco_out/{mag}/eukaryota_odb12/translated_protein.fasta",
            config["outdir"] + "{mag}_input_metadata.tsv",
            #f"{mag}_summary.csv",
            f"{mag}_SuperMatrix.fas",
            f"{mag}_q.aln",
            f"{mag}_ref.aln",
            config["outdir"] + "{mag}_epa_out/{mag}_epa_out.jplace",
            config["outdir"] + "{mag}_epa_out/profile.tsv",
            checkpoint_dir  # This implicitly checks for its existence
        ]

        # Check if the checkpoint has completed by checking the existence of its output directory
        if os.path.exists(checkpoint_dir):
            # Assuming get_genes_from_goneFishing is adjusted to accept a mag parameter
            genes = get_genes_from_goneFishing(mag)
            for gene in genes:
                mag_inputs.extend([
                    config["outdir"] + "{mag}_q_frags/{gene}.fas",
                    config["outdir"] + "{mag}_mafft_out/{gene}.aln"
                ])
        
        # Add the MAG-specific inputs to the overall list
        all_inputs.extend(mag_inputs)
    
    return all_inputs

# checking if config["outdir"] ends with a slash if it doesn't add one
if not config["outdir"].endswith("/"):
    config["outdir"] += "/"
if not config["mag_dir"].endswith("/"):
    config["mag_dir"] += "/"
outdir = config["outdir"]
mag_f = os.listdir(config["mag_dir"])
#print(mag_f)
#### CHANGE THE FILE EXTENSION
#### Swap out this check for file extension
#### with a test to check that the MAG files
#### are actually nucliec acid fastas. 
mag_f = [f for f in mag_f if f.endswith(".fna")]
if mag_f == []:
    raise ValueError("#################\nNo MAGs found in the specified directory.\n#################\n")

# get mag names
genes = []
mags = [os.path.splitext(f)[0] for f in mag_f] 
# for f in mag_f:
#     mag = ""
#     if len(f.split(".")) >= 2:
#         # get mag name
#         mag = '.'.join(f.split(".")[:-1])
#     # append to list
#     else:
#         raise ValueError(
#             "#################\nMake sure MAG file names follow the folowwing format:\n\t[unique id].[file extension] \n#################\n"
#         )
#     print(f"Processing MAG: " + mag)
#     mags.append(mag)

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
        expand(config["mag_dir"] + "{mag}.fna", mag=mags)
    output:
       config["outdir"] + "busco_out/{mag}/summary.txt",
       config["outdir"] + "busco_out/{mag}/eukaryota_odb12/translated_protein.fasta"
    conda:
        "compleasm"
    threads: 22
    log:
        config["outdir"] + "logs/busco/{mag}.log"
    shell:
        """
        echo "Running BUSCO for {wildcards.mag}..."
        echo "Running BUSCO with available threads: {threads}"
        compleasm run -a {input} -t {threads} -l eukaryota -L resources/mb_downloads/ -o {config[outdir]}busco_out/{wildcards.mag} 1> {log} 2> {log}
        """

rule fishing_meta:
    input:
        config["outdir"] + "busco_out/{mag}/eukaryota_odb12/translated_protein.fasta"
    output:
        config["outdir"] + "{mag}_input_metadata.tsv"
    conda:
        "pline_max"
    threads: 1
    priority: 0
    shell:
        "python ./additional_scripts/fishing_meta.py {input} >> {output}"


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
        "sh ./additional_scripts/fishing.sh -t {threads} -i {input}"


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
        "python ./additional_scripts/splitter.py -i {input.tar} -d {input.mag_dir} -o {output} 1> {log} 2> {log}"


rule mafft:
    input:
        query=config["outdir"] + "{mag}_q_frags/{gene}.fas",
        reference=config["outdir"] + "PhyloFishScratch/alignments/{gene}.fas.aln"
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
        lambda wildcards: [
            f"{config['outdir']}{wildcards.mag}_mafft_out/{gene}.trimal"
            for gene in get_superMatrix_targets_for_mag(wildcards.mag)
        ]
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
        python additional_scripts/add_gene_name.py -a ${{i}} -g ${{prot}} -t {wildcards.mag} -o {params.out_dir}{wildcards.mag}_relabeled/${{prot}}.fas
        FIXED_ALNS+=("{params.out_dir}{wildcards.mag}_relabeled/${{prot}}.fas")
        done
        python2 additional_scripts/geneStitcher.py -in ${{FIXED_ALNS[@]}} 
        mv SuperMatrix.fas {output}
        """

rule alignment_splitter:
    input:
        config["outdir"] + "{mag}_SuperMatrix.fas"
    output:
        config["outdir"] + "{mag}_q.aln",
        config["outdir"] + "{mag}_ref.aln"
    conda:
        "pline_max"
    threads: 1
    priority: 0
    log:
        config["outdir"] + "logs/alignment_splitter/{mag}.log"
    shell:
        "python additional_scripts/alignment_splitter.py -a {input} -t {wildcards.mag}"


rule raxml_epa:
    input:
        q_aln= config["outdir"] + "{mag}_q.aln",
        ref_aln= config["outdir"] + "{mag}_ref.aln",
        ref_tree= "resources/ref_concat_PF.tre"
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
        mkdir -p logs/raxml_epa/
        # Run RAxML
        cd {params.out_dir}{wildcards.mag}_epa_out/

        # Execute RAxML with the correct path to the alignment and tree files,
        # specifying only a run name for the output (not a path).
        ulimit -n 65536
        ulimit -s unlimited
        epa-ng --ref-msa {params.out_dir}{mag}_ref.aln \
         --tree resources/ref_concat_PF.tre \
         --query {params.out_dir}{mag}_q.aln \
         --model LG -T {threads}
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
        "pline_max"
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
            --taxon-file {params.out_dir}tax_tree.txt \
            --out-dir {params.out_dir}{wildcards.mag}_epa_out \
            --allow-file-overwriting --best-hit --verbose > {log}
        """
