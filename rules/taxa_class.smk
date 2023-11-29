#!usr/bin/env python

import os
import pandas as pd

# read in all single gene alignments and extract
# gene names for snakemake rule all. Alignments
# are in fasta format and are in the directory
# resources queries/
# files names are formated as gene.q.fas

def get_genes_from_goneFishing():
    # Get the output directory from the checkpoint
    checkpoint_output = checkpoints.goneFishing.get().output[0]

    # List files in the output directory
    files = os.listdir(checkpoint_output)

    # Initialize an empty list to store genes with corresponding trees
    valid_genes = []

    # Check each file to see if a corresponding tree file exists
    for f in files:
        if f.endswith('.fas'):
            gene_name = os.path.splitext(f)[0]
            tree_file = os.path.join("resources/ref_trees", gene_name, gene_name + ".raxml.support")
            if os.path.exists(tree_file):
                valid_genes.append(gene_name)

    return valid_genes


def all_input(wildcards):
    # Get the list of genes after the checkpoint
    genes = get_genes_from_goneFishing()

    # Return the expanded list of inputs
    return expand("resources/busco_out/{mag}/summary.txt", mag=mags) + \
           ["resources/input_metadata.tsv",
            "resources/PhyloFishScratch/",
            "resources/working_dataset"] + \
           expand('resources/q_frags/{gene}.fas' , gene=genes) + \
           expand("resources/mafft_out/{gene}.aln", gene=genes) + \
           expand("resources/epa_out/{gene}/RAxML_portableTree.epa.jplace", gene=genes)


mag_f = os.listdir("resources/mags/")
mag_f = [f for f in mag_f if f.endswith(".fna")]

# get mag names
mags = []
for f in mag_f:
    # get mag name
    mag = f.strip(".fna")
    # append to list
    mags.append(mag)

rule all:
    input:
        all_input

    
rule run_busco:
    input: 
        "resources/mags/{mag}.fna"
    output:
        directory("resources/busco_out/{mag}"),
        "resources/busco_out/{mag}/summary.txt"
    conda:
        "../envs/mb.yaml"
    threads: 22
    log:
        "logs/busco/{mag}.log"
    shell:
        "compleasm run -a {input} -t {threads} -l eukaryota -L resources/mb_downloads/ -o resources/busco_out/{wildcards.mag} 1> {log} 2> {log}"


rule fishing_configurations:
    input:
        "resources/busco_out/"
    output:
        "resources/input_metadata.tsv",
        directory("resources/PhyloFishScratch/")
    conda:
        "../envs/fisher.yaml"
    threads: 22
    shell:
        "python ./additional_scripts/fishing_meta.py {input}"



checkpoint goneFishing:
    input:
        "resources/input_metadata.tsv"
    output:
        directory("resources/working_dataset")
    conda:
        "../envs/fisher.yaml"
    threads: 22
    shell:
        "sh ./additional_scripts/fishing.sh -t {threads} -i {input}"



rule splitter:
    input:
        dir="resources/working_dataset/",
        tar="resources/working_dataset/{gene}.fas"
    output:
       "resources/q_frags/{gene}.fas"
    conda:
        "../envs/raxml-ng.yaml"
    threads: 22
    shell:
        "python ./additional_scripts/splitter.py -i {input.tar} -o {output}"


rule mafft:
    input:
        query="resources/q_frags/{gene}.fas",
        reference="resources/PhyloFishScratch/alignments/{gene}.fas.aln"
    output:
        "resources/mafft_out/{gene}.aln"
    conda:
        "../envs/raxml-ng.yaml"
    threads: 22
    shell:
        "mafft --auto --addfragments {input.query} --keeplength --thread {threads} {input.reference} > {output}"

rule raxml_epa:
    input:
        q_aln="resources/mafft_out/{gene}.aln",
        ref_tree="resources/ref_trees/{gene}/{gene}.raxml.support"
    output:
        "resources/epa_out/{gene}/RAxML_portableTree.epa.jplace"
    conda:
        "../envs/raxml-ng.yaml"
    threads: 22
    shell:
        """
        mkdir -p resources/epa_out/{wildcards.gene}
        cd resources/epa_out/{wildcards.gene}
        raxmlHPC-PTHREADS -f v -T {threads} -s ../../../{input.q_aln} -t ../../../{input.ref_tree} -m PROTGAMMAJTT -n epa
        cp {wildcards.gene}.jplace ../
        """
