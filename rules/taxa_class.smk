#!usr/bin/env python

import os
import pandas as pd

# read in all single gene alignments and extract
# gene names for snakemake rule all. Alignments
# are in fasta format and are in the directory
# resources queries/
# files names are formated as gene.q.fas

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
        expand("resources/busco_out/{mag}/summary.txt", mag=mags),
        "resources/input_metadata.tsv",
        "resources/PhyloFishScratch/",
        "resources/working_dataset"
        #expand("resources/q_alns/{protein}.q.aln", protein=proteins),
        #expand("resources/ref_trees/{protein}/{protein}.raxml.rba", protein=proteins),
        #expand("resources/ref_trees/{protein}/{protein}.raxml.startTree", protein=proteins),
        #expand("resources/ref_trees/{protein}/{protein}.raxml.bestTree", protein=proteins),
        #expand("resources/ref_trees/{protein}/{protein}.raxml.mlTrees", protein=proteins),
        #expand("resources/ref_trees/{protein}/{protein}.raxml.bestModel", protein=proteins),
        #expand("resources/ref_trees/{protein}/{protein}.raxml.support", protein=proteins),
        #expand("resources/epa_out/{protein}/RAxML_portableTree.{protein}.epa.jplace", protein=proteins),
        #expand("resources/epa_out/{protein}/RAxML_info.{protein}.epa", protein=proteins),
        #expand("resources/epa_out/{protein}/RAxML_entropy.{protein}.epa", protein=proteins),
        #expand("resources/epa_out/{protein}/RAxML_classificationLikelihoodWeights.{protein}.epa", protein=proteins),
        #expand("resources/epa_out/{protein}/RAxML_classification.{protein}.epa", protein=proteins),
        #expand("resources/epa_out/{protein}/RAxML_labelledTree.{protein}.epa", protein=proteins),
        #expand("resources/epa_out/{protein}/RAxML_originalLabelledTree.{protein}.epa", protein=proteins)
    
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



rule goneFishing:
    input:
        "resources/input_metadata.tsv"
    output:
        directory("resources/working_dataset")
    conda:
        "../envs/fisher.yaml"
    threads: 22
    shell:
        "sh ./additional_scripts/fishing.sh -t {threads} -i {input}"

#def get_sg_unalns(wildcards):
    # function returns the unaligned single gene
    # fasta files contained in resources/working_dataset/
    # target files end in 
#    return [os.path.]

#rule splitter:
#    input:
#        get_sg_unalns