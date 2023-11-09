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
        "resources/dataset",
        "resources/fisher_out",
        "config.ini"
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
    output:
        "config.ini"
    conda:
        "../envs/fisher.yaml"
    threads: 22
    shell:
        "python ./additional_scripts/fishing_meta.py resources/busco_out"

rule fisher:
    output:
        directory("resources/fisher_out"),
    conda:
        "../envs/fisher.yaml"
    threads: 22
    shell:
        "fisher.py -t {threads} -o {output}"

rule work_dataset:
    input:
        "resources/fisher_out"
    output:
        directory("resources/dataset")
    conda:
        "../envs/fisher.yaml"
    threads: 22
    shell:
        """
        informant.py -i {input}
        shell('working_dataset_constructory.py -i {input} -o {output}
        """




#rule splitter:
#    input:
#        wd=rules.fishing_configurations.output,
#        meta="resources/tax_tree.txt"
#    output:
#        query=directory("resources/queries/")
#    conda:
#        "../envs/raxml-ng.yaml"
#    shell:
#        "python ./additional_scripts/splitter.py -i {input.wd} -t {input.meta} -o {output.query}"

#rule align:
#    input:
#        query="resources/queries/{protein}.q.fas",
#        reference="resources/PhyloFishScratch/"
#    output:
#        "resources/q_alns/{protein}.q.aln"
#    log:
#        "resources/logs/mafft/{protein}.log"
#    conda:
#        "../envs/raxml-ng.yaml"
#    threads: 22
#    shell:
#        "mafft --auto --thread {threads} --addfragments {input.query} --keeplength --globalpair {input.reference} 1> {output} 2> {log}"


#rule ref_trees:
#    input:
#        "resources/ref_alns_reduce/{protein}.final"
#    output:
#        'resources/ref_trees/{protein}/{protein}.raxml.rba',
#        'resources/ref_trees/{protein}/{protein}.raxml.startTree',
#        'resources/ref_trees/{protein}/{protein}.raxml.bestTree',
#        'resources/ref_trees/{protein}/{protein}.raxml.mlTrees',
#        'resources/ref_trees/{protein}/{protein}.raxml.bestModel',
#        'resources/ref_trees/{protein}/{protein}.raxml.support'
#    log:
#        "resources/logs/ref_trees/{protein}.log"
#    params:
#        raxml_out="resources/ref_trees/"
#    conda:
#        "../envs/raxml-ng.yaml"
#    threads: 22
#    shell:
#        "raxml-ng --all --msa {input} --model JTT --prefix {params.raxml_out}/{wildcards.protein}/{wildcards.protein} --tree pars {{10}} --threads {threads} 1> {log} 2> {log}"


#rule raxml_epa:
#    input:
#        query_aln='resources/q_alns/{protein}.q.aln',
#        ref_tree='resources/ref_trees/{protein}/{protein}.raxml.support'
#    output:
#        'resources/epa_out/{protein}/RAxML_portableTree.{protein}.epa.jplace',
#        'resources/epa_out/{protein}/RAxML_info.{protein}.epa',
#        'resources/epa_out/{protein}/RAxML_entropy.{protein}.epa',
#        'resources/epa_out/{protein}/RAxML_classificationLikelihoodWeights.{protein}.epa',
#        'resources/epa_out/{protein}/RAxML_classification.{protein}.epa',
#        'resources/epa_out/{protein}/RAxML_labelledTree.{protein}.epa',
#        'resources/epa_out/{protein}/RAxML_originalLabelledTree.{protein}.epa'
#    log:
#        "resources/logs/raxml-epa/{protein}.log"
#    params:
#        raxml_out="/scratch/jduque2/compact_class/resources/epa_out/"
#    conda:
#        "../envs/raxml-ng.yaml"
#    threads: 22
#    shell:
#        "raxmlHPC-PTHREADS -f v -T {threads} -s {input.query_aln} -t {input.ref_tree} -m PROTGAMMALG4X -n {wildcards.protein}.epa -w {params.raxml_out}{wildcards.protein} 1> {log} 2> {log}"
