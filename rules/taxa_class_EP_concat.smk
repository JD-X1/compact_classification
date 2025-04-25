#!usr/bin/env python

import os
import pandas as pd


# resources mags/

def get_genes_from_goneFishing(mag):
    checkpoint_output = f"resources/{mag}_working_dataset"

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

def get_superMatrix_targets(mag):

    sg_alns = os.listdir(f"resources/{mag}_mafft_out/")
    genes = []
    for f in sg_alns:
        if f.endswith(".aln"):
            gene = f.split(".")[0]
            genes.append(gene)
    return genes


def all_input():
    # Initialize an empty list to collect all inputs
    all_inputs = []
    
    # Iterate over each MAG to generate inputs
    for mag in mags:
        # Dynamically get genes for the current MAG after goneFishing
        # Check if the checkpoint directory for this MAG exists to infer completion
        checkpoint_dir = f"resources/{mag}_working_dataset"
        mag_inputs = [
            f"resources/busco_out/{mag}/summary.txt",
            f"resources/busco_out/{mag}/eukaryota_odb12/translated_protein.fasta",
            f"resources/{mag}_input_metadata.tsv",
            #f"{mag}_summary.csv",
            f"{mag}_SuperMatrix.fas",
            f"{mag}_q.aln",
            f"{mag}_ref.aln",
            f"resources/{mag}_epa_out/{mag}_epa_out.jplace",
            f"resources/{mag}_epa_out/profile.tsv",
            checkpoint_dir  # This implicitly checks for its existence
        ]

        # Check if the checkpoint has completed by checking the existence of its output directory
        if os.path.exists(checkpoint_dir):
            # Assuming get_genes_from_goneFishing is adjusted to accept a mag parameter
            genes = get_genes_from_goneFishing(mag)
            for gene in genes:
                mag_inputs.extend([
                    f"resources/{mag}_q_frags/{gene}.fas",
                    f"resources/{mag}_mafft_out/{gene}.aln"
                ])
        
        # Add the MAG-specific inputs to the overall list
        all_inputs.extend(mag_inputs)
    
    return all_inputs


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
mags = []
genes = []
for f in mag_f:
    # get mag name
    mag = f.strip(".fna")
    # append to list
    mags.append(mag)

rule all:
    input:
        all_input()
    
rule run_busco:
    input: 
        config["mag_dir"] + "{mag}.fna"
    output:
       "resources/busco_out/{mag}/summary.txt",
       "resources/busco_out/{mag}/eukaryota_odb12/translated_protein.fasta"
    conda:
        "../envs/mb.yaml"
    threads: 22
    priority: 0
    log:
        "logs/busco/{mag}.log"
    shell:
        """
        compleasm run -a {input} -t {threads} -l eukaryota -L resources/mb_downloads/ -o resources/busco_out/{wildcards.mag} 1> {log} 2> {log}
        """

rule fishing_meta:
    input:
        "resources/busco_out/{mag}/eukaryota_odb12/translated_protein.fasta"
    output:
        "resources/{mag}_input_metadata.tsv"
    conda:
        "../envs/raxml-ng.yaml"
    threads: 22
    priority: 0
    shell:
        "python ./additional_scripts/fishing_meta.py {input} >> {output}"


checkpoint goneFishing:
    input:
        "resources/{mag}_input_metadata.tsv"
    output:
        directory("resources/{mag}_working_dataset")
    conda:
        "../envs/fisher.yaml"
    threads: 22
    priority: 0
    shell:
        "sh ./additional_scripts/fishing.sh -t {threads} -i {input}"


rule splitter:
    input:
        dir="resources/{mag}_working_dataset/",
        tar="resources/{mag}_working_dataset/{gene}.fas",
        mag_dir=config["mag_dir"]
    output:
       "resources/{mag}_q_frags/{gene}.fas"
    conda:
        "../envs/raxml-ng.yaml"
    threads: 1
    priority: 0
    log:
        "logs/splitter/{mag}_{gene}.log"
    shell:
        "python ./additional_scripts/splitter.py -i {input.tar} -d {input.mag_dir} -o {output} 1> {log} 2> {log}"


rule mafft:
    input:
        query="resources/{mag}_q_frags/{gene}.fas",
        reference="resources/PhyloFishScratch/alignments/{gene}.fas.aln"
    output:
        "resources/{mag}_mafft_out/{gene}.aln"
    conda:
        "../envs/raxml-ng.yaml"
    threads: 22
    priority: 0
    log:
        "logs/mafft/{mag}/{mag}_{gene}_mafft.log"
    shell:
        "mafft --auto --addfragments {input.query} --keeplength --thread {threads} {input.reference} > {output} 2> {log}"

superMatrix_targets = get_superMatrix_targets(mag)
print(len(superMatrix_targets))

rule divvier:
    input:
        "resources/{mag}_mafft_out/{gene}.aln"
    output:
        "resources/{mag}_mafft_out/{gene}.aln.partial.fas",
        "resources/{mag}_mafft_out/{gene}.aln.PP"
    log:
        "logs/divvier/{mag}/{mag}_{gene}_divvier.log"
    conda:
        "../envs/div.yaml"
    threads: 1
    shell:
        """
        divvier -mincol 4 -partial -divvygap {input} > {log} 2> {log}
        """

rule trimal:
    input:
        "resources/{mag}_mafft_out/{gene}.aln.partial.fas"
    output:
        "resources/{mag}_mafft_out/{gene}.trimal"
    log:
        "logs/trimal/{mag}/{mag}_{gene}_trimal.log"
    conda:
        "../envs/trimal.yaml"
    threads: 1
    shell:
        """
        trimal -in {input} -gt 0.8 -out {output} > {log} 2> {log}
        """

rule concat:
    input:
        expand("resources/{mag}_mafft_out/{gene}.trimal", gene=superMatrix_targets, mag=mags)
    output:
        "{mag}_SuperMatrix.fas"
    conda:
        "../envs/pythas_two.yaml"
    threads: 1
    priority: 0
    log:
        "logs/concat/{mag}.log"
    shell:
        """
        FIXED_ALNS=()
        mkdir -p resources/{wildcards.mag}_relabeled
        for i in $(realpath resources/{wildcards.mag}_mafft_out/*trimal)
        do
        prot=$(basename ${{i}} .trimal)
        python additional_scripts/add_gene_name.py -a ${{i}} -g ${{prot}} -t {wildcards.mag} -o resources/{wildcards.mag}_relabeled/${{prot}}.fas
        FIXED_ALNS+=("resources/{wildcards.mag}_relabeled/${{prot}}.fas")
        done
        python2 additional_scripts/geneStitcher.py -in ${{FIXED_ALNS[@]}} 
        mv SuperMatrix.fas {output}
        """

rule alignment_splitter:
    input:
        "{mag}_SuperMatrix.fas"
    output:
        "{mag}_q.aln"
        "{mag}_ref.aln"
    conda:
        "../envs/raxml-ng.yaml"
    threads: 1
    priority: 0
    log:
        "logs/alignment_splitter/{mag}.log"
    shell:
        "python alignment_splitter.py -a {input} -t {wildcards.mag}"


rule raxml_epa:
    input:
        q_aln="{mag}_q.aln",
        ref_aln="{mag}_ref.aln",
        ref_tree="resources/ref_concat.tre"
    output:
        "resources/{mag}_epa_out/{mag}_epa_out.jplace"
    conda:
        "../envs/raxml-ng.yaml"
    threads: 22
    priority: 0
    log:
        "logs/raxml_epa/{mag}_epa.log"
    shell:
        """
        # Ensure the output directory exists
        mkdir -p resources/{wildcards.mag}_epa_out/
        mkdir -p logs/raxml_epa/
        # Run RAxML
        cd resources/{wildcards.mag}_epa_out/

        # Execute RAxML with the correct path to the alignment and tree files,
        # specifying only a run name for the output (not a path).
        ulimit -n 65536
        ulimit -s unlimited
        epa-ng --ref-msa {mag}_ref.aln \
         --tree resources/ref_concat_pruned.tre \
         --query {mag}_q.aln \
         --model LG -T {threads}
        # The actual output file is named according to the RAxML naming convention,
        # incorporating the run name. Ensure this matches your output specification.
        mv epa_result.jplace {output}
        mv epa_info.log {log}
        """

rule gappa:
    input:
        "resources/{mag}_epa_out/{mag}_epa_out.jplace"
    output:
        "resources/{mag}_epa_out/profile.tsv"
    conda:
        "../envs/raxml-ng.yaml"
    threads: 22
    priority: 0
    log:
        "logs/gappa/{mag}.log"
    shell:
        """
        gappa examine assign \
            --jplace-path {input} \
            --taxon-file resources/tax_tree.txt \
            --out-dir resources/{wildcards.mag}_epa_out \
            --allow-file-overwriting --best-hit --verbose > {log}
        """
