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

def get_superMatrix_targets(wildcards):
    ckpt_out = checkpoints.goneFishing.get(mag=wildcards.mag).output[0]
    genes = glob_wildcards(os.path.join(ckpt_out, "{gene}.fas")).gene
    return genes

def sanitize_gene_name(name):
    return name.replace("/", "_").replace(" ", "_").replace("\\", "_")

def get_superMatrix_targets_for_mag(mag):
    ckpt_out = checkpoints.goneFishing.get(mag=mag).output[0]
    gene_files = glob_wildcards(os.path.join(ckpt_out, "{gene}.fas")).gene
    return [sanitize_gene_name(gene) for gene in gene_files]


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
            "resources/PhyloFishScratch/",
            checkpoint_dir  # This implicitly checks for its existence
            ]

        # Check if the checkpoint has completed by checking the existence of its output directory
        if os.path.exists(checkpoint_dir):
            print("################# Checkpoint completed")
            print("################# Checkpoint completed")
            print("################# Checkpoint completed")
            print("################# Checkpoint completed")

            # Assuming get_genes_from_goneFishing is adjusted to accept a mag parameter
            genes = get_genes_from_goneFishing(mag)
            for gene in genes:
                mag_inputs.extend([
                    f"resources/{mag}_q_frags/{gene}.fas",
                    f"resources/{mag}_mafft_out/{gene}.aln",
                    f"resources/{mag}_epa_out/{gene}/RAxML_portableTree.epa.jplace",
                    f"resources/{mag}_epa_out/{gene}/profile.tsv"
                ])
            mag_inputs.extend([
                f"{mag}_summary.csv", 
                f"resources/{mag}_epa_out/profile_summary.tsv",
                ])
        # Add the MAG-specific inputs to the overall list
        all_inputs.extend(mag_inputs)
    
    return all_inputs


mag_f = os.listdir(config["mag_dir"])
print(mag_f)
#### CHANGE THE FILE EXTENSION
#### Swap out this check for file extension
#### with a test to check that the MAG files
#### are actually nucliec acid fastas. 
mag_f = [f for f in mag_f if f.endswith(".fna")]
if mag_f == []:
    raise ValueError("#################\nNo MAGs found in the specified directory.\n#################\n")

# get mag names
mags = []
for f in mag_f:
    # get mag name
    mag = f.split(".")[0]
    # append to list
    mags.append(mag)

rule all:
    input:
        expand("resources/{mag}_working_dataset", mag=mags),
        expand(
            "resources/{mag}_q_frags/{gene}.fas",
            mag=mags,
            gene=lambda wildcards: [
                gene
                for mag in mags
                for gene in get_superMatrix_targets_for_mag(mag)
            ]
        ),
        expand(
            "resources/{mag}_mafft_out/{gene}.aln",
            mag=mags,
            gene=lambda wildcards: [
                gene
                for mag in mags
                for gene in get_superMatrix_targets_for_mag(mag)
            ]
        ),
        expand("resources/{mag}_epa_out/profile_summary.tsv", mag=mags),
        expand("{mag}_summary.csv", mag=mags)
    
rule run_busco:
    input: 
        config["mag_dir"] + "{mag}.fna"
    output:
       "resources/busco_out/{mag}/summary.txt",
       "resources/busco_out/{mag}/eukaryota_odb12/translated_protein.fasta"
    conda:
        "compleasm"
    
    threads: 22
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
        "pline_max"
    shell:
        "python ./additional_scripts/fishing_meta.py {input} >> {output}"


checkpoint goneFishing:
    input:
        "resources/{mag}_input_metadata.tsv"
    output:
        directory("resources/{mag}_working_dataset")
    resources:
        shell_prefix="conda run --no-capture-output -n fisher"
    conda:
        "fisher"
    threads: 22
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
        "pline_max"
    threads: 22
    log:
        "logs/splitter/{mag}_{gene}.log"
    shell:
        "python ./additional_scripts/splitter.py -i {input.tar} -d {input.mag_dir} -o {output} 1> {log} 2> {log}"


rule mafft:
    input:
        query="resources/{mag}_q_frags/{gene}.fas",
        reference="resources/{mag}_PhyloFishScratch/alignments/{gene}.fas.aln"
    output:
        "resources/{mag}_mafft_out/{gene}.aln"
    conda:
        "pline_max"
    threads: 22
    log:
        "logs/mafft/{mag}/{mag}_{gene}_mafft.log"
    shell:
        "mafft --auto --addfragments {input.query} --keeplength --thread {threads} {input.reference} > {output} 2> {log}"

rule split_aln:
    input:
        "resources/{mag}_mafft_out/{gene}.aln"
    output:
        query="resources/{mag}_mafft_out/{gene}/{gene}_q.aln",
        ref="resources/{mag}_mafft_out/{gene}/{gene}_ref.aln"
    conda:
        "pline_max"
    threads: 1
    log:
        "logs/split_aln/{mag}/{gene}.log"
    shell:
        """
        mkdir -p resources/{wildcards.mag}_mafft_out/{wildcards.gene}
        python additional_scripts/alignment_splitter.py -a {input} -t {wildcards.mag} -g {wildcards.gene}> {log}
        mv {wildcards.gene}_q.aln {output.query}
        mv {wildcards.gene}_ref.aln {output.ref}
        """

 
rule raxml_epa:
    input:
        q_aln="resources/{mag}_mafft_out/{gene}/{gene}_q.aln",
        ref_aln="resources/{mag}_mafft_out/{gene}/{gene}_ref.aln",
        ref_tree="resources/ref_trees_PF/{gene}/{gene}.raxml.support" #rename this once final tree files available
    output:
        "resources/{mag}_epa_out/{gene}/{mag}_epa_out.jplace"
    conda:
        "pline_max"
    threads: 22
    log:
        "logs/epa/{mag}/{gene}.log"
    shell:
        """
        mkdir -p resources/{wildcards.mag}_epa_out/{wildcards.gene}
        mkdir -p logs/epa/{wildcards.mag}
        ulimit -n 65536
        ulimit -s unlimited
        epa-ng --redo --ref-msa {input.ref_aln} \
         --tree {input.ref_tree} \
         --query {input.q_aln} \
         --model LG -T {threads} > {log} 2> {log}
        mv epa_result.jplace {output}
        mv epa_info.log {log}
        """

rule gappa:
    input:
        "resources/{mag}_epa_out/{gene}/{mag}_epa_out.jplace"
    output:
        "resources/{mag}_epa_out/{gene}/profile.tsv"
    conda:
        "pline_max"
    threads: 22
    log:
        "logs/gappa/{mag}/{gene}.log"
    shell:
        """
        gappa examine assign \
            --jplace-path {input} \
            --taxon-file resources/tax_tree.txt \
            --out-dir resources/{wildcards.mag}_epa_out/{wildcards.gene} \
            --allow-file-overwriting --best-hit --verbose 1> {log} 2> {log}
        """


rule gappa_summary:
    input:
        lambda wildcards: expand(
            "resources/{mag}_epa_out/{gene}/profile.tsv",
            mag=[wildcards.mag],
            gene=get_superMatrix_targets
        )
    output:
        o1="resources/{mag}_epa_out/profile_summary.tsv",
        o2="{mag}_summary.csv"
    conda:
        "pline_max"
    threads: 22
    shell:
        """
        cat {input} | grep -v "LWR" >> {output.o1}
        python additional_scripts/gappa_parse.py -i {output.o1} -o {output.o2}
        """
