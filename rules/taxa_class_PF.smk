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

# check if config["outdir"] ends with a slash if it doesn't add one
if not config["outdir"].endswith("/"):
    config["outdir"] += "/"

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
        expand(config["outdir"] + "{mag}_epa_out/profile_summary.tsv", mag=mags),
        expand(config["outdir"] + "{mag}_summary.csv", mag=mags)
    
rule run_busco:
    input: 
        config["mag_dir"] + "{mag}.fna"
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
        echo 'Running compleasm with available resources'
        compleasm run -a {input} -t {threads} -l eukaryota -L resources/mb_downloads/ -o {config[outdir]}busco_out/{wildcards.mag} 1> {log} 2> {log}
        """

rule fishing_meta:
    input:
        config["outdir"] + "busco_out/{mag}/eukaryota_odb12/translated_protein.fasta"
    output:
        config["outdir"] + "{mag}_input_metadata.tsv"
    conda:
        "pline_max"
    shell:
        "python ./additional_scripts/fishing_meta.py {input} >> {output}"


checkpoint goneFishing:
    input:
        config["outdir"] + "{mag}_input_metadata.tsv"
    output:
        directory(config["outdir"] + "{mag}_working_dataset")
    conda:
        "fisher"
    threads: 8
    shell:
        "sh ./additional_scripts/fishing.sh -t {threads} -i {wildcards.mag}_input_metadata.tsv -o {config[outdir]}"

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
    log:
        config["outdir"] + "logs/splitter/{mag}_{gene}.log"
    shell:
        "python ./additional_scripts/splitter.py -i {input.tar} -d {input.mag_dir} -o {output} 1> {log} 2> {log}"


rule mafft:
    input:
        query=config["outdir"] + "{mag}_q_frags/{gene}.fas",
        reference=config["outdir"] + "{mag}_PhyloFishScratch/alignments/{gene}.fas.aln"
    output:
        config["outdir"] + "{mag}_mafft_out/{gene}.aln"
    conda:
        "pline_max"
    threads: 20
    log:
        config["outdir"] + "logs/mafft/{mag}/{mag}_{gene}_mafft.log"
    shell:
        "mafft --auto --addfragments {input.query} --keeplength --thread {threads} {input.reference} > {output} 2> {log}"

rule split_aln:
    input:
        config["outdir"] + "{mag}_mafft_out/{gene}.aln"
    output:
        query=config["outdir"] + "{mag}_mafft_out/{gene}/{gene}_q.aln",
        ref=config["outdir"] + "{mag}_mafft_out/{gene}/{gene}_ref.aln"
    conda:
        "pline_max"
    threads: 1
    log:
        config["outdir"] + "logs/split_aln/{mag}/{gene}.log"
    shell:
        """
        mkdir -p {config["outdir"]}{wildcards.mag}_mafft_out/{wildcards.gene}
        python additional_scripts/alignment_splitter.py -a {input} -t {wildcards.mag} -g {wildcards.gene}> {log}
        mv {wildcards.gene}_q.aln {output.query}
        mv {wildcards.gene}_ref.aln {output.ref}
        """

 
rule raxml_epa:
    input:
        q_aln=config["outdir"] + "{mag}_mafft_out/{gene}/{gene}_q.aln",
        ref_aln=config["outdir"] + "{mag}_mafft_out/{gene}/{gene}_ref.aln",
        ref_tree=config["outdir"] + "ref_trees_PF/{gene}/{gene}.raxml.support" #rename this once final tree files available
    output:
        config["outdir"] + "{mag}_epa_out/{gene}/{mag}_epa_out.jplace"
    conda:
        "pline_max"
    threads: 20
    log:
        config["outdir"] + "logs/epa/{mag}/{gene}.log"
    shell:
        """
        mkdir -p {config["outdir"]}{wildcards.mag}_epa_out/{wildcards.gene}
        mkdir -p {config["outdir"]}logs/epa/{wildcards.mag}
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
        config["outdir"] + "{mag}_epa_out/{gene}/{mag}_epa_out.jplace"
    output:
        config["outdir"] + "{mag}_epa_out/{gene}/profile.tsv"
    conda:
        "pline_max"
    threads: 1
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
            config["outdir"] + "{mag}_epa_out/{gene}/profile.tsv",
            mag=[wildcards.mag],
            gene=get_superMatrix_targets
        )
    output:
        o1=config["outdir"] + "{mag}_epa_out/profile_summary.tsv",
        o2=config["outdir"] + "{mag}_summary.csv"
    conda:
        "pline_max"
    threads: 22
    shell:
        """
        cat {input} | grep -v "LWR" >> {output.o1}
        python additional_scripts/gappa_parse.py -i {output.o1} -o {output.o2}
        """
