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
            f"resources/busco_out/{mag}/eukaryota_odb12/translated_protein.fasta", # need to make this persistent across newer verisions of this db
            f"resources/{mag}_input_metadata.tsv",
            "resources/PhyloFishScratch/",
            f"resources/{mag}_epa_out/profile_summary.tsv",
            checkpoint_dir  # This implicitly checks for its existence
        ]

        # Check if the checkpoint has completed by checking the existence of its output directory
        if os.path.exists(checkpoint_dir):
            # Assuming get_genes_from_goneFishing is adjusted to accept a mag parameter
            mag_inputs.append(f"{mag}_summary.csv")
            genes = get_genes_from_goneFishing(mag)
            for gene in genes:
                mag_inputs.extend([
                    f"resources/{mag}_q_frags/{gene}.fas",
                    f"resources/{mag}_mafft_out/{gene}.aln",
                    f"resources/{mag}_epa_out/{gene}/RAxML_portableTree.epa.jplace",
                    f"resources/{mag}_epa_out/{gene}/profile.tsv"
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
genes = []
for f in mag_f:
    # get mag name
    mag = f.split(".")[0]
    # append to list
    mags.append(mag)
print(config["mag_dir"] + "{mag}.fna")

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
        reference="resources/EPDB/alns/{gene}.fas.aln.fixed"
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

rule raxml_epa:
    input:
        q_aln="resources/{mag}_mafft_out/{gene}.aln",
        ref_tree="resources/EPDB/trees/{gene}/{gene}.initial.treefile" #rename this once final tree files available
    output:
        "resources/{mag}_epa_out/{gene}/RAxML_portableTree.epa.jplace"
    conda:
        "../envs/raxml-ng.yaml"
    threads: 22
    priority: 0
    log:
        "logs/raxml_epa/{mag}/{gene}.log"
    shell:
        """
        # Ensure the output directory exists
        mkdir -p resources/{wildcards.mag}_epa_out/{wildcards.gene}
        mkdir -p logs/raxml_epa/{wildcards.mag}
        # Run RAxML
        cd resources/{wildcards.mag}_epa_out/{wildcards.gene}

        # Execute RAxML with the correct path to the alignment and tree files,
        # specifying only a run name for the output (not a path).
        raxmlHPC-PTHREADS -f v -T {threads} \
          -s ../../../{input.q_aln} \
          -t ../../../{input.ref_tree} \
          -m PROTGAMMAJTT -n epa 

        # The actual output file is named according to the RAxML naming convention,
        # incorporating the run name. Ensure this matches your output specification.
        """

rule gappa:
    input:
        "resources/{mag}_epa_out/{gene}/RAxML_portableTree.epa.jplace"
    output:
        "resources/{mag}_epa_out/{gene}/profile.tsv"
    conda:
        "../envs/raxml-ng.yaml"
    threads: 22
    priority: 0
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
        "resources/{mag}_epa_out/{gene}/profile.tsv"
    output:
        o1="resources/{mag}_epa_out/profile_summary.tsv",
        o2="{mag}_summary.csv"
    conda:
        "../envs/raxml-ng.yaml"
    threads: 22
    priority: 0
    shell:
        """
        cat {input} | grep -v "LWR" >> {output.o1}
        python additional_scripts/gappa_parse.py -i {output.o1} -o {output.o2}
        """