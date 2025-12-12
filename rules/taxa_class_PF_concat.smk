#!/usr/bin/env python

import os
import datetime
from pathlib import Path

# ------- Helpers + Constants ------- #

GENOME_EXTS   = [".fna", ".fa", ".fasta", ".fna.gz", ".fa.gz", ".fasta.gz"]
PROTEOME_EXTS = [".faa", ".faa.gz", ".aa.fa", ".aa.fasta"]


def ts():
    return "[{:%Y-%m-%d %H:%M:%S}]".format(datetime.datetime.now())

def log(msg: str) -> None:
    print(f"{ts()}: {msg}")

def sanitize_gene_name(name):
    return name.replace("/", "_").replace(" ", "_").replace("\\", "_")


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

def get_superMatrix_targets_for_mag(mag):
    ckpt_out = checkpoints.goneFishing.get(mag=mag).output[0]
    gene_files = glob_wildcards(os.path.join(ckpt_out, "{gene}.fas")).gene
    # print([gene for gene in gene_files])
    # print([sanitize_gene_name(gene) for gene in gene_files])
    # print("################################################################################")
    return [sanitize_gene_name(gene) for gene in gene_files]

def metaeuk_prefix(mag: str) -> str:
    # Set prefix for metaeuk out files
    return os.path.join(config["outdir"], "metaeuk", mag)

def metaeuk_proteome_path(mag: str) -> str:
    return metaeuk_prefix(mag) + ".faa"

def get_protein_source(wildcards):
    """
    Determine where proteome for this MAG come from.

    - If proteome_input is True, use user-provided proteome + proteome directory.
    - If gene_source == "metaeuk", metaeuk predicted proteome.
    - Else, BUSCO/compleasm predicted proteome.
    """
    mag = wildcards.mag
    if proteome_input:
        return find_mag_file(wildcards)
    if metaeuk_source:
        return metaeuk_proteome_path(mag)
    
    return os.path.join(
        config["outdir"],
        "busco_out",
        mag,
        "eukaryota_odb12",
        "translated_protein.fasta"
    )

def get_plm_proteome_input(wildcards):
    return get_protein_source(wildcards)

# -------------------------------------------- #
# -------    Resources & Pathing   ----------- #
# -------------------------------------------- #


log("Checking for resources directory...")
if os.path.exists("/compact_classification/resources/"):
    RESOURCES_DIR = "/compact_classification/resources/"
    log("Found resources directory at /compact_classification/resources/")
elif os.path.exists("resources/"):
    RESOURCES_DIR = "resources/"
    log("Found resources directory at resources/")
else:
    raise ValueError(
        f"[{ts()}]: If running from source and not Singularity, please ensure that "
        f"the 'resources' directory is present in the current working directory."
    )

log("Checking for additional scripts directory...")
if os.path.exists("/compact_classification/additional_scripts/"):
    ADDITIONAL_SCRIPTS_DIR = "/compact_classification/additional_scripts/"
    log("Found additional scripts directory at /compact_classification/additional_scripts/")
elif os.path.exists("additional_scripts/"):
    ADDITIONAL_SCRIPTS_DIR = "additional_scripts/"
    log("Found additional scripts directory at additional_scripts/")
else:
    raise ValueError(
        f"[{ts()}]: If running from source and not Singularity, please ensure that "
        f"the 'additional_scripts' directory is present and contains the necessary scripts."
    )

# ----------------------------------------------------------------------------------------------------

def get_genes_from_goneFishing(mag):
    checkpoint_output = config["outdir"] + f"{mag}_working_dataset"

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

# Config Normalization
# ---------------------------- #

output_default = os.path.join(os.getcwd(), "output/")
outdir = config.get("outdir", output_default)
outdir = os.path.abspath(outdir)
if not outdir.endswith(os.sep):
    outdir += os.sep
config["outdir"] = outdir

if "mag_dir" not in config:
    raise ValueError(
        "[{:%Y-%m-%d %H:%M:%S}]: No MAG directory specified in config file. Please specify 'mag_dir'.".format(datetime.datetime.now())
    )
mag_dir = os.path.abspath(config["mag_dir"])
if not mag_dir.endswith(os.sep):
    mag_dir += os.sep
config["mag_dir"] = mag_dir

log(f"Command invoked with the following options:")
log(f"Output directory: {config['outdir']}")
log(f"MAG directory: {config['mag_dir']}")

# Bool flags


augustus = bool(config.get("augustus", False))
trim_alignments = bool(config.get("trim", False))
proteome_input = bool(config.get("proteome", False))

# -------------------------------------------------------------------------------------------------

species_tree_flag = bool(config.get("species_tree", False))

log(f"Will use {'Augustus (BUSCO)' if augustus else 'Compleasm'} for BUSCO runs.")
if trim_alignments:
    log("Trimming alignments with trimAl + divvier.")
if proteome_input:
    log("Using proteome input instead of BUSCO Output.")

# Prot Source Logic Handling

gene_source = str(config.get("gene_source", "busco")).strip().lower()
metaeuk_source = (gene_source == "metaeuk")

if metaeuk_source and proteome_input:
    raise ValueError(
        "[{ts()}]: Config conflict: gene_source=metaeuk and proteome=True. "
        "Choose either "
    )

log(f"Using gene source: "
    f"{'pre-predicted external proteome' if proteome_input else gene_source}"
    )

if plmsearch_enabled:
    log("PLMsearch filtering is Enabled (will run PLMembedding + search).")

# -------------------------------------------------------------------------------------------------
# ----- Database Handling & Purging ------- #

DATABASE_TYPE = "PhyloFisher"
PF_DIR = os.path.join(RESOURCES_DIR, "PhyloFisherDatabase_v1.0")
EP_DIR = ""

## Parse database options

purge = False
purge_target = None

if "purge" in config and str(config["purge"]).strip():
    requested = str(config["purge"]).strip()
    purge_target = requested
    purge = True

    if DATABASE_TYPE == "PhyloFisher":
        meta_path = os.path.join(RESOURCES_DIR, "PhyloFisherDatabase_v1.0", "database", "metadata.tsv")
        with open(meta_path, "r") as f:
            header = f.readline().strip().split("\t")
            def col_idx(name: str):
                name = name.lower()
                for i, h in enumerate(header):
                    if h.strip().lower() == name:
                        return i
                return None

            uid_i = col_idx("Unique ID")
            lname_i = col_idx("Long Name")
            if uid_i is None or lname_i is None:
                raise ValueError(
                    f"[{ts()}]: Required columns 'Unique ID' and 'Long Name' not found in {meta_path}."
                )
            
            long2uids = {}
            all_uids = set()
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) <= max(uid_i, lname_i):
                    continue
                uid = parts[uid_i].strip()
                lname = parts[lname_i].strip()
                all_uids.add(uid)
                long2uids.setdefault(lname, set()).add(uid)
        
        if requested in all_uids:
            purge_target = requested
        elif requested in long2uids:
            uids = sorted(long2uids[requested])
            purge_target = ",".join(uids)
        else:
            raise ValueError(
                f"[{ts()}]: Specified purge target '{requested}' not found in database."
            )
    log(f"Will purge the following taxa from the database: {purge_target}")

# ---------------------------------------------------------------------------------------------------- #
# ---------- MAG Tracking & Rule Definitions --------- #
# ---------------------------------------------------------------------------------------------------- #

mag_files = [
    f for f in os.listdir(config["mag_dir"])
    if any(f.endswith(ext) for ext in (GENOME_EXTS + PROTEOME_EXTS))
]

if not mag_files:
    raise ValueError(
        "#################"
        "No MAG files found in the specified directory.\n"
        "#################\n"
        )

mags = []

for f in mag_files:
    parts = f.split(".")
    if len(parts) < 2:
        raise ValueError(
            "#################\n"
            "Make sure MAG file names follow the following format:\n"
            "   [unique id].[file extension] \n"
            "#################\n"
        )
    mag_name = ".".join(parts[:-1])
    if "_" in mag_name or " " in mag_name:
        new_mag_name = mag_name.replace("_", "").replace(" ", "")
        old_path = os.path.join(config["mag_dir"], f)
        new_path = os.path.join(config["mag_dir"], new_mag_name + "." + parts[-1])
        os.rename(old_path, new_path)
        log(f"Renaming MAG file to {f} -> {new_mag_name}.{parts[-1]}")
        mag_name = new_mag_name
        
    log(f"Processing MAG: {mag_name}")
    mags.append(mag_name)


rule all:
    input:
        expand(config["outdir"] + "{mag}_purged_taxa_check.complete", mag=mags),
        expand(config["outdir"] + "{mag}_working_dataset", mag=mags),
        lambda wildcards: [
            f"{config['outdir']}{mag}_q_frags/{gene}.fas"
            for mag in mags
            for gene in get_superMatrix_targets_for_mag(mag)
        ],
        lambda wildcards: [
            f"{config['outdir']}{mag}_mafft_out/{gene}.aln"
            for mag in mags
            for gene in get_superMatrix_targets_for_mag(mag)
        ],
        expand(config["outdir"] + "{mag}_q.aln", mag=mags),
        expand(config["outdir"] + "{mag}_ref.aln", mag=mags),
        expand(config["outdir"] + "{mag}_SuperMatrix.fas", mag=mags),
        expand(config["outdir"] + "{mag}_epa_out/{mag}_epa_out.jplace", mag=mags),
        expand(config["outdir"] + "{mag}_epa_out/profile.tsv", mag=mags),
        expand(config["outdir"] + "{mag}_epa_out/pairwise_qSeqDistance2leaves.tsv", mag=mags),
        expand(config["outdir"] + "species_tree/{mag}_species_tree.treefile", mag=mags) if species_tree_flag else [],
        expand(config["outdir"] + "plm/{mag}_plm_candidates.faa", mag=mags) if plmsearch_enabled else []



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
        set -euo pipefail
        ### compleasm route requires genomic input
        if [[ "{params.busco_mode}" = "proteins" ]]; then
            echo "Error: Proteome input provided but compleasm mode selected. Please provide genomic input for compleasm."
            exit 2
        fi
        echo "Running compleasm for {wildcards.mag}"

        export HOME="{config[outdir]}"
        export XDG_CACHE_HOME="{config[outdir]}.cache"
        mkdir -p "$XDG_CACHE_HOME"
        rm -rf "$HOME/.compleasm" || true
        export HMMSEARCH="/opt/conda/envs/compleasm/bin/hmmsearch"

        # breadcrumbs
        echo "which hmmsearch: $(which hmmsearch)" >> {log}
        echo "HMMSEARCH=$HMMSEARCH"               >> {log}
        
        compleasm run -a {input} -t {threads} \
            -l eukaryota \
            -L {params.resources_dir}/mb_downloads/ \
            -o {config[outdir]}busco_out/{wildcards.mag} \
            1> {log} 2> {log}
        """
        )

rule proc_database:
    input:
        get_protein_source
    output:
        config["outdir"] + "{mag}_purged_taxa_check.complete",
        directory(config["outdir"] + "{mag}_PhyloFishScratch")
    conda:
        "fisher"
    params:
        resources_dir=RESOURCES_DIR,
        target_taxa=purge_target,
        purge_enabled=purge
    threads: 1
    log:
        config["outdir"] + "logs/proc_database/{mag}.log"
    run:
        from pathlib import Path
        from snakemake.shell import shell

        outdir = Path(config["outdir"])
        log_dir = outdir / "logs" / "proc_database"
        log_dir.mkdir(parents=True, exist_ok=True)

        if params.purge_enabled:
            purge_list = outdir / f"{wildcards.mag}_to_purge_list.txt"
            target_taxa = (params.target_taxa or "").strip()

            if not target_taxa:
                raise ValueError(
                    f"[{datetime.datetime.now():%Y-%m-%d %H:%M:%S}]: Purge requested but no target taxa were resolved."
                )

            purge_list.write_text(f"{target_taxa}\n", encoding="utf-8")

            shell(
                r"""
                echo "{purge_list}"
                echo "purge_list: {purge_list}" >> {log}
                echo "purge UID: {params.target_taxa}" >> {log}

                if [ ! -s "{purge_list}" ]; then
                    echo "ERROR: Purge list is empty: {purge_list}" >> {log}
                    exit 2
                fi

                cp -r {params.resources_dir}/PhyloFisherDatabase_v1.0/database {output[1]}

                purge.py \
                    --input "{purge_list}" \
                    --database {output[1]} \
                    1>> {log} 2>&1

                rm -f "{purge_list}"
                touch {output[0]}
                """,
                purge_list=str(purge_list),
            )
        else:
            shell(
                r"""
                cp -r {params.resources_dir}/PhyloFisherDatabase_v1.0/database {config[outdir]}{wildcards.mag}_PhyloFishScratch
                touch {output[0]}
                """
            )

rule metaeuk:
    input:
        find_mag_file
    output:
        proteome = lambda wildcards: metaeuk_proteome_path(wildcards.mag)
    conda:
        "metaeuk"
    threads: workflow.cores
    params:
        metaeuk_db = config.get(
            "metaeuk_db",
            os.path.join(RESOURCES_DIR, "metaeuk_db")
        ),
        out_prefix = lambda wildcards: metaeuk_prefix(wildcards.mag)
    log:
        config["outdir"] + "logs/metaeuk/metaeuk.log"
    shell:
        r"""
        set -euo pipefail

        mkdir -p {config[outdir]}metaeuk/
        mkdir -p {config[outdir]}metaeuk_tmp/{wildcards.mag}/

        # MetaEuk: contigs -> DB -> out_prefix -> tmp_dir
        metaeuk easy-predict \
            {input} \
            {params.metaeuk_db} \
            {params.out_prefix} \
            {config[outdir]}metaeuk_tmp/{wildcards.mag}/ \
            --threads {threads} \
            > {log} 2>&1
        
        if [ -f "{params.out_prefix}.fasta" ]; then
            mv "{params.out_prefix}.fasta" {output.proteome}
        elif [ -f "{params.out_prefix}_predicted_proteins.fasta" ]; then
            mv "{params.out_prefix}_predicted_proteins.fasta" {output.proteome}
        else
            echo "ERROR: MetaEuk output file not found: {params.out_prefix}.fasta" >> {log}
            exit 2
        fi
        """

rule plm_embed_proteome:
    input:
        fasta = get_plm_proteome_input
    output:
        emb   = config["outdir"] + "plm/{mag}_embeds.npy",
        index = config["outdir"] + "plm/{mag}_embed_index.tsv"
    conda:
        "esm_plm"  # your env with esm + torch + biopython + numpy
    threads: int(config.get("plm_threads", 4))
    params:
        ADD_SCRIPTS    = ADDITIONAL_SCRIPTS_DIR,
        model          = config.get("plm_model", "esm2_t33_650M_UR50D"),
        batch_size     = int(config.get("plm_batch_size", 32)),
        max_residues   = int(config.get("plm_max_residues", 1022)),
        chunk_overlap  = int(config.get("plm_chunk_overlap", 256)),
    log:
        config["outdir"] + "logs/plm_embed/{mag}.log"
    shell:
        r"""
        mkdir -p {config[outdir]}plm/ {config[outdir]}logs/plm_embed/

        python {params.ADD_SCRIPTS}/embed_epdb_with_esm.py \
            --fasta {input.fasta} \
            --model {params.model} \
            --batch_size {params.batch_size} \
            --max_residues {params.max_residues} \
            --chunk_overlap {params.chunk_overlap} \
            --out_embeds {output.emb} \
            --out_index {output.index} \
            > {log} 2>&1
        """

rule plmsearch_epdb:
    input:
        emb   = config["outdir"] + "plm/{mag}_embeds.npy",
        index = config["outdir"] + "plm/{mag}_embed_index.tsv"
    output:
        hits       = config["outdir"] + "plm/{mag}_plm_hits.tsv",
        candidates = config["outdir"] + "plm/{mag}_plm_candidates.faa"
    conda:
        "esm_plm"  # same env; needs numpy, pandas, biopython
    threads: 1
    params:
        ADD_SCRIPTS       = ADDITIONAL_SCRIPTS_DIR,
        epdb_embeds       = config.get("plm_epdb_embeds", os.path.join(RESOURCES_DIR, "plm", "epdb_embeds.npy")),
        epdb_meta         = config.get("plm_epdb_meta",   os.path.join(RESOURCES_DIR, "plm", "epdb_plm_meta.tsv")),
        query_fasta       = get_plm_proteome_input,
        top_families      = int(config.get("plm_top_families", 16)),
        family_sim_thresh = float(config.get("plm_family_sim_thresh", 0.15)),
        top_hits          = int(config.get("plm_top_hits", 8)),
        hit_sim_thresh    = float(config.get("plm_hit_sim_thresh", 0.25)),
    log:
        config["outdir"] + "logs/plmsearch/{mag}.log"
    shell:
        r"""
        mkdir -p {config[outdir]}plm/ {config[outdir]}logs/plmsearch/

        python {params.ADD_SCRIPTS}/plmsearch_epdb.py \
          --epdb-embeds  {params.epdb_embeds} \
          --epdb-meta    {params.epdb_meta} \
          --query-embeds {input.emb} \
          --query-index  {input.index} \
          --query-fasta  {params.query_fasta} \
          --out-hits     {output.hits} \
          --out-candidates {output.candidates} \
          --top-families      {params.top_families} \
          --family-sim-thresh {params.family_sim_thresh} \
          --top-hits          {params.top_hits} \
          --hit-sim-thresh    {params.hit_sim_thresh} \
          > {log} 2>&1
        """


rule fishing_meta:
    input:
        lambda wildcards: (
            os.path.join(config["outdir"], "plm", f"{wildcards.mag}_plm_candidates.faa")
            if plmsearch_enabled
            else get_protein_source(wildcards)
        )
    output:
        config["outdir"] + "{mag}_input_metadata.tsv"
    conda:
        "pline_max"
    params:
        ADD_SCRIPTS=ADDITIONAL_SCRIPTS_DIR
    threads: 1
    log:
        config["outdir"] + "logs/fishing_meta/{mag}.log"
    priority: 0
    shell:
        "python {params.ADD_SCRIPTS}/fishing_meta.py -p {input} -o {output} > {log} 2> {log}"


checkpoint goneFishing:
    input:
        meta = config["outdir"] + "{mag}_input_metadata.tsv",
        dbdir = config["outdir"] + "{mag}_PhyloFishScratch/"
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
        "bash {params.ADD_SCRIPTS}/fishing.sh -t {threads} -i {input.meta} -r {params.resources_dir} -o {config[outdir]}"


rule splitter:
    input:
        dir=config["outdir"] + "{mag}_working_dataset/",
        tar=config["outdir"] + "{mag}_working_dataset/{gene}.fas",
        mag_dir=config["mag_dir"]
    output:
       qs=config["outdir"] + "{mag}_q_frags/{gene}.fas",
       ref=config["outdir"] + "{mag}_ref_frags/{gene}.fas"
    conda:
        "pline_max"
    params:
        ADD_SCRIPTS=ADDITIONAL_SCRIPTS_DIR,
    threads: 1
    priority: 0
    log:
        config["outdir"] + "logs/splitter/{mag}_{gene}.log"
    shell:
        "python {params.ADD_SCRIPTS}splitter.py -i {input.tar} -d {input.mag_dir} -o {output.qs} -r {output.ref} 1> {log} 2> {log}"


rule mafft:
    input:
        query=config["outdir"] + "{mag}_q_frags/{gene}.fas",
        reference=config["outdir"] + "{mag}_ref_frags/{gene}.fas"
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
        ADD_SCRIPTS=ADDITIONAL_SCRIPTS_DIR,
        out_dir=config["outdir"],
        ALN_SUFFIX=branch(trim_alignments, ".trimal", ".aln")
    threads: 1
    priority: 0       
    log:
        config["outdir"] + "logs/concat/{mag}.log"
    shell:
        """
        FIXED_ALNS=()
        mkdir -p {params.out_dir}{wildcards.mag}_relabeled
        for i in $(realpath {params.out_dir}{wildcards.mag}_mafft_out/*{params.ALN_SUFFIX});
        do
        prot=$(basename ${{i}} {params.ALN_SUFFIX})
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
        aln=config["outdir"] + "{mag}_ref.aln"
    output:
        config["outdir"] + "{mag}_ref.tre"
    conda:
        "dendropy"
    params:
        ADD_SCRIPTS=ADDITIONAL_SCRIPTS_DIR,
        resources_dir=RESOURCES_DIR,
        target_taxa=purge_target
    threads: 1
    priority: 0
    log:
        config["outdir"] + "logs/sub_tree/{mag}.log"
    shell:
        branch(purge,
        """
        python {params.ADD_SCRIPTS}sub_tree.py -a {input.aln} -t {params.resources_dir}/ref_concat_PF_alt3.tre -p {params.target_taxa} -o {output} > {log} 2> {log}
        """,
        """
        python {params.ADD_SCRIPTS}sub_tree.py -a {input.aln} -t {params.resources_dir}/ref_concat_PF_alt3.tre -o {output} > {log} 2> {log}
        """)

rule epa:
    input:
        q_aln= config["outdir"] + "{mag}_q.aln",
        ref_aln= config["outdir"] + "{mag}_ref.aln",
        ref_tree= config["outdir"] + "{mag}_ref.tre"
    output:
        config["outdir"] + "{mag}_epa_out/{mag}_epa_out.jplace"
    conda:
        "pline_max"
    threads: workflow.cores
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
        #ulimit -n 65536
        #ulimit -s unlimited
        epa-ng --ref-msa {input.ref_aln} \
         --tree {input.ref_tree} \
         --query {input.q_aln} \
         --model LG -T {threads} >{log} 2>&1
        mv epa_result.jplace {output}
        if [ -f epa_info.log ]; then cat epa_info.log >> {log}; rm epa_info.log; fi
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
        out_dir=config["outdir"],
        resources_dir=RESOURCES_DIR,
        tax_tree=os.path.join(RESOURCES_DIR, "tax_tree.txt")
    priority: 0
    log:
        config["outdir"] + "logs/gappa/{mag}.log"
    shell:
        """
        gappa examine assign \
            --jplace-path {input} \
            --taxon-file {params.resources_dir}/tax_tree.txt \
            --out-dir {params.out_dir}{wildcards.mag}_epa_out \
            --allow-file-overwriting --best-hit --verbose > {log}
        if [ ! -s "{input}" ]; then
            echo "ERROR: Missing or empty JPLACE: {input}" >&2
            exit 2
        fi
        
        """

rule jplace_pair_wise_dist_matrix:
    input:
        config["outdir"] + "{mag}_epa_out/{mag}_epa_out.jplace"
    output:
        config["outdir"] + "{mag}_epa_out/pairwise_qSeqDistance2leaves.tsv"
    conda:
        "pline_max"
    threads: 1
    params:
        out_dir=config["outdir"],
        ADD_SCRIPTS=ADDITIONAL_SCRIPTS_DIR
    priority: 0
    log:
        config["outdir"] + "logs/gappa/{mag}_pairwise_distance.log"
    shell:
        """
        python {params.ADD_SCRIPTS}jplace_dist2leaves_csv.py {input} -o {output}
        """

rule species_tree:
    input:
        config["outdir"] + "{mag}_SuperMatrix.fas"
    output:
        config["outdir"] + "species_tree/{mag}_species_tree.treefile"
    conda:
        "raxml-ng"
    threads: workflow.cores
    priority: 0
    log:
        config["outdir"] + "logs/species_tree/{mag}.log"
    shell:
        """
        mkdir -p {config[outdir]}species_tree/
        mkdir -p {config[outdir]}logs/species_tree/
        iqtree -s {input} \
            -m Q.pfam+I+G4 \
            -bb 1000 \
            -nt {threads} \
            -pre {config[outdir]}species_tree/{wildcards.mag}_species_tree \
            1> {log} 2> {log}
        """
