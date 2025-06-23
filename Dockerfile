# syntax=docker/dockerfile:1

FROM continuumio/miniconda3

# Make RUN commands use 'bash --login'

SHELL ["/bin/bash", "--login", "-c"]

# Install Dependencies

### Unincluded system libraries 
RUN apt-get update && apt-get install -y \
    libtiff6 \
    libgl1 \
    libxext6 \
    libsm6 \
    libxrender1 \
    && rm -rf /var/lib/apt/lists/*

### Install mamba to speed up rebuilds
RUN conda install -n base -c conda-forge mamba

### Bring in env files in 
COPY envs/snakemake.yaml envs/
COPY envs/pline_max.yaml envs/
COPY envs/fisher.yaml envs/
COPY envs/mb.yaml envs/
COPY envs/compleasm.yaml envs/
COPY envs/snakemake.yaml envs/

### Run mamba installs
RUN mamba env create -f envs/snakemake.yaml
RUN mamba env create -f envs/pline_max.yaml
RUN mamba env create -f envs/fisher.yaml
RUN mamba env create -f envs/mb.yaml
RUN mamba env create -f envs/compleasm.yaml


# Copy Remaining files
COPY rules/ rules/
COPY additional_scripts/ additional_scripts/
COPY resources/ resources/


# change odb arguement if using other different versions of ODB
# RUN conda activate compleasm \
#     && compleasm download eukaryota \
#     && mv mb_downloads resources/


# RUN wget https://ndownloader.figshare.com/files/29093409 \
#     && tar -xzvf 29093409 \
#     && cp -r PhyloFisherDatabase_v1.0 resources/. \
#     && rm 29093409 \
#     && cd ..



# Make RUN commands use the new environment
SHELL [ "conda", "run", "-n", "snakemake", "/bin/bash", "-c" ]

# Final Configuration
EXPOSE 8000
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "snakemake", "snakemake", "--use-conda", "-p", "--keep-going", "--rerun-incomplete"]
CMD ["-s", "rules/taxa_class_PF.smk", "--cores", "1", "--config", "mag_dir=resources/test", "mode=CONCAT"]