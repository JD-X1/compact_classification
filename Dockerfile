# syntax=docker/dockerfile:1

FROM continuumio/miniconda3

# Make RUN commands use 'bash --login'

SHELL ["/bin/bash", "--login", "-c"]



# Install Dependencies
COPY envs/ envs/
COPY rules/ rules/
COPY additional_scripts/ additional_scripts/
COPY resources/ resources/

RUN conda env create -f envs/snakemake.yaml
RUN conda env create -f envs/pline_max.yaml
RUN conda env create -f envs/fisher.yaml
RUN conda env create -f envs/mb.yaml
RUN conda env create -f envs/compleasm.yaml
# change odb arguement if using other different versions of ODB
RUN conda activate compleasm \
    && compleasm download eukaryota \
    && mv mb_downloads resources/


RUN wget https://ndownloader.figshare.com/files/29093409 \
    && tar -xzvf 29093409 \
    && cp -r PhyloFisherDatabase_v1.0 resources/. \
    && rm 29093409 \
    && cd ..



# Make RUN commands use the new environment
SHELL [ "conda", "run", "-n", "snakemake", "/bin/bash", "-c" ]

# Final Configuration
EXPOSE 8000
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "snakemake", "snakemake", "--use-conda", "-p", "--keep-going", "--rerun-incomplete"]
CMD ["-s", "rules/taxa_class_PF.smk", "--cores", "1", "--config", "mag_dir=resources/test", "mode=CONCAT"]