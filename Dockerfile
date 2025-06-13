# syntax=docker/dockerfile:1

FROM continuumio/miniconda3

# Make RUN commands use 'bash --login'

SHELL ["/bin/bash", "--login", "-c"]

# Install Dependencies
COPY envs/ envs/
RUN conda env create -f envs/snakemake.yaml
RUN conda env 

# Make RUN commands use the new environment
SHELL [ "conda", "run", "-n", "snakemake", "/bin/bash", "-c" ]


# Install App (COPY code)

COPY rules/ rules/
COPY additional_scripts/ additional_scripts/
COPY resources/ resources/


# wget and mv resources



# Final Configuration
EXPOSE 8000
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "snakemake", "snakemake", "-s", "rules/taxa_class_PF_concat.smk", "--use-conda", "-p", "--keep-going", "--rerun-incomplete"]
CMD ["--cores", "1", "--config", "mag_dir=resources/test", "mode=CONCAT"]
