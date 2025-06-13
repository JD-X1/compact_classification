# syntax=docker/dockerfile:1

FROM continuumio/miniconda3

# Make RUN commands use 'bash --login'

SHELL ["/bin/bash", "--login", "-c"]



# Install Dependencies
COPY envs/ envs/
RUN wget https://ndownloader.figshare.com/files/29093409 \
    && mv 29903409 resources/ \
    && cd resources \
    && tar -xzvf 29093409 \
    && rm 29093409 \
    && cd ..

RUN conda env create -f envs/snakemake.yaml
RUN conda env create -f envs/pline_max.yaml
RUN conda env create -f envs/fisher.yaml
RUN conda run -n pline_max /bin/bash -c compleasm download eukaryota \
    && mv mb_downloads resources/
# Make RUN commands use the new environment
SHELL [ "conda", "run", "-n", "snakemake", "/bin/bash", "-c" ]


# Install App (COPY code)

COPY rules/ rules/
COPY additional_scripts/ additional_scripts/
COPY resources/ resources/



# Final Configuration
EXPOSE 8000
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "snakemake", "snakemake", "-s", "rules/taxa_class_PF_concat.smk", "--use-conda", "-p", "--keep-going", "--rerun-incomplete"]
CMD ["--cores", "1", "--config", "mag_dir=resources/test", "mode=CONCAT"]
