# HERACLES — reproducible container.
#
# Build from the repository root (the folder that contains main.py):
#   docker build -f docker/Dockerfile -t heracles:1.0 .
#
# Reference data are NOT baked into the image; mount them at run time:
#   -v /path/to/tRNA_index:/data/tRNA_index:ro
#   -v /path/to/RIsearch2_index:/data/RIsearch2_index:ro
#   -v /path/to/fastq:/data/fastq:ro
#   -v /path/to/results:/data/out
#
# The base image ships conda rather than micromamba on purpose. environment.yml
# ends with a `pip:` section, and micromamba 1.5.x executes that section before
# creating the environment, failing with "Cannot activate, prefix does not
# exist". Conda orders the two correctly, and conda is what the environment was
# verified against.

FROM continuumio/miniconda3:latest

LABEL org.opencontainers.image.title="HERACLES" \
      org.opencontainers.image.description="End-to-end tsRNA analysis pipeline" \
      org.opencontainers.image.source="https://github.com/Heranova-Lifesciences/HERACLES" \
      org.opencontainers.image.licenses="MIT"

# Single source of truth for the environment: the same file a user without
# Docker passes to `conda env create`.
COPY docker/environment.yml /tmp/environment.yml
RUN conda env create -f /tmp/environment.yml \
 && conda clean -a -y \
 && rm -f /tmp/environment.yml

ENV PATH=/opt/conda/envs/heracles/bin:$PATH

WORKDIR /opt/heracles
COPY . /opt/heracles/

# Reference data are mounted here; nothing in the image depends on their content.
ENV HERACLES_DATA=/data
RUN mkdir -p /data/out

# Fail the build early if any component is missing.
# `set -e` is required: chaining with `;` would let a failing check be masked by
# a later passing one.
RUN set -e; \
    bowtie --version > /dev/null; \
    cutadapt --version > /dev/null; \
    trim_galore --version > /dev/null; \
    RIsearch2 2>&1 | grep -q "RIsearch v2"; \
    python -c "import pandas, numpy, Bio, networkx, matplotlib, seaborn, pydeseq2, gseapy"

ENTRYPOINT ["python", "/opt/heracles/main.py"]
CMD ["--help"]
