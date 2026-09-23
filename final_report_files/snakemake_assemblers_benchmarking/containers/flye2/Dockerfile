ARG BASE_IMAGE=quay.io/biocontainers/flye:2.9.6--py313h7fbb527_1
FROM ${BASE_IMAGE}

LABEL org.opencontainers.image.title="Flye 2 assembler"
LABEL org.opencontainers.image.description="Flye 2.9.6 for ONT and PacBio HiFi de novo assembly"
LABEL org.opencontainers.image.version="2.9.6"
LABEL org.opencontainers.image.vendor="Nicolas Ardila F2 benchmarking workflow"
LABEL org.opencontainers.image.base.name="quay.io/biocontainers/flye:2.9.6--py313h7fbb527_1"

RUN command -v flye >/dev/null \
    && flye --version

WORKDIR /work

ENTRYPOINT []
CMD ["/bin/bash"]
