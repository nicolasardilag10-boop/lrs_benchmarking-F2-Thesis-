ARG BASE_IMAGE=quay.io/biocontainers/goldrush:1.2.2--py312h248b54c_1
FROM ${BASE_IMAGE}


LABEL org.opencontainers.image.title="GoldRush assembler"
LABEL org.opencontainers.image.description="GoldRush 1.2.2 long-read assembly pipeline"
LABEL org.opencontainers.image.version="1.2.2"
LABEL org.opencontainers.image.vendor="Nicolas Ardila F2 benchmarking workflow"
LABEL org.opencontainers.image.base.name="quay.io/biocontainers/goldrush:1.2.2--py312h248b54c_1"

RUN command -v goldrush >/dev/null

WORKDIR /work

ENTRYPOINT []
CMD ["/bin/bash"]
