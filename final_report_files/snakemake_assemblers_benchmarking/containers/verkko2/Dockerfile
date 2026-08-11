ARG BASE_IMAGE=quay.io/biocontainers/verkko:2.3.2--hb0edd9e_0
FROM ${BASE_IMAGE}

LABEL org.opencontainers.image.title="Verkko 2 hybrid assembler"
LABEL org.opencontainers.image.description="Verkko 2.3.2 hybrid HiFi and ONT genome assembly"
LABEL org.opencontainers.image.version="2.3.2"
LABEL org.opencontainers.image.vendor="Nicolas Ardila F2 benchmarking workflow"
LABEL org.opencontainers.image.base.name="quay.io/biocontainers/verkko:2.3.2--hb0edd9e_0"

RUN command -v verkko >/dev/null

WORKDIR /work

ENTRYPOINT []
CMD ["/bin/bash"]
