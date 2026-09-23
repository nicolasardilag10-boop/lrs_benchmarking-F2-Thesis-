ARG BASE_IMAGE=quay.io/biocontainers/ntlink:1.3.11--py313h380480a_2
FROM ${BASE_IMAGE}

LABEL org.opencontainers.image.title="ntLink scaffolder"
LABEL org.opencontainers.image.description="ntLink 1.3.11 long-read scaffolding"
LABEL org.opencontainers.image.version="1.3.11"
LABEL org.opencontainers.image.vendor="Nicolas Ardila F2 benchmarking workflow"
LABEL org.opencontainers.image.base.name="quay.io/biocontainers/ntlink:1.3.11--py313h380480a_2"

RUN command -v ntLink >/dev/null

WORKDIR /work

ENTRYPOINT []
CMD ["/bin/bash"]
