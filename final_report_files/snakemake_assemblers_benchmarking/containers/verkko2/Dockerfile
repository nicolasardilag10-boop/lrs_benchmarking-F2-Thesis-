ARG BASE_IMAGE=quay.io/biocontainers/verkko:2.3.2--hb0edd9e_0

# -----------------------------------------------------------------------------
# Stage 1: obtain bc and only the shared libraries bc needs.
# Do not recreate /lib or /lib64 in the final image.
# -----------------------------------------------------------------------------
FROM debian:12-slim AS bc_builder

RUN apt-get update \
    && apt-get install -y --no-install-recommends bc \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir -p /opt/bc-libs \
    && ldd /usr/bin/bc \
       | awk '/=> \// {print $3}' \
       | while read -r lib; do \
             case "$lib" in \
                 */libc.so.*) ;; \
                 *) cp -L "$lib" /opt/bc-libs/ ;; \
             esac; \
         done

# -----------------------------------------------------------------------------
# Stage 2: pinned Verkko BioContainer.
# -----------------------------------------------------------------------------
FROM ${BASE_IMAGE}

LABEL org.opencontainers.image.title="Verkko 2 hybrid assembler"
LABEL org.opencontainers.image.description="Verkko 2.3.2 hybrid HiFi and ONT genome assembly with bc runtime dependency"
LABEL org.opencontainers.image.version="2.3.2"
LABEL org.opencontainers.image.vendor="Nicolas Ardila F2 benchmarking workflow"
LABEL org.opencontainers.image.base.name="quay.io/biocontainers/verkko:2.3.2--hb0edd9e_0"

# Copy only the bc executable and its non-glibc runtime libraries.
COPY --from=bc_builder /usr/bin/bc /usr/local/libexec/bc.real
COPY --from=bc_builder /opt/bc-libs/ /usr/local/lib/bc-libs/

# Provide bc through a small wrapper so the additional libraries affect bc only.
RUN printf '%s\n' \
        '#!/bin/sh' \
        'LD_LIBRARY_PATH="/usr/local/lib/bc-libs${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"' \
        'export LD_LIBRARY_PATH' \
        'exec /usr/local/libexec/bc.real "$@"' \
        > /usr/local/bin/bc \
    && chmod +x /usr/local/bin/bc \
    && command -v verkko >/dev/null \
    && command -v bc >/dev/null \
    && bc --version >/dev/null \
    && test "$(printf '1+1\n' | bc)" = "2"

WORKDIR /work

ENTRYPOINT []
CMD ["/bin/bash"]
