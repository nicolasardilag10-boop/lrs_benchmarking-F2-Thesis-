# ---------------------------------------------------------------------------
# Stage 1: obtain GNU gzip from Debian 12
# ---------------------------------------------------------------------------

FROM debian:12-slim AS gzip_source

RUN apt-get update \
    && apt-get install -y --no-install-recommends gzip \
    && rm -rf /var/lib/apt/lists/*


# ---------------------------------------------------------------------------
# Stage 2: original GoldRush BioContainer
# ---------------------------------------------------------------------------

FROM quay.io/biocontainers/goldrush:1.2.2--py312h248b54c_1

LABEL debug.purpose="Test GNU gzip compatibility with GoldRush/ntLink"

# /usr/local/bin occurs before /usr/bin in PATH.
#
# Therefore this GNU gzip overrides the BusyBox gzip applet for normal
# command resolution without modifying the original BusyBox installation.
COPY --from=gzip_source /usr/bin/gzip /usr/local/bin/gzip


# ---------------------------------------------------------------------------
# Build-time validation
# ---------------------------------------------------------------------------

RUN command -v goldrush >/dev/null \
    && command -v gzip >/dev/null \
    && gzip --version \
    && printf '@test_read\nACGT\n+\nIIII\n' > /tmp/plain.fastq \
    && gzip -f -cd /tmp/plain.fastq > /tmp/passthrough.fastq \
    && cmp /tmp/plain.fastq /tmp/passthrough.fastq \
    && rm -f /tmp/plain.fastq /tmp/passthrough.fastq


WORKDIR /work

ENTRYPOINT []

CMD ["/bin/bash"]
