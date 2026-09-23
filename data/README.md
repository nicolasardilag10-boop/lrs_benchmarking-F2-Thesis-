# Data entry points

This directory provides stable, readable entry points to project data:

- `fastq` → [`../fastq/`](../fastq/)
- `reference` → [`../reference/`](../reference/)

Both are symbolic links. The original directories remain protected because active workflows depend on their root-relative locations. Large sequencing data and indexes are intentionally excluded from version control.

See [`../docs/DATA_MAP.md`](../docs/DATA_MAP.md) before moving or replacing any data file.
