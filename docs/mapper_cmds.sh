# GraphAligner
docker run --rm -it -v $PWD:/data schimar/lrs-graphaligner:latest --help

# ParaHAT — index first, then align
docker run --rm -it -v $PWD:/data --entrypoint ParaHAT-indexer schimar/lrs-parahat:latest <ref.fa> <index_dir>
docker run --rm -it -v $PWD:/data schimar/lrs-parahat:latest -n 4 ParaHAT-aligner -t 8 <index_dir> <reads.fastq> <ref.fa>

# QuickEd
docker run --rm -it -v $PWD:/data schimar/lrs-quicked:latest --help
docker run --rm -it -v $PWD:/data --entrypoint generate_dataset schimar/lrs-quicked:latest --help

# VACmap
docker run --rm -it -v $PWD:/data schimar/lrs-vacmap:latest --help

# VG
docker run --rm -it -v $PWD:/data schimar/lrs-vg:latest help
docker run --rm -it -v $PWD:/data schimar/lrs-vg:latest map --help

# minimap2
docker run --rm -it -v $PWD:/data --entrypoint minimap2 schimar/lrs-minimap2-ntlink:latest --help

# ntLink
docker run --rm -it -v $PWD:/data --entrypoint ntLink schimar/lrs-minimap2-ntlink:latest --help
