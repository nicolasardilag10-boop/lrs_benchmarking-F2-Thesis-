BEGIN {
    quality_characters = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~"
}

NR % 4 == 0 {
    quality_sum = 0

    for (i = 1; i <= length($0); i++) {
        character = substr($0, i, 1)
        phred_score = index(quality_characters, character) - 1
        quality_sum += phred_score
    }

    mean_quality = quality_sum / length($0)

    total_reads++

    if (mean_quality >= 20) {
        q20_reads++
    }

    if (mean_quality >= 30) {
        q30_reads++
    }
}

END {
    printf "%d\t%d\t%.2f\t%d\t%.2f\n",
        total_reads,
        q20_reads,
        100 * q20_reads / total_reads,
        q30_reads,
        100 * q30_reads / total_reads
}
