# ============================================================
# MEAN READ-LENGTH PLOT
#
# This script:
#   1. Reads the validated FASTQ summary table
#   2. Extracts mean read lengths
#   3. Matches ONT and PacBio values by sample
#   4. Creates a paired comparison plot
#   5. Saves the result as a PNG file
# ============================================================


# ------------------------------------------------------------
# 1. Define project, input, and output paths
# ------------------------------------------------------------

project_dir <- path.expand("~/lrs_benchmarking")

input_file <- file.path(
  project_dir,
  "alignment_analysis",
  "tables",
  "fastq_summary.tsv"
)

output_file <- file.path(
  project_dir,
  "alignment_analysis",
  "figures",
  "mean_read_length.png"
)


# ------------------------------------------------------------
# 2. Check that the input table exists
# ------------------------------------------------------------

if (!file.exists(input_file)) {
  stop("Input table not found: ", input_file)
}


# ------------------------------------------------------------
# 3. Read the FASTQ summary table
# ------------------------------------------------------------

fastq <- read.delim(
  input_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 4. Check the required columns
# ------------------------------------------------------------

required_columns <- c(
  "sample",
  "technology",
  "mean_length"
)

missing_columns <- setdiff(
  required_columns,
  names(fastq)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing columns: ",
    paste(missing_columns, collapse = ", ")
  )
}


# ------------------------------------------------------------
# 5. Separate ONT and PacBio values
# ------------------------------------------------------------

ont <- fastq[
  fastq$technology == "ont",
  c("sample", "mean_length")
]

pacbio <- fastq[
  fastq$technology == "pb",
  c("sample", "mean_length")
]


# ------------------------------------------------------------
# 6. Match both technologies by sample
# ------------------------------------------------------------

# merge() ensures that HG002 is compared with HG002,
# HG003 with HG003, and HG004 with HG004.

paired <- merge(
  ont,
  pacbio,
  by = "sample",
  suffixes = c("_ONT", "_PacBio")
)

# Sort samples in a predictable order.
paired <- paired[
  order(paired$sample),
]


# ------------------------------------------------------------
# 7. Validate the paired table
# ------------------------------------------------------------

if (nrow(paired) == 0) {
  stop("No paired ONT and PacBio samples were found")
}

if (any(is.na(paired))) {
  stop("The paired mean-length table contains missing values")
}

cat("Values used in the mean-length plot:\n")
print(paired)


# ------------------------------------------------------------
# 8. Ensure that the output directory exists
# ------------------------------------------------------------

dir.create(
  dirname(output_file),
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 9. Open the PNG graphics device
# ------------------------------------------------------------

png(
  filename = output_file,
  width = 1900,
  height = 1400,
  res = 250
)


# ------------------------------------------------------------
# 10. Adjust plot margins
# ------------------------------------------------------------

par(
  mar = c(5, 6, 4, 7) + 0.1
)


# ------------------------------------------------------------
# 11. Calculate the y-axis range
# ------------------------------------------------------------

all_values <- c(
  paired$mean_length_ONT,
  paired$mean_length_PacBio
)

# Add a small amount of space above and below the data.
y_padding <- diff(range(all_values)) * 0.15

y_limits <- c(
  min(all_values) - y_padding,
  max(all_values) + y_padding
)


# ------------------------------------------------------------
# 12. Create an empty plotting area
# ------------------------------------------------------------

plot(
  NA,
  xlim = c(0.8, 2.2),
  ylim = y_limits,
  xaxt = "n",
  xlab = "Sequencing technology",
  ylab = "Mean read length (bp)",
  main = "Mean read length: ONT versus PacBio"
)

axis(
  side = 1,
  at = c(1, 2),
  labels = c("ONT", "PacBio")
)


# ------------------------------------------------------------
# 13. Draw one paired line per sample
# ------------------------------------------------------------

# Different point shapes distinguish the three samples.

point_shapes <- seq_len(nrow(paired))

for (i in seq_len(nrow(paired))) {

  sample_values <- c(
    paired$mean_length_ONT[i],
    paired$mean_length_PacBio[i]
  )

  lines(
    x = c(1, 2),
    y = sample_values,
    type = "b",
    pch = point_shapes[i],
    lwd = 2,
    cex = 1.2
  )
}


# ------------------------------------------------------------
# 14. Add numeric annotations
# ------------------------------------------------------------

# ONT values are placed slightly to the left.
text(
  x = rep(1, nrow(paired)) - 0.06,
  y = paired$mean_length_ONT,
  labels = format(
    round(paired$mean_length_ONT),
    big.mark = ",",
    scientific = FALSE
  ),
  pos = 2,
  cex = 0.75
)

# PacBio values are placed slightly to the right.
text(
  x = rep(2, nrow(paired)) + 0.06,
  y = paired$mean_length_PacBio,
  labels = format(
    round(paired$mean_length_PacBio),
    big.mark = ",",
    scientific = FALSE
  ),
  pos = 4,
  cex = 0.75
)


# ------------------------------------------------------------
# 15. Add the sample legend
# ------------------------------------------------------------

legend(
  "topright",
  inset = c(-0.22, 0),
  legend = paired$sample,
  pch = point_shapes,
  lty = 1,
  lwd = 2,
  bty = "n",
  xpd = NA
)


# ------------------------------------------------------------
# 16. Close and save the figure
# ------------------------------------------------------------

dev.off()


# ------------------------------------------------------------
# 17. Print confirmation
# ------------------------------------------------------------

cat("\nMean read-length figure created successfully:\n")
cat(output_file, "\n")
