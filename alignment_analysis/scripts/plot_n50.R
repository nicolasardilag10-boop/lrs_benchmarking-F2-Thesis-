# ------------------------------------------------------------
# N50 comparison plot: ONT vs PacBio
# This script reads the fastq_summary.tsv table, extracts N50
# values for both technologies, pairs them by sample, and
# generates a line plot showing per-sample differences.
# ------------------------------------------------------------

# Expand home directory path
project_dir <- path.expand("~/lrs_benchmarking")

# Input TSV containing read-length statistics
input_file <- file.path(
  project_dir,
  "alignment_analysis",
  "tables",
  "fastq_summary.tsv"
)

# Output PNG file for the N50 comparison figure
output_file <- file.path(
  project_dir,
  "alignment_analysis",
  "figures",
  "n50_comparison.png"
)

# Safety check: stop if the input table does not exist
if (!file.exists(input_file)) {
  stop("Input table not found: ", input_file)
}

# Read TSV file (tab-separated, no factor conversion)
fastq <- read.delim(
  input_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# Split data by sequencing technology
# ------------------------------------------------------------

# Extract ONT rows and keep only sample + N50 columns
ont <- fastq[
  fastq$technology == "ont",
  c("sample", "N50")
]

# Extract PacBio rows and keep only sample + N50 columns
pb <- fastq[
  fastq$technology == "pb",
  c("sample", "N50")
]

# ------------------------------------------------------------
# Pair ONT and PacBio values by sample name
# ------------------------------------------------------------

paired <- merge(
  ont,
  pb,
  by = "sample",
  suffixes = c("_ONT", "_PacBio")
)

# ------------------------------------------------------------
# Create PNG output
# ------------------------------------------------------------

png(
  filename = output_file,
  width = 1800,
  height = 1400,
  res = 250
)

# Combine all N50 values to compute y-axis range
y_values <- c(
  paired$N50_ONT,
  paired$N50_PacBio
)

# Empty plot with custom axes and labels
plot(
  NA,
  xlim = c(0.8, 2.2),
  ylim = range(y_values),
  xaxt = "n",
  xlab = "Sequencing technology",
  ylab = "Read-length N50 (bp)",
  main = "Read-length N50: ONT versus PacBio"
)

# Add x-axis labels
axis(
  side = 1,
  at = c(1, 2),
  labels = c("ONT", "PacBio")
)

# ------------------------------------------------------------
# Draw paired lines for each sample
# Each sample gets a unique point shape (pch)
# ------------------------------------------------------------

for (i in seq_len(nrow(paired))) {

  lines(
    x = c(1, 2),
    y = c(
      paired$N50_ONT[i],
      paired$N50_PacBio[i]
    ),
    type = "b",   # both points and lines
    pch = i,      # unique marker per sample
    lwd = 2       # thicker line
  )
}

# Legend showing sample names and their corresponding markers
legend(
  "topright",
  legend = paired$sample,
  pch = seq_len(nrow(paired)),
  lty = 1,
  lwd = 2,
  bty = "n"
)

# Close PNG device
dev.off()

# Print confirmation message
cat("N50 figure created:\n")
cat(output_file, "\n")
