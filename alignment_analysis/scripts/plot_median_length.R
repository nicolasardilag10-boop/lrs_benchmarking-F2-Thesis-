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
  "median_read_length.png"
)

# Check that the table exists
if (!file.exists(input_file)) {
  stop("Input table not found: ", input_file)
}

# Read the summary table
fastq <- read.delim(
  input_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# Select the necessary columns
ont <- fastq[
  fastq$technology == "ont",
  c("sample", "median_length")
]

pb <- fastq[
  fastq$technology == "pb",
  c("sample", "median_length")
]

# Match ONT and PacBio by sample
paired <- merge(
  ont,
  pb,
  by = "sample",
  suffixes = c("_ONT", "_PacBio")
)

# Create the PNG figure
png(
  filename = output_file,
  width = 1800,
  height = 1400,
  res = 250
)

y_values <- c(
  paired$median_length_ONT,
  paired$median_length_PacBio
)

plot(
  NA,
  xlim = c(0.8, 2.2),
  ylim = range(y_values),
  xaxt = "n",
  xlab = "Sequencing technology",
  ylab = "Median read length (bp)",
  main = "Median read length: ONT versus PacBio"
)

axis(
  side = 1,
  at = c(1, 2),
  labels = c("ONT", "PacBio")
)

for (i in seq_len(nrow(paired))) {

  lines(
    x = c(1, 2),
    y = c(
      paired$median_length_ONT[i],
      paired$median_length_PacBio[i]
    ),
    type = "b",
    pch = i,
    lwd = 2
  )
}

legend(
  "topright",
  legend = paired$sample,
  pch = seq_len(nrow(paired)),
  lty = 1,
  lwd = 2,
  bty = "n"
)

dev.off()

cat("Figure successfully created:\n")
cat(output_file, "\n")
