# ============================================================
# FASTQ QUALITY PLOT: Q20 AND Q30
#
# This script:
#   1. Reads the FASTQ summary table
#   2. Validates the Q20 and Q30 values
#   3. Creates a grouped bar plot
#   4. Adds percentage labels without overlapping
#   5. Saves the figure as a PNG file
# ============================================================


# ------------------------------------------------------------
# 1. Define the project directory
# ------------------------------------------------------------

project_dir <- path.expand("~/lrs_benchmarking")


# ------------------------------------------------------------
# 2. Define input and output paths
# ------------------------------------------------------------

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
  "q20_q30_comparison.png"
)


# ------------------------------------------------------------
# 3. Check that the input table exists
# ------------------------------------------------------------

if (!file.exists(input_file)) {
  stop("Input table not found: ", input_file)
}


# ------------------------------------------------------------
# 4. Read the FASTQ summary table
# ------------------------------------------------------------

fastq <- read.delim(
  input_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 5. Check the required columns
# ------------------------------------------------------------

required_columns <- c(
  "sample",
  "technology",
  "Q20_percent",
  "Q30_percent"
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
# 6. Validate the quality values
# ------------------------------------------------------------

if (any(is.na(fastq[, required_columns]))) {
  stop("Required columns contain missing values")
}

if (any(
  fastq$Q20_percent < 0 |
  fastq$Q20_percent > 100
)) {
  stop("Some Q20 percentages are outside 0–100")
}

if (any(
  fastq$Q30_percent < 0 |
  fastq$Q30_percent > 100
)) {
  stop("Some Q30 percentages are outside 0–100")
}

# Q30 must not be higher than Q20
if (any(fastq$Q30_percent > fastq$Q20_percent)) {
  stop("A Q30 percentage is greater than its Q20 percentage")
}


# ------------------------------------------------------------
# 7. Create readable technology labels
# ------------------------------------------------------------

fastq$technology_label <- ifelse(
  fastq$technology == "ont",
  "ONT",
  "PacBio"
)


# ------------------------------------------------------------
# 8. Create sample labels
# ------------------------------------------------------------

sample_labels <- paste(
  fastq$sample,
  fastq$technology_label,
  sep = "-"
)


# ------------------------------------------------------------
# 9. Build the matrix used by barplot()
# ------------------------------------------------------------

quality_matrix <- rbind(
  fastq$Q20_percent,
  fastq$Q30_percent
)

colnames(quality_matrix) <- sample_labels
rownames(quality_matrix) <- c("Q20", "Q30")


# ------------------------------------------------------------
# 10. Show the values in the terminal
# ------------------------------------------------------------

cat("Values used in the plot:\n")
print(quality_matrix)


# ------------------------------------------------------------
# 11. Make sure the output directory exists
# ------------------------------------------------------------

dir.create(
  dirname(output_file),
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 12. Open the PNG graphics device
# ------------------------------------------------------------

png(
  filename = output_file,
  width = 2400,
  height = 1500,
  res = 250
)


# ------------------------------------------------------------
# 13. Adjust plot margins
# ------------------------------------------------------------

# bottom, left, top, right
# We increase the right margin to leave room for the legend.
par(
  mar = c(10, 5, 4, 8) + 0.1
)


# ------------------------------------------------------------
# 14. Define bar colors
# ------------------------------------------------------------

bar_colors <- c("gray35", "gray80")


# ------------------------------------------------------------
# 15. Create the grouped bar plot
# ------------------------------------------------------------

bar_positions <- barplot(
  quality_matrix,
  beside = TRUE,
  ylim = c(0, 112),
  ylab = "Reads meeting quality threshold (%)",
  main = "FASTQ read quality: Q20 and Q30",
  las = 2,
  col = bar_colors,
  border = "black"
)


# ------------------------------------------------------------
# 16. Add horizontal reference lines
# ------------------------------------------------------------

abline(
  h = c(20, 40, 60, 80, 100),
  lty = 3
)


# ------------------------------------------------------------
# 17. Format percentage labels
# ------------------------------------------------------------

format_percentage <- function(values) {
  ifelse(
    values == round(values),
    paste0(round(values), "%"),
    paste0(
      format(
        values,
        nsmall = 1,
        trim = TRUE
      ),
      "%"
    )
  )
}


# ------------------------------------------------------------
# 18. Add Q20 labels
# ------------------------------------------------------------

text(
  x = bar_positions[1, ] - 0.06,
  y = fastq$Q20_percent + 1.2,
  labels = format_percentage(fastq$Q20_percent),
  cex = 0.62,
  xpd = TRUE
)


# ------------------------------------------------------------
# 19. Add Q30 labels
# ------------------------------------------------------------

text(
  x = bar_positions[2, ] + 0.06,
  y = fastq$Q30_percent + 3.0,
  labels = format_percentage(fastq$Q30_percent),
  cex = 0.62,
  xpd = TRUE
)


# ------------------------------------------------------------
# 20. Add legend OUTSIDE the plotting area
# ------------------------------------------------------------

legend(
  "topright",
  inset = c(-0.18, 0),
  legend = c("Q20", "Q30"),
  fill = bar_colors,
  border = "black",
  bty = "n",
  xpd = NA
)