# Filename: plot_h_ratio.R
library(karyoploteR)
library(GenomicRanges)

# --- 1. Read and Prepare Data ---
data_df <- read.csv(snakemake@input[[1]], header = TRUE)

# Ensure chr prefix for karyoploteR
if (!startsWith(as.character(data_df$chr[1]), "chr")) {
  data_df$chr <- paste0("Chr ", data_df$chr)
}

# Define custom genome based on data
genome_sizes <- aggregate(pos ~ chr, data_df, max)
custom_genome <- toGRanges(data.frame(
  chr = genome_sizes$chr,
  start = 1,
  end = genome_sizes$pos
))

# --- 2. Set Up PDF Device ---
pdf(file = snakemake@output[[1]], width = 12, height = 2)

# --- 3. Create Karyotype Plot ---

# Format title based on sample
sample_label <- gsub("_", " x ", snakemake@wildcards[["sample"]])

# Get the default plotting parameters for plot.type 4
pp <- getDefaultPlotParams(plot.type = 4)

# Increase the gap between chromosomes
pp$ideogramgap <- 0.4
# Increase left margin gap
pp$leftmargin <- 0.06

# Using plot.type=4 places all chromosomes in one line
kp <- plotKaryotype(
  genome = custom_genome,
  plot.type = 4,
  plot.params = pp,
  ideogram.plotter = NULL,
  main = sample_label,
  chromosomes = "all"
)

# --- 4. Add Visual Styling (To match the target image) ---

# A. Add a light background grid and shading
kpDataBackground(kp, color = "#FFFFFF")

# B. Add a horizontal threshold line (e.g., at 0.5)
# kpAbline is used for horizontal lines in karyoploteR
kpAbline(kp, h = 0.5, col = "blue", lwd = 1, lty = 2)
# Draw red line at 0.3 and 0.7
kpAbline(kp, h = 0.3, col = "red", lwd = 1, lty = 2)
kpAbline(kp, h = 0.7, col = "red", lwd = 1, lty = 2)
# Draw orange line at 0.4 and 0.6
kpAbline(kp, h = 0.4, col = "orange", lwd = 1, lty = 2)
kpAbline(kp, h = 0.6, col = "orange", lwd = 1, lty = 2)

# C. Add a Y-axis with specific labels
kpAxis(kp, ymin = 0, ymax = 1, numticks = 6, cex = 1)
kpAddLabels(kp, labels = "Frequency", srt = 90, pos = 1, label.margin = 0.05, cex = 1)

# D. Plot the data using a combination of area and lines for clarity
# Using kpArea creates the filled effect often seen in these plots
kpLines(
  karyoplot = kp,
  chr = data_df$chr,
  x = data_df$pos,
  y = data_df$h_ratio,
  col = "#2980B9",
  lwd = 2
)

# --- 5. Finalize ---
dev.off()
