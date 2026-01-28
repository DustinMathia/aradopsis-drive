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
pdf(file = snakemake@output[[1]], width = 12, height = 4)

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

kpDataBackground(kp, color = "#FFFFFF")

# Horizontal threshold lines
# Draw blue line at 0.5
kpAbline(kp, h = 0.5, col = "blue", lwd = 1, lty = 2)
# Draw red line at 0.3 and 0.7
#kpAbline(kp, h = 0.3, col = "red", lwd = 1, lty = 2)
#kpAbline(kp, h = 0.7, col = "red", lwd = 1, lty = 2)
# Draw orange line at 0.4 and 0.6
#kpAbline(kp, h = 0.4, col = "orange", lwd = 1, lty = 2)
#kpAbline(kp, h = 0.6, col = "orange", lwd = 1, lty = 2)

# Vertical centromere lines
# Chromosome 1
v_pos <- snakemake@config[["centromeres"]][["1"]]
kpAbline(kp, v=v_pos, chr="Chr 1", col="black")
# Chromosome 2
v_pos <- snakemake@config[["centromeres"]][["2"]]
kpAbline(kp, v=v_pos, chr="Chr 2", col="black")
# Chromosome 3
v_pos <- snakemake@config[["centromeres"]][["3"]]
kpAbline(kp, v=v_pos, chr="Chr 3", col="black")
# Chromosome 4
v_pos <- snakemake@config[["centromeres"]][["4"]]
kpAbline(kp, v=v_pos, chr="Chr 4", col="black")
# Chromosome 5
v_pos <- snakemake@config[["centromeres"]][["5"]]
kpAbline(kp, v=v_pos, chr="Chr 5", col="black")


# Get dynamic top label
get_non_c_side <- function(filename) {
  # If 'C' is on the left: matches 'C_' and captures (.*), replaces with the capture group
  if (startsWith(filename, "C_")) {
    return(gsub("^C_(.*)$", "\\1", filename))
  }
  # If 'C' is on the right: matches captures (.*) and '_C', replaces with the capture group
  if (endsWith(filename, "_C")) {
    return(gsub("^(.*)_C$", "\\1", filename))
  }
  # If no 'C' is present (e.g., A_B), return the original string
  return(filename)
}

top_label <- sapply(snakemake@wildcards[["sample"]], get_non_c_side)

# C. Add a Y-axis with specific labels
kpAxis(kp, ymin = 0, ymax = 1, numticks = 6, cex = 1)
kpAddLabels(kp, labels = "Segregation Ratio", srt = 90, pos = 1, label.margin = 0.05, cex = 1)

# Bottom label is always "C"
# r0 and r1 define the vertical range for the label placement
kpAddLabels(kp, labels = "C", data.panel = 1, r0=-0.1, r1=0, cex=0.8)

# Add the top label
kpAddLabels(kp, labels = top_label, data.panel = 1, r0=1, r1=1.1, cex=0.8)


# --- Highlight Significant Regions (p <= 0.05) ---
# Filter for significant positions, ensuring no NAs
sig_df <- data_df[!is.na(data_df$p_value) & data_df$p_value <= 0.05, ]

if (nrow(sig_df) > 0) {
  kpSegments(
    karyoplot = kp,
    chr = sig_df$chr,
    x0 = sig_df$pos,
    x1 = sig_df$pos,
    y0 = 0,
    y1 = 1,
    col = adjustcolor("purple", alpha.f = 0.5), # Semi-transparent to prevent occlusion
    lwd = 1
  )
}


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
