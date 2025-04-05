# Load libraries
library(ballgown)
library(RColorBrewer)
library(pheatmap)

# Get command line argument (output directory for ballgown results)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript run_ballgown.R <output/ballgown>")
}
ballgown_output_dir <- args[1]

# Define relative paths (assuming script is run from /Users/cb_ash/Desktop/NGS)
stringtie_output <- file.path("output", "stringtie")          # Location of stringtie results
sample_info_file <- file.path("output", "sample_info.txt")    # Metadata file
plot_file <- file.path(ballgown_output_dir, "sample_correlation_heatmap.png")

# Load sample metadata
sample_info <- read.table(sample_info_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)

# Load ballgown object
sample_paths <- file.path(stringtie_output, sample_info$sample)
bg <- ballgown(samples=sample_paths, meas="all")

# Attach metadata (pData)
pData(bg) <- sample_info

# Quick QC - check sample names and metadata
print(sampleNames(bg))
print(pData(bg))

# Extract expression data
tExpr <- texpr(bg, 'all')     # Full transcript-level expression
fpkmExpr <- texpr(bg, 'FPKM') # Transcript-level FPKM
gExpr <- gexpr(bg)            # Gene-level expression

# Run differential expression analysis (transcript level)
results_transcripts <- stattest(bg, feature='transcript', covariate='condition', getFC=TRUE, meas='FPKM')

# Add gene annotations to results
results_transcripts <- data.frame(
  geneNames = ballgown::geneNames(bg),
  geneIDs = ballgown::geneIDs(bg),
  results_transcripts
)

# Sort and filter significant transcripts
results_transcripts <- results_transcripts[order(results_transcripts$pval), ]
sig_transcripts <- subset(results_transcripts, pval < 0.05)

# Report number of significant transcripts
cat("Number of significant transcripts (pval < 0.05):", nrow(sig_transcripts), "\n")

# Save results to file (inside ballgown output folder)
results_file <- file.path(ballgown_output_dir, "transcript_diffexp_results.tsv")
write.table(results_transcripts, file = results_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Results saved to:", results_file, "\n")

# ---- Visualization: Simple Transcript Clustering ----
# Extract transcript-level expression for all samples
transcript_expression <- texpr(bg, meas="FPKM")

# Compute correlation matrix between samples
transcript_cor <- cor(transcript_expression)

# Assign sample names for heatmap
colnames(transcript_cor) <- sampleNames(bg)
rownames(transcript_cor) <- sampleNames(bg)

# Set conditions as rownames for annotation
rownames(sample_info) <- sample_info$sample
annotation_col <- sample_info["condition"]

# Save heatmap as PNG (inside ballgown output folder)
png(plot_file, width = 1200, height = 800, res = 150)
pheatmap(transcript_cor,
         annotation_col = annotation_col,
         col = colorRampPalette(brewer.pal(9, "GnBu"))(100),
         main = "Transcript Expression Correlation Between Samples",
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         legend = TRUE)
dev.off()

cat("Heatmap saved to:", plot_file, "\n")
