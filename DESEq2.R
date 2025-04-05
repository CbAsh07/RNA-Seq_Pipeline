# Set library path from environment variable
.libPaths(Sys.getenv("DESEQ2_LIB"))

# Load required libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

# Get the arguments from the command line
args <- commandArgs(trailingOnly = TRUE)

# Validate input arguments
if (length(args) != 3) {
  stop("Please provide the featureCounts output directory, DESeq2 output directory, and sample info file as arguments.")
}

# Assign input/output paths from command-line arguments
featurecounts_file <- args[1]
deseq2_output_dir <- args[2]
sample_info_file <- args[3]

# Ensure the provided featureCounts directory is correct
cat("FeatureCounts file provided:", featurecounts_file, "\n")

# Create output directory if it doesn't exist
dir.create(deseq2_output_dir, showWarnings = FALSE, recursive = TRUE)

# Automatically detect the correct counts file
#counts_files <- list.files(featurecounts_file, full.names = TRUE)

#if (length(counts_files) == 0) {
#  stop("Error: No .txt file found in the featureCounts directory.")
#}

# Choose the most relevant file (adjust logic if needed)
#counts_file <- counts_files[1]  # Assuming only one file; refine logic if needed
counts_file <- featurecounts_file

# Debugging: Print the detected file path
cat("Using counts file:", counts_file, "\n")

# Read the counts file
counts <- read.table(counts_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

# Drop metadata columns if they exist
metadata_cols <- c("Chr", "Start", "End", "Strand", "Length")
counts <- counts[, !(colnames(counts) %in% metadata_cols), drop = FALSE]

# Debugging: Print processed counts file structure
cat("Filtered counts file column names:", colnames(counts), "\n")

# Properly clean column names to match the sample names in sample_info
colnames(counts) <- basename(colnames(counts))  # Removes full path (if present)
colnames(counts) <- gsub("_sorted\\.bam$", "", colnames(counts))  # Removes file extensions
colnames(counts) <- make.names(colnames(counts))
cat("Processed counts column names:", colnames(counts), "\n")

# Read the sample information file
sample_info <- read.table(sample_info_file, header = TRUE, row.names = 1, sep = "\t")
rownames(sample_info) <- trimws(rownames(sample_info))  # Trim whitespace
rownames(sample_info) <- make.names(rownames(sample_info))
print(head(sample_info))
print(str(sample_info))

# Ensure that sample names match between counts and sample info
if (!all(colnames(counts) %in% rownames(sample_info))) {
  cat("Error: Mismatch between sample names in counts file and sample information file.\n")
  cat("Counts file column names:", colnames(counts), "\n")
  cat("Sample info file row names:", rownames(sample_info), "\n")
  # Check dimensions of counts and sample info
  cat("Counts file dimensions (genes x samples): ", dim(counts), "\n")
  cat("Sample info file dimensions (samples x metadata columns): ", dim(sample_info), "\n")
  stop("Please check for formatting errors in your input files.")
}

# Reorder sample_info to match counts column order
sample_info <- sample_info[colnames(counts), , drop = FALSE]

# Ensure condition is a factor
sample_info$condition <- as.factor(sample_info$condition)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = counts, colData = sample_info, design = ~ condition)

# Normalize counts
dds <- DESeq(dds)
normalized_counts <- counts(dds, normalized = TRUE)

# Save normalized counts
write.table(normalized_counts, 
            file = file.path(deseq2_output_dir, "normalized_counts.txt"), 
            sep = "\t", quote = FALSE, col.names = NA)

# Extract DEG results
res <- results(dds, contrast = c("condition", "Treated", "Control"))

#unique_conditions <- unique(sample_info$condition)
#if (length(unique_conditions) != 2) {
#  stop("Error: The 'condition' column must have exactly two levels for DE analysis.")
#}
#res <- results(dds, contrast = c("condition", unique_conditions[1], unique_conditions[2]))

# Order results by adjusted p-value
res <- res[order(res$padj, na.last = NA), ]

# Save DEG results
write.table(as.data.frame(res), 
            file = file.path(deseq2_output_dir, "DEG_results.txt"), 
            sep = "\t", quote = FALSE, col.names = NA)

# Summary of significant DEGs
summary(res)

# Filter significant DEGs (p-value < 0.005 and |log2FoldChange| > 1)
sig_res <- subset(res, !is.na(pvalue) & pvalue < 0.005 & abs(log2FoldChange) > 1)

# Save significant DEGs
write.table(as.data.frame(sig_res), 
            file = file.path(deseq2_output_dir, "Significant_DEG_results.txt"), 
            sep = "\t", quote = FALSE, col.names = NA)

# Output the number of significant DEGs
cat("Number of significant DEGs:", nrow(sig_res), "\n")

# MA plot
pdf(file.path(deseq2_output_dir, "MA_plot.pdf"))
plotMA(res, ylim = c(-2, 2), main = "MA Plot")
dev.off()

# Convert results to data frame for volcano plot
res_df <- as.data.frame(res)

# Create and save volcano plot
pdf(file.path(deseq2_output_dir, "Volcano_plot.pdf"))
ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(alpha = 0.5, aes(color = pvalue < 0.005)) +
  theme_minimal() +
  xlab("Log2 Fold Change") +
  ylab("-Log10 P-value") +
  xlim(-10, 10) +
  ylim(0, 50) +
  ggtitle("Volcano Plot") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()
