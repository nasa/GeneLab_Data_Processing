#!/usr/bin/env Rscript

# featurecounts_genes_per_sample.R
# Script to count the number of non-zero genes per sample from FeatureCounts output
# Used in the QUANTIFY_RSEM_GENES process of the RNA-seq workflow

# Load required packages
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# Get the assay suffix from environment if set, otherwise use default
assay_suffix <- Sys.getenv("assay_suffix", "_GLbulkRNAseq")

# Log function
log_message <- function(message) {
  cat(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", message, "\n"))
}

log_message("Starting featurecounts_genes_per_sample.R script")

# Find the FeatureCounts output file
fc_files <- list.files("03-FeatureCounts", pattern = "FeatureCounts.*\\.tsv$", full.names = TRUE)

if (length(fc_files) == 0) {
  stop("No FeatureCounts output file found in 03-FeatureCounts directory")
}

# Use the first file if multiple files are found
fc_file <- fc_files[1]
log_message(paste("Using FeatureCounts file:", fc_file))

# Read sample names
if (file.exists("samples.txt")) {
  samples <- readLines("samples.txt")
  log_message(paste("Found", length(samples), "samples in samples.txt"))
} else {
  log_message("Warning: samples.txt not found. Will use column headers from count table.")
  samples <- NULL
}

# Read FeatureCounts data, skipping the first line which contains command info
log_message("Reading FeatureCounts data...")
tryCatch({
  data <- read.table(fc_file, 
                    header = TRUE, 
                    skip = 1,
                    sep = "\t",
                    comment.char = "",
                    quote = "",
                    stringsAsFactors = FALSE,
                    check.names = FALSE)
  
  # Set row names to Geneid column and remove it
  rownames(data) <- data$Geneid
  data$Geneid <- NULL
  
  log_message(paste("Read data with", nrow(data), "rows and", ncol(data), "columns"))
}, error = function(e) {
  stop(paste("Error reading FeatureCounts file:", e$message))
})

# Remove metadata columns if they exist
metadata_cols <- c("Chr", "Start", "End", "Strand", "Length")
data <- data[, !colnames(data) %in% metadata_cols]

# Clean column names: remove path, .bam extension
colnames(data) <- gsub("^.*/([^/]+)\\.bam$", "\\1", colnames(data))
colnames(data) <- gsub("\\.bam$", "", colnames(data))

log_message("Column names after cleaning:")
log_message(paste(colnames(data), collapse = ", "))

# Verify all samples are present in the data
if (!is.null(samples)) {
  missing_samples <- samples[!samples %in% colnames(data)]
  if (length(missing_samples) > 0) {
    log_message(paste("Warning: The following samples from samples.txt are not found in the count table:", 
                       paste(missing_samples, collapse = ", ")))
  }
}

# Count number of genes with non-zero counts for each sample
log_message("Counting number of genes with non-zero counts for each sample...")
NumNonZeroGenes <- colSums(data > 0)
NumNonZeroGenes <- data.frame(
  "Sample" = names(NumNonZeroGenes),
  "Number_of_genes_with_non_zero_counts" = NumNonZeroGenes,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Set output file name
output_file <- paste0("NumNonZeroGenes", assay_suffix, ".csv")
log_message(paste("Writing results to", output_file))

# Export results
write.csv(NumNonZeroGenes, 
          file = output_file,
          row.names = FALSE,
          quote = TRUE)

log_message("Script completed successfully")
log_message("Session Info:")
print(sessionInfo()) 