#! /usr/bin/env Rscript
library(tximport)
library(tidyverse)

# Ref: https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf
# column 1: gene ID
# column 2: counts for unstranded RNA-seq
# column 3: counts for the 1st read strand aligned with RNA (htseq-count option -s yes)
# column 4: counts for the 2nd read strand aligned with RNA (htseq-count option -s reverse)

# work_dir=""
counts_dir <- "02-STAR_Alignment"
args <- commandArgs(trailingOnly = TRUE)
strandedness <- args[1]
# map to column names
column_target_name <- list(unstranded = "V2", sense = "V3", antisense = "V4")[[strandedness]]
if (is.null(column_target_name)) {
    stop(sprintf("Could not map column based on strandedness provided as CLI argument: '%s'", strandedness))
}
# setwd(file.path(work_dir))

### Pull in sample names ###
samples <- read.csv(Sys.glob("samples.txt"), header = FALSE, row.names = 1, stringsAsFactors = TRUE)

##### Import RSEM Gene Count Data
files <- list.files(file.path(counts_dir), pattern = "_ReadsPerGene.out.tab", full.names = TRUE)
print(sprintf("DEBUG: %s: %s files", Sys.time(), paste(files, collapse = ", ")))

### reorder the genes.results files to match the ordering of the samples in the metadata file
reordering <- sapply(rownames(samples), function(x) grep(paste0("02-STAR_Alignment/", x, "_ReadsPerGene.out.tab$"), files, value = FALSE))
print(sprintf("DEBUG: %s: %s", Sys.time(), toString(reordering)))
files <- files[reordering]
names(files) <- rownames(samples)

# define load function
load_star_table <- function(file, target_column, remove_aggregation_rows = TRUE) {
  #' loads a *_ReadsPerGene.out.tab file. 
  #' Specifically selecting only the 'target_column'
    print(sprintf("DEBUG: %s: Attempting to load column '%s' from '%s'", Sys.time(), target_column, file))
    df <- read.table(file, header = FALSE) %>%
            select(c("V1", target_column)) %>%
            rename(!!names(file) := target_column) %>%
            rename(geneID = V1) %>%
            arrange(geneID)

    if (remove_aggregation_rows) {
        print(sprintf("DEBUG: %s: Removing all rows that start with 'N_', e.g. 'N_unmapped'", Sys.time()))
        df <- df %>% filter(!str_detect(geneID, "^N_"))
    }

    # reformatting to have rownames denote geneID column
    df <- df %>% tibble::column_to_rownames(var = "geneID")
    print(sprintf("DEBUG: %s: Loaded dataframe with column '%s' from file '%s'", Sys.time(), colnames(df), file))
    return(df)
}

# here this is just the initial dataframe to be extended
# by subsequenct sample wise files
df_full <- load_star_table(files[1], column_target_name)

# Note: we skip 1 because that is already loaded
for (i in 2:length(files)) {
    df_new <- load_star_table(files[i], column_target_name)
    print(sprintf("DEBUG: %s: Extending dataframe with column '%s'", Sys.time(), colnames(df_new)))
    df_full <- cbind(df_full, df_new)
}

##### Export unnormalized gene counts table
# setwd(file.path(counts_dir))
write.csv(df_full, file = "STAR_Unnormalized_Counts_GLbulkRNAseq.csv")

##### Count the number of genes with non-zero counts for each sample
num_nonzero_genes <- (as.matrix(colSums(df_full > 0), row.names = 1))
colnames(num_nonzero_genes) <- c("Number of genes with non-zero counts")

##### Export the number of genes with non-zero counts for each sample
# setwd(file.path(counts_dir))
write.csv(num_nonzero_genes, file = "STAR_NumNonZeroGenes_GLbulkRNAseq.csv")

## print session info ##
print("Session Info below: ")
print("")
sessionInfo()