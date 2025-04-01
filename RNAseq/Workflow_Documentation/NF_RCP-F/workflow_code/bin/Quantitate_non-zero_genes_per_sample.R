#! /usr/bin/env Rscript
library(tximport)
library(tidyverse)

#work_dir=""
counts_dir="03-RSEM_Counts"

#setwd(file.path(work_dir))

### Pull in sample names ###
samples <- read.csv(Sys.glob("samples.txt"), header = FALSE, row.names = 1, stringsAsFactors = TRUE)

##### Import RSEM Gene Count Data
files <- list.files(file.path(counts_dir),pattern = ".genes.results", full.names = TRUE)
### reorder the genes.results files to match the ordering of the samples in the metadata file
reordering <- sapply(rownames(samples), function(x)grep(paste0("03-RSEM_Counts/", x,".genes.results$"), files, value=FALSE))
files <- files[reordering]
names(files) <- rownames(samples)
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)

##### Export unnormalized gene counts table
#setwd(file.path(counts_dir))
write.csv(txi.rsem$counts,file='RSEM_Unnormalized_Counts_GLbulkRNAseq.csv')

##### Count the number of genes with non-zero counts for each sample
rawCounts <- txi.rsem$counts
NumNonZeroGenes <- (as.matrix(colSums(rawCounts > 0), row.names = 1))
colnames(NumNonZeroGenes) <- c("Number of genes with non-zero counts")

##### Export the number of genes with non-zero counts for each sample
#setwd(file.path(counts_dir))
write.csv(NumNonZeroGenes,file='RSEM_NumNonZeroGenes_GLbulkRNAseq.csv')

## print session info ##
print("Session Info below: ")
print("")
sessionInfo()
