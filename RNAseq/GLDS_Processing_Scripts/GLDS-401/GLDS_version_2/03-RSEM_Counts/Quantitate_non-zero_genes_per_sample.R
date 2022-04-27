library(tximport)
library(tidyverse)

work_dir="/GLDS-401/processing_scripts/03-RSEM_Counts"
counts_dir="/GLDS-401/03-RSEM_Counts"

setwd(file.path(work_dir))

### Pull in sample names ###
samples <- read.csv(Sys.glob(file.path(work_dir,"samples.txt")), header = FALSE, row.names = 1, stringsAsFactors = TRUE)

##### Import RSEM Gene Count Data
files <- list.files(file.path(counts_dir),pattern = ".genes.results", full.names = TRUE)
### reorder the genes.results files to match the ordering of the samples in the metadata file
files <- files[sapply(rownames(samples), function(x)grep(x, files, value=FALSE, fixed=TRUE))]
names(files) <- rownames(samples)
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)

##### Export unnormalized gene counts table
setwd(file.path(counts_dir))
write.csv(txi.rsem$counts,file='RSEM_Unnormalized_Counts.csv')

##### Count the number of genes with non-zero counts for each sample 
rawCounts <- txi.rsem$counts
NumNonZeroGenes <- (as.matrix(colSums(rawCounts > 0), row.names = 1))
colnames(NumNonZeroGenes) <- c("Number of genes with non-zero counts")

##### Export the number of genes with non-zero counts for each sample
setwd(file.path(counts_dir))
write.csv(NumNonZeroGenes,file='NumNonZeroGenes.csv')

## print session info ##
print("Session Info below: ")
print("")
sessionInfo()
