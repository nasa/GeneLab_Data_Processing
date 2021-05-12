##################################################################################
## R processing script for 16S data of GLDS-280                                 ##
## https://genelab-data.ndc.nasa.gov/genelab/accession/GLDS-280/                ##
##                                                                              ##
## This code as written expects to be run within the Processing_Info/ directory ##
## Processed by Michael D. Lee (Mike.Lee@nasa.gov)                              ##
##################################################################################

# expects to be run as: Rscript full-R-processing.R <left_trunc> <right_trunc> <left_maxEE> <right_maxEE>
  # where <left_trim> and <right_trim> are the values to be passed to the truncLen parameter of dada2's filterAndTrim()
  # and <left_maxEE> and <right_maxEE> are the values to be passed to the maxEE parameter of dada2's filterAndTrim()

# checking 2 arguments provided and they are integers:
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
    stop("Exactly 4 integer arguments are required, see top of script for info.", call.=FALSE)
} else {
    suppressWarnings(left_trunc <- as.integer(args[1]))
    suppressWarnings(right_trunc <- as.integer(args[2]))
    suppressWarnings(left_maxEE <- as.integer(args[3]))
    suppressWarnings(right_maxEE <- as.integer(args[4]))
}

if ( is.na(left_trunc) || is.na(right_trunc) || is.na(left_maxEE) || is.na(right_maxEE) ) {
    stop("All 4 arguments must be integers, see top of script for info.", call.=FALSE)
}

# general procedure comes largely from these sources:
  # https://benjjneb.github.io/dada2/tutorial.html
  # https://astrobiomike.github.io/amplicon/dada2_workflow_ex

BiocManager::version() # 3.9

  # loading libraries
library(dada2); packageVersion("dada2") # 1.12.1
library(DECIPHER); packageVersion("DECIPHER") # 2.12.0
library(biomformat); packageVersion("biomformat") # 1.12.0

    ### general processing ###
  # reading in unique sample names into variable
sample.names <- scan("unique-sample-IDs.txt", what="character")

  # setting variables holding the paths to the forward and reverse primer-trimmed reads
forward_reads <- paste0("../Trimmed_Sequence_Data/", sample.names, "-R1-trimmed.fastq.gz")
reverse_reads <- paste0("../Trimmed_Sequence_Data/", sample.names, "-R2-trimmed.fastq.gz")

  # setting variables holding what will be the output paths of all forward and reverse filtered reads
forward_filtered_reads <- paste0("../Filtered_Sequence_Data/", sample.names, "-R1-filtered.fastq.gz")
reverse_filtered_reads <- paste0("../Filtered_Sequence_Data/", sample.names, "-R2-filtered.fastq.gz")

  # adding sample names to the vectors holding the filtered-reads' paths
names(forward_filtered_reads) <- sample.names
names(reverse_filtered_reads) <- sample.names

  # running filering step
    # reads are written to the files specified in the variables, the "filtered_out" object holds the summary results within R
filtered_out <- filterAndTrim(fwd=forward_reads, forward_filtered_reads, reverse_reads, reverse_filtered_reads, truncLen=c(left_trunc,right_trunc), maxN=0, maxEE=c(left_maxEE,right_maxEE), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)

  # making and writing out summary table that includes counts of filtered reads
    # helper function
getN <- function(x) sum(getUniques(x))

filtered_count_summary_tab <- data.frame(sample=sample.names, cutadapt_trimmed=filtered_out[,1], dada2_filtered=filtered_out[,2])
write.table(filtered_count_summary_tab, "../Filtered_Sequence_Data/filtered-read-counts.tsv", sep="\t", quote=F, row.names=F)

  # learning errors step
forward_errors <- learnErrors(forward_filtered_reads, multithread=TRUE)
reverse_errors <- learnErrors(reverse_filtered_reads, multithread=TRUE)

  # inferring sequences
forward_seqs <- dada(forward_filtered_reads, err=forward_errors, pool="pseudo", multithread=TRUE)
reverse_seqs <- dada(reverse_filtered_reads, err=reverse_errors, pool="pseudo", multithread=TRUE)

  # merging forward and reverse reads
merged_contigs <- mergePairs(forward_seqs, forward_filtered_reads, reverse_seqs, reverse_filtered_reads, verbose=TRUE)

  # generating a sequence table that holds the counts of each sequence per sample
seqtab <- makeSequenceTable(merged_contigs)

  # removing putative chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

  # checking what percentage of sequences were retained after chimera removal
sum(seqtab.nochim)/sum(seqtab) * 100

  # making and writing out a summary table that includes read counts at all steps
raw_and_trimmed_read_counts <- read.table("../Trimmed_Sequence_Data/trimmed-read-counts.tsv", header=T, sep="\t")
    # reading in filtered read counts
filtered_read_counts <- read.table("../Filtered_Sequence_Data/filtered-read-counts.tsv", header=T, sep="\t")

count_summary_tab <- data.frame(raw_and_trimmed_read_counts, dada2_filtered=filtered_read_counts[,3],
                                dada2_denoised_F=sapply(forward_seqs, getN),
                                dada2_denoised_R=sapply(reverse_seqs, getN),
                                dada2_merged=rowSums(seqtab),
                                dada2_chimera_removed=rowSums(seqtab.nochim),
                                final_perc_reads_retained=round(rowSums(seqtab.nochim)/raw_and_trimmed_read_counts$raw_reads * 100, 1),
                                row.names=NULL)

write.table(count_summary_tab, "../Final_Outputs/read-count-tracking.tsv", sep = "\t", quote=F, row.names=F)

    ### assigning taxonomy ###
  # creating a DNAStringSet object from the ASVs
dna <- DNAStringSet(getSequences(seqtab.nochim))

  # downloading reference R taxonomy object (at some point this will be stored somewhere on GeneLab's server and we won't download it, but should leave the code here, just commented out)
download.file("http://www2.decipher.codes/Classification/TrainingSets/SILVA_SSU_r138_2019.RData", "SILVA_SSU_r138_2019.RData")
  # loading reference taxonomy object
load("SILVA_SSU_r138_2019.RData")
  # removing downloaded file
file.remove("SILVA_SSU_r138_2019.RData")

  # classifying
tax_info <- IdTaxa(dna, trainingSet, strand="both", processors=NULL)

    ### generating and writing out standard outputs ###
  # giving our sequences more manageable names (e.g. ASV_1, ASV_2..., rather than the sequence itself)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

  # making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "../Final_Outputs/ASVs.fasta")

  # making and writing out a count table:
asv_tab <- t(seqtab.nochim)
asv_ids <- sub(">", "", asv_headers)
row.names(asv_tab) <- NULL
asv_tab <- data.frame("ASV_ID"=asv_ids, asv_tab, check.names=FALSE)

write.table(asv_tab, "../Final_Outputs/counts.tsv", sep="\t", quote=F, row.names=FALSE)

  # making and writing out a taxonomy table:
    # creating vector of desired ranks
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")

  # creating table of taxonomy and setting any that are unclassified as "NA"
tax_tab <- t(sapply(tax_info, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))

colnames(tax_tab) <- ranks
row.names(tax_tab) <- NULL
tax_tab <- data.frame("ASV_ID"=asv_ids, tax_tab, check.names=FALSE)

write.table(tax_tab, "../Final_Outputs/taxonomy.tsv", sep = "\t", quote=F, row.names=FALSE)

    ### generating and writing out biom file format ###
biom_object <- make_biom(data=asv_tab, observation_metadata=tax_tab)
write_biom(biom_object, "../Final_Outputs/taxonomy-and-counts.biom")

    # making a tsv of combined tax and counts
tax_and_count_tab <- merge(tax_tab, asv_tab)
write.table(tax_and_count_tab, "../Final_Outputs/taxonomy-and-counts.tsv", sep="\t", quote=FALSE, row.names=FALSE)
