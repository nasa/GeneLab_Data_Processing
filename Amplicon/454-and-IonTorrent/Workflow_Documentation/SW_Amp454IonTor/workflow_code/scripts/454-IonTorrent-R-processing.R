##################################################################################
## R processing script for 454/Ion Torrent amplicon data                        ##
## Developed by Michael D. Lee (Mike.Lee@nasa.gov)                              ##
##################################################################################

# as called from the associated Snakefile, this expects to be run as: Rscript full-R-processing.R <otus_fasta_path> <trimmed_dir> <filtered_dir> <final_outputs_directory> <output_prefix> <target_region>

# setting variables used within R:
args <- commandArgs(trailingOnly = TRUE)

suppressWarnings(otus_fasta_path <- args[1])
suppressWarnings(trimmed_dir <- args[2])
suppressWarnings(filtered_dir <- args[3])
suppressWarnings(final_outputs_dir <- args[4])
suppressWarnings(output_prefix <- args[5])
suppressWarnings(target_region <- args[6])

    # loading libraries
library(DECIPHER)
library(biomformat)

    ### assigning taxonomy ###
  # reading OTUs into a DNAStringSet object
dna <- readDNAStringSet(paste0(final_outputs_dir, output_prefix, "OTUs.fasta"))


  # downloading reference R taxonomy object
cat("\n\n  Downloading reference database...\n\n")

if ( target_region == "16S" ) {
    download.file("http://www2.decipher.codes/Classification/TrainingSets/SILVA_SSU_r138_2019.RData", "SILVA_SSU_r138_2019.RData")
    load("SILVA_SSU_r138_2019.RData")
    file.remove("SILVA_SSU_r138_2019.RData")
} else if ( target_region == "ITS" ) {
    download.file("http://www2.decipher.codes/Classification/TrainingSets/UNITE_v2020_February2020.RData", "UNITE_v2020_February2020.RData")
    load("UNITE_v2020_February2020.RData")
    file.remove("UNITE_v2020_February2020.RData")
}


# assigning taxonomy
cat("\n\n  Assigning taxonomy...\n\n")

tax_info <- IdTaxa(dna, trainingSet, strand="both", processors=NULL)

cat("\n\n  Making and writing out tables...\n\n")

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
otu_ids <- names(tax_info)
tax_tab <- data.frame("OTU_ID"=otu_ids, tax_tab, check.names=FALSE)

write.table(tax_tab, paste0(final_outputs_dir, output_prefix, "taxonomy.tsv"), sep = "\t", quote=F, row.names=FALSE)

    # reading in counts table to generate other outputs
otu_tab <- read.table(paste0(final_outputs_dir, output_prefix, "counts.tsv"), sep="\t", header=TRUE, check.names=FALSE)

    # generating and writing out biom file format
biom_object <- make_biom(data=otu_tab, observation_metadata=tax_tab)
write_biom(biom_object, paste0(final_outputs_dir, output_prefix, "taxonomy-and-counts.biom"))

    # making a tsv of combined tax and counts
tax_and_count_tab <- merge(tax_tab, otu_tab)
write.table(tax_and_count_tab, paste0(final_outputs_dir, output_prefix, "taxonomy-and-counts.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

# making final count summary table
cutadapt_tab <- read.table(paste0(trimmed_dir, output_prefix, "trimmed-read-counts.tsv"), sep="\t", header=TRUE)
bbduk_tab <- read.table(paste0(filtered_dir, output_prefix, "filtered-read-counts.tsv"), sep="\t", header=TRUE)[,c(1,3)]
    # re-reading in counts table to this time set first col as rownames (rather than doing it another way)
otu_tab <- read.table(paste0(final_outputs_dir, output_prefix, "counts.tsv"), sep="\t", header=TRUE, check.names=FALSE, row.names = 1)
mapped_sums <- colSums(otu_tab)
mapped_tab <- data.frame(sample=names(mapped_sums), mapped_to_OTUs=mapped_sums, row.names=NULL)

t1 <- merge(cutadapt_tab, bbduk_tab)
count_summary_tab <- merge(t1, mapped_tab)
count_summary_tab$final_perc_reads_retained <- round(count_summary_tab$mapped_to_OTUs / count_summary_tab$raw_reads * 100, 2)

write.table(count_summary_tab, paste0(final_outputs_dir, output_prefix, "read-count-tracking.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

cat("\n\n  Session info:\n\n")
sessionInfo()
