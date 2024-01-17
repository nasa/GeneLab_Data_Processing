##################################################################################
## R processing script for Illumina single-end amplicon data                    ##
## Developed by Michael D. Lee (Mike.Lee@nasa.gov)                              ##
##################################################################################

# as called from the associated Snakefile, this expects to be run as: Rscript full-R-processing.R <left_trunc> <left_maxEE> <TRUE/FALSE - GL trimmed primers or not> <unique-sample-IDs-file> <starting_reads_dir_for_R> <filtered_reads_dir> <input_file_R1_suffix> <filtered_filename_R1_suffix> <final_outputs_directory> <output_prefix> <target_region> <assay_suffix>
    # where <left_trim> is the value to be passed to the truncLen parameter of dada2's filterAndTrim()
    # and <left_maxEE> is the value to be passed to the maxEE parameter of dada2's filterAndTrim()

# checking arguments were provided, first 2 are integers, and setting variables used within R:
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 12) {
    stop("At least 12 positional arguments are required, see top of this R script for more info.", call.=FALSE)
} else {
    suppressWarnings(left_trunc <- as.integer(args[1]))
    suppressWarnings(left_maxEE <- as.integer(args[2]))

    suppressWarnings(GL_trimmed_primers <- args[3])
    suppressWarnings(sample_IDs_file <- args[4])
    suppressWarnings(input_reads_dir <- args[5])
    suppressWarnings(filtered_reads_dir <- args[6])
    suppressWarnings(input_file_R1_suffix <- args[7])
    suppressWarnings(filtered_filename_R1_suffix <- args[8])
    suppressWarnings(final_outputs_dir <- args[9])
    suppressWarnings(output_prefix <- args[10])
    suppressWarnings(target_region <- args[11])
    suppressWarnings(assay_suffix <- args[12])

}

if ( is.na(left_trunc) || is.na(left_maxEE) ) {
    stop("The 2 first arguments must be integers, see top of R script for more info.", call.=FALSE)
}

if ( ! GL_trimmed_primers %in% c("TRUE", "FALSE") ) {
    stop("The third positional argument needs to be 'TRUE' or 'FALSE' for whether or not GL trimmed primers on this dataset, see top of R script for more info.", call.=FALSE)    
} else {
    GL_trimmed_primers <- as.logical(GL_trimmed_primers)
}

# general procedure comes largely from these sources:
    # https://benjjneb.github.io/dada2/tutorial.html
    # https://astrobiomike.github.io/amplicon/dada2_workflow_ex

    # loading libraries
library(dada2); packageVersion("dada2")
library(DECIPHER); packageVersion("DECIPHER")
library(biomformat); packageVersion("biomformat")

    ### general processing ###
    # reading in unique sample names into variable
sample.names <- scan(sample_IDs_file, what="character")

    # setting variables holding the paths to the forward primer-trimmed reads (or "raw" if primers were trimmed prior to submission to GeneLab)
forward_reads <- paste0(input_reads_dir, sample.names, input_file_R1_suffix)

    # setting variables holding what will be the output paths of all forward filtered reads
forward_filtered_reads <- paste0(filtered_reads_dir, sample.names, filtered_filename_R1_suffix)

    # adding sample names to the vectors holding the filtered-reads' paths
names(forward_filtered_reads) <- sample.names

    # running filering step
    # reads are written to the files specified in the variables, the "filtered_out" object holds the summary results within R
filtered_out <- filterAndTrim(fwd=forward_reads, forward_filtered_reads, truncLen=c(left_trunc), maxN=0, maxEE=c(left_maxEE), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=10)

    # making and writing out summary table that includes counts of filtered reads
if ( GL_trimmed_primers ) {

    filtered_count_summary_tab <- data.frame(sample=sample.names, cutadapt_trimmed=filtered_out[,1], dada2_filtered=filtered_out[,2])

} else {

    filtered_count_summary_tab <- data.frame(sample=sample.names, starting_reads=filtered_out[,1], dada2_filtered=filtered_out[,2])

}

write.table(filtered_count_summary_tab, paste0(filtered_reads_dir, output_prefix, "filtered-read-counts_", assay_suffix, ".tsv"), sep="\t", quote=F, row.names=F)

    # learning errors step
forward_errors <- learnErrors(forward_filtered_reads, multithread=10)

    # inferring sequences
forward_seqs <- dada(forward_filtered_reads, err=forward_errors, pool="pseudo", multithread=10)

    # generating a sequence table that holds the counts of each sequence per sample
seqtab <- makeSequenceTable(forward_seqs)

    # removing putative chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=10, verbose=TRUE)

    # checking what percentage of sequences were retained after chimera removal
sum(seqtab.nochim)/sum(seqtab) * 100

    # making and writing out a summary table that includes read counts at all steps
    # helper function
getN <- function(x) sum(getUniques(x))

if ( GL_trimmed_primers ) {
    
    raw_and_trimmed_read_counts <- read.table(paste0(input_reads_dir, output_prefix, "trimmed-read-counts_", assay_suffix, ".tsv"), header=T, sep="\t")
    # reading in filtered read counts
    filtered_read_counts <- read.table(paste0(filtered_reads_dir, output_prefix, "filtered-read-counts_", assay_suffix, ".tsv"), header=T, sep="\t")

    count_summary_tab <- data.frame(raw_and_trimmed_read_counts, dada2_filtered=filtered_read_counts[,3],
                                    dada2_denoised_F=sapply(forward_seqs, getN),
                                    dada2_chimera_removed=rowSums(seqtab.nochim),
                                    final_perc_reads_retained=round(rowSums(seqtab.nochim)/raw_and_trimmed_read_counts$raw_reads * 100, 1),
                                    row.names=NULL)

} else {

    count_summary_tab <- data.frame(filtered_count_summary_tab,
                                dada2_denoised_F=sapply(forward_seqs, getN),
                                dada2_chimera_removed=rowSums(seqtab.nochim),
                                final_perc_reads_retained=round(rowSums(seqtab.nochim)/filtered_count_summary_tab$starting_reads * 100, 1),
                                row.names=NULL)

}

write.table(count_summary_tab, paste0(final_outputs_dir, output_prefix, "read-count-tracking_", assay_suffix, ".tsv"), sep = "\t", quote=F, row.names=F)

    ### assigning taxonomy ###
    # creating a DNAStringSet object from the ASVs
dna <- DNAStringSet(getSequences(seqtab.nochim))

    # downloading reference R taxonomy object (at some point this will be stored somewhere on GeneLab's server and we won't download it, but should leave the code here, just commented out)
cat("\n\n  Downloading reference database...\n\n")
if ( target_region == "16S" ) { 
    download.file("http://www2.decipher.codes/Classification/TrainingSets/SILVA_SSU_r138_2019.RData", "SILVA_SSU_r138_2019.RData")
    # loading reference taxonomy object
    load("SILVA_SSU_r138_2019.RData")
    # removing downloaded file
    file.remove("SILVA_SSU_r138_2019.RData")

    ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")

} else if (target_region == "ITS" ) {

    download.file("http://www2.decipher.codes/Classification/TrainingSets/UNITE_v2020_February2020.RData", "UNITE_v2020_February2020.RData")    
    # loading reference taxonomy object
    load("UNITE_v2020_February2020.RData")
    # removing downloaded file
    file.remove("UNITE_v2020_February2020.RData")

    ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")

} else if (target_region == "18S" ) {

    download.file("http://www2.decipher.codes/Classification/TrainingSets/PR2_v4_13_March2021.RData", "PR2_v4_13_March2021.RData")    
    # loading reference taxonomy object
    load("PR2_v4_13_March2021.RData")
    # removing downloaded file
    file.remove("PR2_v4_13_March2021.RData")

    ranks <- c("kingdom", "division", "phylum", "class", "order", "family", "genus", "species")

} else { 

    cat("\n\n  The requested target_region of ", target_region, " is not accepted.\n\n")
    quit(status = 1)
}

    # classifying
cat("\n\n  Assigning taxonomy...\n\n")
tax_info <- IdTaxa(dna, trainingSet, strand="both", processors=NULL)

    ### generating and writing out standard outputs ###
    # giving our sequences more manageable names (e.g. ASV_1, ASV_2..., rather than the sequence itself)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

# adding additional prefix to ASV headers of target region if one was provided
if ( output_prefix != "" ) {
    for (i in 1:dim(seqtab.nochim)[2]) {
        asv_headers[i] <- paste(">ASV", target_region, i, sep="_")
    }
} else {
    for (i in 1:dim(seqtab.nochim)[2]) {
        asv_headers[i] <- paste(">ASV", i, sep="_")
    }
}

cat("\n\n  Making and writing outputs...\n\n")
    # making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, paste0(final_outputs_dir, output_prefix, "ASVs_", assay_suffix, ".fasta"))

    # making and writing out a count table:
asv_tab <- t(seqtab.nochim)
asv_ids <- sub(">", "", asv_headers)
row.names(asv_tab) <- NULL
asv_tab <- data.frame("ASV_ID"=asv_ids, asv_tab, check.names=FALSE)

write.table(asv_tab, paste0(final_outputs_dir, output_prefix, "counts_", assay_suffix, ".tsv"), sep="\t", quote=F, row.names=FALSE)

    # making and writing out a taxonomy table:
    # vector of desired ranks was created above in ITS/16S target_region if statement

    # creating table of taxonomy and setting any that are unclassified as "NA"
tax_tab <- t(sapply(tax_info, function(x) {
    m <- match(ranks, x$rank)
    taxa <- x$taxon[m]
    taxa[startsWith(taxa, "unclassified_")] <- NA
    taxa
}))

colnames(tax_tab) <- ranks
row.names(tax_tab) <- NULL

# need to add domain values if this is ITS (due to how the reference taxonomy object is structured, which doesn't have a domain entry)
if ( target_region == "ITS" ) {

    # only want to add if at least a kingdom was identified (so we don't add Eukarya if nothing was found)
    new_vec <- ifelse(!is.na(tax_tab[, "kingdom"]) & tax_tab[, "kingdom"] != "NA", "Eukarya", "NA")

    tax_tab <- data.frame("ASV_ID"=asv_ids, "domain" = new_vec, tax_tab, check.names=FALSE)

} else {

    tax_tab <- data.frame("ASV_ID"=asv_ids, tax_tab, check.names=FALSE)

}

# need to change "kingdom" to "domain" if this is 18S (due to how the reference taxonomy object is structured)
if ( target_region == "18S" ) {
    colnames(tax_tab)[colnames(tax_tab) == "kingdom"] <- "domain"
}

write.table(tax_tab, paste0(final_outputs_dir, output_prefix, "taxonomy_", assay_suffix, ".tsv"), sep = "\t", quote=F, row.names=FALSE)

    ### generating and writing out biom file format ###
biom_object <- make_biom(data=asv_tab, observation_metadata=tax_tab)
write_biom(biom_object, paste0(final_outputs_dir, output_prefix, "taxonomy-and-counts_", assay_suffix, ".biom"))

    # making a tsv of combined tax and counts
tax_and_count_tab <- merge(tax_tab, asv_tab)
write.table(tax_and_count_tab, paste0(final_outputs_dir, output_prefix, "taxonomy-and-counts_", assay_suffix, ".tsv"), sep="\t", quote=FALSE, row.names=FALSE)

cat("\n\n  Session info:\n\n")
sessionInfo()
