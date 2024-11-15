#!/usr/bin/env Rscript
###############################################################################
# AUTHOR : OLABIYI ADEREMI OBAYOMI
# DESCRIPTION: A script to generate taxonomy plots at different taxonomy levels
# E-mail: obadbotanist@yahoo.com
# Created: November 2024
# example: Rscript taxonomy_plots.R \
#                  --metadata-table 'mapping/GLDS-487_amplicon_v1_runsheet.csv' \
#                  --feature-table 'data/counts_GLAmpSeq.tsv' \
#                  --taxonomy-table 'data/taxonomy_GLAmpSeq.tsv' \
#                  --group 'groups' \
#                  --samples-column 'Sample Name'
###############################################################################

library(optparse)



######## -------- Get input variables from the command line ----##############

version <- 1.0 

# Input options 
option_list <- list(
  
  make_option(c("-m", "--metadata-table"), type="character", default=NULL, 
              help="path to a comma separated samples metadata file with the 
              group/treatment to be analyzed.",
              metavar="path"),
  
  make_option(c("-f", "--feature-table"), type="character", default=NULL, 
              help="path to a tab separated samples feature table 
              i.e. ASV or OTU table.",
              metavar="path"),
  
  make_option(c("-t", "--taxonomy-table"), type="character", default=NULL, 
              help="path to feature taxonomy table i.e. ASV taxonomy table.",
              metavar="path"),
  
  make_option(c("-g", "--group"), type="character", default="groups", 
              help="Column in metadata to be analyzed",
              metavar="groups"),
  
  make_option(c("-s", "--samples-column"), type="character", default="Sample Name", 
              help="Column in metadata containing the sample names in the feature table. \
                    Deafault: 'Sample Name' ",
              metavar="Sample Name"),
  
  make_option(c("-o", "--output-prefix"), type="character", default="", 
              help="Unique name to tag onto output files. Default: empty string.",
              metavar=""),
  
  make_option(c("-c", "--abundance-cutoff"), type="numeric", default=0.2, 
              help="A fraction defining how abundant features most be to be \
                   analyzes. Default: 1/5. ",
              metavar="0.2"),
  
  make_option(c("-r", "--remove-rare"), type="logical", default=FALSE, 
              help="Should rare features be filtered out?. \
                   Default: FALSE. ", action= "store_true",
              metavar="FALSE"),
  
  make_option(c("-y", "--assay-suffix"), type="character", default="_GLAmpSeq", 
              help="Genelab assay suffix.", metavar="GLAmpSeq"),
  
  make_option(c("--version"), action = "store_true", type="logical", 
              default=FALSE,
              help="Print out version number and exit.", metavar = "boolean")
)


opt_parser <- OptionParser(
  option_list=option_list,
  usage = "Rscript %prog \\
                  --metadata-table 'mapping/GLDS-487_amplicon_v1_runsheet.csv' \\
                  --feature-table 'data/counts_GLAmpSeq.tsv' \\
                  --taxonomy-table 'data/taxonomy_GLAmpSeq.tsv' \\
                  --group 'groups' \\
                  --samples-column 'Sample Name' ",
  description = paste("Author: Olabiyi Aderemi Obayomi",
                      "\nEmail: olabiyi.a.obayomi@nasa.gov",
                      "\nA script to generate taxonomy plots at different taxonomy levels.",
                      "\nIt outputs sample and group taxonomy plots ", sep="")
)




opt <- parse_args(opt_parser)


if (opt$version) {
  cat("taxonomy_plots.R version: ", version, "\n")
  options_tmp <- options(show.error.messages=FALSE)
  on.exit(options(options_tmp))
  stop()
}



if(is.null(opt[["metadata-table"]])) {
  stop("Path to a metadata file must be set.")
}

if(is.null(opt[["feature-table"]])) {
  stop("Path to a feature table e.g. ASV table file must be set.")
}


if(is.null(opt[["taxonomy-table"]])) {
  stop("Path to a metadata file must be set.")
}

if(opt[["group"]] == "groups") {
  message("Alpha diversity will be run on the default 'groups' column \n")
}

if(opt[["samples-column"]] == "Sample Name") {
  message("I will assume that the sample names are in a column named 'Sample Name' \n")
}
library(tidyverse)
library(dendextend)
library(DESeq2)



# ------ Collecting the required input variables ---------- #

# Group in metadata to analyze
group <- opt[["group"]]  # "groups"
samples_column <- opt[["samples-column"]] # "Sample Name"
threads <- opt[["cpus"]] # 8
metadata_file <- opt[["metadata-table"]]
taxonomy_file <-  opt[["taxonomy-table"]]
feature_table_file <- opt[["feature-table"]]
feature <- opt[["feature-type"]]   # "ASV"
output_prefix <- opt[["output-prefix"]]
assay_suffix <- opt[["assay-suffix"]]

# taxon / ASV prevalence cutoff
prevalence_cutoff <- opt[["prevalence-cutoff"]] # 0.15 (15%)
# sample / library read count cutoff
library_cutoff <- opt[["library-cutoff"]]  # 100
diff_abund_out_dir <- "differential_abundance/"
if(!dir.exists(diff_abund_out_dir)) dir.create(diff_abund_out_dir)


de_out_dir <- file.path(plots_dir, "da")
abundance_out_dir <- file.path(de_out_dir, "differential_abundance")
volcano_out_dir <- file.path(de_out_dir, "volcano")



# ------------------------ Read metadata ---------------------------------- #
metadata <- read_csv(metadata_file)  %>% as.data.frame()
rownames(metadata) <- metadata[[samples_column]]



# -------------------------- Read Feature table  -------------------------- #
feature_table <- read_delim(file = feature_table_file) %>% as.data.frame()

# Set the feature id column as the row names of the feature table
# This assumes that the first column contains the feature ids e.g. ASV ID 
rownames(feature_table) <- feature_table[,1]
feature_names <- feature_table[,1]
# Drop the feature column
feature_table <- feature_table[, -1] %>% as.data.frame()
rownames(feature_table) <-  feature_names


# ------------------------ Read Taxonomy table ---------------------------- #
taxonomy <-  read_delim(file = taxonomy_file) %>% as.data.frame()
# Set the feature id column as the row names of the taxonomy table
# This assumes that the first column contains the feature ids e.g. ASV ID 
rownames(taxonomy) <- taxonomy[,1]
taxonomy_table  <- taxonomy[, -1]
feature_names <- rownames(taxonomy_table)
taxonomy_table  <- process_taxonomy(taxonomy_table)
rownames(taxonomy_table) <- feature_names

print(glue("There are {sum(taxonomy_table$domain == 'Other')} features without 
           taxonomy assignments. Dropping them ..."))

# Dropping features that couldn't be assigned taxonomy
taxonomy_table <- taxonomy_table[-which(taxonomy_table$domain == 'Other'),]

# Get long asv taxonomy names and clean
species <- taxonomy_table %>%
  unite(species,domain:species,sep = ";") %>% # Generalize this line -------- 
  pull %>% str_replace_all("Other", "_")

taxonomy_table <- fix_names(taxonomy_table, "Other", ";_")

taxonomy_table[,"species"] <- species


# ---------------------- Subset tables ------------------------------------- #

# Get features common to the taxonomy and feature table 
common_ids <- intersect(rownames(feature_table), rownames(taxonomy_table))

# Subset the feature and taxonomy tables to contain 
# only features found in both table
feature_table <- feature_table[common_ids,]
taxonomy_table <- taxonomy_table[common_ids,]





# 6 Statistically testing for differences

#### pairwise comparisons
unique_groups <- unique(runsheet$groups)
deseq_obj <- phyloseq_to_deseq2(physeq = ASV_physeq, design = ~groups)

# add pseudocount if any 0 count samples are present
if (sum(colSums(counts(deseq_obj)) == 0) > 0) {
  count_data <- counts(deseq_obj) + 1 
  
  count_data <- as.matrix(apply(count_data, 2, as.integer))
  rownames(count_data) <- rownames(counts(deseq_obj))
  colnames(count_data) <- colnames(counts(deseq_obj))
  counts(deseq_obj) <- count_data
}
# https://rdrr.io/bioc/phyloseq/src/inst/doc/phyloseq-mixture-models.R 
deseq_modeled <- tryCatch({
  # Attempt to run DESeq
  DESeq(deseq_obj)
}, error = function(e) {
  message("Error encountered in DESeq, applying alternative method for size factor estimation...")
  
  # Define the geometric mean function
  gm_mean = function(x, na.rm=TRUE) {
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  geoMeans = apply(counts(deseq_obj), 1, gm_mean)
  
  # Apply the alternative size factor estimation method
  deseq_obj <- estimateSizeFactors(deseq_obj, geoMeans=geoMeans)
  
  # Call DESeq again with alternative geom mean size est
  DESeq(deseq_obj)
})

# save final differential abundance counts, individual group comparison results

write.table(counts(deseq_modeled, normalized=TRUE), 
            file = file.path(de_out_dir, paste0(output_prefix,
                                                "normalized_counts", 
                                                assay_suffix, ".tsv")), 
            sep="\t", row.names=TRUE, quote=FALSE)
# make the volcanoplot
plot_comparison <- function(group1, group2) {
  plot_width_inches = 11.1
  plot_height_inches = 8.33
  
  deseq_res <- results(deseq_modeled, contrast = c("groups", group1, group2))
  norm_tab <- counts(deseq_modeled, normalized = TRUE) %>% data.frame()
  
  volcano_data <- as.data.frame(deseq_res)
  
  p_val <- 0.1
  volcano_data <- volcano_data[!is.na(volcano_data$padj), ]
  volcano_data$significant <- volcano_data$padj <= p_val #also logfc cutoff?
  
  ######Long x-axis label adjustments##########
  x_label <- paste("Log2 Fold Change\n(",group1," vs ",group2,")")
  label_length <- nchar(x_label)
  max_allowed_label_length = plot_width_inches * 10
  
  # Construct x-axis label with new line breaks if was too long
  if (label_length > max_allowed_label_length){
    x_label <- paste("Log2 Fold Change\n\n(", group1, "\n vs \n", group2, ")", sep="")
  }
  #######################################

  # ASVs promoted in space on right, reduced on left
  p <- ggplot(volcano_data, aes(x=log2FoldChange, y=-log10(padj), color=significant)) +
    geom_point(alpha=0.7, size=2) +
    scale_color_manual(values=c("black", "red"), 
                       labels=c(paste0("padj > ", p_val), 
                                paste0("padj \u2264 ", p_val))) +
    theme_bw() +
    labs(title="Volcano Plot",
         x=x_label,
         y="-Log10 P-value",
         color=paste0("")) +
    theme(legend.position="top")
  
  # label points and plot
  top_points <- volcano_data %>%
    arrange(padj) %>%
    filter(significant) %>%
    head(10)
  
  volcano_plot <- p + geom_text_repel(data=top_points, aes(label=row.names(top_points)), size=3)
  ggsave(filename=file.path(volcano_out_dir, paste0(output_prefix,
                                                    "volcano_", 
                                                    gsub(" ", "_", group1), 
                                                    "_vs_", 
                                                    gsub(" ", "_", group2), ".png")),
         plot=volcano_plot,
         width = plot_width_inches, height = plot_height_inches, dpi = 300)
  
  write.csv(deseq_res, file = file.path(abundance_out_dir,
                                        paste0(output_prefix, 
                                               gsub(" ", "_", group1),
                                               "_vs_", gsub(" ", "_", group2),
                                               ".csv")))
}


# setting up pairwise comparisons and running
comparisons <- expand.grid(group1 = unique_groups, group2 = unique_groups)
comparisons <- subset(comparisons, group1 != group2)

apply(comparisons, 1, function(pair) plot_comparison(pair['group1'], pair['group2']))
