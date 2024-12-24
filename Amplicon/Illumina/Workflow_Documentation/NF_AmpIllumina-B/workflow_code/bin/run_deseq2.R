#!/usr/bin/env Rscript
###############################################################################
# AUTHOR : OLABIYI ADEREMI OBAYOMI
# DESCRIPTION: A script to perform differential abundance testing using DESeq2
# E-mail: obadbotanist@yahoo.com
# Created: December 2024
# Example: Rscript run_deseq2.R \
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
              group/treatment to be compared.",
              metavar="path"),
  
  make_option(c("-f", "--feature-table"), type="character", default=NULL, 
              help="path to a tab separated samples feature table 
              i.e. ASV or OTU table.",
              metavar="path"),
  
  make_option(c("-a", "--feature-type"), type="character", default="ASV", 
              help="What feature type is in the feature table i.e ASV, OTU etc.
              This will be used to name the feature column in the output table.",
              metavar="ASV"),
  
  make_option(c("-R", "--target-region"), type="character", default="16S", 
              help="Amplicon target region. Options are either 16S, 18S or ITS \
              Default: 16S",
              metavar="16S"),
  
  make_option(c("-t", "--taxonomy-table"), type="character", default=NULL, 
              help="path to feature taxonomy table i.e. ASV taxonomy table.",
              metavar="path"),
  
  make_option(c("-g", "--group"), type="character", default="groups", 
              help="Column in metadata to be analyzed / compared.",
              metavar="groups"),
  
  make_option(c("-s", "--samples-column"), type="character", default="Sample Name", 
              help="Column in metadata containing sample names in the feature table. \
                    Deafault: 'Sample Name' ",
              metavar="Sample Name"),
  
  make_option(c("-o", "--output-prefix"), type="character", default="", 
              help="Unique name to tag onto output files. Default: empty string.",
              metavar=""),
  
  make_option(c("-r", "--remove-rare"), type="logical", default=FALSE, 
              help="Should rare features be filtered out?. If this flag is set \
                   then rare features will be filtered out. Default: FALSE. ",
              action= "store_true", metavar="FALSE"),
  
  make_option(c("-y", "--assay-suffix"), type="character", default="_GLAmpSeq", 
              help="Genelab assay suffix.", metavar="GLAmpSeq"),
  
  make_option(c("-p", "--prevalence-cutoff"), type="numeric", default=0.15, 
              help="If --remove-rare, a numerical fraction between 0 and 1. Taxa with prevalences
              (the proportion of samples in which the taxon is present) less 
              than --prevalence-cutoff will be excluded in the analysis. 
              Default is 0.15, i.e. exclude taxa / features that are not present
              in at least 15% of the samples.",
              metavar="0.15"),
  
  make_option(c("-l", "--library-cutoff"), type="numeric", default=100, 
              help="If --remove-rare, a numerical threshold for filtering samples based on library
              sizes. Samples with library sizes less than lib_cut will be 
              excluded in the analysis. Default is 100. 
              if you do not want to discard any sample then set to 0.",
              metavar="100"),
  
  
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
                      "\nA script to perform differential abundance using DESeq2.",
                      "\nIt outputs a differential abundance statistics table \
                      and volcano plots. ", sep="")
)


opt <- parse_args(opt_parser)


if (opt$version) {
  cat("run_deseq.R version: ", version, "\n")
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
  message("Differential abundance testing will be run on the default 'groups' column \n")
}

if(opt[["samples-column"]] == "Sample Name") {
  message("I will assume that the sample names are in a column named 'Sample Name' \n")
}

library(glue)
library(phyloseq)
library(DESeq2)
library(taxize)
library(ggrepel)
library(tidyverse)

# --------------------------- Functions -------------------------------------#

# Define the geometric mean function
gm_mean <- function(x, na.rm=TRUE) {
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}


remove_rare_features <- function(feature_table, cut_off_percent=3/4){
  
  # feature_table [MATRIX]  feature table matrix with samples as columns and 
  #                         features as rows
  # cut_off_percent [NUMERIC] cut-off fraction  or decimal between 0.001 to 1 
  #                          of the total number of samples to determine the 
  #                          most abundant features. By default it removes 
  #                          features that are not present in 3/4 of the total 
  #                          number of samples
  
  # Filter by occurrence in a fraction of samples
  # Define a cut-off for determining what's rare
  cut_off <- cut_off_percent * ncol(feature_table)
  # Get the occurrence for each feature
  feature_occurence <- rowSums(feature_table > 0)
  # Get names of the abundant features
  abund_features <- names(feature_occurence[feature_occurence >= cut_off])
  # Remove rare features
  abun_features.m <- feature_table[abund_features,]
  return(abun_features.m)
}


process_taxonomy <- function(taxonomy, prefix='\\w__') {
  #function to process a metaphlan2 taxonopmy assigment table
  #1. ~ file_path is a string specifying the taxonomic assignment file name
  #2 prefix ~ is a regular expression specifying the characters to remove
  # from the taxon names  '\\w__'  for greengenes and 'D_\\d__' for SILVA
  
  #taxon_levels <- c("kingdom","phylum","class","order",
  #                  "family","genus","species", "strain")
  
  taxonomy <- apply(X = taxonomy, MARGIN = 2, FUN = as.character) 
  
  #taxonomy[,'species'] <- paste(taxonomy[,'genus'],taxonomy[,'species'])
  # replace NAa with Other and delete the D_num__ prefix from the taxonomy names
  for (rank in colnames(taxonomy)) {
    #delete the taxonomy prefix
    taxonomy[,rank] <- gsub(pattern = prefix, x = taxonomy[, rank],
                            replacement = '')
    indices <- which(is.na(taxonomy[,rank]))
    taxonomy[indices, rank] <- rep(x = "Other", times=length(indices)) 
    #replace empty cell
    indices <- which(taxonomy[,rank] == "")
    taxonomy[indices,rank] <- rep(x = "Other", times=length(indices))
  }
  taxonomy <- apply(X = taxonomy,MARGIN = 2,
                    FUN =  gsub,pattern = "_",replacement = " ") %>% 
    as.data.frame(stringAsfactor=F)
  return(taxonomy)
}



# Function for format a taxonomy assignment table by appending suffix
# to a known name
format_taxonomy_table <- function(taxonomy=taxonomy.m,stringToReplace="Other",
                                  suffix=";Other") {
  
  for (taxa_index in seq_along(taxonomy)) {
    #indices <- which(taxonomy[,taxa_index] == stringToReplace)
    
    indices <- grep(x = taxonomy[,taxa_index], pattern = stringToReplace)
    
    taxonomy[indices,taxa_index] <- 
      paste0(taxonomy[indices,taxa_index-1],
             rep(x = suffix, times=length(indices)))
    
  }
  return(taxonomy)
}


fix_names<- function(taxonomy,stringToReplace,suffix){
  #1~ taxonomy is a taxonomy dataframe with taxonomy ranks as column names
  #2~ stringToReplace is a vector of regex strings specifying what to replace
  #3~ suffix is a string specifying the replacement value
  
  
  for(index in seq_along(stringToReplace)){
    taxonomy <- format_taxonomy_table(taxonomy = taxonomy,
                                      stringToReplace=stringToReplace[index], 
                                      suffix=suffix[index])
  }
  return(taxonomy)
}


taxize_options(ncbi_sleep = 0.8)
# A function to retrieve the NCBI taxonomy id for a given taxonomy name
get_ncbi_ids <- function(taxonomy, target_region){
  
  if(target_region == "ITS"){
    search_string <- "fungi"
  }else if(target_region == "18S"){
    search_string <- "eukaryote"
  }else{
    search_string <- "bacteria"
  }
  
  uid <- get_uid(taxonomy, division_filter = search_string)
  
  tax_ids <- uid[1:length(uid)]
  
  return(tax_ids)
  
}

# ------ Collecting the required input variables ---------- #

# Group in metadata to analyze
group <- opt[["group"]]  # "groups"
samples_column <- opt[["samples-column"]] # "Sample Name"
metadata_file <- opt[["metadata-table"]]
taxonomy_file <-  opt[["taxonomy-table"]]
feature_table_file <- opt[["feature-table"]]
feature <- opt[["feature-type"]]   # "ASV"
target_region <- opt[["target-region"]] # 16S
output_prefix <- opt[["output-prefix"]]
assay_suffix <- opt[["assay-suffix"]]
remove_rare <- opt[["remove-rare"]]
# taxon / ASV prevalence cutoff
prevalence_cutoff <- opt[["prevalence-cutoff"]] # 0.15 (15%)
# sample / library read count cutoff
library_cutoff <- opt[["library-cutoff"]]  # 100
diff_abund_out_dir <- "differential_abundance/deseq2/"
if(!dir.exists(diff_abund_out_dir)) dir.create(diff_abund_out_dir, recursive = TRUE)


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

if(remove_rare){
  
  # Remove samples with less than library-cutoff
  message(glue("Dropping samples with less than {library_cutoff} read counts"))
  feature_table <- feature_table[,colSums(feature_table) >= library_cutoff]
  # Remove rare ASVs
  message(glue("Dropping features with prevalence less than {prevalence_cutoff * 100}%"))
  feature_table <- remove_rare_features(feature_table,
                                        cut_off_percent = prevalence_cutoff)
}

# ------------------------ Read Taxonomy table ---------------------------- #
taxonomy <-  read_delim(file = taxonomy_file) %>% as.data.frame()
# Set the feature id column as the row names of the taxonomy table
# This assumes that the first column contains the feature ids e.g. ASV ID 
rownames(taxonomy) <- taxonomy[,1]
taxonomy_table  <- taxonomy[, -1]
feature_names <- rownames(taxonomy_table)
taxonomy_table  <- process_taxonomy(taxonomy_table)
rownames(taxonomy_table) <- feature_names

message(glue("There are {sum(taxonomy_table$phylum == 'Other')} features without 
           taxonomy assignments. Dropping them ..."))

# Dropping features that couldn't be assigned taxonomy
taxonomy_table <- taxonomy_table[-which(taxonomy_table$phylum == 'Other'),]

# Handle case where no domain was assigned but a phylum wasn't.
if(all(is.na(taxonomy$domain))){
  
  if(target_region == "ITS"){
    
    taxonomy_table$domain <- "Fungi"
    
  }else if(target_region == "18S"){
    
    taxonomy_table$domain <- "Eukaryotes"
    
  }else{
    
    taxonomy_table$domain <- "Bacteria"
  }
  
}

# Removing Chloroplast and Mitochondria Organelle DNA contamination
asvs2drop <- taxonomy_table %>%
  unite(col="taxonomy",domain:species) %>%
  filter(str_detect(taxonomy, "[Cc]hloroplast|[Mn]itochondria")) %>%
  row.names()
taxonomy_table <- taxonomy_table[!(rownames(taxonomy_table) %in% asvs2drop),]

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

#### pairwise comparisons
unique_groups <- unique(metadata[[group]])

# Create phyloseq object
ASV_physeq <- phyloseq(otu_table(feature_table, taxa_are_rows = TRUE),
                       tax_table(as.matrix(taxonomy_table)),
                       sample_data(metadata))

deseq_obj <- phyloseq_to_deseq2(physeq = ASV_physeq,
                                design = reformulate(group))

# Add pseudocount if any 0 count samples are present
if (sum(colSums(counts(deseq_obj)) == 0) > 0) {
  count_data <- counts(deseq_obj) + 1 
  
  count_data <- as.matrix(apply(count_data, 2, as.integer))
  rownames(count_data) <- rownames(counts(deseq_obj))
  colnames(count_data) <- colnames(counts(deseq_obj))
  counts(deseq_obj) <- count_data
}

# Run Deseq
# https://rdrr.io/bioc/phyloseq/src/inst/doc/phyloseq-mixture-models.R 
deseq_modeled <- tryCatch({
  # Attempt to run DESeq
  DESeq(deseq_obj)
}, error = function(e) {
  message("Error encountered in DESeq, applying alternative \
          method for size factor estimation...")
  
  geoMeans <- apply(counts(deseq_obj), 1, gm_mean)
  
  # Apply the alternative size factor estimation method
  deseq_obj <- estimateSizeFactors(deseq_obj, geoMeans=geoMeans)
  
  # Call DESeq again with alternative geom mean size est
  DESeq(deseq_obj)
})


# Get unique group comparison as a matrix
pairwise_comp.m <- utils::combn(metadata[,group] %>% unique, 2)
pairwise_comp_df <- pairwise_comp.m %>% as.data.frame 

colnames(pairwise_comp_df) <- map_chr(pairwise_comp_df,
                                      \(col) str_c(col, collapse = "v"))
comparisons <- colnames(pairwise_comp_df)
names(comparisons) <- comparisons

# Retrieve statistics table
merged_stats_df <-  data.frame(ASV=rownames(feature_table))
colnames(merged_stats_df) <- feature

walk(pairwise_comp_df, function(col){
  
  group1 <- col[1]
  group2 <- col[2]
  
df <- results(deseq_modeled, contrast = c(group, group1, group2)) %>%
  data.frame() %>%
  rownames_to_column(feature) %>% 
  set_names(c(feature ,
              glue("baseMean_({group1})v({group2})"),
              glue("log2FC_({group1})v({group2})"),
              glue("lfcSE_({group1})v({group2})"), 
              glue("stat_({group1})v({group2})"), 
              glue("pvalue_({group1})v({group2})"),
              glue("padj_({group1})v({group2})") 
            ))

            
  merged_stats_df <<- merged_stats_df %>% 
                          dplyr::left_join(df, join_by("ASV"))
})



# Add NCBI id to feature i.e. ASV
tax_names <- map_chr(str_replace_all(taxonomy_table$species, ";_","")  %>%
                       str_split(";"),
                     function(row) row[length(row)])

df <- data.frame(ASV=rownames(taxonomy_table), best_taxonomy=tax_names)

# Pull NCBI IDS for unique taxonomy names
df2 <- data.frame(best_taxonomy = df$best_taxonomy %>%
                    unique()) %>%
  mutate(NCBI_id=get_ncbi_ids(best_taxonomy, target_region),
         .after = best_taxonomy)

df <- df %>%
  left_join(df2, join_by("best_taxonomy")) %>% 
  right_join(merged_stats_df)




group_levels <- metadata[, group] %>% unique() %>% sort()
normalized_table <- counts(deseq_modeled, normalized=TRUE) %>% 
                        as.data.frame() %>%
                        rownames_to_column(feature)


samples <- metadata[[samples_column]]
samplesdropped <- setdiff(x = samples, y = colnames(normalized_table)[-1])
missing_df <- data.frame(ASV=normalized_table[[feature]],
                         matrix(data = NA, 
                                nrow = nrow(normalized_table),
                                ncol = length(samplesdropped)
                         )
)
colnames(missing_df) <- c(feature,samplesdropped)


group_means_df <- normalized_table[feature]
walk(group_levels, function(group_level){
  
  
  mean_col <- glue("Group.Mean_({group_level})")
  std_col <- glue("Group.Stdev_({group_level})")
  
  # Samples that belong to the current group
  Samples <- metadata %>%
    filter(!!sym(group) == group_level) %>%
    pull(!!sym(samples_column))
  # Samples that belong to the current group that are in the normalized table
  Samples <- intersect(colnames(normalized_table), Samples)
  
  temp_df <- normalized_table %>% select(!!feature, all_of(Samples)) %>% 
    rowwise() %>%
    mutate(!!mean_col := mean(c_across(where(is.numeric))),
           !!std_col := sd(c_across(where(is.numeric))) ) %>% 
    select(!!feature,!!sym(mean_col), !!sym(std_col))
  
  group_means_df <<- group_means_df %>% left_join(temp_df)
  
})


# Append Mean and standard deviation
normalized_table <- normalized_table %>%
  rowwise() %>%
  mutate(All.Mean=mean(c_across(where(is.numeric))),
         All.Stdev=sd(c_across(where(is.numeric))) )%>% 
  left_join(missing_df, by = feature) %>% 
  select(!!feature, all_of(samples), All.Mean, All.Stdev)


# Add taxonomy
merged_df <- df  %>%
  left_join(taxonomy_table %>%
              as.data.frame() %>%
              rownames_to_column(feature)) %>% 
  select(!!feature, domain:species,everything()) # Try to generalize

# Merge all prepared tables 
merged_df <- merged_df %>%
  select(!!sym(feature):NCBI_id) %>%
  left_join(normalized_table, by = feature) %>%
  left_join(merged_df) %>% 
  left_join(group_means_df, by = feature) %>% 
  mutate(across(where(is.numeric), ~round(.x, digits=3))) %>% 
  mutate(across(where(is.matrix), as.numeric))


output_file <- glue("{diff_abund_out_dir}/{output_prefix}deseq2_differential_abundance{assay_suffix}.csv")
message("Writing out results of differential abundance using DESeq2...")
write_csv(merged_df,output_file)



# Make volcano plots
walk(pairwise_comp_df, function(col){
  
  group1 <- col[1]
  group2 <- col[2]
  
  plot_width_inches <- 11.1
  plot_height_inches <- 8.33
  p_val <- 0.1 #also logfc cutoff?
  
  deseq_res <- results(deseq_modeled, contrast = c(group, group1, group2))
  volcano_data <- as.data.frame(deseq_res)
  
  
  volcano_data <- volcano_data[!is.na(volcano_data$padj), ]
  volcano_data$significant <- volcano_data$padj <= p_val #also logfc cutoff?
  
  ######Long x-axis label adjustments##########
  x_label <- paste("Log2 Fold Change\n(",group1," vs ",group2,")")
  label_length <- nchar(x_label)
  max_allowed_label_length <- plot_width_inches * 10
  
  # Construct x-axis label with new line breaks if was too long
  if (label_length > max_allowed_label_length){
    x_label <- paste("Log2 Fold Change\n\n(", group1, "\n vs \n", group2, ")", sep="")
  }
  #######################################
  
  # ASVs promoted in space on right, reduced on left
  p <- ggplot(volcano_data, aes(x=log2FoldChange, y=-log10(padj), 
                                color=significant)) +
    geom_point(alpha=0.7, size=2) +
    geom_hline(yintercept = -log10(p_val), linetype = "dashed") +
    scale_color_manual(values=c("black", "red"), 
                       labels=c(paste0("padj > ", p_val), 
                                paste0(" padj \u2264 ", p_val))) +
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
  
  volcano_plot <- p + geom_text_repel(data=top_points, 
                                      aes(label=row.names(top_points)),
                                      size=3)
  
  
  
  
  ggsave(filename=glue("{diff_abund_out_dir}{output_prefix}volcano_{group1}_vs_{group2}.png"),
         plot=volcano_plot,
         width = plot_width_inches, 
         height = plot_height_inches, 
         dpi = 300)
})
