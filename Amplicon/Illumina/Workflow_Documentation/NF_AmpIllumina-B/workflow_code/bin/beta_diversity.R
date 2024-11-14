#!/usr/bin/env Rscript


###############################################################################
# AUTHOR : OLABIYI ADEREMI OBAYOMI
# DESCRIPTION: A script to perform beta diversity analysis.
# E-mail: obadbotanist@yahoo.com
# Created: November 2024
# example: Rscript beta_diversity.R \
#                  --metadata-table 'mapping/GLDS-487_amplicon_v1_runsheet.csv' \
#                  --feature-table 'data/counts_GLAmpSeq.tsv' \
#                  --taxonomy-table 'data/taxonomy_GLAmpSeq.tsv' \
#                  --group 'groups' \
#                  --samples-column 'Sample Name' \
#                  --rarefaction-depth 500
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
  
  make_option(c("-p", "--output-prefix"), type="character", default="", 
              help="Unique name to tag onto output files. Default: empty string.",
              metavar=""),
  
  make_option(c("-d", "--rarefaction-depth"), type="numeric", default=500, 
              help="Minimum rarefaction depth for alpha diversity estimation. \
                   Default: 500. ",
              metavar="500"),
  
  make_option(c("-c", "--abundance-cutoff"), type="numeric", default=0.2, 
              help="A fraction defining how abundant features most be to be \
                   analyzes. Default: 1/5. ",
              metavar="0.2"),
  make_option(c("-r", "--remove-rare"), type="logical", default=FALSE, 
              help="Should rare features be filtered out?. \
                   Default: FALSE. ", action= "store_true",
              metavar="FALSE"),
  
  make_option(c("-l", "--legend-title"), type="character", default="Groups", 
              help="Legend title for alpha diversity plots.",
              metavar="Groups"),
  
  make_option(c("-a", "--assay-suffix"), type="character", default="_GLAmpSeq", 
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
                  --samples-column 'Sample Name' \\
                  --rarefaction-depth 500",
  description = paste("Author: Olabiyi Aderemi Obayomi",
                      "\nEmail: olabiyi.a.obayomi@nasa.gov",
                      "\n A script to perform ASV beta diversity analysis.",
                      "\nIt outputs a dendogram, pcoa and statistics tables. ",
                      sep="")
)


opt <- parse_args(opt_parser)


if (opt$version) {
  cat("beta_diversity.R version: ", version, "\n")
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


library(vegan)
library(phyloseq)
library(here)
library(glue)
library(DESeq2)
library(ggdendro)
library(RColorBrewer)
library(broom)
library(tidyverse)



# ----------------------------- Functions ----------------------------------- #

plot_pcoa <- function(phy, pcoa, addtext=FALSE) {
  
  p <-  plot_ordination(phy, pcoa,
                        color = groups_colname) + 
    geom_point(size = 1) 
  
  if(addtext){
    p <- p + geom_text(aes(label = sample_names), 
                       show.legend = FALSE,
                       hjust = 0.3, vjust = -0.4, size = 4)
  }
  p + 
    labs( 
      x = label_PC1,
      y = label_PC2,
      col = legend_title
    ) +
    coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + 
    scale_color_manual(values = group_colors) +
    theme_bw() + theme(text = element_text(size = 15, face="bold"),
                       legend.direction = "vertical",
                       legend.justification = "center",
                       legend.title = element_text(hjust=0.1)) +
    annotate("text", x = Inf, y = -Inf, 
             label = paste("R2:", toString(round(r2_value, 3))),
             hjust = 1.1, vjust = -2, size = 4)+
    annotate("text", x = Inf, y = -Inf, 
             label = paste("Pr(>F)", toString(round(prf_value,4))),
             hjust = 1.1, vjust = -0.5, size = 4) + ggtitle("PCoA")
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

custom_palette <- c("#1F78B4","#33A02C","#FB9A99","#E31A1C","#6A3D9A",
                    "#FDBF6F", "#FF7F00","#CAB2D6","#FF00FFFF", "#B15928",
                    "#000000","#FFC0CBFF", "#A6CEE3", "#8B864EFF","#F0027F",
                    "#666666","#1B9E77", "#E6AB02","#A6761D","#FFFF00FF",
                    "#00FFFFFF", "#FFFF99", "#B2182B","#FDDBC7","#D1E5F0",
                    "#B2DF8A","#CC0033","#FF00CC","#330033", "#999933",
                    "#FF9933", "#FFFAFAFF",colors()) 

# Remove white colors
pattern_to_filter <- "white|snow|azure|gray|#FFFAFAFF|aliceblue"
custom_palette <- custom_palette[-c(21:23, grep(pattern = pattern_to_filter,
                                                x = custom_palette, 
                                                ignore.case = TRUE))]


###############################################################################
# Things to do:
# 1. add of options to choose normalization method c(vst, relative_abundance, Archison distance)
# 2. Add rarefaction

# Required variables
metadata_file <- here(opt[["metadata-table"]])
features_file <-  here(opt[["feature-table"]]) 
taxonomy_file <-  here(opt[["taxonomy-table"]])
beta_diversity_out_dir <- here("beta_diversity/")
if(!dir.exists(beta_diversity_out_dir)) dir.create(beta_diversity_out_dir)
# Metadata group column name to compare
groups_colname <- opt[["group"]]
sample_colname  <- opt[["samples-column"]]
output_prefix <- opt[["output-prefix"]]
assay_suffix <- opt[["assay-suffix"]]
legend_title <- opt[["legend-title"]] 
remove_rare <- opt[["remove-rare"]]
abundance_cutoff <- opt[["abundance-cutoff"]]


# Read in processed data
metadata <- read_csv(file = metadata_file) %>% as.data.frame()
row.names(metadata) <- metadata[[sample_colname]]
metadata[,sample_colname] <- NULL
group_column_values <- metadata %>% pull(!!sym(groups_colname))
group_levels <- unique(group_column_values)

# Add colors to metadata equals to the number of levels
# in the factor groups column
num_colors <- length(group_levels)
palette <- 'Set1'
number_of_colors_in_palette <- 9
if(num_colors <= number_of_colors_in_palette){
  
  colors <- RColorBrewer::brewer.pal(n = num_colors, name = palette)
  
}else{
  
  colors <- custom_palette[1:num_colors]
}

# Metadata
group_colors <- setNames(colors, group_levels)


# --------------- Process metadata  -------------------------- #
metadata <- metadata %>%
  mutate(color = map_chr(!!sym(groups_colname),
                         function(group) { group_colors[group] }
  ) 
  )
sample_names <- rownames(metadata)
deseq2_sample_names <- make.names(sample_names, unique = TRUE)

short_group_labels <- sprintf("%d: %s", seq_along(group_levels), group_levels)
names(short_group_labels) <- group_levels

sample_info_tab <- metadata %>%
  select(!!groups_colname, color) %>%
  arrange(!!sym(groups_colname))

# Feature or ASV table
feature_table <- read.table(file = features_file, header = TRUE,
                            row.names = 1, sep = "\t")
if(remove_rare){
  # Remove rare ASVs
  feature_table <- remove_rare_features(feature_table,
                                        cut_off_percent=abundance_cutoff)
}


# Preprocess ASV and taxonomy tables

# Taxonomy 
taxonomy_table <-  read.table(file = taxonomy_file, header = TRUE,
                              row.names = 1, sep = "\t")

message(glue("There are {sum(is.na(taxonomy_table$domain))} features without 
           taxonomy assignments. Dropping them..."))


# Dropping features that couldn't be assigned taxonomy
taxonomy_table <- taxonomy_table[-which(is.na(taxonomy_table$domain)),]

# Removing Chloroplast and Mitochondria Organelle DNA contamination
asvs2drop <- taxonomy_table %>%
  unite(col="taxonomy",domain:species) %>%
  filter(str_detect(taxonomy, "[Cc]hloroplast|[Mn]itochondria")) %>%
  row.names()
taxonomy_table <- taxonomy_table[!(rownames(taxonomy_table) %in% asvs2drop),]

# Subset tables 

# Get features common to the taxonomy and feature table 
common_ids <- intersect(rownames(feature_table), rownames(taxonomy_table))

# Subset the feature and taxonomy tables to contain 
# only features found in both table
feature_table <- feature_table[common_ids,]
taxonomy_table <- taxonomy_table[common_ids,]



# Using deseq


# Keep only ASVs with at least 1 count
feature_table <- feature_table[rowSums(feature_table) > 0, ]
# Add +1 pseudocount for VST for vst transformation
feature_table <- feature_table + 1

# Make the order of samples in metadata match the order in feature table
metadata <- metadata[colnames(feature_table),]

# Create VST normalized counts matrix
# ~1 means no design
deseq_counts <- DESeqDataSetFromMatrix(countData = feature_table, 
                                       colData = metadata, 
                                       design = ~1)
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
vst_trans_count_tab <- assay(deseq_counts_vst)




# -----------    Hierarchical Clustering and dendogram plotting
euc_dist <- dist(t(vst_trans_count_tab))
euc_clust <- hclust(d = euc_dist, method = "ward.D2")

# Extract clustering data
hcdata <- dendro_data(euc_clust, type = "rectangle")
segment_data <- segment(hcdata)
label_data <- label(hcdata) %>%
                       left_join(sample_info_tab %>% 
                                    rownames_to_column("label"))

dendogram <- ggplot() +
  geom_segment(data = segment_data, 
               aes(x = x, y = y, xend = xend, yend = yend)
  ) +
  geom_text(data = label_data , 
            aes(x = x, y = y, label = label, 
                color = !!sym(groups_colname) , hjust = 0), 
            size = 4.5, key_glyph = "rect") +
  scale_color_manual(values = group_colors) +
  coord_flip() +
  scale_y_reverse(expand = c(0.2, 0)) +
  labs(color = legend_title) +
  theme_dendro() +
  guides(colour = guide_legend(override.aes = list(size = 5)))+
  theme(legend.key=element_rect(fill=NA),
        text = element_text(face = 'bold'),
        legend.title = element_text(size = 12, face='bold'),
        legend.text = element_text(face = 'bold', size = 11))

ggsave(filename = glue("{beta_diversity_out_dir}/{output_prefix}dendrogram_by_group{assay_suffix}.png"),
       plot = dendogram, width = 14, height = 10, dpi = 300, units = "in")


##### Making PCoA

# Making a phyloseq object with our transformed table
vst_count_phy <- otu_table(object = vst_trans_count_tab, taxa_are_rows = TRUE)
sample_info_tab_phy <- sample_data(sample_info_tab)
vst_physeq <- phyloseq(vst_count_phy, sample_info_tab_phy)

# generating a PCoA with phyloseq
vst_pcoa <- ordinate(physeq = vst_physeq, 
                     method = "PCoA", 
                     distance = "euclidean")
eigen_vals <- vst_pcoa$values$Eigenvalues

# Calculate the percentage of variance
percent_variance <- eigen_vals / sum(eigen_vals) * 100

message("Checking homogeneity of variance across groups")
variance_test <- betadisper(d = euc_dist, group = sample_info_tab[[groups_colname]]) %>%
  anova() %>%
  broom::tidy() %>% 
  mutate(across(where(is.numeric), ~round(.x, digits = 2)))

write_csv(x = variance_test, 
            file = glue("{beta_diversity_out_dir}/{output_prefix}variance_table{assay_suffix}.csv"))



adonis_res <- adonis2(formula = euc_dist ~ sample_info_tab[[groups_colname]])
adonis_test <- adonis_res %>%
                   broom::tidy() %>% 
                   mutate(across(where(is.numeric), ~round(.x, digits = 2)))
write_csv(x = adonis_test, 
            file = glue("{beta_diversity_out_dir}/{output_prefix}adonis_table{assay_suffix}.csv"))


# Retrieving plot labels
r2_value <- adonis_res$R2[1]
prf_value <- adonis_res$`Pr(>F)`[1]
label_PC1 <- sprintf("PC1 [%.1f%%]", percent_variance[1])
label_PC2 <- sprintf("PC2 [%.1f%%]", percent_variance[2])



# Unlabeled PCoA plot
ordination_plot_u <- plot_pcoa(vst_physeq, vst_pcoa) 
ggsave(filename=glue("{beta_diversity_out_dir}/{output_prefix}PCoA_without_labels{assay_suffix}.png"),
       plot=ordination_plot_u, width = 14, height = 8.33, dpi = 300, units = "in")

# Labeled PCoA plot
ordination_plot <- plot_pcoa(vst_physeq, vst_pcoa, addtext=TRUE) 
ggsave(filename=glue("{beta_diversity_out_dir}/{output_prefix}PCoA_w_labels{assay_suffix}.png"), 
       plot=ordination_plot, width = 14, height = 8.33, dpi = 300, units = "in")
