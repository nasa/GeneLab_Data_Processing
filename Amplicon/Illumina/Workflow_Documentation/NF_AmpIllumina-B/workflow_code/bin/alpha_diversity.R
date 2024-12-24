#!/usr/bin/env Rscript


###############################################################################
# AUTHOR : OLABIYI ADEREMI OBAYOMI
# DESCRIPTION: A script to perform alpha diversity analysis.
# E-mail: obadbotanist@yahoo.com
# Created: November 2024
# example: Rscript alpha_diversity.R \
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
  
  make_option(c("-o", "--output-prefix"), type="character", default="", 
              help="Unique name to tag onto output files. Default: empty string.",
              metavar=""),
  
  make_option(c("-d", "--rarefaction-depth"), type="numeric", default=500, 
              help="Minimum rarefaction depth for alpha diversity estimation. \
                   Default: 500. ",
              metavar="500"),

  make_option(c("-r", "--remove-rare"), type="logical", default=FALSE, 
              help="Should rare features be filtered out?. \
                   Default: FALSE. ", action= "store_true",
              metavar="FALSE"),  
    
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
  
  make_option(c("-e", "--legend-title"), type="character", default="Groups", 
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
                      "\n A script to perform ASV alpha diversity analysis.",
                      "\nIt outputs sample and group alpha diversity plots ",
                      "along with group diversity and group statistics table",
                      sep="")
)


opt <- parse_args(opt_parser)

# print(opt)
# stop()


if (opt$version) {
  cat("alpha_diversity.R version: ", version, "\n")
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
library(FSA)
library(multcompView)
library(rstatix)
library(patchwork)
library(RColorBrewer)
library(tidyverse)

# ----------------------------- Functions ----------------------------------- #

calculate_text_size <- function(num_samples, start_samples = 25, min_size = 3) {
  max_size = 11  # Maximum size for up to start_samples
  slope = -0.15
  
  if (num_samples <= start_samples) {
    return(max_size)
  } else {
    # Calculate the current size with the hard coded slope
    current_size = max_size + slope * (num_samples - start_samples)
    
    # Ensure the size doesn't go below the minimum
    return(max(current_size, min_size))
  }
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


# Required variables
metadata_file <- opt[["metadata-table"]] 
features_file <- opt[["feature-table"]]
taxonomy_file <- opt[["taxonomy-table"]]
# Metadata group column name to compare
groups_colname <- opt[["group"]]
sample_colname  <- opt[["samples-column"]]
output_prefix <- opt[["output-prefix"]]
assay_suffix <- opt[["assay-suffix"]]
legend_title <- opt[["legend-title"]]
rarefaction_depth  <- opt[["rarefaction-depth"]]
remove_rare <- opt[["remove-rare"]]
# taxon / ASV prevalence cutoff
prevalence_cutoff <- opt[["prevalence-cutoff"]] # 0.15 (15%)
# sample / library read count cutoff
library_cutoff <- opt[["library-cutoff"]]  # 100
alpha_diversity_out_dir <-"alpha_diversity/"
if(!dir.exists(alpha_diversity_out_dir)) dir.create(alpha_diversity_out_dir)



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
metadata <- metadata %>%
  mutate(color = map_chr(!!sym(groups_colname),
                         function(group) { group_colors[group] }
                         ) 
        )
sample_names <- rownames(metadata)


sample_info_tab <- metadata %>%
  select(!!groups_colname, color) %>%
  arrange(!!sym(groups_colname))


values <- sample_info_tab %>% pull(color) %>% unique()



# Feature or ASV table
feature_table <- read.table(file = features_file, header = TRUE,
                            row.names = 1, sep = "\t")

if(remove_rare){
  
  # Remove samples with less than library-cutoff
  message(glue("Dropping samples with less than {library_cutoff} read counts"))
  feature_table <- feature_table[,colSums(feature_table) >= library_cutoff]
  # Remove rare ASVs
  message(glue("Dropping features with prevalence less than {prevalence_cutoff * 100}%"))
  feature_table <- remove_rare_features(feature_table,
                                        cut_off_percent = prevalence_cutoff)
}


# Taxonomy 
taxonomy_table <-  read.table(file = taxonomy_file, header = TRUE,
                              row.names = 1, sep = "\t")



# Preprocess ASV and taxonomy tables

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



# Create phyloseq object
ASV_physeq <- phyloseq(otu_table(feature_table, taxa_are_rows = TRUE),
                       tax_table(as.matrix(taxonomy_table)), 
                       sample_data(sample_info_tab))

seq_per_sample <- colSums(feature_table) %>% sort()

# minimum value
depth <- min(seq_per_sample)

for (count in seq_per_sample) {
  
  if(count >= rarefaction_depth) {
    depth <- count
    break
    }
  
}

#----- Rarefy sample counts to even depth per sample
ps.rarefied <- rarefy_even_depth(physeq = ASV_physeq, 
                                 sample.size = depth,
                                 rngseed = 1, 
                                 replace = FALSE, 
                                 verbose = FALSE)


# ------------------- Rarefaction curve
# Calculate a rough estimate of the step sample step size for plotting.
# This is meant to keep plotting time constant regardless of sample depth
step <-  (50*depth)/1000

p <- rarecurve(t(otu_table(ps.rarefied)) %>% as.data.frame(),
               step = step,
               col = sample_info_tab[["color"]], 
               lwd = 2, ylab = "ASVs", cex=0.5,
               label = FALSE, tidy = TRUE)


sample_info_tab_names <-  sample_info_tab %>% rownames_to_column("Site")

p <- p %>% left_join(sample_info_tab_names, by = "Site")



# Sample rarefaction curves

rareplot <- ggplot(p, aes(x = Sample, y = Species,
                          group = Site, color = !!sym(groups_colname))) + 
  geom_line() + 
  scale_color_manual(values = group_colors) +
  labs(x = "Number of Sequences", y = "Number of ASVs", color = legend_title) + 
  theme_bw() +
  theme(legend.position = "right",
        text = element_text(face = 'bold', size = 15),
        legend.text = element_text(face = 'bold', size = 14),
        legend.direction = "vertical",
        legend.justification = "center",
        legend.box.just = "center",
        legend.title = element_text(size = 15, face='bold'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(t = 10, r = 20, b = 10, l = 10, unit = "pt"))

ggsave(filename = glue("{alpha_diversity_out_dir}/{output_prefix}rarefaction_curves{assay_suffix}.png"),
       plot=rareplot, width = 14, height = 8.33, dpi = 300)

# ------------------  Richness and diversity estimates  ------------------#

# Statistics table
diversity_metrics <- c("Observed", "Chao1", "Shannon", "Simpson")
names(diversity_metrics) <- diversity_metrics
diversity.df <- estimate_richness(ps.rarefied, 
                                  measures = diversity_metrics) %>% 
               select(-se.chao1) %>%
               rownames_to_column("samples")
               

merged_table  <- metadata %>%
  rownames_to_column("samples") %>%
  inner_join(diversity.df)

diversity_stats <- map_dfr(.x = diversity_metrics, function(metric){
  
  res <- dunnTest(merged_table[,metric],merged_table[,groups_colname])
  
  df <- res$res %>%
    separate(col = Comparison, into = c("group1", "group2"), sep = " - ") %>% 
    mutate(Metric=metric)  %>% 
    rename(p=P.unadj, p.adj=P.adj) %>% 
    mutate(p.format=round(p,digits = 2))
  
  add_significance(df, p.col='p', output.col = 'p.signif') %>% 
    select(Metric,group1, group2, Z, p, p.adj, p.format, p.signif) %>% 
    mutate(across(where(is.numeric), ~round(.x, digits = 2)))
  
})



write_csv(x = diversity_stats, 
            file = glue("{alpha_diversity_out_dir}/{output_prefix}statistics_table{assay_suffix}.csv"))


comp_letters <- data.frame(group = group_levels)
colnames(comp_letters) <- groups_colname

walk(.x = diversity_metrics, function(metric=.x){
  
  sub_comp <- diversity_stats %>% filter(Metric == metric)
  p_values <- sub_comp$p
  names(p_values) <- paste(sub_comp$group1,sub_comp$group2, sep = "-")
  
  letters_df <-  enframe(multcompView::multcompLetters(p_values)$Letters,
                      name = groups_colname,
                      value = glue("{metric}_letter"))
  comp_letters <<- comp_letters %>% left_join(letters_df)
   
})


# Summary table

diversity_table <- metadata %>%
  rownames_to_column("samples") %>%
  inner_join(diversity.df) %>%
  group_by(!!sym(groups_colname))  %>% 
  summarise(N=n(), across(Observed:Simpson,
                   .fns=list(mean=mean, se=se),
                   .names = "{.col}_{.fn}")) %>% 
  mutate(across(where(is.numeric), ~round(.x, digits = 2))) %>%
  left_join(comp_letters) %>% 
  mutate(Observed = glue("{Observed_mean} ± {Observed_se}{Observed_letter}"),
         Chao1 = glue("{Chao1_mean} ± {Chao1_se}{Chao1_letter}"),
         Shannon = glue("{Shannon_mean} ± {Shannon_se}{Shannon_letter}"),
         Simpson = glue("{Simpson_mean} ± {Simpson_se}{Simpson_letter}")
         ) %>%
  select (-contains("_"))

write_csv(x = diversity_table, 
            file = glue("{alpha_diversity_out_dir}/{output_prefix}summary_table{assay_suffix}.csv"))


# ------------------ Make richness by sample dot plots ---------------------- #

number_of_samples <- length(rownames(sample_info_tab))
richness_sample_label_size <- calculate_text_size(number_of_samples)
metrics2plot <- c("Observed", "Shannon")
names(metrics2plot) <- metrics2plot

samples_order <- metadata %>% arrange(!!sym(groups_colname)) %>% rownames()

richness_by_sample <- plot_richness(ps.rarefied, color = groups_colname,
                                    measures = metrics2plot)

richness_by_sample <-  ggplot(richness_by_sample$data %>% 
                                mutate(samples = factor(samples, 
                                                        levels=samples_order)),
                              aes(x=samples, y=value, colour = !!sym(groups_colname))) +
  geom_point() +
  geom_errorbar(aes(ymin=value-se, ymax = value+se),
                width=0.2, position=position_dodge(0.9)) +
  facet_wrap(~variable, scales = "free_y") +
  scale_color_manual(values = group_colors) + 
  theme_bw() +labs(x = NULL, color = legend_title, y="Alpha Diversity Measure") +
  theme(
    text = element_text(face = 'bold', size = 15),
    legend.text = element_text(face = 'bold', size = 14),
    legend.position = "bottom",
    legend.direction = "vertical",
    legend.justification = "center",
    legend.box.just = "center",
    legend.title = element_text(face = 'bold', size = 15, hjust = 0.09),
    axis.text.x = element_text(angle = 90,
                               size = richness_sample_label_size,
                               vjust = 0.5,  # Vertically center the text
                               hjust = 1),
    axis.ticks.length=unit(-0.15, "cm"),
    strip.text =  element_text(size = 14,face ='bold')
  )


ggsave(filename = glue("{alpha_diversity_out_dir}/{output_prefix}richness_and_diversity_estimates_by_sample{assay_suffix}.png"),
       plot=richness_by_sample, width = 14, height = 8.33, dpi = 300, units = "in")


# ------------------- Make richness by group box plots ----------------------- #


richness_by_group <- plot_richness(ps.rarefied, x = groups_colname, 
                                   color = groups_colname,
                                   measures = metrics2plot)


p <- map(.x = metrics2plot, .f = function(metric){
  
  
  p <-  ggplot(richness_by_group$data %>%  filter(variable == metric), 
              aes(x=!!sym(groups_colname), y=value, fill=!!sym(groups_colname)) 
              ) +
    geom_point() +
    geom_boxplot() +
    scale_fill_manual(values = group_colors) + 
    theme_bw() + labs(fill = legend_title, x = NULL, y= metric) +
    theme(
      text = element_text(size = 15, face = 'bold'),
      legend.text = element_text(face = 'bold', size = 14),
      legend.position = "right",
      legend.direction = "vertical",
      legend.justification = "center",
      legend.box.just = "center",
      legend.title = element_text(face = 'bold', size = 15),
      axis.text.x = element_blank(),
      axis.ticks.length=unit(-0.15, "cm"),
      strip.text =  element_text(size = 14,face ='bold')
    ) 
  
  
  summary_table <- p$data %>%
    select(!!sym(groups_colname), value) %>% 
    group_by(!!sym(groups_colname)) %>%
    summarise(max=max(value), range=max(value)-min(value)) %>%
    left_join(comp_letters %>%
                select(groups, label= !!sym( glue("{metric}_letter") ) 
                       ) 
              )
  text_size <- 6
  
  # Calculate a constant to add to the max value of each group
  # to determine where each group text will be added 
  toAdd <- if_else(condition = max(summary_table$range) <= 5,
                   true = min(summary_table$range),
                   false = (median(summary_table$range) - min(summary_table$range)) / 20
                    )
  
  # Add text to plot
  p + geom_text(data=summary_table,
                                    mapping = aes(y=max+toAdd, 
                                                  label=label,
                                                  fontface = "bold"),
                                  size = text_size)


})

richness_by_group <- wrap_plots(p, ncol = 2, guides =  'collect')

width <- 3.6 * length(group_levels)
ggsave(filename = glue("{alpha_diversity_out_dir}/{output_prefix}richness_and_diversity_estimates_by_group{assay_suffix}.png"),
       plot=richness_by_group, width = width, height = 8.33, dpi = 300, units = "in")

