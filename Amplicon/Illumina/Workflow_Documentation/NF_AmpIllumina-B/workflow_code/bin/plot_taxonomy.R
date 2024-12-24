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
                  --samples-column 'Sample Name' ",
  description = paste("Author: Olabiyi Aderemi Obayomi",
                      "\nEmail: olabiyi.a.obayomi@nasa.gov",
                      "\nA script to generate taxonomy plots at different taxonomy levels.",
                      "\nIt outputs sample and group taxonomy plots ", sep="")
)




opt <- parse_args(opt_parser)


if (opt$version) {
  cat("plot_taxonomy.R version: ", version, "\n")
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
  message("Taxonomy plots will be grouped by the default 'groups' column \n")
}

if(opt[["samples-column"]] == "Sample Name") {
  message("I will assume that the sample names are in a column named 'Sample Name' \n")
}



library(glue)
library(tools)
library(here)
library(tidyverse)


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
  # Function to process a taxonopmy assignment table
  #1. ~ taxonomy is a string specifying the taxonomic assignment file name
  #2 prefix ~ is a regular expression specifying the characters to remove
  # from the taxon names  '\\w__'  for greengenes and 'D_\\d__' for SILVA
  
  
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

# Function to format a taxonomy assignment table by appending a suffix
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

# A function to generate taxon level count matrix based on a taxonomy table and
# an existing feature table
make_feature_table <- function(count_matrix,taxonomy,
                               taxon_level, samples2keep=NULL){

  # EAMPLE:
  # make_feature_table(count_matrix = feature_counts_matrix,
  #                    taxonomy = taxonomy_table, taxon_level = "Phylum")

  feature_counts_df <- data.frame(taxon_level=taxonomy[,taxon_level],
                                  count_matrix, check.names = FALSE,
                                  stringsAsFactors = FALSE)

  feature_counts_df <- aggregate(.~taxon_level,data = feature_counts_df,
                                 FUN = sum)
  rownames(feature_counts_df) <- feature_counts_df[,"taxon_level"]
  feature_table <- feature_counts_df[,-1]
  # Retain only taxa found in at least one sample
  taxa2keep <- rowSums(feature_table) > 0
  feature_table <- feature_table[taxa2keep,]

  if(!is.null(samples2keep)){
    feature_table <- feature_table[,samples2keep]
    # Retain only taxa found in at least one sample
    taxa2keep <- rowSums(feature_table) > 0
    feature_table <- feature_table[taxa2keep,]
  }

  return(feature_table)
}


# Function to group rare taxa or return a table with the rare taxa
group_low_abund_taxa <- function(abund_table, threshold=0.05,
                                 rare_taxa=FALSE) {
  # abund_table is a relative abundance matrix with taxa as columns and  samples as rows
  #rare_taxa is a boolean specifying if only rare taxa should be returned
  #If set to TRU then a table with only the rare taxa will be returned 
  #intialize an empty vector that will contain the indices for the
  #low abundance columns/ taxa to group
  taxa_to_group <- c()
  #intialize the index variable of species with low abundance (taxa/columns)
  index <- 1
  
  #loop over every column or taxa check to see if the max abundance is less than the set threshold
  #if true save the index in the taxa_to_group vector variable
  while(TRUE){
    
    for (column in seq_along(abund_table)){
      if(max(abund_table[,column]) < threshold ){
        taxa_to_group[index] <- column
        index = index + 1
      }
    }
    if(is.null(taxa_to_group)){
      threshold <- readline("please provide a higher threshold for grouping rare taxa, only numbers are allowed   ")
      threshold <- as.numeric(threshold)
    }else{
      break
    }
    
  }
  
  
  if(rare_taxa){
    abund_table <- abund_table[,taxa_to_group,drop=FALSE]
  }else{
    #remove the low abundant taxa or columns
    abundant_taxa <-abund_table[,-(taxa_to_group), drop=FALSE]
    #get the rare taxa
    # rare_taxa <-abund_table[,taxa_to_group]
    rare_taxa <- subset(x = abund_table, select = taxa_to_group)
    #get the proportion of each sample that makes up the rare taxa
    rare <- rowSums(rare_taxa)
    #bind the abundant taxa to the rae taxa
    abund_table <- cbind(abundant_taxa,rare)
    #rename the columns i.e the taxa
    colnames(abund_table) <- c(colnames(abundant_taxa),"Rare")
  }
  
  return(abund_table)
}


# Function to collapse the samples in an oTU table with a defined function(fun)
# based on a group in metadata 
collapse_samples <- function(taxon_table,metadata,group,fun=sum, 
                             convertToRelativeAbundance=FALSE){
  # function to collapse the samples in an oTU table with a defined function(fun)  based on a group in metadata 
  # taxon_table - a matrix count table with samples as rows and features/OTUs as columns
  # metadata - a dataframe to containing the group to collapse samples by. Sample names must be the rownames of the metadata
  # group - an independent factor variable within the metadata to collapse the samples by
  # fun - a function without brackets to apply in order to collapse the samples
  # convertToRelativeAbundance - a boolean set to TRUE OR FALSE if the taxon_table shout be converted to relative abundance
  # default is FALSE
  common.ids <- intersect(rownames(taxon_table),rownames(metadata))
  metadata <- droplevels(metadata[common.ids,,drop=FALSE])
  taxon_table <- taxon_table[common.ids,,drop=FALSE]
  taxon_table <- cbind(subset(x = metadata, select=group),taxon_table)
  
  taxon_table <- aggregate(reformulate(termlabels = group, response = '.'),
                           data = taxon_table, FUN = fun)
  rownames(taxon_table) <- taxon_table[,1]
  taxon_table <- taxon_table[,-1]
  if(convertToRelativeAbundance){
    taxon_table <- t(apply(X = taxon_table, MARGIN = 1, FUN = function(x) x/sum(x)))
  }
  
  final <- list(taxon_table,metadata)
  names(final) <- c("taxon_table","metadata")
  return(final)
}



publication_format <- theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.ticks.length=unit(-0.15, "cm"),
        axis.text.x=element_text(margin=ggplot2::margin(t=0.5,r=0,b=0,l=0,unit ="cm")),
        axis.text.y=element_text(margin=ggplot2::margin(t=0,r=0.5,b=0,l=0,unit ="cm")), 
        axis.title = element_text(size = 18,face ='bold.italic', color = 'black'), 
        axis.text = element_text(size = 16,face ='bold', color = 'black'),
        legend.position = 'right', 
        legend.title = element_text(size = 15,face ='bold', color = 'black'),
        legend.text = element_text(size = 14,face ='bold', color = 'black'),
        strip.text =  element_text(size = 14,face ='bold', color = 'black'))





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
remove_rare <- opt[["remove-rare"]]
# taxon / ASV prevalence cutoff
prevalence_cutoff <- opt[["prevalence-cutoff"]] # 0.15 (15%)
# sample / library read count cutoff
library_cutoff <- opt[["library-cutoff"]]  # 100
taxonomy_plots_out_dir <- "taxonomy_plots/"
if(!dir.exists(taxonomy_plots_out_dir)) dir.create(taxonomy_plots_out_dir)


# Read in processed data
metadata <- read_csv(file = metadata_file) %>% as.data.frame()
row.names(metadata) <- metadata[[sample_colname]]
metadata[,sample_colname] <- NULL
group_column_values <- metadata %>% pull(!!sym(groups_colname))
group_levels <- unique(group_column_values)

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


feature_names <- rownames(taxonomy_table)
taxonomy_table  <- process_taxonomy(taxonomy_table)
rownames(taxonomy_table) <- feature_names
taxonomy_table <- fix_names(taxonomy_table, "Other", ";_")



# Subset tables 

# Get features common to the taxonomy and feature table 
common_ids <- intersect(rownames(feature_table), rownames(taxonomy_table))

# Subset the feature and taxonomy tables to contain 
# only features found in both table
feature_table <- feature_table[common_ids,]
taxonomy_table <- taxonomy_table[common_ids,]


# -------------------------Prepare feature tables -------------------------- #
taxon_levels <- colnames(taxonomy_table)
names(taxon_levels) <- taxon_levels
taxon_tables <- map(.x = taxon_levels,
                    .f = make_feature_table,
                    count_matrix = feature_table,
                    taxonomy = taxonomy_table)


# ----------------------- Sample abundance plots -------------------------- #
group_rare <- TRUE
samples_order <- metadata %>% arrange(!!sym(groups_colname)) %>% rownames()
dont_group <- c("phylum")
# In percentage
thresholds <- c(phylum=1,class=3, order=3, family=8, genus=8, species=9)
# Convert from wide to long format
# -1 drops the kingdom level since all the microbes are bacteria
relAbundace_tbs_rare_grouped <- map2(.x = taxon_levels[-1],
                                     .y = taxon_tables[-1], 
                                     .f = function(taxon_level=.x,
                                                   taxon_table=.y){
                                      
                                       
                                       taxon_table <- apply(X = taxon_table, MARGIN = 2,
                                                            FUN = function(x) x/sum(x)) * 100
                                       
                                       
                                       taxon_table <- as.data.frame(taxon_table %>% t())
                                       if(group_rare && !(taxon_level %in% dont_group)){
                                         
                                         taxon_table <- group_low_abund_taxa(taxon_table %>%  
                                                                               as.data.frame(check.names=FALSE,
                                                                                             stringAsFactor=FALSE),
                                                                             threshold =  thresholds[taxon_level])
                                         
                                       }
                                       taxon_table$samples <- rownames(taxon_table)
                                       
                                       
                                       # Change data frame from wide to long format
                                       taxon_table <- taxon_table %>% 
                                         pivot_longer(cols = -samples, names_to = taxon_level, values_to = "relativeAbundance")
                                       taxon_table$samples <- factor(x = taxon_table$samples, 
                                                                     levels = samples_order)
                                       return(taxon_table)
                                     })


x_lab <- "Samples"
y_lab <- "Relative abundance (%)"  
x <-  'samples'  
y <- "relativeAbundance"
facet_by <- reformulate(groups_colname)
number_of_samples <- length(samples_order)
plot_width <- 0.6 * number_of_samples

# Make sample plots
walk2(.x = relAbundace_tbs_rare_grouped, .y = taxon_levels[-1], 
                           .f = function(relAbundace_tb, taxon_level){
                             
                             df <- relAbundace_tb %>%
                               left_join(metadata %>% rownames_to_column("samples"))
                             
                          p <- ggplot(data =  df, mapping = aes(x= !!sym(x), y=!!sym(y) )) +
                               geom_col(aes(fill = !!sym(taxon_level) )) + 
                               facet_wrap(facet_by, scales = "free", nrow = 1) +
                               publication_format +
                               labs(x = x_lab , y = y_lab, fill= tools::toTitleCase(taxon_level)) + 
                               scale_fill_manual(values = custom_palette) +
                               theme(axis.text.x=element_text(
                                 margin=margin(t=0.5,r=0,b=0,l=0,unit ="cm"),
                                 angle = 90, 
                                 hjust = 0.5, vjust = 0.5)) + 
                               labs(x=NULL)
                          
                          ggsave(filename = glue("{taxonomy_plots_out_dir}/{output_prefix}samples_{taxon_level}{assay_suffix}.png"),
                                 plot=p, width = plot_width, height = 8.5, dpi = 300)
                          
                           })



# ------------------------ Group abundance plots ----------------------------- #


# In percentage
thresholds <- c(phylum=1,class=2, order=2, family=2, genus=2, species=2)

# Convert from wide to long format for every treatment group of interest
group_rare <- TRUE
maximum_number_of_taxa <- 500

group_relAbundace_tbs <- map2(.x = taxon_levels[-1], .y = taxon_tables[-1], 
                                     .f = function(taxon_level=.x, taxon_table=.y){
                                       
                                       taxon_table <- as.data.frame(taxon_table %>% t()) 
                                       taxon_table <- (collapse_samples(taxon_table = taxon_table,
                                                                        metadata = metadata, group = groups_colname, 
                                                                        convertToRelativeAbundance = TRUE)$taxon_table * 100 )  %>%  
                                         as.data.frame(check.names=FALSE)
                                       
                                       if(ncol(taxon_table) > maximum_number_of_taxa){
                                         group_rare <- TRUE
                                       }
                                       
                                       if(group_rare){
                                         taxon_table <- group_low_abund_taxa(taxon_table %>%  
                                                                               as.data.frame(check.names=FALSE,
                                                                                             stringAsFactor=FALSE),
                                                                             threshold =  thresholds[taxon_level])
                                         group_rare <- FALSE
                                       }
                                       
                                       taxon_table[,groups_colname] <- rownames(taxon_table)
                                       
                                       
                                       # Change from wide to long format
                                       taxon_table <- taxon_table %>% 
                                         pivot_longer(cols =  -!!sym(groups_colname),
                                                      names_to = taxon_level,
                                                      values_to = "relativeAbundance")
                              
                                       return(taxon_table)
                                       
                                     })


# Make bar plots
y_lab <- "Relative abundance (%)"  
y <- "relativeAbundance"
number_of_groups <- length(group_levels)
plot_width <- 2 * number_of_groups
walk2(.x = group_relAbundace_tbs, .y = taxon_levels[-1], 
                           .f = function(relAbundace_tb=.x, taxon_level=.y){
                             
                             p <- ggplot(data =  relAbundace_tb, mapping = aes(x= !!sym(groups_colname), y = !!sym(y)   )) +
                               geom_col(aes(fill = !!sym(taxon_level))) + 
                               publication_format +
                               theme(axis.text.x=element_text(
                                 margin=margin(t=0.5,r=0,b=0,l=0,unit ="cm"),
                                 angle = 90, 
                                 hjust = 0.5, vjust = 0.5)) + 
                               labs(x = NULL , y = y_lab, fill = tools::toTitleCase(taxon_level)) + 
                               scale_fill_manual(values = custom_palette)
                             ggsave(filename = glue("{taxonomy_plots_out_dir}/{output_prefix}groups_{taxon_level}{assay_suffix}.png"),
                                    plot=p, width = plot_width, height = 10, dpi = 300)
                           })


