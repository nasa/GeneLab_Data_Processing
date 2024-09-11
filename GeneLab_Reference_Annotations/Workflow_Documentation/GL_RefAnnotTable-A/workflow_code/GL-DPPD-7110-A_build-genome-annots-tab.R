#!/usr/bin/env Rscript
# Written by Mike Lee
# GeneLab script for generating organism-specific gene annotation tables
# Example usage: Rscript GL-DPPD-7110-A_build-genome-annots-tab.R 'Mus musculus'

# Define variables associated with current pipeline and annotation table versions
GL_DPPD_ID <- "GL-DPPD-7110-A"
ref_tab_path <- "https://raw.githubusercontent.com/nasa/GeneLab_Data_Processing/master/GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv"
readme_path <- "https://github.com/nasa/GeneLab_Data_Processing/tree/master/GeneLab_Reference_Annotations/Workflow_Documentation/GL_RefAnnotTable-A/README.md"

# List currently supported organisms 
currently_accepted_orgs <- c("Arabidopsis thaliana", "Bacillus subtilis", "Brachypodium distachyon", 
                             "Caenorhabditis elegans", "Danio rerio", "Drosophila melanogaster", 
                             "Escherichia coli", "Homo sapiens", "Lactobacillus acidophilus", 
                             "Mus musculus", "Mycobacterium marinum", "Oryza sativa", 
                             "Oryzias latipes", "Pseudomonas aeruginosa", "Rattus norvegicus", 
                             "Saccharomyces cerevisiae", "Salmonella enterica", "Serratia liquefaciens", 
                             "Staphylococcus aureus", "Streptococcus mutans", "Vibrio fischeri")


#########################################################################
############### Pull in and check command line arguments ################
#########################################################################

# Pull in command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Get the target organism (CLI argument 1) and check that it is listed in currently_accepted_orgs
validate_arguments <- function(args, supported_orgs) {
  if (length(args) < 1) {
    stop("One positional argument is required that specifies the target organism. Available options are:\n", paste(supported_orgs, collapse = ", "))
  }
  
  # Convert the first argument to uppercase
  target_organism <- toupper(args[1])
  
  # Check if the uppercased target organism is in the uppercased supported_orgs
  if (!target_organism %in% sapply(supported_orgs, toupper)) {
    stop(paste0("'", target_organism, "' is not currently supported."))
  }
  
  return(args[1])
}

target_organism <- validate_arguments(args, currently_accepted_orgs)

# If provided, get the reference table URL from CLI arguments (CLI argument 2) and update ref_tab_path
ref_tab_path <- if (length(args) >= 2) args[2] else ref_tab_path


#########################################################################
######################## Set up environment #############################
#########################################################################

required_packages <- c("tidyverse", "STRINGdb", "PANTHER.db", "rtracklayer")
# Check for required packages other than the org-specific db #
report_package_needed <- function(package_name) {
  cat(paste0("\n  The package '", package_name, "' is required. Please see:\n"))
  cat("    https://github.com/nasa/GeneLab_Data_Processing/tree/master/GeneLab_Reference_Annotations/Workflow_Documentation/GL_RefAnnotTable-A/README.md\n\n")
  quit()
}

# Check and report missing packages other than the org-specific db
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    report_package_needed(pkg)
  }
}

# Import libraries
library(tidyverse)
library(STRINGdb)
library(PANTHER.db)
library(rtracklayer)


#########################################################################
############## Define variables and output file names ###################
#########################################################################

# Set timeout time to ensure annotation file downloads will complete
options(timeout = 600)

ref_table <- tryCatch(
  read.csv(ref_tab_path),
  error = function(e) {
    message <- paste("Error: Unable to read the reference table from the path provided. Please check the path and try again.\nPath:", ref_tab_path)
    stop(message)
  }
)

# Get target organism information
target_info <- ref_table %>%
  filter(species == target_organism)

# Extract the relevant columns from the reference table
target_taxid <- target_info$taxon # Taxonomic identifier
target_org_db <- target_info$annotations # org.eg.db R package
gtf_link <- target_info$gtf # Path to reference assembly GTF
target_short_name <- target_info$name # PANTHER / UNIPROT short name; blank if not available
ref_source <- target_info$ref_source # Reference files source  

# Error handling for missing values
if (is.na(target_taxid) || is.na(target_org_db) || is.na(target_organism) || is.na(gtf_link)) {
  stop(paste("Error: Missing data for target organism", target_organism, "in reference table."))
}

# Create output filenames
base_gtf_filename <- basename(gtf_link)
base_output_name <- str_replace(base_gtf_filename, ".gtf.gz", "")

# Add the species name to base_output_name if the reference source is not ENSEMBL
if (!(ref_source %in% c("ensembl_plants", "ensembl_bacteria", "ensembl"))) {
  base_output_name <- paste(str_replace(target_organism, " ", "_"), base_output_name, sep = "_")
}

out_table_filename <- paste0(base_output_name, "-GL-annotations.tsv")
out_log_filename <- paste0(base_output_name, "-GL-build-info.txt")

# Check if output file already exists and if it does, exit without overwriting
if ( file.exists(out_table_filename) ) {
  cat("\n-------------------------------------------------------------------------------------------------\n")
  cat(paste0("\n  The file that would be created, '", out_table_filename, "', exists already.\n"))
  cat(paste0("  We don't want to overwrite it accidentally. Move it and run this again if wanting to proceed.\n"))
  cat("\n-------------------------------------------------------------------------------------------------\n")
  quit()
}


#############################################
######## Load annotation databases  #########
#############################################

# Set timeout time to ensure annotation file downloads will complete
options(timeout = 600)

####### GTF ##########

# Create the GTF dataframe from its path, unique gene identities in the reference assembly are under 'gene_id'
GTF <- rtracklayer::import(gtf_link)
GTF <- data.frame(GTF)

###### org.db ########

# Define a function to load the specified org.db package for a given target organism
install_and_load_org_db <- function(target_organism, target_org_db, ref_tab_path) {
  if (!is.na(target_org_db) && target_org_db != "") {
    # Attempt to install the package from Bioconductor
    BiocManager::install(target_org_db, ask = FALSE)
    
    # Check if the package was successfully loaded
    if (!requireNamespace(target_org_db, quietly = TRUE)) {
      # If not, attempt to create it locally using a helper script
      source("install-org-db.R")
      target_org_db <- install_annotations(target_organism, ref_tab_path)
    }
  } else {
    # If target_org_db is NA or empty, create it locally using the helper script
    source("install-org-db.R")
    target_org_db <- install_annotations(target_organism, ref_tab_path)
  }
  
  # Load the package into the R session
  library(target_org_db, character.only = TRUE)
}

# Define list of supported organisms which do not use annotations from an org.db
no_org_db <- c("Lactobacillus acidophilus", "Mycobacterium marinum", "Oryza sativa", "Pseudomonas aeruginosa",
              "Serratia liquefaciens", "Staphylococcus aureus", "Streptococcus mutans", "Vibrio fischeri")

# Run the function unless the target_organism is in no_org_db
if (!(target_organism %in% no_org_db) && (target_organism %in% currently_accepted_orgs)) {
  install_and_load_org_db(target_organism, target_org_db, ref_tab_path)
}


############################################
######## Build annotation table ############
############################################

# Initialize table from GTF

# Define GTF keys based on the target organism; gene_id conrains unique gene IDs in the reference assembly. Defaults to ENSEMBL

gtf_keytype_mappings <- list(
  "Arabidopsis thaliana" = c(gene_id = "TAIR"),
  "Bacillus subtilis" = c(gene_id = "ENSEMBL", gene_name = "SYMBOL"),
  "Brachypodium distachyon" = c(gene_id = "ENSEMBL", transcript_id = "ACCNUM"),
  "Caenorhabditis elegans" = c(gene_id = "ENSEMBL"),  
  "Escherichia coli" = c(gene_id = "ENSEMBL", gene_name = "SYMBOL"),
  "Lactobacillus acidophilus" = c(gene_id = "LOCUS", old_locus_tag = "OLD_LOCUS", gene = "SYMBOL", product = "GENENAME", Ontology_term = "GO"),
  "Mycobacterium marinum" = c(gene_id = "LOCUS", old_locus_tag = "OLD_LOCUS", gene = "SYMBOL", product = "GENENAME", Ontology_term = "GO"),
  "Pseudomonas aeruginosa" = c(gene_id = "LOCUS", gene = "SYMBOL", product = "GENENAME", Ontology_term = "GO"),
  "Salmonella enterica" = c(gene_id = "ENSEMBL", db_xref = "ENTREZID"),
  "Serratia liquefaciens" = c(gene_id = "LOCUS", old_locus_tag = "OLD_LOCUS", gene = "SYMBOL", product = "GENENAME", Ontology_term = "GO"),
  "Staphylococcus aureus" = c(gene_id = "LOCUS", gene = "SYMBOL", product = "GENENAME", Ontology_term = "GO"),
  "Streptococcus mutans" = c(gene_id = "LOCUS", old_locus_tag = "OLD_LOCUS", gene = "SYMBOL", product = "GENENAME", Ontology_term = "GO"),
  "Vibrio fischeri" = c(gene_id = "LOCUS", old_locus_tag = "OLD_LOCUS", gene = "SYMBOL", product = "GENENAME", Ontology_term = "GO"),
  "default" = c(gene_id = "ENSEMBL")
)

# Get the key types for the target organism or use the default
wanted_gtf_keytypes <- if (!is.null(gtf_keytype_mappings[[target_organism]])) {
  gtf_keytype_mappings[[target_organism]]
} else {
  c(gene_id = "ENSEMBL")
}

# Initialize the annotation table from the GTF, keeping only the wanted_gtf_keytypes
annot_gtf <- GTF[, names(wanted_gtf_keytypes), drop = FALSE]
annot_gtf <- annot_gtf %>% distinct()

# Rename the columns in the annot_gtf dataframe according to the key types
colnames(annot_gtf) <- wanted_gtf_keytypes

# Save the name of the primary key type (gene_id) being used 
primary_keytype <- wanted_gtf_keytypes[1]

# Filter out unwanted genes from the GTF

# Define filtering criteria for specific organisms
filter_criteria <- list(
  "Bacillus subtilis" = "^BSU",
  "Drosophila melanogaster" = "^RR",
  "Saccharomyces cerevisiae" = "^Y[A-Z0-9]{6}-?[A-Z]?$",
  "Escherichia coli" = "^b[0-9]{4}$"
)

# Apply the filter if there's a specific criterion for the target organism
filter_pattern <- filter_criteria[[target_organism]]

if (!is.null(filter_pattern)) {
  if (target_organism == "Drosophila melanogaster") {
    annot_gtf <- annot_gtf %>% filter(!grepl(filter_pattern, !!sym(primary_keytype)))
  } else {
    annot_gtf <- annot_gtf %>% filter(grepl(filter_pattern, !!sym(primary_keytype)))
  }
}

# Remove "Gene:" labels on ENTREZ IDs 
if (target_organism == "Salmonella enterica") { 
  annot_gtf <- annot_gtf %>% dplyr::mutate(ENTREZID = gsub("^GeneID:", "", ENTREZID)) %>% as.data.frame
}

#########################################################################
########################### Add org.db keys #############################
#########################################################################

annot_orgdb <- annot_gtf

# Define the initial keys to pull from the organism-specific database
orgdb_keytypes_list <- list(
  "Brachypodium distachyon" = c("GENENAME", "REFSEQ", "ENTREZID"),
  "Escherichia coli" = c("GENENAME", "REFSEQ", "ENTREZID"),
  "Caenorhabditis elegans" = c("SYMBOL", "GENENAME", "REFSEQ", "ENTREZID", "GO"),
  "Salmonella enterica" = c("SYMBOL", "GENENAME", "REFSEQ"),
  "Saccharomyces cerevisiae" = c("GENENAME", "ALIAS", "REFSEQ", "ENTREZID"),
  "default" = c("SYMBOL", "GENENAME", "REFSEQ", "ENTREZID")
)

# Add entries for organisms in no_org_db as character(0) (no keys wanted from the org.db)
for (organism in no_org_db) {
  orgdb_keytypes_list[[organism]] <- character(0)
}

wanted_org_db_keytypes <- if (target_organism %in% names(orgdb_keytypes_list)) {
  orgdb_keytypes_list[[target_organism]]
} else {
  orgdb_keytypes_list[["default"]]
}

# Define mappings for query and keytype based on target organism
orgdb_keytype_mappings <- list(
  "Bacillus subtilis" = list(query = "SYMBOL", keytype = "SYMBOL"),
  "Brachypodium distachyon" = list(query = "ACCNUM", keytype = "ACCNUM"),
  "Caenorhabditis elegans" = list(query = primary_keytype, keytype = "ENSEMBL"),
  "Escherichia coli" = list(query = "SYMBOL", keytype = "SYMBOL"),
  "Salmonella enterica" = list(query = "ENTREZID", keytype = "ENTREZID"),
  "default" = list(query = primary_keytype, keytype = primary_keytype)
)

# Define the orgdb_query, this is the key type that will be used to map to the org.db
orgdb_query <- if (!is.null(orgdb_keytype_mappings[[target_organism]])) {
  orgdb_keytype_mappings[[target_organism]][["query"]]
} else {
  orgdb_keytype_mappings[["default"]][["query"]]
}

# Define the orgdb_keytype, this is the name of the key type in the org.db 
orgdb_keytype <- if (!is.null(orgdb_keytype_mappings[[target_organism]])) {
  orgdb_keytype_mappings[[target_organism]][["keytype"]]
} else {
  orgdb_keytype_mappings[["default"]][["keytype"]]
}

# Function to remove version numbers from ACCNUM keys and match them for BRADI
match_accnum <- function(annot_table, org_db, query_col, keytype_col, target_column) {
  # Remove version numbers from the ACCNUM keys in the GTF annotations
  cleaned_annot_keys <- sub("\\..*", "", annot_table[[query_col]])
  
  # Retrieve and remove version numbers from the org.db keys
  orgdb_keys <- keys(org_db, keytype = keytype_col)
  cleaned_orgdb_keys <- sub("\\..*", "", orgdb_keys)
  
  # Create a lookup table for matching cleaned keys to original keys
  lookup_table <- setNames(orgdb_keys, cleaned_orgdb_keys)
  
  # Match cleaned GTF keys to original org.db keys
  matched_keys <- lookup_table[cleaned_annot_keys]
  
  # Use the matched keys to retrieve the target annotations from org.db
  mapIds(org_db, keys = matched_keys, keytype = keytype_col, column = target_column, multiVals = "list")
}

# Loop through the desired key types and add annotations to the GTF table
for (keytype in wanted_org_db_keytypes) {
  # Check if keytype is a valid column in the target org.db
  if (keytype %in% columns(get(target_org_db, envir = .GlobalEnv))) {
    if (target_organism == "Brachypodium distachyon" && orgdb_query == "ACCNUM") {
      # For BRADI: use the match_accnum function to map to org.db ACCNUM entries
      org_matches <- match_accnum(annot_orgdb, get(target_org_db, envir = .GlobalEnv), query_col = orgdb_query, keytype_col = orgdb_keytype, target_column = keytype)
    } else {
      # Default mapping for other organisms
      org_matches <- mapIds(get(target_org_db, envir = .GlobalEnv), keys = annot_orgdb[[orgdb_query]], keytype = orgdb_keytype, column = keytype, multiVals = "list")
    }
    # Add the mapped annotations to the GTF table
    annot_orgdb[[keytype]] <- sapply(org_matches, function(x) paste(x, collapse = "|"))
  } else {
    # Set column to NA if keytype is not present in org.db
    annot_orgdb[[keytype]] <- NA
  }
}

# For SALTY, reorder columns to mtach other tables
if (target_organism == "Salmonella enterica") { # Reorder columns to match others; was mismatched since ENTREZ came from GTF
  annot_orgdb <- annot_orgdb[, c("ENSEMBL", "SYMBOL", "GENENAME", "REFSEQ", "ENTREZID")]
}

# For YEAST, Rename ALIAS to GENENAME 
if (target_organism == "Saccharomyces cerevisiae") {
  colnames(annot_orgdb) <- c("ENSEMBL", "SYMBOL", "GENENAME", "REFSEQ", "ENTREZID")
}

#########################################################################
########################### Add STRING IDs ##############################
#########################################################################

# Define organisms that do not use STRING annotations
no_stringdb <- c("Pseudomonas aeruginosa", "Staphylococcus aureus")

# Define the key type used for mapping to STRING
stringdb_query_list <- list(
  "Lactobacillus acidophilus" = "OLD_LOCUS",
  "Mycobacterium marinum" = "OLD_LOCUS",
  "Serratia liquefaciens" = "OLD_LOCUS",
  "Streptococcus mutans" = "OLD_LOCUS",
  "Vibrio fischeri" = "OLD_LOCUS",
  "default" = primary_keytype
)

# Define the key type for mapping in STRING, using the default if necessary
stringdb_query <- if (!is.null(stringdb_query_list[[target_organism]])) {
  stringdb_query_list[[target_organism]]
} else {
  stringdb_query_list[["default"]]
}

# Handle organisms which do not use the GTF's gene_id keys to map to STRING 
# These are microbial species for which NCBI references were used rather than ENSEMBL,
# for which the STRING accessions match the GTF's gene_name keys, but not the gene_id keys.
uses_old_locus <- c("Lactobacillus acidophilus", "Mycobacterium marinum", "Serratia liquefaciens", "Streptococcus mutans", "Vibrio fischeri")
# Handle STRING annotation processing based on the target organism
if (target_organism %in% uses_old_locus) {
  # If the target organism is one of the NOENTRY organisms, handle the OLD_LOCUS splitting
  annot_stringdb <- annot_orgdb %>%
    separate_rows(!!sym(stringdb_query), sep = ",", convert = TRUE) %>%
    distinct() %>%
    as.data.frame()
} else {
  # For other organisms, collapse on the primary key
  annot_stringdb <- annot_orgdb %>% distinct()
  annot_stringdb <- annot_stringdb %>%
    group_by(!!sym(primary_keytype)) %>%
    summarise(across(everything(), ~paste(unique(na.omit(.))[unique(na.omit(.)) != ""], collapse = "|")), .groups = 'drop') %>%
    as.data.frame()
}

# Replace "BSU_" with "BSU" in the primary_keytype column for BACSU before STRING mapping
if (target_organism == "Bacillus subtilis") {
  annot_stringdb[[stringdb_query]] <- gsub("^BSU_", "BSU", annot_stringdb[[stringdb_query]])
}

# Map alternative taxonomy IDs for organisms not directly supported by STRING
taxid_map <- list(
  "Saccharomyces cerevisiae" = 4932,
  "Brassica rapa" = 51351,
  "Serratia liquefaciens" = 614
)

# Assign the alternative taxonomy identifier if applicable
target_taxid <- if (!is.null(taxid_map[[target_organism]])) {
  taxid_map[[target_organism]]
} else {
  target_taxid
}

# Initialize string_map
string_map <- NULL

# If the target organism is supported by STRING, get STRING annotations
if (!(target_organism %in% no_stringdb)) {
  string_db <- STRINGdb$new(version = "12.0", species = target_taxid, score_threshold = 0)
  string_map <- string_db$map(annot_stringdb, stringdb_query, removeUnmappedRows = FALSE, takeFirst = FALSE)
}
if (!is.null(string_map)) {
  annot_stringdb <- annot_stringdb %>%
    group_by(!!sym(primary_keytype)) %>%
    summarise(across(everything(), ~paste(unique(na.omit(.))[unique(na.omit(.)) != ""], collapse = "|")), .groups = 'drop')
  
  string_map <- string_map %>%
    group_by(!!sym(primary_keytype)) %>%
    summarise(across(everything(), ~paste(unique(na.omit(.))[unique(na.omit(.)) != ""], collapse = "|")), .groups = 'drop')
}

if (!is.null(string_map)) {
  # Determine the appropriate join key
  join_key <- if (target_organism %in% c("Lactobacillus acidophilus", "Mycobacterium marinum", "Serratia liquefaciens", "Streptococcus mutans", "Vibrio fischeri")) {
    primary_keytype
  } else {
    stringdb_query
  }
  
  # Add temporary column to add string IDs to annotation table
  annot_stringdb <- annot_stringdb %>%
    mutate(join_key = toupper(!!sym(join_key)))
  
  string_map <- string_map %>%
    mutate(join_key = toupper(!!sym(join_key)))
  
  # Join STRING IDs to the annotation table
  annot_stringdb <- left_join(annot_stringdb, string_map %>% dplyr::select(join_key, STRING_id), by = "join_key") %>%
    dplyr::select(-join_key)
}

# Undo the "BSU_" to "BSU" replacement for BACSU after STRING mapping
if (target_organism == "Bacillus subtilis") {
  annot_stringdb[[stringdb_query]] <- gsub("^BSU", "BSU_", annot_stringdb[[stringdb_query]])
}

annot_stringdb <- as.data.frame(annot_stringdb)

#########################################################################
################ Add Gene Ontology (GO) slim IDs ########################
#########################################################################

# Define organisms that do not use PANTHER annotations 
no_panther_db <- c("Caenorhabditis elegans", "Mycobacterium marinum", "Oryza sativa", "Staphylococcus aureus", "Lactobacillus acidophilus", "Serratia liquefaciens", "Streptococcus mutans", "Vibrio fischeri", "Pseudomonas aeruginosa")

annot_pantherdb <- annot_stringdb

if (!(target_organism %in% no_panther_db)) {
  
  # Define the key type in the annotation table used to map to PANTHER DB
  pantherdb_query = "ENTREZID"
  pantherdb_keytype = "ENTREZ"
  
  # Retrieve target organism PANTHER GO slim annotations database using the UNIPROT / PANTHER short name
  pthOrganisms(PANTHER.db) <- target_short_name
  
  # Define a function to retrieve GO slim IDs for a given gene's ENTREZIDs, which may include entries separated by a "|"
  get_go_slim_ids <- function(entrez_id) {
    if (is.na(entrez_id) || entrez_id == "NA") {
      return("NA")
    }
    
    entrez_ids <- unlist(strsplit(entrez_id, "|", fixed = TRUE))
    go_ids <- lapply(entrez_ids, function(id) {
      mapIds(PANTHER.db, keys = id, keytype = pantherdb_keytype, column = "GOSLIM_ID", multiVals = "list")
    })
    
    # Flatten the list and remove duplicates
    go_ids <- unique(unlist(go_ids))
    
    if (length(go_ids) == 0) {
      return("NA")
    } else {
      return(paste(go_ids, collapse = "|"))
    }
  }
  
  # Apply the GO slim ID mapping function to all valid rows
  annot_pantherdb <- annot_pantherdb %>%
    mutate(GOSLIM_IDS = sapply(get(pantherdb_query), get_go_slim_ids))
}


#########################################################################
############# Export annotation table and build info ####################
#########################################################################

# Group by primary key to remove any remaining unjoined or duplicate rows
annot <- annot_pantherdb %>%
  group_by(!!sym(primary_keytype)) %>%
  summarise(across(everything(), ~paste(unique(na.omit(.))[unique(na.omit(.)) != ""], collapse = "|")), .groups = 'drop')

# If "GO" column exists, move it to the end to keep columns in consistent order across organisms
if ("GO" %in% names(annot)) {
  go_column <- annot$GO
  annot$GO <- NULL
  annot$GO <- go_column
}

# Sort the annotation table based on primary keytype gene IDs
annot <- annot %>% arrange(.[[1]])

# Replace any blank cells with NA 
annot[annot == "" | annot == "NA"] <- NA

# Export the annotation table
write.table(annot, out_table_filename, sep = "\t", quote = FALSE, row.names = FALSE)

# Define the date when the annotation table was generated
date_generated <- format(Sys.time(), "%d-%B-%Y")

# Export annotation build information
writeLines(paste(c("Based on:\n    ", GL_DPPD_ID), collapse = ""), out_log_filename)
write(paste(c("\nBuild done on:\n    ", date_generated), collapse = ""), out_log_filename, append = TRUE)
write(paste(c("\nUsed gtf file:\n    ", gtf_link), collapse = ""), out_log_filename, append = TRUE)
if (!(target_organism %in% no_org_db)) {
  write(paste(c("\nUsed ", target_org_db, " version:\n    ", packageVersion(target_org_db) %>% as.character()), collapse = ""), out_log_filename, append = TRUE)
}
write(paste(c("\nUsed STRINGdb version:\n    ", packageVersion("STRINGdb") %>% as.character()), collapse = ""), out_log_filename, append = TRUE)
write(paste(c("\nUsed PANTHER.db version:\n    ", packageVersion("PANTHER.db") %>% as.character()), collapse = ""), out_log_filename, append = TRUE)

write("\n\nAll session info:\n", out_log_filename, append = TRUE)
write(capture.output(sessionInfo()), out_log_filename, append = TRUE)
