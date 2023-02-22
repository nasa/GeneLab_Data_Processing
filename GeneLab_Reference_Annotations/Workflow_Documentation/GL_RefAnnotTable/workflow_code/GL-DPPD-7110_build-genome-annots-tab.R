#!/usr/bin/env Rscript

# Maintained by Mike Lee 
# GeneLab script for generating organism ENSEMBL annotation tables
# Example usage: Rscript GL-DPPD-7110_build-genome-annots-tab.R MOUSE

GL_DPPD_ID <- "GL-DPPD-7110"

#########################################################################
############### Pull In and Check Command Line Arguments ################
#########################################################################


## Import command line arguments ##

args <- commandArgs(trailingOnly = TRUE)

## Define currently acceptable input organisms (matching names in ref organisms.csv table) ##

currently_accepted_orgs <- c("ARABIDOPSIS",
                             "FLY",
                             "HUMAN",
                             "MOUSE",
                             "RAT",
                             "WORM",
                             "YEAST",
                             "ZEBRAFISH",
                             "BACSU")
                             # "ECOLI") # need to work out links to ref fasta and gtf files

## Check that at least one positional command line argument was provided ##

if ( length(args) < 1 ) {
    cat("\n  One positional argument is required that specifies the target organism. Currently available include:\n")

    for ( item in currently_accepted_orgs ) {

        cat(paste0("\n        ", item))
    }

    cat("\n\n")

    quit()

} else {

    suppressWarnings(target_organism <- toupper(args[1]))

}


## Check that the positional argument provided is acceptable ##

if ( ! target_organism %in% currently_accepted_orgs ) {

    cat(paste0("\n  '", args[1], "' isn't a valid entry.\n"))

    cat("\n  The currently available organisms include:\n")

    for ( item in currently_accepted_orgs ) {

        cat(paste0("\n        ", item))
    }

    cat("\n\n")

    quit()

}


## checking for required packages other than the org-specific db ##

# helper function for pointing to GL setup page if missing a package
report_package_needed <- function(package_name) {
    cat(paste0("\n  The package '", package_name, "' is required. Please see:\n"))
    cat("    https://github.com/nasa/GeneLab_Data_Processing/tree/master/GeneLab_Reference_Annotations/Workflow_Documentation/GL_RefAnnotTable/README.md\n\n")
    quit()
}

# checking and reporting
if (!requireNamespace("tidyverse", quietly = TRUE))
    report_package_needed("tidyverse")

if (!requireNamespace("BiocManager", quietly = TRUE))
    report_package_needed("BiocManager")

if (!requireNamespace("STRINGdb", quietly = TRUE))
    report_package_needed("STRINGdb")

if (!requireNamespace("PANTHER.db", quietly = TRUE))
    report_package_needed("PANTHER.db")

if (!requireNamespace("rtracklayer", quietly = TRUE))
    report_package_needed("rtracklayer")

#########################################################################
######################## Set Up Environment #############################
#########################################################################

## Import libraries ##

library(tidyverse)
library(STRINGdb)
library(PANTHER.db)
library(rtracklayer)


## Set the primary annotation keytype, TAIR for Arabidopsis, ENSEMBL for all other organisms ##

if ( target_organism == "ARABIDOPSIS" ) {

    primary_keytype <- "TAIR"

} else {

    primary_keytype <- "ENSEMBL"

}

## Define annotation keys to retrieve ##

wanted_keys_vec <- c("SYMBOL", "GENENAME", "REFSEQ", "ENTREZID")

## Define links to tables containing species-specific annotation info ##

ref_tab_link <-
    
"https://raw.githubusercontent.com/nasa/GeneLab_Data_Processing/master/GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110/GL-DPPD-7110_annotations.csv"


#########################################################################
############## Define Variables and Output File Names ###################
#########################################################################


## Set timeout time to ensure annotation file downloads will complete ##

options(timeout = 600)

## Read in tables containing species-specific annotation info ##

ref_table <- read.csv(ref_tab_link)

## Retrieve and define target organism taxid, annotation database name, and scientific name ##

target_taxid <- ref_table %>%
    filter(name == target_organism) %>%
    pull(taxon)

target_org_db <- ref_table %>%
    filter(name == target_organism) %>%
    pull(annotations)

target_species_designation <- ref_table %>%
    filter(name == target_organism) %>%
    pull(species)

## Define link to Ensembl annotation gtf file for the target organism ##

gtf_link <- ref_table %>%
    filter(species == target_species_designation) %>%
    pull(gtf)

## Create output files names ##

base_gtf_filename <- basename(gtf_link)
base_output_name <- str_replace(base_gtf_filename, ".gtf.gz", "")

out_table_filename <- paste0(base_output_name, "-GL-annotations.tsv")
out_log_filename <- paste0(base_output_name, "-GL-build-info.txt")

## Check if output file already exists and if it does, exit without overwriting ##

if ( file.exists(out_table_filename) ) {

    cat("\n-------------------------------------------------------------------------------------------------\n")
    cat(paste0("\n  The file that would be created, '", out_table_filename, "', exists already.\n"))
    cat(paste0("  We don't want to overwrite it accidentally. Move it and run this again if wanting to proceed.\n"))
    cat("\n-------------------------------------------------------------------------------------------------\n")

    quit()

}


#########################################################################
######## Load Annotation Databases and Retrieve Unique Gene IDs #########
#########################################################################


## Import Ensembl annotation gtf file for the target organism ##

gtf_obj <- import(gtf_link)

## Define unique Ensembl IDs ##

unique_IDs <- gtf_obj$gene_id %>% unique()

## Remove gtf object to conserve RAM, since it is no longer needed ##

rm(gtf_obj)

## Redefine target organism annotation database ##

ann.dbi <- target_org_db

## Install target organism annotation database if not already installed, then load the annotation database library ##

if ( ! require(ann.dbi, character.only = TRUE)) {

    BiocManager::install(ann.dbi, ask = FALSE)

}

library(ann.dbi, character.only = TRUE)


#########################################################################
######################## Build Annotation Table #########################
#########################################################################

## Begin annotation table using unique IDs of the primary keytype ##

annot <- data.frame(unique_IDs)
colnames(annot) <- primary_keytype

## Retrieve and add additional annotation keys as table columns ##

for ( key in wanted_keys_vec ) {

    if ( key %in% columns(eval(parse(text = ann.dbi), env = .GlobalEnv))) {

        new_list <- mapIds(eval(parse(text = ann.dbi), env = .GlobalEnv), keys = unique_IDs, keytype = primary_keytype, column = key, multiVals = "list")

        # they come as lists when we accept the multiple hits, so converting to character strings here
        annot[[key]] <- sapply(new_list, paste, collapse = "|")

    } else {

        # if the annotation DB didn't have any of the wanted key types, that column will be missing
        # adding in here as an empty column
        annot[key] <- NA
    }
}


#########################################################################
########################### Add STRING IDs ##############################
#########################################################################

## Retrieve target organism STRING protein-protein interaction database and create STRING ID map to the primary keytype ##

# for YEAST, the only one in STRINGdb is the primary taxid 4932, so switching to that here
if ( target_organism == "YEAST" ) {
    target_taxid <- 4932
}

string_db <- STRINGdb$new(version = "11.5", species = target_taxid, score_threshold = 0)
string_map <- string_db$map(annot, primary_keytype, removeUnmappedRows = FALSE, takeFirst = FALSE)


## Adding some blank lines just for spacing on print-out ##
cat("\n\n")

## Create a table using the gene IDs of the primary keytype as row names and a column containing STRING IDs. For genes containig multiple STRING IDs, combine all STRING IDs for each gene into one row and separate each ID with a '|' ##

tab_with_multiple_STRINGids_combined <-
    data.frame(row.names = annot[[primary_keytype]])

for ( curr_gene_ID in row.names(tab_with_multiple_STRINGids_combined) ) {

    curr_STRING_ids <- string_map %>%
        filter(!!rlang::sym(primary_keytype) == curr_gene_ID) %>%
        pull(STRING_id) %>% paste(collapse = "|")

    tab_with_multiple_STRINGids_combined[curr_gene_ID, "STRING_id"] <- curr_STRING_ids

}

## Move the primary keytype gene IDs back to being a column in the STRING ID table (since they were switched to row names above) ##

tab_with_multiple_STRINGids_combined <-
    tab_with_multiple_STRINGids_combined %>%
    rownames_to_column(primary_keytype)

## Add the STRING ID column to the annotation table ##

annot <- dplyr::left_join(annot,
                          tab_with_multiple_STRINGids_combined,
                          by = primary_keytype)



#########################################################################
################ Add Gene Ontology (GO) slim IDs ########################
#########################################################################


## Retrieve target organism PANTHER GO slim annotations database ##

pthOrganisms(PANTHER.db) <- target_organism

## Use ENTREZ IDs to map genes to respective PANTHER GO slim annotation(s) ##

## Note: Since there can be none (indicated in the annotation table as "NA"), one, or multiple ENTREZ IDs for a gene, this section contains 3 distinct parts to handle each of those scenarios and create a new column in the annotation table containg the GO slim IDs ## 

for ( curr_row in 1:dim(annot)[1] ) {

    curr_entry <- annot[curr_row, "ENTREZID"]

    ## For genes without an ENTREZ ID ##
    if ( curr_entry == "NA" ) {

        annot[curr_row, "GOSLIM_IDS"] <- "NA"

    } else if ( ! grepl("|", curr_entry, fixed = TRUE) ) {

        ## For genes with one ENTREZ ID ##
        curr_GO_IDs <- mapIds(PANTHER.db, keys = curr_entry, keytype = "ENTREZ", column = "GOSLIM_ID", multiVals = "list") %>% unlist() %>% as.vector()

        ## Add "NA" to the GO slim column for ENTREZ IDs that do not contain a respective GO slim ID ##
        if ( is.null(curr_GO_IDs) ) {

            curr_GO_IDs <- "NA"
        }

        annot[curr_row, "GOSLIM_IDS"] <- paste(curr_GO_IDs, collapse = "|")

    } else {

        ## For genes with multiple ENTREZ ID ##
        ## Note: In this scenario, the ENTREZ IDs for each gene are first split with a '|' to separate the IDs, then the GO slim ID(s) for each ENTREZ ID are collected and combined, then duplicates are removed, and the final list of GO slim IDs for each gene are added in a single row, separated with a '|' ## 

        ## Split the ENTREZ IDs ##
        curr_entry_vec <- strsplit(curr_entry, "|", fixed = TRUE)

        ## Start a vector of current GO slim IDs ##
        curr_GO_IDs <- vector()

        ## Collect and combine GO slim ID(s) for each ENTREZ ID ##
        for ( curr_entry in curr_entry_vec ) {

            new_GO_IDs <- mapIds(PANTHER.db, keys = curr_entry, keytype = "ENTREZ", column = "GOSLIM_ID", multiVals = "list") %>% unlist() %>% as.vector()

            ## Add new GO slim IDs to the GO slim IDs vector ##
            curr_GO_IDs <- c(curr_GO_IDs, new_GO_IDs)

        }

        ## Remove duplicate GO slim IDs ##
        curr_GO_IDs <- unique(curr_GO_IDs)

        ## Add "NA" to the GO slim vector for ENTREZ IDs that do not contain a respective GO slim ID ##
        if ( length(curr_GO_IDs) == 0 ) {

            curr_GO_IDs <- "NA"
        }

        ## Add additional GO slim IDs to the GOSLIM ID column in the annotation table ## 
        annot[curr_row, "GOSLIM_IDS"] <- paste(curr_GO_IDs, collapse = "|")

    }

}


#########################################################################
############# Export Annotation Table and Build Info ####################
#########################################################################

## Sort the annotation table based on primary keytype gene IDs ##

annot <- annot %>% arrange(.[[1]])

## Replacing any blank cells with NA ##
annot[annot == ""] <- NA

## Export the annotation table ##

write.table(annot, out_table_filename, sep = "\t", quote = FALSE, row.names = FALSE)

## Define the date the annotation table was generated ##

date_generated <- format(Sys.time(), "%d-%B-%Y")

## Export annotation table build info ##

writeLines(paste(c("Based on:\n    ", GL_DPPD_ID), collapse = ""), out_log_filename)
write(paste(c("\nBuild done on:\n    ", date_generated), collapse = ""), out_log_filename, append = TRUE)
write(paste(c("\nUsed gtf file:\n    ", gtf_link), collapse = ""), out_log_filename, append = TRUE)
write(paste(c("\nUsed ", ann.dbi, " version:\n    ", packageVersion(ann.dbi) %>% as.character()), collapse = ""), out_log_filename, append = TRUE)
write(paste(c("\nUsed STRINGdb version:\n    ", packageVersion("STRINGdb") %>% as.character()), collapse = ""), out_log_filename, append = TRUE)
write(paste(c("\nUsed PANTHER.db version:\n    ", packageVersion("PANTHER.db") %>% as.character()), collapse = ""), out_log_filename, append = TRUE)

write("\n\nAll session info:\n", out_log_filename, append = TRUE)
write(capture.output(sessionInfo()), out_log_filename, append = TRUE)
