# GeneLab pipeline for generating reference annotation tables

> **This page holds an overview and instructions for how GeneLab generates reference annotation tables. The GeneLab reference annotation table used to add annotations to processed data files are indicated in the exact processing scripts provided for each GLDS dataset under the respective omics datatype subdirectory.**  

---

**Date:** July 26, 2022  
**Revision:** -  
**Document Number:** GL-DPPD-7110  

**Submitted by:**  
Mike Lee (GeneLab Data Processing Team)

**Approved by:**  
Amanda Saravia-Butler (GeneLab Data Processing Lead)  
Sylvain Costes (GeneLab Project Manager)  
Samrawit Gebre (GeneLab Deputy Project Manager and Interim GeneLab Configuration Manager)  
Jonathan Galazka (GeneLab Project Scientist)

---

# Table of contents  

- [GeneLab pipeline for generating reference annotation tables](#genelab-pipeline-for-generating-reference-annotation-tables)
- [Table of contents](#table-of-contents)
- [Software used](#software-used)
- [Annotation table build overview with example commands](#annotation-table-build-overview-with-example-commands)
  - [0. Set Up Environment](#0-set-up-environment)
  - [1. Define Variables and Output File Names](#1-define-variables-and-output-file-names)
  - [2. Load Annotation Databases and Retrieve Unique Gene IDs](#2-load-annotation-databases-and-retrieve-unique-gene-ids)
  - [3. Build Initial Annotation Table](#3-build-initial-annotation-table)
  - [4. Add STRING IDs](#4-add-string-ids)
  - [5. Add Gene Ontology (GO) slim IDs](#5-add-gene-ontology-go-slim-ids)
  - [6. Export Annotation Table and Build Info](#6-export-annotation-table-and-build-info)


---

# Software used  

|Program|Version|Relevant Links|
|:------|:------:|:-------------|
|R|4.2.1|[https://www.r-project.org/](https://www.r-project.org/)|
|Bioconductor|3.15|[https://bioconductor.org](https://bioconductor.org)|
|tidyverse|1.3.2|[https://www.tidyverse.org](https://www.tidyverse.org)|
|STRINGdb|2.8.4|[https://www.bioconductor.org/packages/release/bioc/html/STRINGdb.html](https://www.bioconductor.org/packages/release/bioc/html/STRINGdb.html)|
|PANTHER.db|1.0.11|[https://bioconductor.org/packages/release/data/annotation/html/PANTHER.db.html](https://bioconductor.org/packages/release/data/annotation/html/PANTHER.db.html)|
|rtracklayer|1.56.1|[https://bioconductor.org/packages/release/bioc/html/rtracklayer.html](https://bioconductor.org/packages/release/bioc/html/rtracklayer.html)
|org.Hs.eg.db|3.15.0|[https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html](https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html)|
|org.Mm.eg.db|3.15.0|[https://bioconductor.org/packages/release/data/annotation/html/org.Mm.eg.db.html](https://bioconductor.org/packages/release/data/annotation/html/org.Mm.eg.db.html)|
|org.Rn.eg.db|3.15.0|[https://bioconductor.org/packages/release/data/annotation/html/org.Rn.eg.db.html](https://bioconductor.org/packages/release/data/annotation/html/org.Rn.eg.db.html)
|org.Dm.eg.db|3.15.0|[https://bioconductor.org/packages/release/data/annotation/html/org.Dm.eg.db.html](https://bioconductor.org/packages/release/data/annotation/html/org.Dm.eg.db.html)|
|org.Ce.eg.db|3.15.0|[https://bioconductor.org/packages/release/data/annotation/html/org.Ce.eg.db.html](https://bioconductor.org/packages/release/data/annotation/html/org.Ce.eg.db.html)|
|org.At.tair.db|3.15.0|[https://bioconductor.org/packages/release/data/annotation/html/org.At.tair.db.html](https://bioconductor.org/packages/release/data/annotation/html/org.At.tair.db.html)|
|org.EcK12.eg.db|3.15.0|[https://bioconductor.org/packages/release/data/annotation/html/org.EcK12.eg.db.html](https://bioconductor.org/packages/release/data/annotation/html/org.EcK12.eg.db.html)|
|org.Sc.sgd.db|3.15.0|[https://bioconductor.org/packages/release/data/annotation/html/org.Sc.sgd.db.html](https://bioconductor.org/packages/release/data/annotation/html/org.Sc.sgd.db.html)|

---

# Annotation table build overview with example commands  

> Current GeneLab annotation tables are available on [figshare](https://figshare.com/), exact links for each reference organism are provided in the [GL-DPPD-7110_annotations.csv](GL-DPPD-7110_annotations.csv) file.  
> 
> **[Ensembl Reference Files](https://www.ensembl.org/index.html) Used:**
> - Animals: Ensembl release 107
> - Plants: Ensembl plants release 54
> - Bacteria: Ensembl bacteria release 54


---

This example below is done for *Mus musculus*. All code is executed in R.

## 0. Set Up Environment

```R
target_organism == "MOUSE"

GL_DPPD_ID <- "GL-DPPD-7110"

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


## Set timeout time to allow more time for annotation file downloads to complete ##
options(timeout = 600)
```

---

## 1. Define Variables and Output File Names

```R
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

## Create output file names ##
base_gtf_filename <- basename(gtf_link)
base_output_name <- str_replace(base_gtf_filename, ".gtf.gz", "")

out_table_filename <- paste0(base_output_name, "-GL-annotations.tsv")
out_log_filename <- paste0(base_output_name, "-GL-build-info.txt")
```

<br>

---

## 2. Load Annotation Databases and Retrieve Unique Gene IDs

```R
## Import Ensembl annotation gtf file for the target organism ##
gtf_obj <- import(gtf_link)

## Define unique Ensembl IDs ##
unique_IDs <- gtf_obj$gene_id %>% unique()

## Remove gtf object to conserve RAM, since it is no longer needed ##
rm(gtf_obj)

## Define target organism annotation database ##
ann.dbi <- target_org_db

## Install target organism annotation database if not already installed, then load the annotation database library ##
if ( ! require(ann.dbi, character.only = TRUE)) {

    BiocManager::install(ann.dbi, ask = FALSE)

}

library(ann.dbi, character.only = TRUE)
```

<br>

---

## 3. Build Initial Annotation Table

```R
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

```

<br>

---

## 4. Add STRING IDs

```R
## Retrieve target organism STRING protein-protein interaction database and create STRING ID map to the primary keytype ##
string_db <- STRINGdb$new(version = "11.5", species = target_taxid, score_threshold = 0)
string_map <- string_db$map(annot, primary_keytype, removeUnmappedRows = FALSE, takeFirst = FALSE)

## Create a table using the gene IDs of the primary keytype as row names and a column containing STRING IDs. ##
## For genes containig multiple STRING IDs, combine all STRING IDs for each gene into one row and separate each ID with a '|' ##
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
```

<br>

---

## 5. Add Gene Ontology (GO) slim IDs

```R
## Retrieve target organism PANTHER GO slim annotations database ##
pthOrganisms(PANTHER.db) <- target_organism

## Use ENTREZ IDs to map genes to respective PANTHER GO slim annotation(s) ##
# Note: Since there can be none (indicated in the annotation table as "NA"), one, or 
# multiple ENTREZ IDs for a gene, this section contains 3 distinct parts to handle
# each of those scenarios and create a new column in the annotation table containg the GO slim IDs

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
        # Note: In this scenario, the ENTREZ IDs for each gene are first split with a '|' to
        # separate the IDs, then the GO slim ID(s) for each ENTREZ ID are collected and
        # combined, then duplicates are removed, and the final list of GO slim IDs for
        # each gene are added in a single row, separated with a '|'

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
```

<br>

---

## 6. Export Annotation Table and Build Info

```R
## Sort the annotation table based on primary keytype gene IDs ##
annot <- annot %>% arrange(.[[1]])

## Replacing any blank cells with NA ##
annot[annot == ""] <- NA

## Export the annotation table using the file name defined in Step 1 ##
write.table(annot, out_table_filename, sep = "\t", quote = FALSE, row.names = FALSE)

## Define the date the annotation table was generated ## 
date_generated <- format(Sys.time(), "%d-%B-%Y")

## Export annotation table build info using the file name defined in Step 1 ##
writeLines(paste(c("Based on:\n    ", GL_DPPD_ID), collapse = ""), out_log_filename)
write(paste(c("Build done on:\n    ", date_generated), collapse = ""), out_log_filename, append = TRUE)
write(paste(c("\nUsed gtf file:\n    ", gtf_link), collapse = ""), out_log_filename, append = TRUE)
write(paste(c("\nUsed ", ann.dbi, " version:\n    ", packageVersion(ann.dbi) %>% as.character()), collapse = ""), out_log_filename, append = TRUE)
write(paste(c("\nUsed STRINGdb version:\n    ", packageVersion("STRINGdb") %>% as.character()), collapse = ""), out_log_filename, append = TRUE)
write(paste(c("\nUsed PANTHER.db version:\n    ", packageVersion("PANTHER.db") %>% as.character()), collapse = ""), out_log_filename, append = TRUE)

write("\n\nAll session info:\n", out_log_filename, append = TRUE)
write(capture.output(sessionInfo()), out_log_filename, append = TRUE)
```

<br>

---

**Pipeline Input data:**

- No input files required, but a target organism must be specified as a positional command line argument

**Pipeline Output data:**

- *-GL-annotations.tsv (Tab delineated table of gene annotations, used to add gene annotations in other GeneLab processing pipelines)
- *-GL-build-info.txt (Text file containing information used to create the annotation table, including tool and tool versions and date of creation)
