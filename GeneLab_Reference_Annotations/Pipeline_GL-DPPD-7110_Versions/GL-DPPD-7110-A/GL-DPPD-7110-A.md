# GeneLab Pipeline for Generating Reference Annotation Tables  

> **This page provides an overview and instructions for how GeneLab generates reference annotation tables. The GeneLab reference annotation table used to add annotations to processed data files is indicated in the exact processing scripts provided for each GLDS dataset under the respective omics datatype subdirectory.**  
 
---

**Date:** August 12, 2024  
**Revision:** -A  
**Document Number:** GL-DPPD-7110-A  

**Submitted by:**  
Alexis Torres and Crystal Han (GeneLab Data Processing Team)  

**Approved by:**  
Sylvain Costes (OSDR Project Manager)  
Samrawit Gebre (GeneLab Deputy Project Manager and Acting GeneLab Configuration Manager)  
Lauren Sanders (OSDR Project Scientist)  
Amanda Saravia-Butler (GeneLab Science Lead)  
Barbara Novak (GeneLab Data Processing Lead)  

---

## Updates from Previous Version

- **Updated Software:**
  - R version updated from 4.1.3 to 4.4.0.  
  - Bioconductor version updated from 3.15.1 to 3.19.1.  

- **Ensembl Releases:**
  - Animals: Updated from release 107 to 112  
  - Plants: Updated from release 54 to 59   
  - Bacteria: Updated from release 54 to 59  

- **New Organism Support:**
  1. Bacillus subtilis, subsp. subtilis 168   
  2. Brachypodium distachyon  
  3. Escherichia coli, str. K-12 substr. MG1655   
  4. Oryzias latipes   
  5. Lactobacillus acidophilus NCFM  
  6. Mycobacterium marinum M  
  7. Oryza sativa Japonica   
  8. Pseudomonas aeruginosa UCBPP-PA14  
  9. Salmonella enterica subsp. enterica serovar Typhimurium str. LT2   
  10. Serratia liquefaciens ATCC 27592   
  11. Staphylococcus aureus MRSA252   
  12. Streptococcus mutans UA159  
  13. Vibrio fischeri ES114   

- **Added NCBI as a Reference Source:** 
  FASTA and GTF files were sourced from NCBI for the following organisms: 
  1. Lactobacillus acidophilus NCFM   
  2. Mycobacterium marinum M   
  3. Pseudomonas aeruginosa UCBPP-PA14   
  4. Serratia liquefaciens ATCC 27592   
  5. Staphylococcus aureus MRSA252   
  6. Streptococcus mutans UA159   
  7. Vibrio fischeri ES114    

- **org.db Creation:**  
  Added functionality to create an annotation database using `AnnotationForge`. This is applicable to organisms without a maintained annotation database package in Bioconductor (e.g., `org.Hs.eg.db`). Currently, this approach is in use for the following organisms:  
  1. Bacillus subtilis, subsp. subtilis 168   
  2. Brachypodium distachyon   
  3. Escherichia coli, str. K-12 substr. MG1655   
  4. Oryzias latipes   
  5. Salmonella enterica subsp. enterica serovar Typhimurium str. LT2   

The pipeline is designed to annotate unique gene IDs in a reference assembly, map them to organism-specific `org.db` databases for additional annotations, integrate STRING DB IDs, and use PANTHER to obtain GO slim IDs based on ENTREZ IDs.

The default columns in the annotation table are:  
- ENSEMBL (or TAIR), SYMBOL, GENENAME, REFSEQ, ENTREZID, STRING_id, GOSLIM_IDS

- For organisms with FASTA and GTF files sourced from NCBI, the LOCUS, OLD_LOCUS, SYMBOL, GENENAME, and GO annotations were directly derived from the GTF file. The `GO` column contains GO terms. `OLD_LOCUS`, or `old_locus_tag` in the GTF was retained when needed to map to STRING IDs.  
- Missing columns indicate the absence of corresponding data for that organism

1. **Brachypodium distachyon (BRADI)**:   
   - Columns: ENSEMBL, ACCNUM, SYMBOL, GENENAME, REFSEQ, ENTREZID, STRING_id, GOSLIM_IDS    
     > Note: GTF `transcript_id` entries were matched with `ACCNUM` keys in the `org.db` and saved as `ACCNUM`

2. **Caenorhabditis elegans (WORM)**:   
   - Columns: ENSEMBL, SYMBOL, GENENAME, REFSEQ, ENTREZID, STRING_id   
     > Note: org.db ENTREZ keys did not match PANTHER ENTREZ keys so the empty `GOSLIM_IDS` column was ommitted

3. **Lactobacillus acidophilus (NCFM)**:   
   - Columns: LOCUS, OLD_LOCUS, SYMBOL, GENENAME, GO, STRING_id   

4. **Mycobacterium marinum (MMARINUMM)**:  
   - Columns: LOCUS, OLD_LOCUS, SYMBOL, GENENAME, GO, STRING_id   

5. **Oryza sativa Japonica (ORYSJ)**:  
   - Columns: ENSEMBL, STRING_id   

6. **Pseudomonas aeruginosa UCBPP-PA14 (PA14)**:  
   - Columns: LOCUS, SYMBOL, GENENAME, GO    

7. **Serratia liquefaciens ATCC 27592 (ATCC27592)**:  
   - Columns: LOCUS, OLD_LOCUS, SYMBOL, GENENAME, GO, STRING_id   

8. **Staphylococcus aureus MRSA252 (MRSA252)**:  
   - Columns: LOCUS, SYMBOL, GENENAME, GO  

9. **Streptococcus mutans UA159 (UA159)**:  
   - Columns: LOCUS, OLD_LOCUS, SYMBOL, GENENAME, GO, STRING_id  

10. **Vibrio fischeri ES114 (ES114)**:  
   - Columns: LOCUS, OLD_LOCUS, SYMBOL, GENENAME, GO, STRING_id   

---

# Table of Contents

- [GeneLab Pipeline for Generating Reference Annotation Tables](#genelab-pipeline-for-generating-reference-annotation-tables)
- [Table of Contents](#table-of-contents)
- [Software Used](#software-used)
- [Annotation Table Build Overview with Example Commands](#annotation-table-build-overview-with-example-commands)
  - [0. Set Up Environment](#0-set-up-environment)
  - [1. Define Variables and Output File Names](#1-define-variables-and-output-file-names)
  - [2. Load Annotation Databases](#2-load-annotation-databases)
  - [3. Build Initial Annotation Table](#3-build-initial-annotation-table)
  - [4. Add org.db Keys](#4-add-orgdb-keys)
  - [5. Add STRING IDs](#5-add-string-ids)
  - [6. Add Gene Ontology (GO) Slim IDs](#6-add-gene-ontology-go-slim-ids)
  - [7. Export Annotation Table and Build Info](#7-export-annotation-table-and-build-info)



---

# Software Used  

| Program       | Version | Relevant Links |
|:--------------|:-------:|:---------------|
| R             |  4.4.0  | [https://www.r-project.org/](https://www.r-project.org/) |
| Bioconductor  | 3.19.1  | [https://bioconductor.org](https://bioconductor.org) |
| tidyverse     |  2.0.0  | [https://www.tidyverse.org](https://www.tidyverse.org) |
| STRINGdb      | 2.16.0  | [https://www.bioconductor.org/packages/release/bioc/html/STRINGdb.html](https://www.bioconductor.org/packages/release/bioc/html/STRINGdb.html) |
| PANTHER.db    | 1.0.12  | [https://bioconductor.org/packages/release/data/annotation/html/PANTHER.db.html](https://www.bioconductor.org/packages/release/data/annotation/html/PANTHER.db.html) |
| rtracklayer   | 1.64.0  | [https://bioconductor.org/packages/release/bioc/html/rtracklayer.html](https://www.bioconductor.org/packages/release/bioc/html/rtracklayer.html) |
| org.At.tair.db| 3.19.1  | [https://bioconductor.org/packages/release/data/annotation/html/org.At.tair.db.html](https://www.bioconductor.org/packages/release/data/annotation/html/org.At.tair.db.html) |
| org.Ce.eg.db  | 3.19.1  | [https://bioconductor.org/packages/release/data/annotation/html/org.Ce.eg.db.html](https://www.bioconductor.org/packages/release/data/annotation/html/org.Ce.eg.db.html) |
| org.Dm.eg.db  | 3.19.1  | [https://bioconductor.org/packages/release/data/annotation/html/org.Dm.eg.db.html](https://www.bioconductor.org/packages/release/data/annotation/html/org.Dm.eg.db.html) |
| org.Dr.eg.db  | 3.19.1  | [https://bioconductor.org/packages/release/data/annotation/html/org.Dr.eg.db.html](https://www.bioconductor.org/packages/release/data/annotation/html/org.Dr.eg.db.html) |
| org.Hs.eg.db  | 3.19.1  | [https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html](https://www.bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html) |
| org.Mm.eg.db  | 3.19.1  | [https://bioconductor.org/packages/release/data/annotation/html/org.Mm.eg.db.html](https://www.bioconductor.org/packages/release/data/annotation/html/org.Mm.eg.db.html) |
| org.Rn.eg.db  | 3.19.1  | [https://bioconductor.org/packages/release/data/annotation/html/org.Rn.eg.db.html](https://www.bioconductor.org/packages/release/data/annotation/html/org.Rn.eg.db.html) |
| org.Sc.sgd.db | 3.19.1  | [https://bioconductor.org/packages/release/data/annotation/html/org.Sc.sgd.db.html](https://www.bioconductor.org/packages/release/data/annotation/html/org.Sc.sgd.db.html) |

---

# Annotation table build overview with example commands  

> Current GeneLab annotation tables are available on [figshare](https://figshare.com/), exact links for each reference organism are provided in the [GL-DPPD-7110-A_annotations.csv](GL-DPPD-7110-A_annotations.csv) file.  
> 
> **[Ensembl Reference Files](https://www.ensembl.org/index.html) Used:**
> - Animals: Ensembl release 112
> - Plants: Ensembl plants release 59
> - Bacteria: Ensembl bacteria release 59


---

This example below is done for *Mus musculus*. All code is executed in R.

## 0. Set Up Environment

```R
# Define variables associated with current pipeline and annotation table versions
GL_DPPD_ID <- "GL-DPPD-7110-A"
ref_tab_path <- "https://raw.githubusercontent.com/nasa/GeneLab_Data_Processing/master/GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv"
readme_path <- "https://github.com/nasa/GeneLab_Data_Processing/tree/master/GeneLab_Reference_Annotations/Workflow_Documentation/GL_RefAnnotTable-A/README.md"

# List currently supported organisms 
currently_accepted_orgs <- c("ARABIDOPSIS", "BACSU", "BRADI", "WORM", "ZEBRAFISH",
                             "FLY", "ECOLI", "HUMAN", "NCFM", "MOUSE",
                             "MMARINUMM", "ORYSJ", "ORYLA", "PA14", "RAT",
                             "YEAST", "SALTY", "ATCC27592", "MRSA252", "UA159",
                             "ES114")

# Import libraries
library(tidyverse)
library(STRINGdb)
library(PANTHER.db)
library(rtracklayer)
```

---

## 1. Define Variables and Output File Names

```R
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
  filter(name == target_organism)

# Extract the relevant columns from the reference table
target_taxid <- target_info$taxon # Taxonomic identifier
target_org_db <- target_info$annotations # org.eg.db R package
target_species_designation <- target_info$species # Full species name
gtf_link <- target_info$gtf # Path to reference assembly GTF

# Error handling for missing values
if (is.na(target_taxid) || is.na(target_org_db) || is.na(target_species_designation) || is.na(gtf_link)) {
  stop(paste("Error: Missing data for target organism", target_organism, "in reference table."))
}

# Create output filenames
base_gtf_filename <- basename(gtf_link)
base_output_name <- str_replace(base_gtf_filename, ".gtf.gz", "")

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
```

<br>

---

## 2. Load Annotation Databases

```R
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
no_org_db <- c("NCFM", "MMARINUMM", "ORYSJ", "PA14", "ATCC27592", "MRSA252", "UA159", "ES114")

# Run the function unless the target_organism is in no_org_db
if (!(target_organism %in% no_org_db) && (target_organism %in% currently_accepted_orgs)) {
  install_and_load_org_db(target_organism, target_org_db, ref_tab_path)
}
```

<br>

---

## 3. Build Initial Annotation Table

```R
# Initialize table from GTF

# Define GTF keys based on the target organism; gene_id conrains unique gene IDs in the reference assembly. Defaults to ENSEMBL

gtf_keytype_mappings <- list(
  ARABIDOPSIS = c(gene_id = "TAIR"),
  BACSU = c(gene_id = "ENSEMBL", gene_name = "SYMBOL"),
  BRADI = c(gene_id = "ENSEMBL", transcript_id = "ACCNUM"),
  WORM = c(gene_id = "ENSEMBL"),  
  ECOLI = c(gene_id = "ENSEMBL", gene_name = "SYMBOL"),
  NCFM = c(gene_id = "LOCUS", old_locus_tag = "OLD_LOCUS", gene = "SYMBOL", product = "GENENAME", Ontology_term = "GO"),
  MMARINUMM = c(gene_id = "LOCUS", old_locus_tag = "OLD_LOCUS", gene = "SYMBOL", product = "GENENAME", Ontology_term = "GO"),
  PA14 = c(gene_id = "LOCUS", gene = "SYMBOL", product = "GENENAME", Ontology_term = "GO"),
  SALTY = c(gene_id = "ENSEMBL", db_xref = "ENTREZID"),
  ATCC27592 = c(gene_id = "LOCUS", old_locus_tag = "OLD_LOCUS", gene = "SYMBOL", product = "GENENAME", Ontology_term = "GO"),
  MRSA252 = c(gene_id = "LOCUS", gene = "SYMBOL", product = "GENENAME", Ontology_term = "GO"),
  UA159 = c(gene_id = "LOCUS", old_locus_tag = "OLD_LOCUS", gene = "SYMBOL", product = "GENENAME", Ontology_term = "GO"),
  ES114 = c(gene_id = "LOCUS", old_locus_tag = "OLD_LOCUS", gene = "SYMBOL", product = "GENENAME", Ontology_term = "GO"),
  default = c(gene_id = "ENSEMBL")
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
  BACSU = "^BSU",
  FLY = "^RR",
  YEAST = "^Y[A-Z0-9]{6}-?[A-Z]?$",
  ECOLI = "^b[0-9]{4}$"
)

# Apply the filter if there's a specific criterion for the target organism
filter_pattern <- filter_criteria[[target_organism]]

if (!is.null(filter_pattern)) {
  if (target_organism == "FLY") {
    annot_gtf <- annot_gtf %>% filter(!grepl(filter_pattern, !!sym(primary_keytype)))
  } else {
    annot_gtf <- annot_gtf %>% filter(grepl(filter_pattern, !!sym(primary_keytype)))
  }
}

# Remove "Gene:" labels on ENTREZ IDs 
if (target_organism == "SALTY") { 
  annot_gtf <- annot_gtf %>% dplyr::mutate(ENTREZID = gsub("^GeneID:", "", ENTREZID)) %>% as.data.frame
}
```

<br>

---

## 4. Add org.db Keys

```R
annot_orgdb <- annot_gtf

# Define the initial keys to pull from the organism-specific database
orgdb_keytypes_list <- list(
  BRADI = c("GENENAME", "REFSEQ", "ENTREZID"),
  ECOLI = c("GENENAME", "REFSEQ", "ENTREZID"),
  WORM = c("SYMBOL", "GENENAME", "REFSEQ", "ENTREZID", "GO"),
  SALTY = c("SYMBOL", "GENENAME", "REFSEQ"),
  YEAST = c("GENENAME", "ALIAS", "REFSEQ", "ENTREZID"),
  default = c("SYMBOL", "GENENAME", "REFSEQ", "ENTREZID")
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
  BACSU = list(query = "SYMBOL", keytype = "SYMBOL"),
  BRADI = list(query = "ACCNUM", keytype = "ACCNUM"),
  WORM = list(query = primary_keytype, keytype = "ENSEMBL"),
  ECOLI = list(query = "SYMBOL", keytype = "SYMBOL"),
  SALTY = list(query = "ENTREZID", keytype = "ENTREZID"),
  default = list(query = primary_keytype, keytype = primary_keytype)
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

# Function to clean and match ACCNUM keys for BRADI
clean_and_match_accnum <- function(annot_table, org_db, query_col, keytype_col, target_column) {
  # Clean the ACCNUM keys in the GTF annotations
  cleaned_annot_keys <- sub("\\..*", "", annot_table[[query_col]])
  
  # Retrieve and clean the org.db keys
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
    if (target_organism == "BRADI" && orgdb_query == "ACCNUM") {
      # For BRADI: use the clean_and_match_accnum function to map to org.db ACCNUM entries
      org_matches <- clean_and_match_accnum(annot_orgdb, get(target_org_db, envir = .GlobalEnv), query_col = orgdb_query, keytype_col = orgdb_keytype, target_column = keytype)
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
if (target_organism == "SALTY") { # Reorder columns to match others; was mismatched since ENTREZ came from GTF
  annot_orgdb <- annot_orgdb[, c("ENSEMBL", "SYMBOL", "GENENAME", "REFSEQ", "ENTREZID")]
}

# For YEAST, Rename ALIAS to GENENAME 
if (target_organism == "YEAST") {
  colnames(annot_orgdb) <- c("ENSEMBL", "SYMBOL", "GENENAME", "REFSEQ", "ENTREZID")
}
```

<br>

---

## 5. Add STRING IDs

```R
# Define organisms that do not use STRING annotations
no_stringdb <- c("PA14", "MRSA252")

# Define the key type used for mapping to STRING
stringdb_query_list <- list(
  NCFM = "OLD_LOCUS",
  MMARINUMM = "OLD_LOCUS",
  ATCC27592 = "OLD_LOCUS",
  UA159 = "OLD_LOCUS",
  ES114 = "OLD_LOCUS",
  default = primary_keytype
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
uses_old_locus <- c("NCFM", "MMARINUMM", "ATCC27592", "UA159", "ES114")
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
if (target_organism == "BACSU") {
  annot_stringdb[[stringdb_query]] <- gsub("^BSU_", "BSU", annot_stringdb[[stringdb_query]])
}

# Map alternative taxonomy IDs for organisms not directly supported by STRING
taxid_map <- list(
  YEAST = 4932,
  BRARP = 51351,
  ATCC27592 = 614
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
  join_key <- if (target_organism %in% c("NCFM", "MMARINUMM", "ATCC27592", "UA159", "ES114")) {
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
if (target_organism == "BACSU") {
  annot_stringdb[[stringdb_query]] <- gsub("^BSU", "BSU_", annot_stringdb[[stringdb_query]])
}

annot_stringdb <- as.data.frame(annot_stringdb)
```

<br>

---

## 6. Add Gene Ontology (GO) slim IDs

```R
# Define organisms that do not use PANTHER annotations 
no_panther_db <- c("WORM", "MMARINUMM", "ORYSJ", "MRSA252", "NCFM", "ATCC27592", "UA159", "ES114", "PA14")

annot_pantherdb <- annot_stringdb

if (!(target_organism %in% no_panther_db)) {
  
  # Define the key type in the annotation table used to map to PANTHER DB
  pantherdb_query = "ENTREZID"
  pantherdb_keytype = "ENTREZ"
  
  # Retrieve target organism PANTHER GO slim annotations database
  pthOrganisms(PANTHER.db) <- target_organism
  
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
```

<br>

---

## 7. Export Annotation Table and Build Info

```R
# Group by primary key to remove any remaining unjoined or duplicate rows
annot <- annot_pantherdb %>%
  group_by(!!sym(primary_keytype)) %>%
  summarise(across(everything(), ~paste(unique(na.omit(.))[unique(na.omit(.)) != ""], collapse = "|")), .groups = 'drop')

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
```

<br>

---

**Pipeline Input data:**

- No input files required, but a target organism must be specified as a positional command line argument

**Pipeline Output data:**

- *-GL-annotations.tsv (Tab delineated table of gene annotations, used to add gene annotations in other GeneLab processing pipelines)
- *-GL-build-info.txt (Text file containing information used to create the annotation table, including tool and tool versions and date of creation)
