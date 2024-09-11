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
  - tidyverse version updated from 1.3.2 to 2.0.0.  
  - STRINGdb version updated from 2.8.4 to 2.16.0 (DB version: 12.0).   
  - PANTHER.db version updated from 1.0.11 to 1.0.12 (DB version: 18.0).  
  - rtracklayer version updated from 1.56.1 to 1.64.0.

- **Added Software:**
  - AnnotationForge version 1.46.0.
  - biomaRt version 2.60.1.
  - GO.db version 3.19.1 (DB schema version 2.1)    

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
  Added functionality to create an annotation database using `AnnotationForge`. This is applicable to organisms without a maintained annotation database package in Bioconductor (e.g., `org.Hs.eg.db`). This approach was used for the following organisms:  
  1. Bacillus subtilis, subsp. subtilis 168   
  2. Brachypodium distachyon   
  3. Escherichia coli, str. K-12 substr. MG1655   
  4. Oryzias latipes   
  5. Salmonella enterica subsp. enterica serovar Typhimurium str. LT2   

The pipeline is designed to annotate unique gene IDs in a reference assembly, map them to organism-specific `org.db` databases for additional annotations, integrate STRING DB IDs, and use PANTHER to obtain GO slim IDs based on ENTREZ IDs.

The default columns in the annotation table are:  
- ENSEMBL (or TAIR), SYMBOL, GENENAME, REFSEQ, ENTREZID, STRING_id, GOSLIM_IDS

- For organisms with FASTA and GTF files sourced from NCBI, the LOCUS, OLD_LOCUS, SYMBOL, GENENAME, and GO annotations were directly derived from the GTF file. The `GO` column contains GO terms. `OLD_LOCUS`, or `old_locus_tag` in the GTF was retained when needed to map to STRING IDs.  
- Missing columns indicate the absence of corresponding data for that organism.

1. **Brachypodium distachyon**:   
   - Columns: ENSEMBL, ACCNUM, SYMBOL, GENENAME, REFSEQ, ENTREZID, STRING_id, GOSLIM_IDS    
     > Note: GTF `transcript_id` entries were matched with `ACCNUM` keys in the `org.db` and saved as `ACCNUM`

2. **Caenorhabditis elegans**:   
   - Columns: ENSEMBL, SYMBOL, GENENAME, REFSEQ, ENTREZID, STRING_id   
     > Note: org.db ENTREZ keys did not match PANTHER ENTREZ keys so the empty `GOSLIM_IDS` column was ommitted

3. **Lactobacillus acidophilus**:   
   - Columns: LOCUS, OLD_LOCUS, SYMBOL, GENENAME, STRING_id, GO   

4. **Mycobacterium marinum**:  
   - Columns: LOCUS, OLD_LOCUS, SYMBOL, GENENAME, STRING_id, GO   

5. **Oryza sativa Japonica**:  
   - Columns: ENSEMBL, STRING_id   

6. **Pseudomonas aeruginosa UCBPP-PA14**:  
   - Columns: LOCUS, SYMBOL, GENENAME, GO    

7. **Serratia liquefaciens ATCC 27592**:  
   - Columns: LOCUS, OLD_LOCUS, SYMBOL, GENENAME, STRING_id, GO   

8. **Staphylococcus aureus MRSA252**:  
   - Columns: LOCUS, SYMBOL, GENENAME, GO  

9. **Streptococcus mutans UA159**:  
   - Columns: LOCUS, OLD_LOCUS, SYMBOL, GENENAME, STRING_id, GO  

10. **Vibrio fischeri ES114**:  
   - Columns: LOCUS, OLD_LOCUS, SYMBOL, GENENAME, STRING_id, GO   

---

# Table of Contents

- [Software Used](#software-used)
- [Annotation Table Build Overview with Example Commands](#annotation-table-build-overview-with-example-commands)
  - [0. Set Up Environment](#0-set-up-environment)
  - [1. Define Variables and Output File Names](#1-define-variables-and-output-file-names)
  - [2. Create the Organism Package if it is Not Hosted by Bioconductor](#2-create-the-organism-package-if-it-is-not-hosted-by-bioconductor)
  - [3. Load Annotation Databases](#3-load-annotation-databases)
  - [4. Build Initial Annotation Table](#4-build-initial-annotation-table)
  - [5. Add org.db Keys](#5-add-orgdb-keys)
  - [6. Add STRING IDs](#6-add-string-ids)
  - [7. Add Gene Ontology (GO) Slim IDs](#7-add-gene-ontology-go-slim-ids)
  - [8. Export Annotation Table and Build Info](#8-export-annotation-table-and-build-info)

---

# Software Used  

| Program         | Version | Relevant Links |
|:----------------|:-------:|:---------------|
| R               |  4.4.0  | [https://www.r-project.org/](https://www.r-project.org/) |
| Bioconductor    | 3.19.1  | [https://bioconductor.org](https://bioconductor.org) |
| tidyverse       |  2.0.0  | [https://www.tidyverse.org](https://www.tidyverse.org) |
| STRINGdb        | 2.16.0  | [https://www.bioconductor.org/packages/release/bioc/html/STRINGdb.html](https://www.bioconductor.org/packages/release/bioc/html/STRINGdb.html) |
| PANTHER.db      | 1.0.12  | [https://bioconductor.org/packages/release/data/annotation/html/PANTHER.db.html](https://www.bioconductor.org/packages/release/data/annotation/html/PANTHER.db.html) |
| rtracklayer     | 1.64.0  | [https://bioconductor.org/packages/release/bioc/html/rtracklayer.html](https://www.bioconductor.org/packages/release/bioc/html/rtracklayer.html) |
| org.At.tair.db  | 3.19.1  | [https://bioconductor.org/packages/release/data/annotation/html/org.At.tair.db.html](https://www.bioconductor.org/packages/release/data/annotation/html/org.At.tair.db.html) |
| org.Ce.eg.db    | 3.19.1  | [https://bioconductor.org/packages/release/data/annotation/html/org.Ce.eg.db.html](https://www.bioconductor.org/packages/release/data/annotation/html/org.Ce.eg.db.html) |
| org.Dm.eg.db    | 3.19.1  | [https://bioconductor.org/packages/release/data/annotation/html/org.Dm.eg.db.html](https://www.bioconductor.org/packages/release/data/annotation/html/org.Dm.eg.db.html) |
| org.Dr.eg.db    | 3.19.1  | [https://bioconductor.org/packages/release/data/annotation/html/org.Dr.eg.db.html](https://www.bioconductor.org/packages/release/data/annotation/html/org.Dr.eg.db.html) |
| org.Hs.eg.db    | 3.19.1  | [https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html](https://www.bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html) |
| org.Mm.eg.db    | 3.19.1  | [https://bioconductor.org/packages/release/data/annotation/html/org.Mm.eg.db.html](https://www.bioconductor.org/packages/release/data/annotation/html/org.Mm.eg.db.html) |
| org.Rn.eg.db    | 3.19.1  | [https://bioconductor.org/packages/release/data/annotation/html/org.Rn.eg.db.html](https://www.bioconductor.org/packages/release/data/annotation/html/org.Rn.eg.db.html) |
| org.Sc.sgd.db   | 3.19.1  | [https://bioconductor.org/packages/release/data/annotation/html/org.Sc.sgd.db.html](https://www.bioconductor.org/packages/release/data/annotation/html/org.Sc.sgd.db.html) |
| AnnotationForge | 1.46.0  | [https://bioconductor.org/packages/AnnotationForge](https://bioconductor.org/packages/AnnotationForge) |
| biomaRt         |  2.60.1  | [https://bioconductor.org/packages/biomaRt](https://bioconductor.org/packages/biomaRt) |
| GO.db           |  3.19.1  | [https://bioconductor.org/packages/GO.db](https://bioconductor.org/packages/GO.db) |

---

# Annotation table build overview with example commands  

Current GeneLab annotation tables are available on [figshare](https://figshare.com/), exact links for each reference organism are provided in the [GL-DPPD-7110-A_annotations.csv](GL-DPPD-7110-A_annotations.csv) file.  

**[Ensembl Reference Versions](https://www.ensembl.org/index.html):**
- Animals: Ensembl release 112
- Plants: Ensembl plants release 59
- Bacteria: Ensembl bacteria release 59  

**Database Versions:**
- STRINGdb: 12.0  
- PANTHERdb: 18.0  
  > Note: The values in the 'name' column of [GL-DPPD-7110-A_annotations.csv](GL-DPPD-7110-A_annotations.csv) (e.g., HUMAN, MOUSE, RAT) are derived from the short names used in PANTHER. These short names are subject to change.  
- GO.db:
  - GO ontology file updated on 2024-01-17
  - Entrez gene data updated on 2024-03-12
  - DB schema version 2.1



---

This example below is done for *Mus musculus*. All code is executed in R.

## 0. Set Up Environment

```R
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

# Import libraries
library(tidyverse)
library(STRINGdb)
library(PANTHER.db)
library(rtracklayer)
```
**Input Data:**

- None (This is an initial setup step using predefined variables)

**Output Data:**

- GL_DPPD_ID (GeneLab Data Processing Pipeline Document ID)
- ref_tab_path (path to the reference table CSV file)
- readme_path (path to the README file)
- currently_accepted_orgs (list of currently supported organisms)

<br>

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
  filter(species == target_organism)

# Extract the relevant columns from the reference table
target_taxid <- target_info$taxon # Taxonomic identifier
target_org_db <- target_info$annotations # org.eg.db R package
target_species_designation <- target_info$species # Full species name
gtf_link <- target_info$gtf # Path to reference assembly GTF
target_short_name <- target_info$name # PANTHER / UNIPROT short name; blank if not available
ref_source <- target_info$ref_source # Reference files source  

# Error handling for missing values
if (is.na(target_taxid) || is.na(target_org_db) || is.na(target_species_designation) || is.na(gtf_link)) {
  stop(paste("Error: Missing data for target organism", target_organism, "in reference table."))
}

# Create output filenames
base_gtf_filename <- basename(gtf_link)
base_output_name <- str_replace(base_gtf_filename, ".gtf.gz", "")

# Add the species name to base_output_name if the reference source is not ENSEMBL
if (!(ref_source %in% c("ensembl_plants", "ensembl_bacteria", "ensembl"))) {
  base_output_name <- paste(str_replace(target_species_designation, " ", "_"), base_output_name, sep = "_")
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
```
**Input Data:**

- ref_tab_path (path to the reference table CSV file, output from [step 0](#0-set-up-environment))
- target_organism (name of the target organism for which annotations are being generated)

**Output Data:**

- target_taxid (taxonomic identifier for the target organism)
- target_org_db (name of the org.db R package for the target organism)
- target_species_designation (full species name of the target organism)
- gtf_link (URL to the GTF file for the target organism)
- target_short_name (PANTHER/UNIPROT short name for the target organism)
- ref_source (source of the reference files, e.g., "ensembl", "ensembl_plants", "ensembl_bacteria", "ncbi")
- out_table_filename (name of the output annotation table file)
- out_log_filename (name of the output log file)

<br>

---

## 2. Create the Organism Package if it is Not Hosted by Bioconductor

```R
# Use AnnotationForge's makeOrgPackageFromNCBI function with default settings to create the organism-specific org.db R package from available NCBI annotations

# Try to download the org.db from Bioconductor, build it locally if installation fails
BiocManager::install(target_org_db, ask = FALSE) 
if (!requireNamespace(target_org_db, quietly = TRUE)) { 
  tryCatch({
    # Parse organism's name in the reference table to create the org.db name (target_org_db)
    genus_species <- strsplit(target_species_designation, " ")[[1]]
    if (length(genus_species) < 1) {
        stop("Species designation is not correctly formatted: ", target_species_designation)
    }
    genus <- genus_species[1]
    species <- ifelse(length(genus_species) > 1, genus_species[2], "")
    strain <- ref_table %>%
      filter(species == target_organism) %>%
      pull(strain) %>%
      gsub("[^A-Za-z0-9]", "", .)
    if (!is.na(strain) && strain != "") {
        species <- paste0(species, strain)
    }
    target_org_db <- paste0("org.", substr(genus, 1, 1), species, ".eg.db")
    
    BiocManager::install(c("AnnotationForge", "biomaRt", "GO.db"), ask = FALSE) 
    library(AnnotationForge)
    makeOrgPackageFromNCBI( 
        version = "0.1",
        author = "Your Name <your.email@example.com>",
        maintainer = "Your Name <your.email@example.com>",
        outputDir = "./",
        tax_id = target_taxid,
        genus = genus,
        species = species
    )
    install.packages(file.path("./", target_org_db), repos = NULL, type = "source", quiet = TRUE)
    cat(paste0("'", target_org_db, "' has been successfully built and installed.\n"))
  }, error = function(e) {
      stop("Failed to build and load the package: ", target_org_db, "\nError: ", e$message)
  })
  target_org_db <- install_annotations(target_organism, ref_tab_path)
}
```

**Input Data:**

- target_org_db (name of the org.db R package for the target organism, output from [step 1](#1-define-variables-and-output-file-names))
- target_species_designation (full species name of the target organism, output from [step 1](#1-define-variables-and-output-file-names))
- ref_table (reference table containing organism-specific information, output from [step 1](#1-define-variables-and-output-file-names))
- target_organism (name of the target organism, output from [step 1](#1-define-variables-and-output-file-names))
- target_taxid (taxonomic identifier for the target organism, output from [step 1](#1-define-variables-and-output-file-names))

**Output Data:**

- target_org_db (updated name of the org.db R package, if it was created locally)
- Locally installed org.db package (if the package is not available on Bioconductor, a new package is created and installed)

<br>

---

## 3. Load Annotation Databases

```R
# Set timeout time to ensure annotation file downloads will complete
options(timeout = 600)

####### GTF ##########

# Create the GTF dataframe from its path, unique gene identities in the reference assembly are under 'gene_id'
GTF <- rtracklayer::import(gtf_link)
GTF <- data.frame(GTF)

###### org.db ########

# Load the package into the R session
library(target_org_db, character.only = TRUE)

# Define list of supported organisms which do not use annotations from an org.db
no_org_db <- c("Lactobacillus acidophilus", "Mycobacterium marinum", "Oryza sativa", "Pseudomonas aeruginosa",
              "Serratia liquefaciens", "Staphylococcus aureus", "Streptococcus mutans", "Vibrio fischeri")

# Run the function unless the target_organism is in no_org_db
if (!(target_organism %in% no_org_db) && (target_organism %in% currently_accepted_orgs)) {
  install_and_load_org_db(target_organism, target_org_db, ref_tab_path)
}
```

**Input Data:**

- gtf_link (URL to the GTF file for the target organism, output from [step 1](#1-define-variables-and-output-file-names))
- target_org_db (name of the org.eg.db R package for the target organism, output from [steps 1](#1-define-variables-and-output-file-names) and [2](#2-create-the-organism-package-if-it-is-not-hosted-by-bioconductor))
- target_organism (name of the target organism, output from [step 1](#1-define-variables-and-output-file-names))
- currently_accepted_orgs (list of currently supported organisms, output from [step 0](#0-set-up-environment))
- ref_tab_path (path to the reference table CSV file, output from [step 0](#0-set-up-environment))

**Output Data:**

- GTF (data frame containing the GTF file for the target organism)
- no_org_db (list of organisms that do not use org.db annotations due to inconsistent gene names across GTF and org.db)

<br>

---

## 4. Build Initial Annotation Table

```R
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
```

**Input Data:**

- GTF (data frame containing the parsed GTF file for the target organism, output from [step 3](#3-load-annotation-databases))
- target_organism (target organism's full species name, output from [step 1](#1-define-variables-and-output-file-names))
- gtf_keytype_mappings (list of keys to extract from the GTF, for each organism)

**Output Data:**

- annot_gtf (initial annotation table derived from the GTF file, containing only the relevant columns for the target organism)
- primary_keytype (the name of the primary key type being used, e.g., "ENSEMBL", "TAIR", "LOCUS", based on the GTF gene_id entries)

<br>

---

## 5. Add org.db Keys

```R
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
```

**Input Data:**

- annot_gtf (initial annotation table derived from the GTF file, output from [step 4](#4-build-initial-annotation-table))
- target_organism (target organism's full species name, output from [step 1](#1-define-variables-and-output-file-names))
- no_org_db (list of organisms that do not use annotations from an org.db, output from [step 3](#3-load-annotation-databases))
- primary_keytype (the name of the primary key type being used, output from [step 4](#4-build-initial-annotation-table))
- target_org_db (name of the org.eg.db R package for the target organism, output from [steps 1](#1-define-variables-and-output-file-names) and [2](#2-create-the-organism-package-if-it-is-not-hosted-by-bioconductor))

**Output Data:**

- annot_orgdb (updated annotation table with additional keys from the organism-specific org.db)
- orgdb_query (the key type used to map to the org.db)
- orgdb_keytype (the name of the key type in the org.db)

<br>

---

## 6. Add STRING IDs

```R
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
```

**Input Data:**

- annot_orgdb (annotation table with GTF and org.db annotations, output from [step 5](#5-add-orgdb-keys))
- target_organism (target organism's full species name, output from [step 1](#1-define-variables-and-output-file-names))
- primary_keytype (the name of the primary key type being used, output from [step 4](#4-build-initial-annotation-table))
- target_taxid (taxonomic identifier for the target organism, output from [step 1](#1-define-variables-and-output-file-names))

**Output Data:**

- annot_stringdb (updated annotation table with added STRING IDs)
- no_stringdb (list of organisms that do not use STRING annotations)
- stringdb_query (the key type used for mapping to STRING database)
- uses_old_locus (list of organisms where GTF gene_id entries do not match those in STRING, so entries in OLD_LOCUS are used to query STRING)

<br>

---

## 7. Add Gene Ontology (GO) slim IDs

```R
# Define organisms that do not use PANTHER annotations 
no_panther_db <- c("Caenorhabditis elegans", "Mycobacterium marinum", "Oryza sativa", "Staphylococcus aureus", "Lactobacillus acidophilus", "Serratia liquefaciens", "Streptococcus mutans", "Vibrio fischeri", "Pseudomonas aeruginosa")

annot_pantherdb <- annot_stringdb

if (!(target_organism %in% no_panther_db)) {
  
  # Define the key type in the annotation table used to map to PANTHER DB
  pantherdb_query = "ENTREZID"
  pantherdb_keytype = "ENTREZ"
  
  # Retrieve target organism PANTHER GO slim annotations database using the UNIPROT / PANTHER short name
  target_short_name <- target_species_designation
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
```

**Input Data:**

- annot_orgdb (annotation table with GTF and org.db annotations, output from [step 5](#5-add-orgdb-keys))
- target_organism (target organism's full species name, output from [step 1](#1-define-variables-and-output-file-names))
- primary_keytype (the name of the primary key type being used, output from [step 4](#4-build-initial-annotation-table))
- target_taxid (taxonomic identifier for the target organism, output from [step 1](#1-define-variables-and-output-file-names))

**Output Data:**

- annot_stringdb (updated annotation table with added STRING IDs)
- no_stringdb (list of organisms that do not use STRING annotations)
- stringdb_query (the key type used for mapping to STRING database)
- uses_old_locus (list of organisms where the 'gene_id' column in the GTF dataframe does not match STRING identifiers, so the 'old_locus_tag' column from the GTF dataframe is used to query STRING instead)

<br>

---

## 8. Export Annotation Table and Build Info

```R
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
```

**Input Data:**

- annot_pantherdb (annotation table with GTF, org.db, STRING, and PANTHER annotations, output from [step 7](#7-add-gene-ontology-go-slim-ids))
- primary_keytype (the name of the primary key type being used, output from [step 4](#4-build-initial-annotation-table))
- out_table_filename (name of the output annotation table file, output from [step 1](#1-define-variables-and-output-file-names))
- out_log_filename (name of the output log file, output from [step 1](#1-define-variables-and-output-file-names))
- GL_DPPD_ID (GeneLab Data Processing Pipeline Document ID, output from [step 0](#0-set-up-environment))
- gtf_link (URL to the GTF file for the target organism, output from [step 1](#1-define-variables-and-output-file-names))
- target_org_db (name of the org.eg.db R package for the target organism, output from [steps 1](#1-define-variables-and-output-file-names) and [2](#2-create-the-organism-package-if-it-is-not-hosted-by-bioconductor))
- no_org_db (list of organisms that do not use org.db annotations, output from [step 3](#3-load-annotation-databases))

**Output Data:**

- annot (final annotation table with annotations from the GTF, org.db, STRING, and PANTHER)
- ***-GL-annotations.tsv** (annot saved as a tab-delimited table file)
- ***-GL-build-info.txt** (annotation table build information log file)

<br>

---

**Pipeline Input data:**

- No input files required, but a target organism must be specified as a positional command line argument

**Pipeline Output data:**

- ***-GL-annotations.tsv** (Tab-delineated table of gene annotations, used to add gene annotations in other GeneLab processing pipelines)
- ***-GL-build-info.txt** (Text file containing information used to create the annotation table, including tool and tool versions and date of creation)
