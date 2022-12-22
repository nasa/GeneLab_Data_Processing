# GeneLab bioinformatics processing pipeline for Agilent 1-channel microarray data <!-- omit in toc -->

> **This page holds an overview and instructions for how GeneLab processes Agilent 1-channel microarray datasets. Exact processing commands and GL-DPPD-XXXX version used for specific datasets are provided with their processed data in the [GeneLab Data Systems 
(GLDS) repository](https://genelab-data.ndc.nasa.gov/genelab/projects).**  
---

**Date:** TBD MONTH NN, 2023  
**Revision:** -  
**Document Number:** GL-DPPD-XXXX   

**Submitted by:**  
Jonathan Oribello (GeneLab Data Processing Team)

**Approved by:**  
Sylvain Costes (GeneLab Project Manager) 
Samrawit Gebre (GeneLab Deputy Project Manager) 
Amanda Saravia-Butler (GeneLab Data Processing Lead) 
Lauren Sanders (acting GeneLab Project Scientist) 

---

# Table of contents  <!-- omit in toc -->

- [Software used](#software-used)
- [General processing overview with example commands](#general-processing-overview-with-example-commands)
    - [1. Create Sample RunSheet](#1-create-sample-runsheet)
    - [2. Load Metadata and Raw Data](#2-load-metadata-and-raw-data)
    - [3. Raw Data Quality Assessment](#3-raw-data-quality-assessment)
      - [3a. Density Plot](#3a-density-plot)
      - [3b. Pseudo Image Plots](#3b-pseudo-image-plots)
      - [3c. MA Plots](#3c-ma-plots)
      - [3d. Foreground-Background Plots](#3d-foreground-background-plots)
      - [3e. Boxplots](#3e-boxplots)
    - [4. Background Correction](#4-background-correction)
    - [5. Between Array Normalization](#5-between-array-normalization)
    - [6. Normalized Data Quality Assessment](#6-normalized-data-quality-assessment)
      - [6a. Density Plot](#6a-density-plot)
      - [6b. Pseudo Image Plots](#6b-pseudo-image-plots)
      - [6c. MA Plots](#6c-ma-plots)
      - [6d. Boxplots](#6d-boxplots)
    - [7. Probeset Differential Expression (DE)](#7-probeset-differential-expression)
      - [7a. Add Probeset Annotations](#7a-add-probeset-annotations)
      - [7b. Summarize Biomart Mapping vs. Manufacturer Mapping](#7b-summarize-biomart-mapping-vs-manufacturer-mapping)
      - [7c. Generate Design Matrix](#7c-generate-design-matrix)
      - [7d. Perform Individual Probe Level DE](#7d-perform-individual-probe-level-de)
      - [7e. Add Additional Columns and Format DE Table](#7e-add-additional-columns-and-format-de-table)

---

# Software used  

|Program|Version|Relevant Links|
|:------|:------:|:-------------|
|R|4.1.3|[https://www.bioinformatics.babraham.ac.uk/projects/fastqc/](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)|
|DT|0.26|[https://github.com/rstudio/DT](https://github.com/rstudio/DT)|
|dplyr|1.0.10|[https://dplyr.tidyverse.org/](https://dplyr.tidyverse.org/)|
|stringr|1.5.0|[https://stringr.tidyverse.org](https://stringr.tidyverse.org)|
|R.utils|2.12.2|[https://github.com/HenrikBengtsson/R.utils](https://github.com/HenrikBengtsson/R.utils)|
|limma|3.50.3|[https://bioconductor.org/packages/3.14/bioc/html/limma.html](https://bioconductor.org/packages/3.14/bioc/html/limma.html)|
|glue|R1.6.2|[https://glue.tidyverse.org/](https://glue.tidyverse.org/)|
|biomaRt|2.50.0|[https://bioconductor.org/packages/3.14/bioc/html/biomaRt.html](https://bioconductor.org/packages/3.14/bioc/html/biomaRt.html)|
|scales|1.2.1|[https://scales.r-lib.org/](https://scales.r-lib.org/)|
|dp_tools|1.2.0|[https://github.com/J-81/dp_tools](https://github.com/J-81/dp_tools)|
|singularity|3.9|[https://sylabs.io/](https://sylabs.io/)|
|Quarto|1.1.251|[https://quarto.org/](https://quarto.org/)|


---

# General processing overview with example commands  

> Exact processing commands for specific datasets that have been released are provided with their processed data in the [GLDS repository](https://genelab-data.ndc.nasa.gov/genelab/projects).
> 
> All output files in **bold** are published with the microarray processed data in the [GLDS repository](https://genelab-data.ndc.nasa.gov/genelab/projects). 

---

## 1. Create Sample RunSheet

> Notes: 
> - Rather than running the commands below to create the runsheet needed for processing, the runsheet may also be created manually by following the [file specification](../Workflow_Documentation/NF_Agile1CMP-A/examples/README.md).
> 
> - These command line tools are part of the [dp_tools](https://github.com/J-81/dp_tools) program.

```bash
### Download the *ISA.zip file from the GeneLab Repository ###

dpt-get-isa-archive \
 --accession OSD-###

### Parse the metadata from the *ISA.zip file to create a sample runsheet ###

dpt-isa-to-runsheet --accession OSD-### \
 --config-type microarray \
 --config-version Latest \
 --isa-archive *ISA.zip
```

**Parameter Definitions:**

- `--accession OSD-###` - OSD accession ID (replace ### with the OSD number being processed), used to retrieve the urls for the ISA archive and raw expression files hosted on the GeneLab Repository
- `--config-type` - Instructs the script to extract the metadata required for `microarray` processing from the ISA archive
- `--config-version` - Specifies the `dp-tools` configuration version to use, a value of `Latest` will specify the most recent version
- `--isa-archive` - Specifies the *ISA.zip file for the respective GLDS dataset, downloaded in the `dpt-get-isa-archive` command


**Input Data:**

- No input data required but the OSD accession ID needs to be indicated, which is used to download the respective ISA archive 

**Output Data:**

- *ISA.zip (compressed ISA directory containing Investigation, Study, and Assay (ISA) metadata files for the respective OSD dataset, used to define sample groups - the *ISA.zip file is located in the [OSD repository](https://osdr.nasa.gov/bio/repo/search?q=&data_source=cgene,alsda&data_type=study) under 'Study Files' -> 'metadata')

- {OSD-Accession-ID}_microarray_v{version}_runsheet.csv\# (table containing metadata required for processing, version denotes the dp_tools schema used to specify the metadata to extract from the ISA archive)

<br>

---

## 2. Load Metadata and Raw Data

> Note: Steps 2 - 7 are done in R

```R
# fileEncoding removes strange characters from the column names
df_rs <- read.csv(params$runsheet, check.names = FALSE, fileEncoding = 'UTF-8-BOM') 


# Define paths to raw data files
runsheetPathsAreURIs <- function(df_runsheet) {
  allTrue(stringr::str_starts(df_runsheet$`Array Data File Path`, "https"))
}


# Download raw data files
downloadFilesFromRunsheet <- function(df_runsheet) {
  urls <- df_runsheet$`Array Data File Path`
  destinationFiles <- df_runsheet$`Array Data File Name`

  mapply(function(url, destinationFile) {
    print(paste0("Downloading from '", url, "' TO '", destinationFile, "'"))
    if ( file.exists(destinationFile ) ) {
      warning(paste( "Using Existing File:", destinationFile ))
    } else {
      download.file(url, destinationFile)
    }
  }, urls, destinationFiles)

  destinationFiles # Return these paths
}

if ( runsheetPathsAreURIs(df_rs) ) {
  print("Determined Raw Data Locations are URIS")
  local_paths <- downloadFilesFromRunsheet(df_rs)
} else {
  print("Or Determined Raw Data Locations are local paths")
  local_paths <- df_rs$`Array Data File Path`
}


# Decompress files if needed
if ( allTrue(stringr::str_ends(local_paths, ".gz")) ) {
  print("Determined these files are gzip compressed... Decompressing now")
  # This does the decompression
  lapply(local_paths, R.utils::gunzip, remove = FALSE, overwrite = TRUE)
  # This removes the .gz extension to get the decompressed filenames
  local_paths <- vapply(local_paths, 
                        stringr::str_replace, # Run this function against each item in 'local_paths'
                        FUN.VALUE = character(1),  # Execpt an character vector as a return
                        USE.NAMES = FALSE,  # Don't use the input to assign names for the returned list
                        pattern = ".gz$", # first argument for applied function
                        replacement = ""  # second argument for applied function
                        )
}

df_local_paths <- data.frame(`Sample Name` = df_rs$`Sample Name`, `Local Paths` = local_paths, check.names = FALSE)


# Load raw data into R object
raw_data <- limma::read.maimages(df_local_paths$`Local Paths`, 
                                 source = "agilent",  # Specify platform
                                 green.only = TRUE, # Specify one-channel design
                                 names = df_local_paths$`Sample Name` # Map column names as Sample Names (instead of default filenames)
                                 )


# Summarize raw data
print(paste0("Number of Arrays: ", dim(raw_data)[2]))
print(paste0("Number of Probes: ", dim(raw_data)[1]))
```

**Input Data:**

- `params$runsheet` (Path to runsheet, output from [Step 1](#1-create-sample-runsheet))

**Output Data:**

- `raw_data` (R object containing raw microarray data)

    > Note: The raw data R object will be used to generate quality assessment (QA) plots in the next step.

<br>

---

## 3. Raw Data Quality Assessment

<br>

### 3a. Density Plot

```R
limma::plotDensities(raw_data, 
                     log = TRUE, 
                     legend = "topright")
```

**Input Data:**

- `raw_data` (raw data R object created in [Step 2](#2-load-metadata-and-raw-data) above)

**Output Data:**

- Plot containing the density of raw intensities for each array (raw intensity values with background intensity values subtracted; lack of overlap indicates a need for normalization)

<br>

### 3b. Pseudo Image Plots

```R
agilentImagePlot <- function(eListRaw) {
  # Adapted from this discussion: https://support.bioconductor.org/p/15523/
  copy_raw_data <- eListRaw
  copy_raw_data$genes$Block <- 1 # Agilent arrays only have one block
  names(copy_raw_data$genes)[2] <- "Column"
  copy_raw_data$printer <- limma::getLayout(copy_raw_data$genes)

  r <- copy_raw_data$genes$Row
  c <- copy_raw_data$genes$Column
  nr <- max(r)
  nc <- max(c)
  y <- rep(NA,nr*nc)
  i <- (r-1)*nc+c
  for ( array_i in seq(colnames(copy_raw_data$E)) ) {
    y[i] <- log2(copy_raw_data$E[,array_i])
    limma::imageplot(y,copy_raw_data$printer, main = rownames(copy_raw_data$targets)[array_i])
  }
}

agilentImagePlot(raw_data)
```

**Input Data:**

- `raw_data` (raw data R object created in [Step 2](#2-load-metadata-and-raw-data) above)

**Output Data:**

- Psuedo images of each array before background correction and normalization

<br>

### 3c. MA Plots

```R
for ( array_i in seq(colnames(raw_data$E)) ) {
  sample_name <- rownames(raw_data$targets)[array_i]
  limma::plotMA(raw_data,array=array_i,xlab="Average log-expression",ylab="Expression log-ratio (this sample vs. others)", main = sample_name, status=raw_data$genes$ControlType)
}
```

**Input Data:**

- `raw_data` (raw data R object created in [Step 2](#2-load-metadata-and-raw-data) above)

**Output Data:**

- M (log ratio) vs. A (average log expression) plot for each array before background correction and normalization (negative and positive control probes are in green and red, respectively)

<br>

### 3d. Foreground-Background Plots

```R
for ( array_i in seq(colnames(raw_data$E)) ) {
  sample_name <- rownames(raw_data$targets)[array_i]
  limma::plotFB(raw_data, array = array_i, xlab = "log2 Background", ylab = "log2 Foreground", main = sample_name) 
}
```

**Input Data:**

- `raw_data` (raw data R object created in [Step 2](#2-load-metadata-and-raw-data) above)

**Output Data:**

- Foreground vs. background expression plot for each array before background correction and normalization 

<br>

### 3e. Boxplots

```R
boxplotExpressionSafeMargin <- function(data) {
  longest_sample_name_length <- max(nchar(rownames(data$targets))) * 1
  bottom_margin <- min(35, longest_sample_name_length)
  par(mar=c(bottom_margin,2,1,1))
  boxplot(log2(data$E), las=2)
}

boxplotExpressionSafeMargin(raw_data)
```

**Input Data:**

- `raw_data` (raw data R object created in [Step 2](#2-load-metadata-and-raw-data) above)

**Output Data:**

- Boxplot of raw expression data for each array before background correction and normalization 

<br>

---

## 4. Background Correction

```R
norm_data <- limma::backgroundCorrect(raw_data, method = "normexp")
```

**Input Data:**

- `raw_data` (raw data R object created in [Step 2](#2-load-metadata-and-raw-data) above)

**Output Data:**

- `norm_data` (R object containing background-corrected microarray data)

  >   
  > Note: Background correction was performed using the `normexp` method, TODO: insert description of what this does

<br>

---

## 5. Between Array Normalization

```R
# Normalize background-corrected data using the quantile method
norm_data <- limma::normalizeBetweenArrays(norm_data, method = "quantile")

# Summarize background-corrected and normalized data
print(paste0("Number of Arrays: ", dim(norm_data)[2]))
print(paste0("Number of Probes: ", dim(norm_data)[1]))
```

**Input Data:**

- `norm_data` (R object containing background-corrected microarray data created in [Step 4](#4-background-correction) above)

**Output Data:**

- `norm_data` (R object containing background-corrected and normalized microarray data)
 
  >   
  > Note: Background correction was performed using the `quantile` method, TODO: insert description of what this does

<br>

---

## 6. Normalized Data Quality Assessment

<br>

### 6a. Density Plot

```R
limma::plotDensities(norm_data, 
                     log = TRUE, 
                     legend = "topright")
```

**Input Data:**

- `norm_data` (R object containing background-corrected and normalized microarray data created in [Step 5](#5-between-array-normalization) above)

**Output Data:**

- Plot containing the density of background-corrected and normalized intensities for each array (near complete overlap is expected after normalization)

<br>

### 6b. Pseudo Image Plots

```R
agilentImagePlot(norm_data)
```

**Input Data:**

- `norm_data` (R object containing background-corrected and normalized microarray data created in [Step 5](#5-between-array-normalization) above)

**Output Data:**

- Psuedo images of each array after background correction and normalization

<br>

### 6c. MA Plots

```R
for ( array_i in seq(colnames(norm_data$E)) ) {
  sample_name <- rownames(norm_data$targets)[array_i]
  limma::plotMA(norm_data,array=array_i,xlab="Average log-expression",ylab="Expression log-ratio (this sample vs. others)", main = sample_name, status=norm_data$genes$ControlType)
}
```

**Input Data:**

- `norm_data` (R object containing background-corrected and normalized microarray data created in [Step 5](#5-between-array-normalization) above)

**Output Data:**

- M (log ratio) vs. A (average log expression) plot for each array after background correction and normalization (negative and positive control probes are in green and red, respectively)

<br>

### 6d. Boxplots

```R
boxplotExpressionSafeMargin(norm_data)
```

**Input Data:**

- `norm_data` (R object containing background-corrected and normalized microarray data created in [Step 5](#5-between-array-normalization) above)

**Output Data:**

- Boxplot of expression data for each array after background correction and normalization 

<br>

---

## 7. Probeset Differential Expression

<br>

### 7a. Add Probeset Annotations

```R
shortenedOrganismName <- function(long_name) {
  #' Convert organism names like 'Homo Sapiens' into 'hsapiens'
  tokens <- long_name %>% stringr::str_split(" ", simplify = TRUE)
  genus_name <- tokens[1]

  species_name <- tokens[2]

  short_name <- stringr::str_to_lower(paste0(substr(genus_name, start = 1, stop = 1), species_name))

  return(short_name)
}


# locate dataset
expected_dataset_name <- shortenedOrganismName(unique(df_rs$organism)) %>% stringr::str_c("_gene_ensembl")
print(paste0("Expected dataset name: '", expected_dataset_name, "'"))


# Specify Ensembl version used in current GeneLab reference annotations
ENSEMBL_VERSION <- '107'

ensembl <- biomaRt::useEnsembl(biomart = "genes", 
                               dataset = expected_dataset_name,
                               version = ENSEMBL_VERSION)
print(ensembl)


getBioMartAttribute <- function(df_rs, params) {
  #' Returns resolved biomart attribute

  # check if runsheet has Array Design REF
  if ( !is.null(df_rs$`Array Design REF`) ) {
    print("Using attribute name sourced from runsheet")
    return(unique(df_rs$`Array Design REF`))
  } else {
    print("Could not find 'Array Design REF' in runsheet, falling back to parameters")
  }

  # check if a fallback has been given via params
  if ( !is.null(params$biomart_attribute ) ) {
    print("Using attribute name sourced from parameters")
    return(params$biomart_attribute)
  }

  # finally throw an error if neither guard condition was true
  stop("No valid biomart attribute identified")
}

expected_attribute_name <- getBioMartAttribute(df_rs, params)
print(paste0("Expected attribute name: '", expected_attribute_name, "'"))

probe_ids <- norm_data$genes$ProbeName

# DEBUG clause
if ( is.integer(params$DEBUG_limit_biomart_query) ) {
  warning(paste("DEBUG MODE: Limiting query to", params$DEBUG_limit_biomart_query, "entries"))
  probe_ids <- probe_ids[1:params$DEBUG_limit_biomart_query]
}


# Create probe map
df_mapping <- biomaRt::getBM(
    attributes = c(
        expected_attribute_name,
        "ensembl_gene_id",
        "uniprot_gn_symbol",
        #"go_id", # Full GO IDS
        "goslim_goa_accession" # GO SLIM IDS, we use these to be consistent with bulkRNASeq and other transcriptomics assays
        ), filters = expected_attribute_name, 
        values = c(probe_ids), 
        mart = ensembl)


# Convert list of multi-mapped genes to string
listToUniquePipedString <- function(str_list) {
  #! convert lists into strings denoting unique elements separated by '|' characters
  #! e.g. c("GO1","GO2","GO2","G03") -> "GO1|GO2|GO3"
  return(toString(unique(str_list)) %>% stringr::str_replace_all(pattern = stringr::fixed(", "), replacement = "|"))
}

unique_probe_ids <- df_mapping %>% 
                      dplyr::group_by(!!sym(expected_attribute_name)) %>% 
                      dplyr::summarise(
                        SYMBOL = listToUniquePipedString(uniprot_gn_symbol),
                        ENSEMBL = listToUniquePipedString(ensembl_gene_id),
                        GOSLIM_IDS = listToUniquePipedString(goslim_goa_accession)
                        ) %>%
                      # Count number of ensembl IDS mapped
                      dplyr::mutate( 
                        count_ENSEMBL_mappings = 1 + stringr::str_count(ENSEMBL, ",")
                      )

norm_data$genes <- norm_data$genes %>% 
  dplyr::left_join(unique_probe_ids, by = c("ProbeName" = expected_attribute_name ) ) %>%
  dplyr::mutate( count_ENSEMBL_mappings = ifelse(is.na(ENSEMBL), 0, count_ENSEMBL_mappings) )
```

**Input Data:**

- `df_rs$organism` (organism specified in the runsheet created in [Step 1](#1-create-sample-runsheet))
- `df_rs$'Array Design REF'` (array design reference specified in the runsheet created in [Step 1](#1-create-sample-runsheet))
- Ensembl version (reference organism Ensembl version indicated in the `ensemblVersion` column of the [GL-DPPD-7110_annotations.csv](../../GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110/GL-DPPD-7110_annotations.csv) GeneLab Annotations file)
- `probe_ids` (probe IDs for each probe associated with the `norm_data` R object containing background-corrected and normalized microarray data created in [Step 5](#5-between-array-normalization))

**Output Data:**

- `df_mapping` (R object containing a map of each probe ID to the respective Biomart Ensembl ID(s), gene symbol(s) and GO slim ID(s)) 
- `unique_probe_ids` (R object containing probe IDs that map uniquely to an Ensembl ID and feature)
- Should this just be `norm_data$genes`???

<br>

### 7b. Summarize Biomart Mapping vs. Manufacturer Mapping

```R
describeMapping <- function(stats) {
  my_label <- scales::label_percent(accuracy = 0.01, prefix="(", suffix = "%)")
  percentOfUniqueProbes <- function(num) {
    #' Gives number as a percent of unique probes 
    
    return(my_label(num / stats$count_unique_probe))
  }

  print(paste("  Mapped (original): ", stats$count_total_original_mapped, percentOfUniqueProbes(stats$count_total_original_mapped)))
  print(paste("  Mapped (biomart): ", stats$count_total_biomart_mapped, percentOfUniqueProbes(stats$count_total_biomart_mapped)))
  print(paste("  Ensembl One-To-One Mapped (biomart): ", stats$count_total_biomart_1to1_mapped, percentOfUniqueProbes(stats$count_total_biomart_1to1_mapped)))
  print(paste("  Ensembl Multi-Mapped (biomart): ", stats$count_total_biomart_multi_mapped, percentOfUniqueProbes(stats$count_total_biomart_multi_mapped)))
  print(paste("  Not mapped (unique original): ", stats$count_unique_original_unmapped, percentOfUniqueProbes(stats$count_unique_original_unmapped)))
  print(paste("  Not mapped (unique biomart): ", stats$count_unique_biomart_unmapped, percentOfUniqueProbes(stats$count_unique_biomart_unmapped)))
  print(paste("  Not mapped (shared biomart & original): ", stats$count_both_original_and_biomart_unmapped, percentOfUniqueProbes(stats$count_both_original_and_biomart_unmapped)))

  print("More Info")
  print(paste("  Not mapped (total original mapping): ", stats$count_total_original_unmapped, percentOfUniqueProbes(stats$count_total_original_unmapped)))
  print(paste("  Not mapped (total biomart mapping): ", stats$count_total_biomart_unmapped, percentOfUniqueProbes(stats$count_total_biomart_unmapped)))
  print(paste("  Features (genes) (stats$count_biomart_unique_mapped_features): ", stats$count_biomart_unique_mapped_features))
  print(paste("  Features (genes) (stats$count_original_unique_mapped_features): ", stats$count_original_unique_mapped_features))
  print(paste("  Percent fewer unique features: ", my_label(-(stats$count_biomart_unique_mapped_features - stats$count_original_unique_mapped_features) / stats$count_original_unique_mapped_features)))
}

calculateMappingStats <- function(df_genes) {

  stats <- list()

  df_genes.no.control.probes <- df_genes %>% dplyr::filter( ControlType == 0 )
  print(paste("Original Probe Count: ", length(df_genes$ProbeName)))
  print(paste("Original Non-Control Probe Count: ", length(df_genes.no.control.probes$ProbeName)))
  stats$count_unique_probe <- length(unique(df_genes.no.control.probes$ProbeName))

  stats$count_total_original_unmapped <- length(
    unique(
      df_genes.no.control.probes %>%
      dplyr::filter( (ProbeName == SystematicName) ) %>%
      dplyr::pull( ProbeName )
    )
  )
  stats$count_total_original_mapped <- length(
    unique(
      df_genes.no.control.probes %>%
      dplyr::filter( (ProbeName != SystematicName) ) %>%
      dplyr::pull( ProbeName )
    )
  )
  stats$count_total_biomart_unmapped <- length(
    unique(
      df_genes.no.control.probes %>%
      dplyr::filter( is.na(SYMBOL) ) %>%
      dplyr::pull( ProbeName )
    )
  )
  stats$count_total_biomart_mapped <- length(
    unique(
      df_genes.no.control.probes %>%
      dplyr::filter( !is.na(SYMBOL) ) %>%
      dplyr::pull( ProbeName )
    )
  )
  stats$count_total_biomart_1to1_mapped <- length(
    unique(
      df_genes.no.control.probes %>%
      dplyr::filter( !is.na(SYMBOL) & (count_ENSEMBL_mappings == 1) ) %>%
      dplyr::pull( ProbeName )
    )
  )
  stats$count_total_biomart_multi_mapped <- length(
    unique(
      df_genes.no.control.probes %>%
      dplyr::filter( !is.na(SYMBOL) & (count_ENSEMBL_mappings > 1) ) %>%
      dplyr::pull( ProbeName )
    )
  )
  stats$count_unique_original_unmapped <- length(
    unique(
      df_genes.no.control.probes %>%
      dplyr::filter( (ProbeName == SystematicName) & !is.na(SYMBOL) ) %>%
      dplyr::pull( ProbeName )
    )
  )
  stats$count_unique_biomart_unmapped <- length(
    unique(
      df_genes.no.control.probes %>%
      dplyr::filter( ( count_ENSEMBL_mappings == 0 ) & (ProbeName != SystematicName) ) %>%
      dplyr::pull( ProbeName )
    )
  )
  stats$count_both_original_and_biomart_unmapped <- length(
    unique(
      df_genes.no.control.probes %>%
      dplyr::filter( is.na(SYMBOL) & (ProbeName == SystematicName) ) %>%
      dplyr::pull( ProbeName )
    )
  )

  stats$count_biomart_unique_mapped_features <- df_genes.no.control.probes %>%
      dplyr::filter( count_ENSEMBL_mappings == 1  ) %>%
      dplyr::pull( ENSEMBL ) %>%
      unique() %>%
      length()

  stats$count_original_unique_mapped_features <- df_genes.no.control.probes %>%
      dplyr::filter( ProbeName != GeneName  ) %>%
      dplyr::pull( GeneName ) %>%
      unique() %>%
      length()

  return(stats)
}


calculateMappingStats(norm_data$genes) %>% describeMapping()
```

**Input Data:**

- ???

**Output Data:**

- Indicate what's included in the mapping summary - let's go over each of these.

<br>

### 7c. Generate Design Matrix

```R
# Pull all factors for each sample in the study from the runsheet created in Step 1
runsheetToDesignMatrix <- function(runsheet_path) {
    df = read.csv(runsheet_path)
    # get only Factor Value columns
    factors = as.data.frame(df[,grep("Factor.Value", colnames(df), ignore.case=TRUE)])
    colnames(factors) = paste("factor",1:dim(factors)[2], sep= "_")
    
    # Load metadata from runsheet csv file
    compare_csv = data.frame(sample_id = df[,c("Sample.Name")], factors)

    # Create data frame containing all samples and respective factors
    study <- as.data.frame(compare_csv[,2:dim(compare_csv)[2]])
    colnames(study) <- colnames(compare_csv)[2:dim(compare_csv)[2]]
    rownames(study) <- compare_csv[,1] 
    
    # Format groups and indicate the group that each sample belongs to
    if (dim(study)[2] >= 2){
        group<-apply(study,1,paste,collapse = " & ") # concatenate multiple factors into one condition per sample
    } else{
        group<-study[,1]
    }
    group_names <- paste0("(",group,")",sep = "") # human readable group names
    group <- sub("^BLOCKER_", "",  make.names(paste0("BLOCKER_", group))) # group naming compatible with R models, this maintains the default behaviour of make.names with the exception that 'X' is never prepended to group namesnames(group) <- group_names
    names(group) <- group_names

    # Format contrasts table, defining pairwise comparisons for all groups
    contrast.names <- combn(levels(factor(names(group))),2) # generate matrix of pairwise group combinations for comparison
    contrasts <- apply(contrast.names, MARGIN=2, function(col) sub("^BLOCKER_", "",  make.names(paste0("BLOCKER_", stringr::str_sub(col, 2, -2)))))
    contrast.names <- c(paste(contrast.names[1,],contrast.names[2,],sep = "v"),paste(contrast.names[2,],contrast.names[1,],sep = "v")) # format combinations for output table files names
    contrasts <- cbind(contrasts,contrasts[c(2,1),])
    colnames(contrasts) <- contrast.names
    sampleTable <- data.frame(condition=factor(group))
    rownames(sampleTable) <- df[,c("Sample.Name")]

    condition <- sampleTable[,'condition']
    names_mapping <- as.data.frame(cbind(safe_name = as.character(condition), original_name = group_names))

    design <- model.matrix(~ 0 + condition)
    design_data <- list( matrix = design, mapping = names_mapping, groups = as.data.frame( cbind(sample = df[,c("Sample.Name")], group = group_names) ) )
    return(design_data)
}


# Loading metadata from runsheet csv file
design_data <- runsheetToDesignMatrix(params$runsheet)
design <- design_data$matrix
```

**Input Data:**

- `params$runsheet` (Path to runsheet, output from [Step 1](#1-create-sample-runsheet))

**Output Data:**

- `design` (R object containing the study design matrix, indicating the group that each sample belongs to)

<br>

### 7d. Perform Individual Probe Level DE

```R
lmFitPairwise <- function(norm_data, design) {
    #' Perform all pairwise comparisons

    #' Approach based on limma manual section 17.4 (version 3.52.4)

    fit <- limma::lmFit(norm_data, design)

    # Create Contrast Model
    fit.groups <- colnames(fit$design)[which(fit$assign == 1)]
    combos <- combn(fit.groups,2)
    contrasts<-c(paste(combos[1,],combos[2,],sep = "-"),paste(combos[2,],combos[1,],sep = "-")) # format combinations for limma:makeContrasts
    cont.matrix <- limma::makeContrasts(contrasts=contrasts,levels=design)
    contrast.fit <- limma::contrasts.fit(fit, cont.matrix)

    contrast.fit <- limma::eBayes(contrast.fit,trend=TRUE,robust=TRUE)
    return(contrast.fit)
}

# Calculate results
res <- lmFitPairwise(norm_data, design)

# Print DE table, without filtering
NO_FILTER_DGE = paste0(params$id, "_INTERIM_no_filtering.csv")
limma::write.fit(res, adjust = 'BH', 
                file = NO_FILTER_DGE,
                row.names = FALSE,
                quote = TRUE,
                sep = ",")
```

**Input Data:**

- `norm_data` (R object containing background-corrected and normalized microarray data created in [Step 5](#5-between-array-normalization))
- `design` (R object containing the study design matrix, indicating the group that each sample belongs to, created in [Step 7c](#7c-generate-design-matrix) above)

**Output Data:**

- `*_INTERIM_no_filtering.csv` (???)

<br>

### 7e. Add Additional Columns and Format DE Table 

```R
## Reformat Table for consistency across DE analyses tables within GeneLab ##

# Read in DE table 
df_interim <- read.csv(NO_FILTER_DGE)

# Reformat column names
reformat_names <- function(colname, group_name_mapping) {
  new_colname <- colname  %>% 
                  stringr::str_replace(pattern = "^P.value.adj.condition", replacement = "Adj.p.value_") %>%
                  stringr::str_replace(pattern = "^P.value.condition", replacement = "P.value_") %>%
                  stringr::str_replace(pattern = "^Coef.condition", replacement = "Log2fc_") %>% # This is the Log2FC as per: https://rdrr.io/bioc/limma/man/writefit.html
                  stringr::str_replace(pattern = "^t.condition", replacement = "T.stat_") %>%
                  stringr::str_replace(pattern = stringr::fixed("Genes.ProbeUID"), replacement = "PROBEID") %>% 
                  stringr::str_replace(pattern = stringr::fixed("Genes.SYMBOL"), replacement = "SYMBOL") %>% 
                  stringr::str_replace(pattern = stringr::fixed("Genes.ENSEMBL"), replacement = "ENSEMBL") %>% 
                  stringr::str_replace(pattern = stringr::fixed("Genes.GOSLIM_IDS"), replacement = "GOSLIM_IDS") %>% 
                  stringr::str_replace(pattern = ".condition", replacement = "v")
  
  # remap to group names before make.names was applied
  for ( i in seq(nrow(group_name_mapping)) ) {
    safe_name <- group_name_mapping[i,]$safe_name
    original_name <- group_name_mapping[i,]$original_name
    new_colname <- new_colname %>% stringr::str_replace(pattern = stringr::fixed(safe_name), replacement = original_name)
  }

  return(new_colname)
}

df_interim <- df_interim %>% dplyr::rename_with( reformat_names, group_name_mapping = design_data$mapping )


# Concatenate expression values for each sample
df_interim <- df_interim %>% dplyr::bind_cols(norm_data$E)


## Add Group Wise Statistics ##

# Group mean and standard deviations for normalized expression values are computed and added to the table

unique_groups <- unique(design_data$group$group)
for ( i in seq_along(unique_groups) ) {
  current_group <- unique_groups[i]
  current_samples <- design_data$group %>% 
                      dplyr::group_by(group) %>%
                      dplyr::summarize(
                        samples = sort(unique(sample))
                      ) %>%
                      dplyr::filter(
                        group == current_group
                      ) %>% 
                      dplyr::pull()
  #### THIS IS DIFFERENT FROM HTML ####                    
  print(glue::glue("Computing mean and standard deviation for Group {i} of {length(unique_groups)}"))
  print(glue::glue("Group: {current_group}"))
  print(glue::glue("Samples in Group: '{toString(current_samples)}'"))
  
  df_interim <- df_interim %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate( 
      "Group.Mean_{current_group}" := mean(c_across(current_samples)),
      "Group.Stdev_{current_group}" := sd(c_across(current_samples)),
      ) %>% 
    dplyr::ungroup() %>%
    as.data.frame()
}

print("Remove extra columns from final table")

# These columns are data mapped to column PROBEID as per the original Manufacturer and can be linked as needed
colnames_to_remove = c(
  "Genes.Row",
  "Genes.Col",
  "Genes.Start",
  "Genes.Sequence",
  "Genes.ControlType",
  "Genes.ProbeName",
  "Genes.GeneName",
  "Genes.SystematicName",
  "Genes.Description"
  # "Genes.count_ENSEMBL_mappings", Keep this
)

df_interim <- df_interim %>% dplyr::select(-any_of(colnames_to_remove))

# Save to file
write.csv(df_interim, paste0(params$id, "_differential_expression.csv"), row.names = FALSE)
```

**Input Data:**

- `*_INTERIM_no_filtering.csv` (???, output from [Step 7d](#7d-perform-individual-probe-level-de) above)

**Output Data:**

- SampleTable.csv\# (table containing samples and their respective groups - ### I don't see this outputted here and shoudn't it be outputted in 7c? ###)
- contrasts.csv\# (table containing all pairwise comparisons - ### I don't see this outputted here and shoudn't it be outputted in 7c? ###)
- visualization_output_table.csv (file used to generate GeneLab DGE visualizations - ### I don't see where this outputted here ###)
- visualization_PCA_table.csv (file used to generate GeneLab PCA plots - ### I don't see where this calculated or outputted here ###)
- differential_expression.csv\# (table containing normalized counts for each sample, group statistics, Limma probe DE results for each pairwise comparison, and gene annotations)
