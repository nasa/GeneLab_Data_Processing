# GeneLab bioinformatics processing pipeline for Agilent 1-channel microarray data <!-- omit in toc -->

> **This page holds an overview and instructions for how GeneLab processes Agilent 1-channel microarray datasets\*. Exact processing commands and GL-DPPD-7112 version used for specific datasets are provided with their processed data in the [GeneLab Data Systems 
(GLDS) repository](https://genelab-data.ndc.nasa.gov/genelab/projects).**  
>
> \* The pipeline detailed below is currently used for animal studies only, it will be updated soon for processing plants and microbe microarray data.
---

**Date:** February 13, 2023  
**Revision:** -  
**Document Number:** GL-DPPD-7112   

**Submitted by:**  
Jonathan Oribello (GeneLab Data Processing Team)

**Approved by:**  
Sylvain Costes (GeneLab Project Manager)  
Samrawit Gebre (GeneLab Deputy Project Manager)  
Amanda Saravia-Butler (GeneLab Data Processing Lead)  
Lauren Sanders (acting GeneLab Project Scientist)  

---

# Table of contents <!-- omit in toc -->

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
    - [7. Probe Differential Expression (DE)](#7-probe-differential-expression-de)
      - [7a. Add Probe Annotations](#7a-add-probe-annotations)
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
|dplyr|1.0.10|[https://dplyr.tidyverse.org](https://dplyr.tidyverse.org)|
|stringr|1.5.0|[https://stringr.tidyverse.org](https://stringr.tidyverse.org)|
|R.utils|2.12.2|[https://github.com/HenrikBengtsson/R.utils](https://github.com/HenrikBengtsson/R.utils)|
|limma|3.50.3|[https://bioconductor.org/packages/3.14/bioc/html/limma.html](https://bioconductor.org/packages/3.14/bioc/html/limma.html)|
|glue|1.6.2|[https://glue.tidyverse.org](https://glue.tidyverse.org)|
|ggplot2|3.4.0|[https://ggplot2.tidyverse.org](https://ggplot2.tidyverse.org)|
|biomaRt|2.50.0|[https://bioconductor.org/packages/3.14/bioc/html/biomaRt.html](https://bioconductor.org/packages/3.14/bioc/html/biomaRt.html)|
|matrixStats|0.63.0|[https://github.com/HenrikBengtsson/matrixStats](https://github.com/HenrikBengtsson/matrixStats)|
|statmod|1.5.0|[https://github.com/cran/statmod](https://github.com/cran/statmod)|
|dp_tools|1.3.0|[https://github.com/J-81/dp_tools](https://github.com/J-81/dp_tools)|
|singularity|3.9|[https://sylabs.io](https://sylabs.io)|
|Quarto|1.2.313|[https://quarto.org](https://quarto.org)|

---

# General processing overview with example commands  

> Exact processing commands for specific datasets that have been released are provided with their processed data in the [GLDS repository](https://genelab-data.ndc.nasa.gov/genelab/projects).
> 
> All output files in **bold** are published with the Agilent 1-channel microarray processed data in the [GLDS repository](https://genelab-data.ndc.nasa.gov/genelab/projects). 

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

- **{OSD-Accession-ID}_microarray_v{version}_runsheet.csv** (table containing metadata required for processing, version denotes the dp_tools schema used to specify the metadata to extract from the ISA archive)

<br>

---

## 2. Load Metadata and Raw Data

> Note: Steps 2 - 7 are done in R

```R
### Install R packages if not already installed ###

install.packages("tidyverse")
install.packages("R.utils")
install.packages("glue")
install.packages("matrixStats")
install.packages("statmod")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.14")
BiocManager::install("limma")
BiocManager::install("biomaRt")


### Import libraries ###

library(tidyverse)
library(dplyr)
library(stringr)
library(R.utils)
library(glue)
library(matrixStats)
library(limma)
library(biomaRt)
library(statmod)



# Define path to runsheet
runsheet <- "/path/to/runsheet/{OSD-Accession-ID}_microarray_v{version}_runsheet.csv"

## Set up output structure

# Output Constants
DIR_RAW_DATA <- "00-RawData"
DIR_NORMALIZED_EXPRESSION <- "01-limma_NormExp"
DIR_DGE <- "02-limma_DGE"

dir.create(DIR_RAW_DATA)
dir.create(DIR_NORMALIZED_EXPRESSION)
dir.create(DIR_DGE)

# fileEncoding removes strange characters from the column names
df_rs <- read.csv(runsheet, check.names = FALSE, fileEncoding = 'UTF-8-BOM') 

## Determines the organism specific annotation file to use based on the organism in the runsheet
fetch_organism_specific_annotation_file_path <- function(organism) {
  # Uses the GeneLab GL-DPPD-7110_annotations.csv file to find the organism specific annotation file path
  # Raises an exception if the organism does not have an associated annotation file yet
  

  all_organism_table <- read.csv("https://raw.githubusercontent.com/nasa/GeneLab_Data_Processing/GL_RefAnnotTable_1.0.0/GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110/GL-DPPD-7110_annotations.csv")

  annotation_file_path <- all_organism_table %>% dplyr::filter(species == organism) %>% dplyr::pull(genelab_annots_link)

  # Guard clause: Ensure annotation_file_path populated
  # Else: raise exception for unsupported organism
  if (length(annotation_file_path) == 0) {
    stop(glue::glue("Organism supplied '{organism}' is not supported. See the following url for supported organisms: https://github.com/nasa/GeneLab_Data_Processing/blob/GL_RefAnnotTable_1.0.0/GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110/GL-DPPD-7110_annotations.csv.  Supported organisms will correspond to a row based on the 'species' column and include a url in the 'genelab_annots_link' column of that row"))
  }

  return(annotation_file_path)
}
annotation_file_path <- fetch_organism_specific_annotation_file_path(unique(df_rs$organism))

allTrue <- function(i_vector) {
  if ( length(i_vector) == 0 ) {
    stop(paste("Input vector is length zero"))
  }
  all(i_vector)
}

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


# uncompress files if needed
if ( allTrue(stringr::str_ends(local_paths, ".gz")) ) {
  print("Determined these files are gzip compressed... uncompressing now")
  # This does the uncompression
  lapply(local_paths, R.utils::gunzip, remove = FALSE, overwrite = TRUE)
  # This removes the .gz extension to get the uncompressed filenames
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

- `runsheet` (Path to runsheet, output from [Step 1](#1-create-sample-runsheet))

**Output Data:**

- `df_rs` (R dataframe containing information from the runsheet)
- `raw_data` (R object containing raw microarray data)

    > Note: The raw data R object will be used to generate quality assessment (QA) plots in the next step.

<br>

---

## 3. Raw Data Quality Assessment

<br>

### 3a. Density Plot

```R
# Plot settings
par(
  xpd = TRUE # Ensure legend can extend past plot area
)

number_of_sets = ceiling(dim(raw_data)[2] / 30) # Set of 30 samples, used to scale plot

limma::plotDensities(raw_data, 
                     log = TRUE, 
                     legend = FALSE)
legend("topright", legend = colnames(raw_data),
        lty = 1, # Solid line
        col = 1:ncol(raw_data), # Ensure legend color is in sync with plot
        ncol = number_of_sets, # Set number of columns by number of sets
        cex = max(0.5, 1 + 0.2 - (number_of_sets*0.2)) # Reduce scale by 20% for each column beyond 1, minimum of 0.5
      )
```

**Input Data:**

- `raw_data` (raw data R object created in [Step 2](#2-load-metadata-and-raw-data) above)

**Output Data:**

- Plot containing the density of raw intensities for each array (raw intensity values with background intensity values subtracted; lack of overlap indicates a need for normalization)

<br>

### 3b. Pseudo Image Plots

```R
agilentImagePlot <- function(eListRaw, transform_func = identity) {
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
    y[i] <- transform_func(copy_raw_data$E[,array_i])
    limma::imageplot(y,copy_raw_data$printer, main = rownames(copy_raw_data$targets)[array_i])
  }
}

agilentImagePlot(raw_data, transform_func = function(expression_matrix) log2(expression_matrix + 1))
```

**Input Data:**

- `raw_data` (raw data R object created in [Step 2](#2-load-metadata-and-raw-data) above)

**Output Data:**

- Pseudo images of each array before background correction and normalization 

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

- M (log ratio of the subject array vs a pseudo-reference, the mean of all other arrays) vs. A (average log expression) plot for each array before background correction and normalization (negative and positive control probes are in green and red, respectively)

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
boxplotExpressionSafeMargin <- function(data, transform_func = identity, xlab = "Log2 Intensity") {
  # Basic box plot
  df_data <- as.data.frame(transform_func(data$E))
  ggplot2::ggplot(stack(df_data), ggplot2::aes(x=values, y=ind)) + 
    ggplot2::geom_boxplot() + 
    ggplot2::scale_y_discrete(limits=rev) +
    ggplot2::labs(y= "Sample Name", x = xlab)
}

boxplotExpressionSafeMargin(raw_data, transform_func = log2)
```

**Input Data:**

- `raw_data` (raw data R object created in [Step 2](#2-load-metadata-and-raw-data) above)

**Output Data:**

- Boxplot of raw expression data for each array before background correction and normalization 

<br>

---

## 4. Background Correction

```R
background_corrected_data <- limma::backgroundCorrect(raw_data, method = "normexp")
```

**Input Data:**

- `raw_data` (raw data R object created in [Step 2](#2-load-metadata-and-raw-data) above)

**Output Data:**

- `background_corrected_data` (R object containing background-corrected microarray data)

  >   
  > Note: Background correction was performed using the `normexp` method as recommended by [Ritchie, M.E., et al.](http://bioinformatics.oxfordjournals.org/content/23/20/2700), which performs background correction and quantile normalization using the control probes by utilizing the `normexp.fit.control` function to estimate the parameters required by normal+exponential(normexp) convolution model with the help of negative control probes, followed by the `normexp.signal` function to perform the background correction.

<br>

---

## 5. Between Array Normalization

```R
# Normalize background-corrected data using the quantile method
norm_data <- limma::normalizeBetweenArrays(background_corrected_data, method = "quantile")

# Summarize background-corrected and normalized data
print(paste0("Number of Arrays: ", dim(norm_data)[2]))
print(paste0("Number of Probes: ", dim(norm_data)[1]))
```

**Input Data:**

- `background_corrected_data` (R object containing background-corrected microarray data created in [Step 4](#4-background-correction) above)

**Output Data:**

- `norm_data` (R object containing background-corrected and normalized microarray data)
 
  >   
  > Note: Normalization was performed using the `quantile` method, which forces the entire empirical distribution of all arrays to be identical

<br>

---

## 6. Normalized Data Quality Assessment

<br>

### 6a. Density Plot

```R
# Plot settings
par(
  xpd = TRUE # Ensure legend can extend past plot area
)

number_of_sets = ceiling(dim(norm_data)[2] / 30) # Set of 30 samples, used to scale plot

limma::plotDensities(norm_data, 
                     log = TRUE, 
                     legend = FALSE)
legend("topright", legend = colnames(norm_data),
        lty = 1, # Solid line
        col = 1:ncol(norm_data), # Ensure legend color is in sync with plot
        ncol = number_of_sets, # Set number of columns by number of sets
        cex = max(0.5, 1 + 0.2 - (number_of_sets*0.2)) # Reduce scale by 20% for each column beyond 1, minimum of 0.5
      )
```

**Input Data:**

- `norm_data` (R object containing background-corrected and normalized microarray data created in [Step 5](#5-between-array-normalization) above)

**Output Data:**

- Plot containing the density of background-corrected and normalized intensities for each array (near complete overlap is expected after normalization)

<br>

### 6b. Pseudo Image Plots

```R
agilentImagePlot(norm_data, 
                 transform_func = function(expression_matrix) log2(2**expression_matrix + 1) # Compute as log2 of normalized expression after adding a +1 offset to prevent negative values in the pseudoimage
                 )
```

**Input Data:**

- `norm_data` (R object containing background-corrected and normalized microarray data created in [Step 5](#5-between-array-normalization) above)

**Output Data:**

- Pseudo images of each array after background correction and normalization 

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

- M (log ratio of the subject array vs a pseudo-reference, the mean of all other arrays) vs. A (average log expression) plot for each array after background correction and normalization (negative and positive control probes are in green and red, respectively)

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

## 7. Probe Differential Expression (DE)

<br>

### 7a. Add Probe Annotations

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


getBioMartAttribute <- function(df_rs) {
  #' Returns resolved biomart attribute source from runsheet

  # check if runsheet has 'biomart_attribute' column
  if ( !is.null(df_rs$`biomart_attribute`) ) {
    print("Using attribute name sourced from runsheet")
    # Format according to biomart needs
    formatted_value <- unique(df_rs$`biomart_attribute`) %>% 
                        stringr::str_replace_all(" ","_") %>% # Replace all spaces with underscore
                        stringr::str_to_lower() # Lower casing only
    return(formatted_value)
  } else {
    stop("ERROR: Could not find 'biomart_attribute' in runsheet")
  }
}

expected_attribute_name <- getBioMartAttribute(df_rs)
print(paste0("Expected attribute name: '", expected_attribute_name, "'"))

probe_ids <- unique(norm_data$genes$ProbeName)


# Create probe map
# Run Biomart Queries in chunks to prevent request timeouts
#   Note: If timeout is occuring (possibly due to larger load on biomart), reduce chunk size
CHUNK_SIZE= 8000
probe_id_chunks <- split(probe_ids, ceiling(seq_along(probe_ids) / CHUNK_SIZE))
df_mapping <- data.frame()
for (i in seq_along(probe_id_chunks)) {
  probe_id_chunk <- probe_id_chunks[[i]]
  print(glue::glue("Running biomart query chunk {i} of {length(probe_id_chunks)}. Total probes IDS in query ({length(probe_id_chunk)})"))
  chunk_results <- biomaRt::getBM(
      attributes = c(
          expected_attribute_name,
          "ensembl_gene_id"
          ), 
          filters = expected_attribute_name, 
          values = probe_id_chunk, 
          mart = ensembl)

  df_mapping <- df_mapping %>% dplyr::bind_rows(chunk_results)
  Sys.sleep(10) # Slight break between requests to prevent back-to-back requests
}

# Convert list of multi-mapped genes to string
listToUniquePipedString <- function(str_list) {
  #! convert lists into strings denoting unique elements separated by '|' characters
  #! e.g. c("GO1","GO2","GO2","G03") -> "GO1|GO2|GO3"
  return(toString(unique(str_list)) %>% stringr::str_replace_all(pattern = stringr::fixed(", "), replacement = "|"))
}

unique_probe_ids <- df_mapping %>% 
                      dplyr::group_by(!!sym(expected_attribute_name)) %>% 
                      dplyr::summarise(
                        ENSEMBL = listToUniquePipedString(ensembl_gene_id)
                        ) %>%
                      # Count number of ensembl IDS mapped
                      dplyr::mutate( 
                        count_ENSEMBL_mappings = 1 + stringr::str_count(ENSEMBL, stringr::fixed("|"))
                      )

norm_data$genes <- norm_data$genes %>% 
  dplyr::left_join(unique_probe_ids, by = c("ProbeName" = expected_attribute_name ) ) %>%
  dplyr::mutate( count_ENSEMBL_mappings = ifelse(is.na(ENSEMBL), 0, count_ENSEMBL_mappings) )
```

**Input Data:**

- `df_rs$organism` (organism specified in the runsheet created in [Step 1](#1-create-sample-runsheet))
- `df_rs$'Array Design REF'` (array design reference specified in the runsheet created in [Step 1](#1-create-sample-runsheet))
- ENSEMBL_VERSION (reference organism Ensembl version indicated in the `ensemblVersion` column of the [GL-DPPD-7110_annotations.csv](../../GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110/GL-DPPD-7110_annotations.csv) GeneLab Annotations file)
- `norm_data$genes` (Manufacturer's probe metadata, including probe IDs and sequence position gene annotations associated with the `norm_data` R object containing background-corrected and normalized microarray data created in [Step 5](#5-between-array-normalization))

**Output Data:**

- `norm_data$genes` (Probe metadata, updated to include gene annotations specified by [Biomart](https://bioconductor.org/packages/3.14/bioc/html/biomaRt.html))

<br>

### 7b. Summarize Biomart Mapping vs. Manufacturer Mapping

```R
# Pie Chart with Percentages
slices <- c(
    'Control probes' = nrow(norm_data$gene %>% dplyr::filter(ControlType != 0) %>% dplyr::distinct(ProbeName)), 
    'Unique Mapping' = nrow(norm_data$gene %>% dplyr::filter(ControlType == 0) %>% dplyr::filter(count_ENSEMBL_mappings == 1) %>% dplyr::distinct(ProbeName)), 
    'Multi Mapping' = nrow(norm_data$gene %>% dplyr::filter(ControlType == 0) %>% dplyr::filter(count_ENSEMBL_mappings > 1) %>% dplyr::distinct(ProbeName)), 
    'No Mapping' = nrow(norm_data$gene %>% dplyr::filter(ControlType == 0) %>% dplyr::filter(count_ENSEMBL_mappings == 0) %>% dplyr::distinct(ProbeName))
)
pct <- round(slices/sum(slices)*100)
chart_names <- names(slices)
chart_names <- glue::glue("{names(slices)} ({slices})") # add count to labels
chart_names <- paste(chart_names, pct) # add percents to labels
chart_names <- paste(chart_names,"%",sep="") # ad % to labels
pie(slices,labels = chart_names, col=rainbow(length(slices)),
    main=glue::glue("Biomart Mapping to Ensembl Primary Keytype\n {nrow(norm_data$gene %>% dplyr::distinct(ProbeName))} Total Unique Probes")
    )

original_mapping_rate = nrow(norm_data$gene %>% dplyr::filter(ControlType == 0) %>% dplyr::filter(ProbeName != SystematicName) %>% dplyr::distinct(ProbeName))
print(glue::glue("Original Manufacturer Reported Mapping Count: {original_mapping_rate}"))
print(glue::glue("Biomart Unique Mapping Count: {slices[['Unique Mapping']]}"))
```

**Input Data:**

- `norm_data$genes` (Probe metadata, updated to include gene annotations specified by [Biomart](https://bioconductor.org/packages/3.14/bioc/html/biomaRt.html), output from [Step 7a](#7a-add-probe-annotations) above)

**Output Data:**

- A pie chart denoting the biomart mapping rates for each unique probe ID
- A printout denoting the count of unique mappings for both the original manufacturer mapping and the biomart mapping

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
    design_data <- list( matrix = design, mapping = names_mapping, groups = as.data.frame( cbind(sample = df[,c("Sample.Name")], group = group_names) ), contrasts = contrasts )
    return(design_data)
}


# Loading metadata from runsheet csv file
design_data <- runsheetToDesignMatrix(runsheet)
design <- design_data$matrix

# Write SampleTable.csv and contrasts.csv file
write.csv(design_data$groups, file.path(DIR_DGE, "SampleTable.csv"), row.names = FALSE)
write.csv(design_data$contrasts, file.path(DIR_DGE, "contrasts.csv"))
```

**Input Data:**

- `runsheet` (Path to runsheet, output from [Step 1](#1-create-sample-runsheet))

**Output Data:**

- `design` (R object containing the limma study design matrix, indicating the group that each sample belongs to)
- **SampleTable.csv** (table containing samples and their respective groups)
- **contrasts.csv** (table containing all pairwise comparisons)

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
limma::write.fit(res, adjust = 'BH', 
                file = "INTERIM.csv",
                row.names = FALSE,
                quote = TRUE,
                sep = ",")
```

**Input Data:**

- `norm_data` (R object containing background-corrected and normalized microarray data created in [Step 5](#5-between-array-normalization))
- `design` (R object containing the limma study design matrix, indicating the group that each sample belongs to, created in [Step 7c](#7c-generate-design-matrix) above)

**Output Data:**

- INTERIM.csv (Statistical values from individual probe level DE analysis, including:
  - Log2fc between all pairwise comparisons
  - T statistic for all pairwise comparison tests
  - P value for all pairwise comparison tests)

<br>

### 7e. Add Additional Columns and Format DE Table 

```R
## Reformat Table for consistency across DE analyses tables within GeneLab ##

# Read in DE table 
df_interim <- read.csv("INTERIM.csv")

# Reformat column names
reformat_names <- function(colname, group_name_mapping) {
  new_colname <- colname  %>% 
                  stringr::str_replace(pattern = "^P.value.adj.condition", replacement = "Adj.p.value_") %>%
                  stringr::str_replace(pattern = "^P.value.condition", replacement = "P.value_") %>%
                  stringr::str_replace(pattern = "^Coef.condition", replacement = "Log2fc_") %>% # This is the Log2FC as per: https://rdrr.io/bioc/limma/man/writefit.html
                  stringr::str_replace(pattern = "^t.condition", replacement = "T.stat_") %>%
                  stringr::str_replace(pattern = stringr::fixed("Genes.ProbeName"), replacement = "ProbeName") %>% 
                  stringr::str_replace(pattern = stringr::fixed("Genes.count_ENSEMBL_mappings"), replacement = "count_ENSEMBL_mappings") %>% 
                  stringr::str_replace(pattern = stringr::fixed("Genes.ProbeUID"), replacement = "ProbeUID") %>% 
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
                    
  print(glue::glue("Computing mean and standard deviation for Group {i} of {length(unique_groups)}"))
  print(glue::glue("Group: {current_group}"))
  print(glue::glue("Samples in Group: '{toString(current_samples)}'"))
  
  df_interim <- df_interim %>% 
    dplyr::mutate( 
      "Group.Mean_{current_group}" := rowMeans(dplyr::select(., all_of(current_samples))),
      "Group.Stdev_{current_group}" := matrixStats::rowSds(as.matrix(dplyr::select(., all_of(current_samples)))),
      ) %>% 
    dplyr::ungroup() %>%
    as.data.frame()
}

all_samples <- design_data$group %>% dplyr::pull(sample)
df_interim <- df_interim %>% 
  dplyr::mutate( 
    "All.mean" := rowMeans(dplyr::select(., all_of(all_samples))),
    "All.stdev" := matrixStats::rowSds(as.matrix(dplyr::select(., all_of(all_samples)))),
    ) %>% 
  dplyr::ungroup() %>%
  as.data.frame()

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
  "Genes.Description",
  "AveExpr" # Replaced by 'All.mean' column
  # "Genes.count_ENSEMBL_mappings", Keep this
)

df_interim <- df_interim %>% dplyr::select(-any_of(colnames_to_remove))

## Concatenate annotations for genes (for uniquely mapped probes) ##
### Read in annotation table for the appropriate organism ###
annot <- read.table(
            annotation_file_path,
            sep = "\t",
            header = TRUE,
            quote = "",
            comment.char = "",
        )

# Join annotation table and uniquely mapped data

# Determine appropriate keytype
map_primary_keytypes <- c(
  'Caenorhabditis elegans' = 'ENSEMBL',
  'Danio rerio' = 'ENSEMBL',
  'Drosophila melanogaster' = 'ENSEMBL',
  'Rattus norvegicus' = 'ENSEMBL',
  'Saccharomyces cerevisiae' = 'ENSEMBL',
  'Homo sapiens' = 'ENSEMBL',
  'Mus musculus' = 'ENSEMBL',
  'Arabidopsis thaliana' = 'TAIR'
)

df_interim <- merge(
                annot,
                df_interim,
                by = map_primary_keytypes[[unique(df_rs$organism)]],
                # ensure all original dge rows are kept.
                # If unmatched in the annotation database, then fill missing with NAN
                all.y = TRUE
            )

## Reorder columns before saving to file
ANNOTATIONS_COLUMN_ORDER = c(
  map_primary_keytypes[[unique(df_rs$organism)]],
  "SYMBOL",
  "GENENAME",
  "REFSEQ",
  "ENTREZID",
  "STRING_id",
  "GOSLIM_IDS"
)

PROBE_INFO_COLUMN_ORDER = c(
  "ProbeUID",
  "ProbeName",
  "count_ENSEMBL_mappings"
)
SAMPLE_COLUMN_ORDER <- all_samples
generate_prefixed_column_order <- function(subjects, prefixes) {
  #' Return a vector of columns based on subject and given prefixes
  #'  Used for both contrasts and groups column name generation
  
  # Track order of columns
  final_order = c()

  # For each contrast
  for (subject in subjects) {
    # Generate column names for each prefix and append to final_order
    for (prefix in prefixes) {
      final_order <- append(final_order, glue::glue("{prefix}{subject}"))
    }
  }
  return(final_order)
}
STAT_COLUMNS_ORDER <- generate_prefixed_column_order(
  subjects = colnames(design_data$contrasts),
  prefixes = c(
    "Log2fc_",
    "T.stat_",
    "P.value_",
    "Adj.p.value_"
    )
  )
ALL_SAMPLE_STATS_COLUMNS_ORDER <- c(
  "All.mean",
  "All.stdev",
  "F",
  "F.p.value"
)

GROUP_MEAN_COLUMNS_ORDER <- generate_prefixed_column_order(
  subjects = unique(design_data$groups$group),
  prefixes = c(
    "Group.Mean_"
    )
  )
GROUP_STDEV_COLUMNS_ORDER <- generate_prefixed_column_order(
  subjects = unique(design_data$groups$group),
  prefixes = c(
    "Group.Stdev_"
    )
  )
FINAL_COLUMN_ORDER <- c(
  ANNOTATIONS_COLUMN_ORDER, 
  PROBE_INFO_COLUMN_ORDER, 
  SAMPLE_COLUMN_ORDER, 
  STAT_COLUMNS_ORDER, 
  ALL_SAMPLE_STATS_COLUMNS_ORDER, 
  GROUP_MEAN_COLUMNS_ORDER,
  GROUP_STDEV_COLUMNS_ORDER
  )

## Assert final column order includes all columns from original table
if (!setequal(FINAL_COLUMN_ORDER, colnames(df_interim))) {
  FINAL_COLUMN_ORDER_STRING <- paste(FINAL_COLUMN_ORDER, collapse = ":::::")
  stop(glue::glue("Column reordering attempt resulted in different sets of columns than orignal. Order attempted: {FINAL_COLUMN_ORDER_STRING}"))
}

## Perform reordering
df_interim <- df_interim %>% dplyr::relocate(dplyr::all_of(FINAL_COLUMN_ORDER))

# Save to file
write.csv(df_interim, file.path(DIR_DGE, "differential_expression.csv"), row.names = FALSE)

### Generate and export PCA table for GeneLab visualization plots
## Only use positive expression values, negative values can make up a small portion ( < 0.5% ) of normalized expression values and cannot be log transformed
exp_raw <- log2(norm_data$E) # negatives get converted to NA
exp_raw <- na.omit(norm_data$E)
PCA_raw <- prcomp(t(exp_raw), scale = FALSE)
write.csv(PCA_raw$x,
          file.path(DIR_DGE, "visualization_PCA_table.csv")
          )

## Generate raw intensity matrix that includes annotations
raw_data_matrix <- background_corrected_data$genes %>% 
                    dplyr::select(ProbeUID, ProbeName) %>%
                    dplyr::bind_cols(background_corrected_data$E) %>% 
                    dplyr::left_join(unique_probe_ids, by = c("ProbeName" = expected_attribute_name ) ) %>%
                    dplyr::mutate( count_ENSEMBL_mappings = ifelse(is.na(ENSEMBL), 0, count_ENSEMBL_mappings) )

raw_data_matrix_annotated <- merge(
                annot,
                raw_data_matrix,
                by = map_primary_keytypes[[unique(df_rs$organism)]],
                # ensure all original dge rows are kept.
                # If unmatched in the annotation database, then fill missing with NAN
                all.y = TRUE
            )

## Reorder columns to match DGE output
FINAL_COLUMN_ORDER <- c(
  ANNOTATIONS_COLUMN_ORDER, 
  PROBE_INFO_COLUMN_ORDER, 
  SAMPLE_COLUMN_ORDER
  )

raw_data_matrix_annotated <- raw_data_matrix_annotated %>% 
  dplyr::relocate(dplyr::all_of(FINAL_COLUMN_ORDER))

write.csv(raw_data_matrix_annotated, file.path(DIR_RAW_DATA, "raw_intensities.csv"), row.names = FALSE)


## Generate normalized expression matrix that includes annotations
norm_data_matrix <- norm_data$genes %>% 
                    dplyr::select(ProbeUID, ProbeName) %>%
                    dplyr::bind_cols(norm_data$E) %>% 
                    dplyr::left_join(unique_probe_ids, by = c("ProbeName" = expected_attribute_name ) ) %>%
                    dplyr::mutate( count_ENSEMBL_mappings = ifelse(is.na(ENSEMBL), 0, count_ENSEMBL_mappings) )

norm_data_matrix_annotated <- merge(
                annot,
                norm_data_matrix,
                by = map_primary_keytypes[[unique(df_rs$organism)]],
                # ensure all original dge rows are kept.
                # If unmatched in the annotation database, then fill missing with NAN
                all.y = TRUE
            )

norm_data_matrix_annotated <- norm_data_matrix_annotated %>% 
  dplyr::relocate(dplyr::all_of(FINAL_COLUMN_ORDER))

write.csv(norm_data_matrix_annotated, file.path(DIR_NORMALIZED_EXPRESSION, "normalized_expression.csv"), row.names = FALSE)
```

**Input Data:**

- INTERIM.csv (Statistical values from individual probe level DE analysis, output from [Step 7d](#7d-perform-individual-probe-level-de) above)
- `annotation_file_path` (Annotation file url from 'genelab_annots_link' column of [GL-DPPD-7110_annotations.csv](https://github.com/nasa/GeneLab_Data_Processing/blob/GL_RefAnnotTable_1.0.0/GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110/GL-DPPD-7110_annotations.csv) corresponding to the subject organism)
- `primary_keytype` (Keytype to join annotation table and microarray probes, dependent on organism, e.g. mus musculus uses 'ENSEMBL')
- `background_corrected_data` (R object containing background-corrected microarray data)
- `norm_data` (R object containing background-corrected and normalized microarray data created in [Step 5](#5-between-array-normalization))

**Output Data:**

- **differential_expression.csv** (table containing normalized counts for each sample, group statistics, Limma probe DE results for each pairwise comparison, and gene annotations)
- visualization_PCA_table.csv (file used to generate GeneLab PCA plots)
- **raw_intensities.csv** (table containing the background corrected unnormalized intensity values for each sample including gene annotations)
- **normalized_expression.csv** (table containing the background corrected, normalized intensity values for each sample including gene annotations)