# GeneLab bioinformatics processing pipeline for Agilent 1-channel microarray data <!-- omit in toc -->

> **This page holds an overview and instructions for how GeneLab processes Agilent 1-channel microarray datasets\*. Exact processing commands and GL-DPPD-7112 version used for specific datasets are provided with their processed data in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/).**  
>

---

**Date:** March XX, 2025
**Revision:** -A
**Document Number:** GL-DPPD-7112-A

**Submitted by:**  
Crystal Han (GeneLab Data Processing Team)

**Approved by:**  
Samrawit Gebre (OSDR Project Manager)
Lauren Sanders (OSDR Project Scientist)
Amanda Saravia-Butler (GeneLab Science Lead)
Barbara Novak (GeneLab Data Processing Lead)

---

## Updates from previous version

Updated [Ensembl Reference Files](../../../GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv) to the following releases:
- Animals: Ensembl release 112
- Plants: Ensembl plants release 59
- Bacteria: Ensembl bacteria release 59

Software Updates:

| Program | Previous Version | New Version    |
|:--------|:-----------------|:---------------|
|R|4.1.3|4.4.2|
|DT|0.26|0.33|
|dplyr|1.0.10|1.1.4|
|stringr|1.5.0|1.5.1|
|R.utils|2.12.2|2.12.3|
|limma|3.50.3|3.62.2|
|glue|1.6.2|1.8.0|
|ggplot2|3.4.0|3.5.1|
|biomaRt|2.50.0|2.62.0|
|matrixStats|0.63.0|1.5.0|
|dp_tools|1.3.4|1.3.5|
|Quarto|1.2.313|1.6.40|

Custom Annotations

- Added ability to use custom gene annotations when annotations are not available in Biomart or Ensembl FTP for *Arabidopsis thaliana*, see [Step 7](#7-probe-annotations)

---

# Table of contents <!-- omit in toc -->

- [Software used](#software-used)
- [General processing overview with example commands](#general-processing-overview-with-example-commands)
    - [1. Create Sample RunSheet](#1-create-sample-runsheet)
    - [2. Load Data](#2-load-data)
      - [2a. Load Libraries and Define Input Parameters](#2a-load-libraries-and-define-input-parameters)
      - [2b. Define Custom Functions](#2b-define-custom-functions)
      - [2c. Load Metadata and Raw Data](#2c-load-metadata-and-raw-data)
      - [2d. Load Annotation Metadata](#2d-load-annotation-metadata)
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
    - [7. Probe Annotations](#7-probe-annotations)
      - [7a. Get Probe Annotations](#7a-get-probe-annotations)
      - [7b. Summarize Gene Mapping](#7b-summarize-gene-mapping)
      - [7c. Generate Annotated Raw and Normalized Expression Tables](#7c-generate-annotated-raw-and-normalized-expression-tables)
    - [8. Perform Probe Differential Expression (DE)](#8-perform-probe-differential-expression-de)
      - [8a. Generate Design Matrix](#8a-generate-design-matrix)
      - [8b. Perform Individual Probe Level DE](#8b-perform-individual-probe-level-de)
      - [8c. Add Annotation and Stats Columns and Format DE Table](#8c-add-annotation-and-stats-columns-and-format-de-table)

---

# Software used  

|Program|Version|Relevant Links|
|:------|:------:|:-------------|
|R|4.4.2|[https://www.r-project.org/](https://www.r-project.org/)|
|DT|0.33|[https://github.com/rstudio/DT](https://github.com/rstudio/DT)|
|dplyr|1.1.4|[https://dplyr.tidyverse.org](https://dplyr.tidyverse.org)|
|stringr|1.5.1|[https://stringr.tidyverse.org](https://stringr.tidyverse.org)|
|R.utils|2.12.3|[https://github.com/HenrikBengtsson/R.utils](https://github.com/HenrikBengtsson/R.utils)|
|limma|3.62.2|[https://bioconductor.org/packages/3.14/bioc/html/limma.html](https://bioconductor.org/packages/3.14/bioc/html/limma.html)|
|glue|1.8.0|[https://glue.tidyverse.org](https://glue.tidyverse.org)|
|ggplot2|3.5.1|[https://ggplot2.tidyverse.org](https://ggplot2.tidyverse.org)|
|biomaRt|2.62.0|[https://bioconductor.org/packages/3.14/bioc/html/biomaRt.html](https://bioconductor.org/packages/3.14/bioc/html/biomaRt.html)|
|matrixStats|1.5.0|[https://github.com/HenrikBengtsson/matrixStats](https://github.com/HenrikBengtsson/matrixStats)|
|statmod|1.5.0|[https://github.com/cran/statmod](https://github.com/cran/statmod)|
|dp_tools|1.3.5|[https://github.com/J-81/dp_tools](https://github.com/J-81/dp_tools)|
|singularity|3.9|[https://sylabs.io](https://sylabs.io)|
|Quarto|1.6.40|[https://quarto.org](https://quarto.org)|

---

# General processing overview with example commands  

> Exact processing commands and output files listed in **bold** below are included with each Microarray processed dataset in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/).

---

## 1. Create Sample RunSheet

> Notes: 
> - Rather than running the commands below to create the runsheet needed for processing, the runsheet may also be created manually by following the [file specification](../Workflow_Documentation/NF_MAAgilent1ch/examples/runsheet/README.md).
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

## 2. Load Data

> Note: Steps 2 - 8 are done in R

<br>

### 2a. Load Libraries and Define Input Parameters

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
```

<br>

### 2b. Define Custom Functions

#### all_true()
<details>
  <summary>wraps R <code>base::all()</code> function; overrides default behavior for empty input vector</summary>

  ```R
  all_true <- function(i_vector) {
    if ( length(i_vector) == 0 ) {
      stop(paste("Input vector is length zero"))
    }
    all(i_vector)
  }
  ```

  **Function Parameter Definitions:**
  - `i_vector=` - a vector of logical values

  **Returns:** a logical of length 1; `TRUE` if all values are true, `FALSE` otherwise; stops and returns an error if input vector is empty
</details>

#### runsheet_paths_are_URIs()
<details>
  <summary>tests if paths provided in runsheet dataframe are URIs</summary>

  ```R
  runsheet_paths_are_URIs <- function(df_runsheet) {
    all_true(stringr::str_starts(df_runsheet$`Array Data File Path`, "https"))
  }
  ```

  **Custom Functions Used:**
  - [all_true()](#all_true)

  **Function Parameter Definitions:**
  - `df_runsheet=` - a dataframe containing the sample runsheet information

  **Returns:** a logical of length 1; `TRUE` if all values in the `Array Data File Path` of the runsheet start with "https", `FALSE` otherwise; stops and returns an error if input vector is empty
</details>

#### download_files_from_runsheet()
<details>
  <summary>downloads the raw data files</summary>

  ```R
  download_files_from_runsheet <- function(df_runsheet) {
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
  ```

  **Function Parameter Definitions:**
  - `df_runsheet=` - a dataframe containing the sample runsheet information

  **Returns:** a list of filenames that were downloaded; same as the `Array Data File Name` in the sample runsheet
</details>

#### fetch_organism_specific_annotation_table()
<details>
  <summary>determines the organism specific annotation file to use based on the provided organism name</summary>

  ```R
  fetch_organism_specific_annotation_table <- function(organism, annotation_table_link) {
    # Uses the latest GeneLab annotations table to find the organism specific annotation file path and ensembl version
    # Raises an exception if the organism does not have an associated annotation file or ensembl version yet
    
    all_organism_table <- read.csv(annotation_table_link)

    annotation_table <- all_organism_table %>% dplyr::filter(species == organism)

    # Guard clause: Ensure annotation_table populated
    # Else: raise exception for unsupported organism
    if (nrow(annotation_table) == 0 || annotation_table$genelab_annots_link == "" || is.na(annotation_table$ensemblVersion)) {
      stop(glue::glue("Organism supplied '{organism}' is not supported. See the following url for supported organisms: {annotation_table_link}.  Supported organisms will correspond to a row based on the 'species' column and include a url in the 'genelab_annots_link' column of that row and a version number in the 'ensemblVersion' column."))
    }

    return(annotation_table)
  }
  ```

  **Function Parameter Definitions:**
  - `organism=` - a string containing the name of the organism (as found in the species column of the GeneLab annotation table)
  - `annotation_table_link=` - a string specifying the URL or path to latest GeneLab Annotations file, see [GL-DPPD-7110-A_annotations.csv](../../../GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv)

  **Returns:** a dataframe containing all rows in the GeneLab annotations file that match the specified organism
</details>

#### agilent_image_plot()
<details>
  <summary>plots pseudo images of each array</summary>

  ```R
  agilent_image_plot <- function(eListRaw, transform_func = identity) {
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
  ```

  **Function Parameter Definitions:**
  - `eListRaw=` - R object containing expression matrix to be plotted
  - `transform_func=identity` - function used to transform expression matrix before plotting

  **Returns:** pseudo images of each array
</details>

#### boxplot_expression_safe_margin()
<details>
  <summary>plots boxplot of expression data for each array</summary>

  ```R
  boxplot_expression_safe_margin <- function(data, transform_func = identity, xlab = "Log2 Intensity") {
    # Basic box plot
    df_data <- as.data.frame(transform_func(data$E))
    ggplot2::ggplot(stack(df_data), ggplot2::aes(x=values, y=ind)) + 
      ggplot2::geom_boxplot() + 
      ggplot2::scale_y_discrete(limits=rev) +
      ggplot2::labs(y= "Sample Name", x = xlab)
  }
  ```

  **Function Parameter Definitions:**
  - `data=` - R object containing expression matrix to be plotted
  - `transform_func=identity` - function used to transform expression matrix before plotting
  - `xlab="Log2 Intensity"` - string containing x-axis label for the plot

  **Returns:** boxplot of expression data for each array
</details>

#### shortened_organism_name()
<details>
  <summary>shortens organism names, for example 'Homo Sapiens' to 'hsapiens'</summary>

  ```R
  shortened_organism_name <- function(long_name) {
    #' Convert organism names like 'Homo Sapiens' into 'hsapiens'
    tokens <- long_name %>% stringr::str_split(" ", simplify = TRUE)
    genus_name <- tokens[1]

    species_name <- tokens[2]

    short_name <- stringr::str_to_lower(paste0(substr(genus_name, start = 1, stop = 1), species_name))

    return(short_name)
  }
  ```

  **Function Parameter Definitions:**
  - `long_name=` - a string containing the long name of the organism

  **Returns:** a string containing the short name of the organism
</details>

#### get_biomart_attribute()
<details>
  <summary>retrieves resolved biomart attribute source from runsheet dataframe</summary>

  ```R
  get_biomart_attribute <- function(df_rs) {
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
  ```

  **Function Parameter Definitions:**
  - `df_rs=` - a dataframe containing the sample runsheet information

  **Returns:** a string containing the formatted value from the `biomart_attribute` column of the runsheet, with all spaces converted to underscores and uppercase converted to lowercase; if no `biomart_attribute` exists in the runsheet, stop and return an error
</details>

#### get_ensembl_genomes_mappings_from_ftp()
<details>
  <summary>obtains mapping table directly from ftp; useful when biomart live service no longer exists for desired version</summary>

  ```R
  get_ensembl_genomes_mappings_from_ftp <- function(organism, ensembl_genomes_portal, ensembl_genomes_version, biomart_attribute) {
    request_url <- glue::glue("https://ftp.ebi.ac.uk/ensemblgenomes/pub/{ensembl_genomes_portal}/release-{ensembl_genomes_version}/mysql/{ensembl_genomes_portal}_mart_{ensembl_genomes_version}/{organism}_eg_gene__efg_{biomart_attribute}__dm.txt.gz")

    print(glue::glue("Mappings file URL: {request_url}"))

    # Create a temporary file name
    temp_file <- tempfile(fileext = ".gz")

    # Download the gzipped table file using the download.file function
    download.file(url = request_url, destfile = temp_file, method = "libcurl") # Use 'libcurl' to support ftps

    # Uncompress the file
    uncompressed_temp_file <- tempfile()
    gzcon <- gzfile(temp_file, "rt")
    content <- readLines(gzcon)
    writeLines(content, uncompressed_temp_file)
    close(gzcon)


    # Load the data into a dataframe
    mapping <- read.table(uncompressed_temp_file, # Read the uncompressed file
                          # Add column names as follows: MAPID, TAIR, PROBEID
                          col.names = c("MAPID", "ensembl_gene_id", biomart_attribute),
                          header = FALSE, # No header in original table
                          sep = "\t") # Tab separated

    # Clean up temporary files
    unlink(temp_file)
    unlink(uncompressed_temp_file)

    return(mapping)
  }
  ```

  **Function Parameter Definitions:**
  - `organism=` - a string containing the name of the organism (formatted using `shortened_organism_name()`)
  - `ensembl_genomes_portal=` - a string containing the name of the genomes portal, for example 'plants'
  - `ensembl_genomes_version=` - a string containing the version of Ensembl to use
  - `biomart_attribute=` - a string containing the biomart attribute (formatted using `get_biomart_attribute()`)

  **Returns:** a dataframe containing the mapping between Ensembl ID and probe ID, as obtained via FTP
</details>

#### list_to_unique_piped_string()
<details>
  <summary>converts character vector into string denoting unique elements separated by '|' characters</summary>

  ```R
  list_to_unique_piped_string <- function(str_list) {
    #! Convert vector of multi-mapped genes to string separated by '|' characters
    #! e.g. c("GO1","GO2","GO2","G03") -> "GO1|GO2|GO3"
    return(toString(unique(str_list)) %>% stringr::str_replace_all(pattern = stringr::fixed(", "), replacement = "|"))
  }
  ```

  **Function Parameter Definitions:**
  - `str_list=` - vector of character elements

  **Returns:** a string containing the unique elements from `str_list` concatenated together, separated by '|' characters
</details>

#### runsheet_to_design_matrix()
<details>
  <summary>loads the GeneLab runsheet into a list of dataframes</summary>

  ```R
  runsheet_to_design_matrix <- function(runsheet_path) {
      # Pull all factors for each sample in the study from the runsheet created in Step 1
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
  ```

  **Function Parameter Definitions:**
  - `runsheet_path=` - a string containing the path to the runsheet generated in [Step 1](#1-create-sample-runsheet)

  **Returns:** a list of R objects containing the sample information and metadata
  - `design_data$matrix` - a design (or model) matrix describing the conditions in the dataset
  - `design_data$mapping` - a dataframe mapping the human-readable group names to the names of the conditions modified for use in R
  - `design_data$groups` - a dataframe of group names and contrasts for each sample
  - `design_data$contrasts` - a matrix of all pairwise comparisons of the groups
</details>

#### lm_fit_pairwise()
<details>
  <summary>performs all pairwise comparisons using <code>limma::lmFit()</code></summary>

  ```R
  lm_fit_pairwise <- function(norm_data, design) {
      # Approach based on limma manual section 17.4 (version 3.52.4)
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
  ```

  **Function Parameter Definitions:**
  - `norm_data=` - an R object containing log-ratios or log-expression values for a series of arrays, with rows corresponding to genes and columns to samples
  - `design=` - the design matrix of the microarray experiment, with rows corresponding to samples and columns to coefficients to be estimated

  **Returns:** an R object of class `MArrayLM`
</details>

#### reformat_names()
<details>
  <summary>reformats column names for consistency across DE analyses tables within GeneLab</summary>

  ```R
  reformat_names <- function(colname, group_name_mapping) {
    new_colname <- colname  %>% 
                    stringr::str_replace(pattern = "^P.value.adj.condition", replacement = "Adj.p.value_") %>%
                    stringr::str_replace(pattern = "^P.value.condition", replacement = "P.value_") %>%
                    stringr::str_replace(pattern = "^Coef.condition", replacement = "Log2fc_") %>% # This is the Log2FC as per: https://rdrr.io/bioc/limma/man/writefit.html
                    stringr::str_replace(pattern = "^t.condition", replacement = "T.stat_") %>%
                    stringr::str_replace(pattern = "^Genes\\.", replacement = "") %>%
                    stringr::str_replace(pattern = ".condition", replacement = "v")
    
    # remap to group names before make.names was applied
    unique_group_name_mapping <- unique(group_name_mapping) %>% arrange(-nchar(safe_name))
    for ( i in seq(nrow(unique_group_name_mapping)) ) {
      safe_name <- unique_group_name_mapping[i,]$safe_name
      original_name <- unique_group_name_mapping[i,]$original_name
      new_colname <- new_colname %>% stringr::str_replace(pattern = stringr::fixed(safe_name), replacement = original_name)
    }

    return(new_colname)
  }
  ```

  **Function Parameter Definitions:**
  - `colnames=` - a character vector containing the column names to reformat
  - `group_name_mapping=` - a dataframe mapping the original human-readable group names to the R modified safe names

  **Returns:** a character vector containing the formatted column names
</details>

#### generate_prefixed_column_order()
<details>
  <summary>creates a vector of column names based on subject and given prefixes; used for both contrasts and groups column name generation</summary>

  ```R
  generate_prefixed_column_order <- function(subjects, prefixes) {
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
  ```

  **Function Parameter Definitions:**
  - `subjects` - a character vector containing subject strings to add prefixes to 
  - `prefixes` - a character vector of prefixes to add to the beginning of each subject string

  **Returns:** a character vector with all possible combinations of prefix + subject
</details>

<br>

### 2c. Load Metadata and Raw Data

```R
# fileEncoding removes strange characters from the column names
df_rs <- read.csv(runsheet, check.names = FALSE, fileEncoding = 'UTF-8-BOM') 

if ( runsheet_paths_are_URIs(df_rs) ) {
  print("Determined Raw Data Locations are URIS")
  local_paths <- download_files_from_runsheet(df_rs)
} else {
  print("Or Determined Raw Data Locations are local paths")
  local_paths <- df_rs$`Array Data File Path`
}

# uncompress files if needed
if ( all_true(stringr::str_ends(local_paths, ".gz")) ) {
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

# Handle raw data which lacks certain replaceable column data

## This likely arises as Agilent Feature Extraction (the process that generates the raw data files on OSDR) 
##   gives some user flexibilty in what probe column to ouput

## Missing ProbeUID "Unique integer for each unique probe in a design"
### Source: https://www.agilent.com/cs/library/usermanuals/public/GEN-MAN-G4460-90057.pdf Page 178
### Remedy: Assign unique integers for each probe

if ( !("ProbeUID" %in% colnames(raw_data$genes)) ) {
  # Assign unique integers for each probe
  print("Assigning `ProbeUID` as original files did not include them")
  raw_data$genes$ProbeUID <- seq_len(nrow(raw_data$genes))
}

# Summarize raw data
print(paste0("Number of Arrays: ", dim(raw_data)[2]))
print(paste0("Number of Probes: ", dim(raw_data)[1]))
```

**Custom Functions Used:**

- [all_true()](#all_true)
- [runsheet_paths_are_URIs()](#runsheet_paths_are_URIs)
- [download_files_from_runsheet()](#download_files_from_runsheet)

**Input Data:**

- `runsheet` (Path to runsheet, output from [Step 1](#1-create-sample-runsheet))

**Output Data:**

- `df_rs` (R dataframe containing information from the runsheet)
- `raw_data` (R object containing raw microarray data)

    > Note: The raw data R object will be used to generate quality assessment (QA) plots in the next step.

<br>

### 2d. Load Annotation Metadata

```R
# If using custom annotation, local_annotation_dir is path to directory containing annotation file and annotation_config_path is path/url to config file
local_annotation_dir <- NULL # <path/to/custom_annotation>
annotation_config_path <- NULL # <path/to/config_file>

annotation_table_link <- "https://raw.githubusercontent.com/nasa/GeneLab_Data_Processing/GL_RefAnnotTable-A_1.1.0/GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv"

annotation_table <- fetch_organism_specific_annotation_table(unique(df_rs$organism))

annotation_file_path <- annotation_table$genelab_annots_link
ensembl_version <- as.character(annotation_table$ensemblVersion)
```

**Custom Functions Used:**

- [fetch_organism_specific_annotation_table()](#fetch_organism_specific_annotation_table)

**Input Data:**

- `local_annotation_dir` (Path to local annotation directory if using custom annotations, see [Step 7a](#7a-get-probe-annotations))

    > Note: If not using custom annotations, leave `local_annotation_dir` as `NULL`.

- `annotation_config_path` (URL or path to annotation config file if using custom annotations, see [Step 7a](#7a-get-probe-annotations))

    > Note: If not using custom annotations, leave `annotation_config_path` as `NULL`.

- `df_rs$organism` (organism specified in the runsheet created in [Step 1](#1-create-sample-runsheet))
- `annotation_table_link` (URL or path to latest GeneLab Annotations file, see [GL-DPPD-7110-A_annotations.csv](../../../GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv))

**Output Data:**

- `annotation_file_path` (reference organism annotation file url indicated in the 'genelab_annots_link' column of the GeneLab Annotations file provided in `annotation_table_link`)
- `ensembl_version` (reference organism Ensembl version indicated in the 'ensemblVersion' column of the GeneLab Annotations file provided in `annotation_table_link`)

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

- `raw_data` (raw data R object created in [Step 2c](#2c-load-metadata-and-raw-data) above)

**Output Data:**

- Plot containing the density of raw intensities for each array (raw intensity values with background intensity values subtracted; lack of overlap indicates a need for normalization)

<br>

### 3b. Pseudo Image Plots

```R
agilent_image_plot(raw_data, transform_func = function(expression_matrix) log2(expression_matrix + 1))
```

**Custom Functions Used:**

- [agilent_image_plot()](#agilent_image_plot)

**Input Data:**

- `raw_data` (raw data R object created in [Step 2c](#2c-load-metadata-and-raw-data) above)

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

- `raw_data` (raw data R object created in [Step 2c](#2c-load-metadata-and-raw-data) above)

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

- `raw_data` (raw data R object created in [Step 2c](#2c-load-metadata-and-raw-data) above)

**Output Data:**

- Foreground vs. background expression plot for each array before background correction and normalization 

<br>

### 3e. Boxplots

```R
boxplot_expression_safe_margin(raw_data, transform_func = log2)
```

**Custom Functions Used:**

- [boxplot_expression_safe_margin()](#boxplot_expression_safe_margin)

**Input Data:**

- `raw_data` (raw data R object created in [Step 2c](#2c-load-metadata-and-raw-data) above)

**Output Data:**

- Boxplot of raw expression data for each array before background correction and normalization 

<br>

---

## 4. Background Correction

```R
background_corrected_data <- limma::backgroundCorrect(raw_data, method = "normexp")
```

**Input Data:**

- `raw_data` (raw data R object created in [Step 2c](#2c-load-metadata-and-raw-data) above)

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
agilent_image_plot(norm_data, 
                 transform_func = function(expression_matrix) log2(2**expression_matrix + 1) # Compute as log2 of normalized expression after adding a +1 offset to prevent negative values in the pseudoimage
                 )
```

**Custom Functions Used:**

- [agilent_image_plot()](#agilent_image_plot)

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
boxplot_expression_safe_margin(norm_data)
```

**Custom Functions Used:**

- [boxplot_expression_safe_margin()](#boxplot_expression_safe_margin)

**Input Data:**

- `norm_data` (R object containing background-corrected and normalized microarray data created in [Step 5](#5-between-array-normalization) above)

**Output Data:**

- Boxplot of expression data for each array after background correction and normalization 

<br>

---

## 7. Probe Annotations

<br>

### 7a. Get Probe Annotations

```R
organism <- shortened_organism_name(unique(df_rs$organism))
annot_key <- ifelse(organism %in% c("athaliana"), 'TAIR', 'ENSEMBL')

if (organism %in% c("athaliana")) {
  ENSEMBL_VERSION = ensembl_version
  ensembl_genomes_portal = "plants"
  print(glue::glue("Using ensembl genomes ftp to get specific version of probe id mapping table. Ensembl genomes portal: {ensembl_genomes_portal}, version: {ENSEMBL_VERSION}"))
  expected_attribute_name <- get_biomart_attribute(df_rs)
  df_mapping <- get_ensembl_genomes_mappings_from_ftp(
    organism = organism,
    ensembl_genomes_portal = ensembl_genomes_portal,
    ensembl_genomes_version = ENSEMBL_VERSION,
    biomart_attribute = expected_attribute_name
    )

  # TAIR from the mapping tables tend to be in the format 'AT1G01010.1' but the raw data has 'AT1G01010'
  # So here we remove the '.NNN' from the mapping table where .NNN is any number
  df_mapping$ensembl_gene_id <- stringr::str_replace_all(df_mapping$ensembl_gene_id, "\\.\\d+$", "")

  use_custom_annot <- FALSE
} else {
  # Use biomart from main Ensembl website which archives keep each release on the live service
  # locate dataset
  expected_dataset_name <- shortened_organism_name(unique(df_rs$organism)) %>% stringr::str_c("_gene_ensembl")
  print(paste0("Expected dataset name: '", expected_dataset_name, "'"))

  expected_attribute_name <- get_biomart_attribute(df_rs)
  print(paste0("Expected attribute name: '", expected_attribute_name, "'"))

  # Specify Ensembl version used in current GeneLab reference annotations
  ENSEMBL_VERSION <- ensembl_version

  print(glue::glue("Using Ensembl biomart to get specific version of mapping table. Ensembl version: {ENSEMBL_VERSION}"))

  # Check if organism/array design is supported in biomart
  use_custom_annot <- TRUE

  ensembl <- biomaRt::useEnsembl(biomart = "genes", version = ENSEMBL_VERSION)
  ensembl_datasets <- biomaRt::listDatasets(ensembl)
  if (expected_dataset_name %in% ensembl_datasets$dataset) {
    ensembl <- biomaRt::useEnsembl(biomart = "genes", dataset = expected_dataset_name, version = ENSEMBL_VERSION)
    ensembl_attributes <- biomaRt::listAttributes(ensembl)
    if (expected_attribute_name %in% ensembl_attributes$name) {
      use_custom_annot <- FALSE
    }
  }

  if (use_custom_annot) {
    unloadNamespace("biomaRt")
  } else {
    print(ensembl)

    probe_ids <- unique(norm_data$genes$ProbeName)

    # Create probe map
    # Run Biomart Queries in chunks to prevent request timeouts
    #   Note: If timeout is occuring (possibly due to larger load on biomart), reduce chunk size
    CHUNK_SIZE= 1500
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

      if (nrow(chunk_results) > 0) {
        df_mapping <- df_mapping %>% dplyr::bind_rows(chunk_results)
      }
      
      Sys.sleep(10) # Slight break between requests to prevent back-to-back requests
    }
  }
}

# At this point, we have df_mapping from either the biomart live service or the ensembl genomes ftp archive depending on the organism
# If no df_mapping obtained (e.g., not supported in biomart), use custom annotations; otherwise, merge in-house annotations to df_mapping

if (use_custom_annot) {
  expected_attribute_name <- 'ProbeName'

  annot_type <- 'NO_CUSTOM_ANNOT'
  if (!is.null(local_annotation_dir) && !is.null(annotation_config_path)) {
    config_df <- read.csv(annotation_config_path, row.names=1)
    if (unique(df_rs$`biomart_attribute`) %in% row.names(config_df)) {
      annot_config <- config_df[unique(df_rs$`biomart_attribute`), ]
      annot_type <- annot_config$annot_type[[1]]
    } else {
      warning(paste0("No entry for '", unique(df_rs$`biomart_attribute`), "' in provided config file: ", annotation_config_path))
    }
  } else {
    warning("Need to provide both local_annotation_dir and annotation_config_path to use custom annotation.")
  }

  if (annot_type == 'agilent') {
    unique_probe_ids <- read.delim(
      file.path(local_annotation_dir, annot_config$annot_filename[[1]]),
      header = TRUE, na.strings = c('NA', '')
    )[c('ProbeID', 'EnsemblID', 'GeneSymbol', 'GeneName', 'RefSeqAccession', 'EntrezGeneID', 'GO')]

    stopifnot(nrow(unique_probe_ids) == length(unique(unique_probe_ids$ProbeID)))

    # Clean columns
    unique_probe_ids$GO <- purrr::map_chr(stringr::str_extract_all(unique_probe_ids$GO, 'GO:\\d{7}'), ~paste0(unique(.), collapse = '|')) %>% stringr::str_replace('^$', NA_character_)

    names(unique_probe_ids) <- c('ProbeName', 'ENSEMBL', 'SYMBOL', 'GENENAME', 'REFSEQ', 'ENTREZID', 'GOSLIM_IDS')

    unique_probe_ids$STRING_id <- NA_character_

    gene_col <- 'ENSEMBL'
    if (sum(!is.na(unique_probe_ids$ENTREZID)) > sum(!is.na(unique_probe_ids$ENSEMBL))) {
      gene_col <- 'ENTREZID'
    }
    if (sum(!is.na(unique_probe_ids$SYMBOL)) > max(sum(!is.na(unique_probe_ids$ENTREZID)), sum(!is.na(unique_probe_ids$ENSEMBL)))) {
      gene_col <- 'SYMBOL'
    }

    unique_probe_ids <- unique_probe_ids %>%
                        dplyr::mutate( 
                          count_gene_mappings = 1 + stringr::str_count(get(gene_col), stringr::fixed("|")),
                          gene_mapping_source = gene_col
                        )
  } else if (annot_type == 'custom') {
    unique_probe_ids <- read.csv(
      file.path(local_annotation_dir, annot_config$annot_filename[[1]]),
      header = TRUE, na.strings = c('NA', '')
    )
  } else {
    annot_cols <- c('ProbeName', 'ENTREZID', 'SYMBOL', 'GENENAME', 'ENSEMBL', 'REFSEQ', 'GOSLIM_IDS', 'STRING_id', 'count_gene_mappings', 'gene_mapping_source')
    unique_probe_ids <- setNames(data.frame(matrix(NA_character_, nrow = 1, ncol = length(annot_cols))), annot_cols)
  }
} else {
  annot <- read.table(
      as.character(annotation_file_path),
      sep = "\t",
      header = TRUE,
      quote = "",
      comment.char = ""
  )

  unique_probe_ids <- df_mapping %>% 
                        dplyr::mutate(dplyr::across(!!sym(expected_attribute_name), as.character)) %>% # Ensure probe ids treated as character type
                        dplyr::group_by(!!sym(expected_attribute_name)) %>% 
                        dplyr::summarise(
                          ENSEMBL = list_to_unique_piped_string(ensembl_gene_id)
                          ) %>%
                        # Count number of ensembl IDS mapped
                        dplyr::mutate( 
                          count_gene_mappings = 1 + stringr::str_count(ENSEMBL, stringr::fixed("|")),
                          gene_mapping_source = annot_key
                        ) %>%
                        dplyr::left_join(annot, by = c("ENSEMBL" = annot_key))
}

norm_data$genes <- norm_data$genes %>% 
  dplyr::left_join(unique_probe_ids, by = c("ProbeName" = expected_attribute_name ) ) %>%
  dplyr::mutate( count_gene_mappings := ifelse(is.na(count_gene_mappings), 0, count_gene_mappings) ) %>%
  dplyr::mutate( gene_mapping_source := unique(unique_probe_ids$gene_mapping_source) )
```

**Custom Functions Used:**

- [shortened_organism_name()](#shortened_organism_name)
- [get_biomart_attribute()](#get_biomart_attribute)
- [get_ensembl_genomes_mappings_from_ftp()](#get_ensembl_genomes_mappings_from_ftp)
- [list_to_unique_piped_string()](#list_to_unique_piped_string)

**Input Data:**

- `df_rs$organism` (organism specified in the runsheet created in [Step 1](#1-create-sample-runsheet))
- `df_rs$biomart_attribute` (array design biomart identifier specified in the runsheet created in [Step 1](#1-create-sample-runsheet))
- `annotation_file_path` (reference organism annotation file url indicated in the 'genelab_annots_link' column of the GeneLab Annotations file provided in `annotation_table_link`, output from [Step 2d](#2d-load-annotation-metadata))
- `ensembl_version` (reference organism Ensembl version indicated in the 'ensemblVersion' column of the GeneLab Annotations file provided in `annotation_table_link`, output from [Step 2d](#2d-load-annotation-metadata))
- `annot_key` (keytype to join annotation table and microarray probes, dependent on organism, e.g. mus musculus uses 'ENSEMBL')
- `local_annotation_dir` (path to local annotation directory if using custom annotations, output from [Step 2d](#2d-load-annotation-metadata))
- `annotation_config_path` (URL or path to annotation config file if using custom annotations, output from [Step 2d](#2d-load-annotation-metadata))

  > Note: See [Agilent_array_annotations.csv](../Array_Annotations/Agilent_array_annotations.csv) for the latest config file used at GeneLab. This file can also be created manually by following the [file specification](../Workflow_Documentation/NF_MAAgilent1ch/examples/annotations/README.md).

- `norm_data$genes` (Manufacturer's probe metadata, including probe IDs and sequence position gene annotations associated with the `norm_data` R object containing background-corrected and normalized microarray data created in [Step 5](#5-between-array-normalization))

**Output Data:**

- `unique_probe_ids` (R object containing probe ID to gene annotation mappings)
- `norm_data$genes` (Probe metadata, updated to include gene annotations specified by [Biomart](https://bioconductor.org/packages/3.14/bioc/html/biomaRt.html) or custom annotations)

<br>

### 7b. Summarize Gene Mapping

```R
# Pie Chart with Percentages
slices <- c(
    'Control probes' = nrow(norm_data$gene %>% dplyr::filter(ControlType != 0) %>% dplyr::distinct(ProbeName)), 
    'Unique Mapping' = nrow(norm_data$gene %>% dplyr::filter(ControlType == 0) %>% dplyr::filter(count_gene_mappings == 1) %>% dplyr::distinct(ProbeName)), 
    'Multi Mapping' = nrow(norm_data$gene %>% dplyr::filter(ControlType == 0) %>% dplyr::filter(count_gene_mappings > 1) %>% dplyr::distinct(ProbeName)), 
    'No Mapping' = nrow(norm_data$gene %>% dplyr::filter(ControlType == 0) %>% dplyr::filter(count_gene_mappings == 0) %>% dplyr::distinct(ProbeName))
)
pct <- round(slices/sum(slices)*100)
chart_names <- names(slices)
chart_names <- glue::glue("{names(slices)} ({slices})") # add count to labels
chart_names <- paste(chart_names, pct) # add percents to labels
chart_names <- paste(chart_names,"%",sep="") # ad % to labels
pie(slices,labels = chart_names, col=rainbow(length(slices)),
    main=glue::glue("Mapping to Primary Keytype\n {nrow(norm_data$gene %>% dplyr::distinct(ProbeName))} Total Unique Probes")
    )

original_mapping_rate = nrow(norm_data$gene %>% dplyr::filter(ControlType == 0) %>% dplyr::filter(ProbeName != SystematicName) %>% dplyr::distinct(ProbeName))
print(glue::glue("Original Manufacturer Reported Mapping Count: {original_mapping_rate}"))
print(glue::glue("Unique Mapping Count: {slices[['Unique Mapping']]}"))
```

**Input Data:**

- `norm_data$genes` (Probe metadata, updated to include gene annotations specified by [Biomart](https://bioconductor.org/packages/3.14/bioc/html/biomaRt.html) or custom annotations, output from [Step 7a](#7a-get-probe-annotations) above)

**Output Data:**

- A pie chart denoting the gene mapping rates for each unique probe ID
- A printout denoting the count of unique mappings for both the original manufacturer mapping and the gene mapping

<br>

### 7c. Generate Annotated Raw and Normalized Expression Tables

```R
## Reorder columns before saving to file
ANNOTATIONS_COLUMN_ORDER = c(
  annot_key,
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
  "count_gene_mappings",
  "gene_mapping_source"
)
SAMPLE_COLUMN_ORDER <- df_rs$`Sample Name`

## Generate raw intensity matrix that includes annotations
raw_data_matrix <- background_corrected_data$genes %>% 
                    dplyr::select(ProbeUID, ProbeName) %>%
                    dplyr::bind_cols(background_corrected_data$E) %>% 
                    dplyr::left_join(unique_probe_ids, by = c("ProbeName" = expected_attribute_name ) ) %>%
                    dplyr::mutate( count_gene_mappings = ifelse(is.na(count_gene_mappings), 0, count_gene_mappings) ) %>%
                    dplyr::mutate( gene_mapping_source = unique(unique_probe_ids$gene_mapping_source) )

## Perform reordering
FINAL_COLUMN_ORDER <- c(
  ANNOTATIONS_COLUMN_ORDER, 
  PROBE_INFO_COLUMN_ORDER, 
  SAMPLE_COLUMN_ORDER
  )

raw_data_matrix <- raw_data_matrix %>% 
  dplyr::relocate(dplyr::all_of(FINAL_COLUMN_ORDER))

write.csv(raw_data_matrix, file.path(DIR_RAW_DATA, "raw_intensities_GLmicroarray.csv"), row.names = FALSE)

## Generate normalized expression matrix that includes annotations
norm_data_matrix <- norm_data$genes %>% 
                    dplyr::select(ProbeUID, ProbeName) %>%
                    dplyr::bind_cols(norm_data$E) %>% 
                    dplyr::left_join(unique_probe_ids, by = c("ProbeName" = expected_attribute_name ) ) %>%
                    dplyr::mutate( count_gene_mappings = ifelse(is.na(count_gene_mappings), 0, count_gene_mappings) ) %>%
                    dplyr::mutate( gene_mapping_source = unique(unique_probe_ids$gene_mapping_source) )

norm_data_matrix <- norm_data_matrix %>% 
  dplyr::relocate(dplyr::all_of(FINAL_COLUMN_ORDER))

write.csv(norm_data_matrix, file.path(DIR_NORMALIZED_EXPRESSION, "normalized_expression_GLmicroarray.csv"), row.names = FALSE)
```

**Input Data:**

- `df_rs` (R dataframe containing information from the runsheet, output from [Step 2c](#2c-load-metadata-and-raw-data))
- `annot_key` (keytype to join annotation table and microarray probes, dependent on organism, e.g. mus musculus uses 'ENSEMBL', defined in [Step 7a](#7a-get-probe-annotations))
- `background_corrected_data` (R object containing background-corrected microarray data, output from [Step 4](#4-background-correction))
- `norm_data` (R object containing background-corrected and normalized microarray data, output from [Step 5](#5-between-array-normalization))
- `unique_probe_ids` (R object containing probe ID to gene annotation mappings, output from [Step 7a](#7a-get-probe-annotations))

**Output Data:**

- **raw_intensities_GLmicroarray.csv** (table containing the background corrected unnormalized intensity values for each sample including gene annotations)
- **normalized_expression_GLmicroarray.csv** (table containing the background corrected, normalized intensity values for each sample including gene annotations)

## 8. Perform Probe Differential Expression (DE)

> Note: Run differential expression analysis only if there are at least 2 replicates per factor group.

<br>

### 8a. Generate Design Matrix

```R
# Loading metadata from runsheet csv file
design_data <- runsheet_to_design_matrix(runsheet)
design <- design_data$matrix

# Write SampleTable.csv and contrasts.csv file
write.csv(design_data$groups, file.path(DIR_DGE, "SampleTable_GLmicroarray.csv"), row.names = FALSE)
write.csv(design_data$contrasts, file.path(DIR_DGE, "contrasts_GLmicroarray.csv"))
```

**Custom Functions Used:**

- [runsheet_to_design_matrix()](#runsheet_to_design_matrix)

**Input Data:**

- `runsheet` (Path to runsheet, output from [Step 1](#1-create-sample-runsheet))

**Output Data:**

- `design_data` (a list of R objects containing the sample information and metadata
  - `design_data$matrix` - a design (or model) matrix describing the conditions in the dataset
  - `design_data$mapping` - a dataframe mapping the human-readable group names to the names of the conditions modified for use in R
  - `design_data$groups` - a dataframe of group names and contrasts for each sample
  - `design_data$contrasts` - a matrix of all pairwise comparisons of the groups)
- `design` (R object containing the limma study design matrix, indicating the group that each sample belongs to)
- **SampleTable_GLmicroarray.csv** (table containing samples and their respective groups)
- **contrasts_GLmicroarray.csv** (table containing all pairwise comparisons)

<br>

### 8b. Perform Individual Probe Level DE

```R
# Calculate results
res <- lm_fit_pairwise(norm_data, design)

# Print DE table, without filtering
limma::write.fit(res, adjust = 'BH', 
                file = "INTERIM.csv",
                row.names = FALSE,
                quote = TRUE,
                sep = ",")
```

**Custom Functions Used:**

- [lm_fit_pairwise()](#lm_fit_pairwise)

**Input Data:**

- `norm_data` (R object containing background-corrected and normalized microarray data created in [Step 5](#5-between-array-normalization))
- `design` (R object containing the limma study design matrix, indicating the group that each sample belongs to, created in [Step 8a](#8a-generate-design-matrix) above)

**Output Data:**

- INTERIM.csv (Statistical values from individual probe level DE analysis, including:
  - Log2fc between all pairwise comparisons
  - T statistic for all pairwise comparison tests
  - P value for all pairwise comparison tests)

<br>

### 8c. Add Annotation and Stats Columns and Format DE Table

```R
## Reformat Table for consistency across DE analyses tables within GeneLab ##

# Read in DE table 
df_interim <- read.csv("INTERIM.csv")

print("Remove extra columns from final table")

# These columns are data mapped to column PROBEID as per the original Manufacturer and can be linked as needed
colnames_to_remove = c(
  "Genes.Row",
  "Genes.Col",
  "Genes.Start",
  "Genes.Sequence",
  "Genes.ControlType",
  "Genes.GeneName",
  "Genes.SystematicName",
  "Genes.Description",
  "AveExpr" # Replaced by 'All.mean' column
)

df_interim <- df_interim %>% dplyr::select(-any_of(colnames_to_remove))

df_interim <- df_interim %>% dplyr::rename_with(reformat_names, .cols = matches('\\.condition|^Genes\\.'), group_name_mapping = design_data$mapping)


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

df_interim <- df_interim %>% 
  dplyr::mutate( 
    "All.mean" := rowMeans(dplyr::select(., all_of(SAMPLE_COLUMN_ORDER))),
    "All.stdev" := matrixStats::rowSds(as.matrix(dplyr::select(., all_of(SAMPLE_COLUMN_ORDER)))),
    ) %>% 
  dplyr::ungroup() %>%
  as.data.frame()

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

GROUP_MEAN_STDEV_COLUMNS_ORDER <- generate_prefixed_column_order(
  subjects = unique(design_data$groups$group),
  prefixes = c(
    "Group.Mean_",
    "Group.Stdev_"
  )
)

FINAL_COLUMN_ORDER <- c(
  ANNOTATIONS_COLUMN_ORDER, 
  PROBE_INFO_COLUMN_ORDER, 
  SAMPLE_COLUMN_ORDER, 
  STAT_COLUMNS_ORDER, 
  ALL_SAMPLE_STATS_COLUMNS_ORDER, 
  GROUP_MEAN_STDEV_COLUMNS_ORDER
  )

## Assert final column order includes all columns from original table
if (!setequal(FINAL_COLUMN_ORDER, colnames(df_interim))) {
  write.csv(FINAL_COLUMN_ORDER, "FINAL_COLUMN_ORDER.csv")
  NOT_IN_DF_INTERIM <- paste(setdiff(FINAL_COLUMN_ORDER, colnames(df_interim)), collapse = ":::")
  NOT_IN_FINAL_COLUMN_ORDER <- paste(setdiff(colnames(df_interim), FINAL_COLUMN_ORDER), collapse = ":::")
  stop(glue::glue("Column reordering attempt resulted in different sets of columns than original. Names unique to 'df_interim': {NOT_IN_FINAL_COLUMN_ORDER}. Names unique to 'FINAL_COLUMN_ORDER': {NOT_IN_DF_INTERIM}."))
}

## Perform reordering
df_interim <- df_interim %>% dplyr::relocate(dplyr::all_of(FINAL_COLUMN_ORDER))

# Save to file
write.csv(df_interim, file.path(DIR_DGE, "differential_expression_GLmicroarray.csv"), row.names = FALSE)
```

**Custom Functions Used:**

- [reformat_names()](#reformat_names)
- [generate_prefixed_column_order()](#generate_prefixed_column_order)

**Input Data:**

- `design_data` (a list of R objects containing the sample information and metadata, output from [Step 8a](#8a-generate-design-matrix) above)
- INTERIM.csv (Statistical values from individual probe level DE analysis, output from [Step 8b](#8b-perform-individual-probe-level-de) above)
- `norm_data` (R object containing background-corrected and normalized microarray data created in [Step 5](#5-between-array-normalization))

**Output Data:**

- **differential_expression_GLmicroarray.csv** (table containing normalized expression for each sample, group statistics, Limma probe DE results for each pairwise comparison, and gene annotations)


> All steps of the Microarray pipeline are performed using R markdown and the completed R markdown is rendered (via Quarto) as an html file (**NF_MAAgilent1ch_v\*_GLmicroarray.html**) and published in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/) for the respective dataset.
