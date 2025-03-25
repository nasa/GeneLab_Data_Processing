# GeneLab bioinformatics processing pipeline for Affymetrix microarray data <!-- omit in toc -->

> **This page holds an overview and instructions for how GeneLab processes Affymetrix microarray datasets. Exact processing commands and GL-DPPD-7114 version used for specific GeneLab datasets (GLDS) are provided with their processed data in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo).**  
> 
> \* The pipeline detailed below is currently used for animal and *Arabidopsis thaliana* studies only, it will be updated soon for processing microbe microarray data and other plant data.

---

**Date:** February XX, 2025  
**Revision:** -A  
**Document Number:** GL-DPPD-7114-A   

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
|tibble|3.1.8|3.2.1|
|stringr|1.5.0|1.5.1|
|R.utils|2.12.2|2.12.3|
|oligo|1.58.0|1.70.0|
|limma|3.50.3|3.62.2|
|glue|1.6.2|1.8.0|
|biomaRt|2.50.0|2.62.0|
|matrixStats|0.63.0|1.5.0|
|dp_tools|1.3.4|1.3.5|
|Quarto|1.2.313|1.6.40|

MA Plots

- Added support for plotting HTAFeatureSet data

Custom Annotations

- Added ability to use custom gene annotations when annotations are not available in Biomart or Ensembl FTP for *Arabidopsis thaliana*, see [Step 8](#8-probeset-annotations)

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
    - [3d. Boxplots](#3d-boxplots)
  - [4. Background Correction](#4-background-correction)
  - [5. Between Array Normalization](#5-between-array-normalization)
  - [6. Normalized Data Quality Assessment](#6-normalized-data-quality-assessment)
    - [6a. Density Plot](#6a-density-plot)
    - [6b. Pseudo Image Plots](#6b-pseudo-image-plots)
    - [6c. MA Plots](#6c-ma-plots)
    - [6d. Boxplots](#6d-boxplots)
  - [7. Probeset Summarization](#7-probeset-summarization)
  - [8. Probeset Annotations](#8-probeset-annotations)
    - [8a. Get Probeset Annotations](#8a-get-probeset-annotations)
    - [8b. Summarize Gene Mapping](#8b-summarize-gene-mapping)
    - [8c. Generate Annotated Raw and Normalized Expression Tables](#8c-generate-annotated-raw-and-normalized-expression-tables)
  - [9. Perform Probeset Differential Expression (DE)](#9-perform-probeset-differential-expression-de)
    - [9a. Generate Design Matrix](#9a-generate-design-matrix)
    - [9b. Perform Individual Probeset Level DE](#9b-perform-individual-probeset-level-de)
    - [9c. Add Annotation and Stats Columns and Format DE Table](#9c-add-annotation-and-stats-columns-and-format-de-table)

---

# Software used  

|Program|Version|Relevant Links|
|:------|:------:|:-------------|
|R|4.4.2|[https://www.r-project.org/](https://www.r-project.org/)|
|DT|0.33|[https://github.com/rstudio/DT](https://github.com/rstudio/DT)|
|dplyr|1.1.4|[https://dplyr.tidyverse.org](https://dplyr.tidyverse.org)|
|tibble|3.2.1|[https://tibble.tidyverse.org](https://tibble.tidyverse.org)|
|stringr|1.5.1|[https://stringr.tidyverse.org](https://stringr.tidyverse.org)|
|R.utils|2.12.3|[https://github.com/HenrikBengtsson/R.utils](https://github.com/HenrikBengtsson/R.utils)|
|oligo|1.70.0|[https://bioconductor.org/packages/3.14/bioc/html/oligo.html](https://bioconductor.org/packages/3.14/bioc/html/oligo.html)|
|limma|3.62.2|[https://bioconductor.org/packages/3.14/bioc/html/limma.html](https://bioconductor.org/packages/3.14/bioc/html/limma.html)|
|glue|1.8.0|[https://glue.tidyverse.org](https://glue.tidyverse.org)|
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
> - Rather than running the commands below to create the runsheet needed for processing, the runsheet may also be created manually by following the [file specification](../Workflow_Documentation/NF_MAAffymetrix/examples/runsheet/README.md).
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

> Note: Steps 2 - 9 are done in R

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
BiocManager::install("oligo")


## Note: Only dplyr is explicitly loaded. Other library functions are called with explicit namespace (e.g. LIBRARYNAME::FUNCTION)
library(dplyr) # Ensure infix operator is available, methods should still reference dplyr namespace otherwise
options(dplyr.summarise.inform = FALSE) # Don't print out '`summarise()` has grouped output by 'group'. You can override using the `.groups` argument.'

# Define path to runsheet
runsheet <- "/path/to/runsheet/{OSD-Accession-ID}_microarray_v{version}_runsheet.csv"

## Set up output structure

# Output Constants
DIR_RAW_DATA <- "00-RawData"
DIR_NORMALIZED_EXPRESSION <- "01-oligo_NormExp"
DIR_DGE <- "02-limma_DGE"

dir.create(DIR_RAW_DATA)
dir.create(DIR_NORMALIZED_EXPRESSION)
dir.create(DIR_DGE)

## Save original par settings
##   Par may be temporarily changed for plotting purposes and reset once the plotting is done

original_par <- par()
options(preferRaster=TRUE) # use Raster when possible to avoid antialiasing artifacts in images

options(timeout=1000) # ensure enough time for data downloads
```

<br>

### 2b. Define Custom Functions

#### retry_with_delay()
<details>
  <summary>utility function to improve robustness of function calls; used to remedy intermittent internet issues during runtime</summary>

  ```R
  retry_with_delay <- function(func, ...) {
    max_attempts = 5
    initial_delay = 10
    delay_increase = 30
    attempt <- 1
    current_delay <- initial_delay
    while (attempt <= max_attempts) {
      result <- tryCatch(
        expr = func(...),
        error = function(e) e
      )

      if (!inherits(result, "error")) {
        return(result)
      } else {
        if (attempt < max_attempts) {
          message(paste("Retry attempt", attempt, "failed for function with name <", deparse(substitute(func)) ,">. Retrying in", current_delay, "second(s)..."))
          Sys.sleep(current_delay)
          current_delay <- current_delay + delay_increase
        } else {
          stop(paste("Max retry attempts reached. Last error:", result$message))
        }
      }

      attempt <- attempt + 1
    }
  }
  ```

  **Function Parameter Definitions:**
  - `func=` - specifies the function to wrap
  - `...` - other arguments passed on to the `func`

  **Returns:** the output of the wrapped function
</details>

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

#### shortened_organism_name()
<details>
  <summary>shortens organism names, for example 'Homo Sapiens' to 'hsapiens'</summary>

  ```R
  shortened_organism_name <- function(long_name) {
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
                          # Add column names as follows: MAPID, TAIR, PROBESETID
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

  **Returns:** a dataframe containing the mapping between Ensembl ID and probeset ID, as obtained via FTP
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
      df <- read.csv(runsheet, check.names = FALSE) %>% 
                dplyr::mutate_all(function(x) iconv(x, "latin1", "ASCII", sub="")) # Convert all characters to ascii, when not possible, remove the character    # get only Factor Value columns
      factors = as.data.frame(df[,grep("Factor.Value", colnames(df), ignore.case=TRUE)])
      colnames(factors) = paste("factor",1:dim(factors)[2], sep= "_")
      
      # Load metadata from runsheet csv file
      compare_csv = data.frame(sample_id = df[,c("Sample Name")], factors)

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
      rownames(sampleTable) <- df[,c("Sample Name")]

      condition <- sampleTable[,'condition']
      names_mapping <- as.data.frame(cbind(safe_name = as.character(condition), original_name = group_names))

      design <- model.matrix(~ 0 + condition)
      design_data <- list( matrix = design, mapping = names_mapping, groups = as.data.frame( cbind(sample = df[,c("Sample Name")], group = group_names) ), contrasts = contrasts )
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
df_rs <- read.csv(runsheet, check.names = FALSE) %>% 
          dplyr::mutate_all(function(x) iconv(x, "latin1", "ASCII", sub="")) # Convert all characters to ascii, when not possible, remove the character

if ( runsheet_paths_are_URIs(df_rs) ) {
  print("Determined Raw Data Locations are URIS")
  local_paths <- retry_with_delay(download_files_from_runsheet, df_rs)
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
# Retry with delay here to accomodate oligo's automatic loading of annotation packages and occasional internet related failures to load
raw_data <- retry_with_delay(
              oligo::read.celfiles,
              df_local_paths$`Local Paths`,
              sampleNames = df_local_paths$`Sample Name`# Map column names as Sample Names (instead of default filenames)
            )

# Summarize raw data
print(paste0("Number of Arrays: ", dim(raw_data)[2]))
print(paste0("Number of Probes: ", dim(raw_data)[1]))
```

**Custom Functions Used:**

- [retry_with_delay()](#retry_with_delay)
- [all_true()](#all_true)
- [runsheet_paths_are_URIs()](#runsheet_paths_are_uris)
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

annotation_table <- retry_with_delay(fetch_organism_specific_annotation_table, unique(df_rs$organism), annotation_table_link)

annotation_file_path <- annotation_table$genelab_annots_link
ensembl_version <- as.character(annotation_table$ensemblVersion)
```

**Custom Functions Used:**

- [retry_with_delay()](#retry_with_delay)
- [fetch_organism_specific_annotation_table()](#fetch_organism_specific_annotation_table)

**Input Data:**

- `local_annotation_dir` (Path to local annotation directory if using custom annotations, see [Step 8a](#8a-get-probeset-annotations))

    > Note: If not using custom annotations, leave `local_annotation_dir` as `NULL`.

- `annotation_config_path` (URL or path to annotation config file if using custom annotations, see [Step 8a](#8a-get-probeset-annotations))

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
scale_factor = 0.2 # Default scale factor

if (max(nchar(colnames(raw_data@assayData$exprs))) > 35 & number_of_sets > 1) { # Scale more if sample names are long
  scale_factor = if_else(number_of_sets == 2, 0.4, 0.25)
}

oligo::hist(raw_data, 
            transfo=log2, # Log2 transform raw intensity values
            which=c("all"),
            nsample=10000, # Number of probes to plot
            main = "Density of raw intensities for multiple arrays")
legend("topright", legend = colnames(raw_data@assayData$exprs),
        lty = c(1,2,3,4,5), # Seems like oligo::hist cycles through these first five line types
        col = oligo::darkColors(n = ncol(raw_data)), # Ensure legend color is in sync with plot
        ncol = number_of_sets, # Set number of columns by number of sets
        cex = max(0.35, 1 + scale_factor - (number_of_sets*scale_factor)) # Reduce for each column beyond 1 with minimum of 35%
      )

# Reset par
par(original_par)
```

**Input Data:**

- `raw_data` (raw data R object created in [Step 2c](#2c-load-metadata-and-raw-data) above)

**Output Data:**

- Plot containing the density of raw intensities for each array (lack of overlap indicates a need for normalization)

<br>

### 3b. Pseudo Image Plots

```R
for ( i in seq_along(1:ncol(raw_data))) {
  oligo::image(raw_data[,i], 
    transfo = log2,
    main = colnames(raw_data)[i]
    )
}
```

**Input Data:**

- `raw_data` (raw data R object created in [Step 2c](#2c-load-metadata-and-raw-data) above)

**Output Data:**

- Pseudo images of each array before background correction and normalization

<br>

### 3c. MA Plots

```R
if (inherits(raw_data, "GeneFeatureSet")) {
  print("Raw data is a GeneFeatureSet, using exprs() to access expression values and adding 0.0001 to avoid log(0)")
} else if (inherits(raw_data, "ExpressionSet") || inherits(raw_data, "ExpressionFeatureSet") || inherits(raw_data, "HTAFeatureSet")) { 
  print(paste0("Raw data is ", class(raw_data), ". Using default approach for this class for MA Plot"))
}

if (inherits(raw_data, "GeneFeatureSet")) {
  MA_plot <- oligo::MAplot(
    exprs(raw_data) + 0.0001,
    transfo=log2,
    ylim=c(-2, 4),
    main="" # This function uses 'main' as a suffix to the sample name. Here we want just the sample name, thus here main is an empty string
  )
} else if (inherits(raw_data, "ExpressionSet") || inherits(raw_data, "ExpressionFeatureSet") || inherits(raw_data, "HTAFeatureSet")) { 
  MA_plot <- oligo::MAplot(
    raw_data,
    ylim=c(-2, 4),
    main="" # This function uses 'main' as a suffix to the sample name. Here we want just the sample name, thus here main is an empty string
  )
} else {
  stop(glue::glue("No strategy for MA plots for {class(raw_data)}"))
}
```

**Input Data:**

- `raw_data` (raw data R object created in [Step 2c](#2c-load-metadata-and-raw-data) above)

**Output Data:**

- `MA_plot` (M (log ratio of the subject array vs a pseudo-reference, the median of all other arrays) vs. A (average log expression) plot for each array before background correction and normalization)

<br>


### 3d. Boxplots

```R
max_samplename_length <- max(nchar(colnames(raw_data)))
dynamic_lefthand_margin <- max(max_samplename_length * 0.7, 10)
par(
  mar = c(8, dynamic_lefthand_margin, 8, 2) + 0.1, # mar is the margin around the plot. c(bottom, left, top, right)
  xpd = TRUE
  ) 
boxplot <- oligo::boxplot(raw_data[, rev(colnames(raw_data))], # Here we reverse column order to ensure descending order for samples in horizontal boxplot 
                          transfo=log2, # Log2 transform raw intensity values
                          which=c("all"),
                          nsample=10000, # Number of probes to plot
                          las = 1, # las specifies the orientation of the axis labels. 1 = always horizontal
                          ylab="",
                          xlab="log2 Intensity",
                          main = "Boxplot of raw intensities \nfor perfect match and mismatch probes",
                          horizontal = TRUE
                          )
title(ylab = "Sample Name", mgp = c(dynamic_lefthand_margin-2, 1, 0))
# Reset par
par(original_par)
```

**Input Data:**

- `raw_data` (raw data R object created in [Step 2c](#2c-load-metadata-and-raw-data) above)

**Output Data:**

- `boxplot` (Boxplot of raw expression data for each array before background correction and normalization)

<br>

---

## 4. Background Correction

```R
background_corrected_data <- raw_data %>% oligo::backgroundCorrect(method="rma")
```

**Input Data:**

- `raw_data` (raw data R object created in [Step 2c](#2c-load-metadata-and-raw-data) above)

**Output Data:**

- `background_corrected_data` (R object containing background-corrected microarray data)

  >   
  > Note: Background correction was performed using the oligo `rma` method, specifically "Convolution Background Correction"

<br>

---

## 5. Between Array Normalization

```R
# Normalize background-corrected data using the quantile method
norm_data <- oligo::normalize(background_corrected_data, 
                              method = "quantile",
                              target = "core" # Use oligo default: core metaprobeset mappings
                              )
                              
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

oligo::hist(norm_data,
            transfo=log2, # Log2 transform normalized intensity values
            which=c("all"),
            nsample=10000, # Number of probes to plot
            main = "Density of normalized intensities for multiple arrays")
legend("topright", legend = colnames(norm_data@assayData$exprs),
        lty = c(1,2,3,4,5), # Seems like oligo::hist cycles through these first five line types
        col = oligo::darkColors(n = ncol(norm_data)), # Ensure legend color is in sync with plot
        ncol = number_of_sets, # Set number of columns by number of sets
        cex = max(0.35, 1 + scale_factor - (number_of_sets*scale_factor)) # Reduce for each column beyond 1 with minimum of 35%
      )

# Reset par
par(original_par)
```

**Input Data:**

- `norm_data` (R object containing background-corrected and normalized microarray data created in [Step 5](#5-between-array-normalization) above)

**Output Data:**

- Plot containing the density of background-corrected and normalized intensities for each array (near complete overlap is expected after normalization)

<br>

### 6b. Pseudo Image Plots

```R
for ( i in seq_along(1:ncol(norm_data))) {
  oligo::image(norm_data[,i], 
    transfo = log2,
    main = colnames(norm_data)[i]
    )
}
```

**Input Data:**

- `norm_data` (R object containing background-corrected and normalized microarray data created in [Step 5](#5-between-array-normalization) above)

**Output Data:**

- Pseudo images of each array after background correction and normalization 

<br>

### 6c. MA Plots

```R
MA_plot <- oligo::MAplot(
    norm_data, 
    ylim=c(-2, 4),
    main="" # This function uses 'main' as a suffix to the sample name. Here we want just the sample name, thus here main is an empty string
)
```

**Input Data:**

- `norm_data` (R object containing background-corrected and normalized microarray data created in [Step 5](#5-between-array-normalization) above)

**Output Data:**

- `MA_plot` (M (log ratio of the subject array vs a pseudo-reference, the median of all other arrays) vs. A (average log expression) plot for each array after background correction and normalization)

<br>

### 6d. Boxplots

```R
max_samplename_length <- max(nchar(colnames(norm_data)))
dynamic_lefthand_margin <- max(max_samplename_length * 0.7, 10)
par(
  mar = c(8, dynamic_lefthand_margin, 8, 2) + 0.1, # mar is the margin around the plot. c(bottom, left, top, right)
  xpd = TRUE
  ) 
boxplot <- oligo::boxplot(norm_data[, rev(colnames(norm_data))], # Here we reverse column order to ensure descending order for samples in horizontal boxplot
                          transfo=log2, # Log2 transform normalized intensity values
                          which=c("all"),
                          nsample=10000, # Number of probes to plot
                          las = 1, # las specifies the orientation of the axis labels. 1 = always horizontal
                          ylab="",
                          xlab="log2 Intensity",
                          main = "Boxplot of normalized intensities \nfor perfect match and mismatch probes",
                          horizontal = TRUE
                          )
title(ylab = "Sample Name", mgp = c(dynamic_lefthand_margin-2, 1, 0))
# Reset par
par(original_par)
```

**Input Data:**

- `norm_data` (R object containing background-corrected and normalized microarray data created in [Step 5](#5-between-array-normalization) above)

**Output Data:**

- `boxplot` (Boxplot of expression data for each array after background correction and normalization)

<br>

---

## 7. Probeset Summarization

```R
probeset_level_data <- oligo::rma(norm_data, 
                                  normalize=FALSE, 
                                  background=FALSE
                                  )

# Summarize background-corrected and normalized data
print("Summarized Probeset Level Data Below")
print(paste0("Number of Arrays: ", dim(probeset_level_data)[2]))
print(paste0("Total Number of Probes Assigned To A Probeset: ", dim(oligo::getProbeInfo(probeset_level_data, target="core")['man_fsetid'])[1])) # man_fsetid means 'Manufacturer Probeset ID'. Ref: https://support.bioconductor.org/p/57191/
print(paste0("Number of Probesets: ", dim(unique(oligo::getProbeInfo(probeset_level_data, target="core")['man_fsetid']))[1])) # man_fsetid means 'Manufacturer Probeset ID'. Ref: https://support.bioconductor.org/p/57191/
```

**Input Data:**

- `norm_data` (R object containing background-corrected and normalized microarray data created in [Step 5](#5-between-array-normalization) above)

**Output Data:**

- `probeset_level_data` (R object containing probeset level expression values after summarization of normalized probeset level data)

<br>

---

## 8. Probeset Annotations

<br>

### 8a. Get Probeset Annotations

```R
organism <- shortened_organism_name(unique(df_rs$organism))
annot_key <- ifelse(organism %in% c("athaliana"), 'TAIR', 'ENSEMBL')

if (organism %in% c("athaliana")) {
  ENSEMBL_VERSION = ensembl_version
  ensembl_genomes_portal = "plants"
  print(glue::glue("Using ensembl genomes ftp to get specific version of probeset id mapping table. Ensembl genomes portal: {ensembl_genomes_portal}, version: {ENSEMBL_VERSION}"))
  expected_attribute_name <- get_biomart_attribute(df_rs)
  df_mapping <- retry_with_delay(
      get_ensembl_genomes_mappings_from_ftp,
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

    ensembl <- biomaRt::useEnsembl(biomart = "genes", 
                                  dataset = expected_dataset_name,
                                  version = ENSEMBL_VERSION)
    print(ensembl)

    # Some probe_ids for affy_hta_2_0 may end in .hg.1 instead of .hg (how it is in biomaRt), leading to 0 results returned
    if (expected_attribute_name == 'affy_hta_2_0') {
      rownames(probeset_level_data) <- stringr::str_replace(rownames(probeset_level_data), '\\.hg\\.1$', '.hg')
    }

    probe_ids <- rownames(probeset_level_data)

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
# If no df_mapping obtained (e.g., organism not supported in biomart), use custom annotations; otherwise, merge in-house annotations to df_mapping

if (use_custom_annot) {
  expected_attribute_name <- 'ProbesetID'

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

  if (annot_type == '3prime-IVT') {
    unique_probe_ids <- read.csv(
      file.path(local_annotation_dir, annot_config$annot_filename[[1]]),
      skip = 13, header = TRUE, na.strings = c('NA', '---')
    )[c('Probe.Set.ID', 'Entrez.Gene', 'Gene.Symbol', 'Gene.Title', 'Ensembl', 'RefSeq.Transcript.ID', 'RefSeq.Protein.ID', 'Gene.Ontology.Biological.Process', 'Gene.Ontology.Cellular.Component', 'Gene.Ontology.Molecular.Function')]

    # Clean columns
    unique_probe_ids$Gene.Symbol <- purrr::map_chr(stringr::str_split(unique_probe_ids$Gene.Symbol, stringr::fixed(' /// ')), ~paste0(unique(.), collapse = "|")) %>% stringr::str_replace('NA', NA_character_)
    unique_probe_ids$Gene.Title <- purrr::map_chr(stringr::str_split(unique_probe_ids$Gene.Title, stringr::fixed(' /// ')), ~paste0(unique(.), collapse = "|")) %>% stringr::str_replace('NA', NA_character_)
    unique_probe_ids$Entrez.Gene <- purrr::map_chr(stringr::str_split(unique_probe_ids$Entrez.Gene, stringr::fixed(' /// ')), ~paste0(unique(.), collapse = "|")) %>% stringr::str_replace('NA', NA_character_)
    unique_probe_ids$Ensembl <- purrr::map_chr(stringr::str_split(unique_probe_ids$Ensembl, stringr::fixed(' /// ')), ~paste0(unique(.), collapse = "|")) %>% stringr::str_replace('NA', NA_character_)

    unique_probe_ids$RefSeq <- paste(unique_probe_ids$RefSeq.Transcript.ID, unique_probe_ids$RefSeq.Protein.ID)
    unique_probe_ids$RefSeq <- purrr::map_chr(stringr::str_extract_all(unique_probe_ids$RefSeq, '[A-Z]+_[\\d.]+'), ~paste0(unique(.), collapse = "|")) %>% stringr::str_replace('^$', NA_character_)

    unique_probe_ids$GO <- paste(unique_probe_ids$Gene.Ontology.Biological.Process, unique_probe_ids$Gene.Ontology.Cellular.Component, unique_probe_ids$Gene.Ontology.Molecular.Function)
    unique_probe_ids$GO <- purrr::map_chr(stringr::str_extract_all(unique_probe_ids$GO, '\\d{7}'), ~paste0('GO:', unique(.), collapse = "|")) %>% stringr::str_replace('^GO:$', NA_character_)

    unique_probe_ids <- unique_probe_ids[c('Probe.Set.ID', 'Entrez.Gene', 'Gene.Symbol', 'Gene.Title', 'Ensembl', 'RefSeq', 'GO')]
    names(unique_probe_ids) <- c('ProbesetID', 'ENTREZID', 'SYMBOL', 'GENENAME', 'ENSEMBL', 'REFSEQ', 'GOSLIM_IDS')

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
    annot_cols <- c('ProbesetID', 'ENTREZID', 'SYMBOL', 'GENENAME', 'ENSEMBL', 'REFSEQ', 'GOSLIM_IDS', 'STRING_id', 'count_gene_mappings', 'gene_mapping_source')
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
                        dplyr::mutate(dplyr::across(!!sym(expected_attribute_name), as.character)) %>% # Ensure probeset ids treated as character type
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

probeset_expression_matrix <- oligo::exprs(probeset_level_data)

probeset_expression_matrix.gene_mapped <- probeset_expression_matrix %>% 
  as.data.frame() %>%
  tibble::rownames_to_column(var = "ProbesetID") %>% # Ensure rownames (probeset IDs) can be used as join key
  dplyr::left_join(unique_probe_ids, by = c("ProbesetID" = expected_attribute_name ) ) %>%
  dplyr::mutate( count_gene_mappings := ifelse(is.na(count_gene_mappings), 0, count_gene_mappings) ) %>%
  dplyr::mutate( gene_mapping_source := unique(unique_probe_ids$gene_mapping_source) )
```

**Custom Functions Used:**

- [retry_with_delay()](#retry_with_delay)
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

  > Note: See [Affymetrix_array_annotations.csv](../Array_Annotations/Affymetrix_array_annotations.csv) for the latest config file used at GeneLab. This file can also be created manually by following the [file specification](../Workflow_Documentation/NF_MAAffymetrix/examples/annotations/README.md).

- `probeset_level_data` (R object containing probeset level expression values after summarization of normalized probeset level data, output from [Step 7](#7-probeset-summarization))

**Output Data:**

- `unique_probe_ids` (R object containing probeset ID to gene annotation mappings)
- `probeset_expression_matrix.gene_mapped` (R object containing probeset level expression values after summarization of normalized probeset level data combined with gene annotations specified by [Biomart](https://bioconductor.org/packages/3.14/bioc/html/biomaRt.html) or custom annotations)

<br>

### 8b. Summarize Gene Mapping

```R
# Pie Chart with Percentages
slices <- c(
    'Unique Mapping' = nrow(probeset_expression_matrix.gene_mapped %>% dplyr::filter(count_gene_mappings == 1) %>% dplyr::distinct(ProbesetID)), 
    'Multi Mapping' = nrow(probeset_expression_matrix.gene_mapped %>% dplyr::filter(count_gene_mappings > 1) %>% dplyr::distinct(ProbesetID)), 
    'No Mapping' = nrow(probeset_expression_matrix.gene_mapped %>% dplyr::filter(count_gene_mappings == 0) %>% dplyr::distinct(ProbesetID))
)
pct <- round(slices/sum(slices)*100)
chart_names <- names(slices)
chart_names <- glue::glue("{names(slices)} ({slices})") # add count to labels
chart_names <- paste(chart_names, pct) # add percents to labels
chart_names <- paste(chart_names,"%",sep="") # ad % to labels
pie(slices,labels = chart_names, col=rainbow(length(slices)),
    main=glue::glue("Mapping to Primary Keytype\n {nrow(probeset_expression_matrix.gene_mapped %>% dplyr::distinct(ProbesetID))} Total Unique Probesets")
    )

print(glue::glue("Unique Mapping Count: {slices[['Unique Mapping']]}"))
```

**Input Data:**

- `probeset_expression_matrix.gene_mapped` (R object containing probeset level expression values after summarization of normalized probeset level data combined with gene annotations specified by [Biomart](https://bioconductor.org/packages/3.14/bioc/html/biomaRt.html) or custom annotations, output from [Step 8a](#8a-get-probeset-annotations) above)

**Output Data:**

- A pie chart denoting the gene mapping rates for each unique probeset ID
- A printout denoting the count of unique mappings for gene mapping

<br>

### 8c. Generate Annotated Raw and Normalized Expression Tables

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

SAMPLE_COLUMN_ORDER <- df_rs$`Sample Name`

probeset_expression_matrix.gene_mapped <- probeset_expression_matrix.gene_mapped %>% dplyr::rename( !!annot_key := ENSEMBL )

## Output column subset file with just normalized probeset level expression values
write.csv(
  probeset_expression_matrix.gene_mapped[c(
  ANNOTATIONS_COLUMN_ORDER,
  "ProbesetID",
  "count_gene_mappings",
  "gene_mapping_source",
  SAMPLE_COLUMN_ORDER)
  ], file.path(DIR_NORMALIZED_EXPRESSION, "normalized_expression_probeset_GLmicroarray.csv"), row.names = FALSE)

## Determine column order for probe level tables

PROBE_INFO_COLUMN_ORDER = c(
  "ProbesetID",
  "ProbeID",
  "count_gene_mappings",
  "gene_mapping_source"
)

FINAL_COLUMN_ORDER <- c(
  ANNOTATIONS_COLUMN_ORDER, 
  PROBE_INFO_COLUMN_ORDER, 
  SAMPLE_COLUMN_ORDER
)

## Generate raw intensity matrix that includes annotations

background_corrected_data_annotated <- oligo::exprs(background_corrected_data) %>% 
  as.data.frame() %>%
  tibble::rownames_to_column(var = "fid") %>% # Ensure rownames (probeset IDs) can be used as join key
  dplyr::mutate(dplyr::across(fid, as.integer)) %>% # Ensure fid is integer type, consistent with getProbeInfo typing
  dplyr::right_join(oligo::getProbeInfo(background_corrected_data), by = "fid") %>% # Add 'man_fsetid' via mapping based on fid
  dplyr::rename( ProbesetID = man_fsetid ) %>% # Rename from getProbeInfo name to ProbesetID
  dplyr::rename( ProbeID = fid ) %>% # Rename from getProbeInfo name to ProbeID
  dplyr::left_join(unique_probe_ids, by = c("ProbesetID" = expected_attribute_name ) ) %>% # Join with biomaRt ENSEMBL mappings
  dplyr::mutate( count_gene_mappings := ifelse(is.na(count_gene_mappings), 0, count_gene_mappings) ) %>% # Convert NA mapping to 0
  dplyr::mutate( gene_mapping_source := unique(unique_probe_ids$gene_mapping_source) ) %>%
  dplyr::rename( !!annot_key := ENSEMBL )

## Perform reordering
background_corrected_data_annotated <- background_corrected_data_annotated %>% 
  dplyr::relocate(dplyr::all_of(FINAL_COLUMN_ORDER))

write.csv(background_corrected_data_annotated, file.path(DIR_RAW_DATA, "raw_intensities_probe_GLmicroarray.csv"), row.names = FALSE)

## Generate normalized expression matrix that includes annotations
norm_data_matrix_annotated <- oligo::exprs(norm_data) %>% 
  as.data.frame() %>%
  tibble::rownames_to_column(var = "fid") %>% # Ensure rownames (probeset IDs) can be used as join key
  dplyr::mutate(dplyr::across(fid, as.integer)) %>% # Ensure fid is integer type, consistent with getProbeInfo typing
  dplyr::right_join(oligo::getProbeInfo(norm_data), by = "fid") %>% # Add 'man_fsetid' via mapping based on fid
  dplyr::rename( ProbesetID = man_fsetid ) %>% # Rename from getProbeInfo name to ProbesetID
  dplyr::rename( ProbeID = fid ) %>% # Rename from getProbeInfo name to ProbeID
  dplyr::left_join(unique_probe_ids, by = c("ProbesetID" = expected_attribute_name ) ) %>%
  dplyr::mutate( count_gene_mappings := ifelse(is.na(count_gene_mappings), 0, count_gene_mappings) ) %>% # Convert NA mapping to 0
  dplyr::mutate( gene_mapping_source := unique(unique_probe_ids$gene_mapping_source) ) %>%
  dplyr::rename( !!annot_key := ENSEMBL ) 

norm_data_matrix_annotated <- norm_data_matrix_annotated %>% 
  dplyr::relocate(dplyr::all_of(FINAL_COLUMN_ORDER))

write.csv(norm_data_matrix_annotated, file.path(DIR_NORMALIZED_EXPRESSION, "normalized_intensities_probe_GLmicroarray.csv"), row.names = FALSE)
```

**Input Data:**

- `df_rs` (R dataframe containing information from the runsheet, output from [Step 2c](#2c-load-metadata-and-raw-data))
- `annot_key` (keytype to join annotation table and microarray probes, dependent on organism, e.g. mus musculus uses 'ENSEMBL', defined in [Step 8a](#8a-get-probeset-annotations))
- `probeset_expression_matrix.gene_mapped` (R object containing probeset level expression values after summarization of normalized probeset level data combined with gene annotations specified by [Biomart](https://bioconductor.org/packages/3.14/bioc/html/biomaRt.html) or custom annotations, output from [Step 8a](#8a-get-probeset-annotations) above)
- `background_corrected_data` (R object containing background-corrected microarray data, output from [Step 4](#4-background-correction))
- `norm_data` (R object containing background-corrected and normalized microarray data, output from [Step 5](#5-between-array-normalization))
- `unique_probe_ids` (R object containing probeset ID to gene annotation mappings, output from [Step 8a](#8a-get-probeset-annotations))

**Output Data:**

- **normalized_expression_probeset_GLmicroarray.csv** (table containing the background corrected, normalized probeset expression values for each sample. The ProbesetID is the unique index column.)
- **raw_intensities_probe_GLmicroarray.csv** (table containing the background corrected, unnormalized probe intensity values for each sample including gene annotations. The ProbeID is the unique index column.)
- **normalized_intensities_probe_GLmicroarray.csv** (table containing the background corrected, normalized probe intensity values for each sample including gene annotations.  The ProbeID is the unique index column.)

## 9. Perform Probeset Differential Expression (DE)

> Note: Run differential expression analysis only if there are at least 2 replicates per factor group.

<br>

### 9a. Generate Design Matrix

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

- `runsheet` (path to runsheet, output from [Step 1](#1-create-sample-runsheet))

**Output Data:**

- `design_data` (a list of R objects containing the sample information and metadata
  - `design_data$matrix` - the limma study design matrix, indicating the group that each sample belongs to
  - `design_data$mapping` - a dataframe of conditions and group names
  - `design_data$groups` - a dataframe of samples and group names
  - `design_data$contrasts` - a matrix of all pairwise comparisons of the groups)
- `design` (R object containing the limma study design matrix, indicating the group that each sample belongs to)
- **SampleTable_GLmicroarray.csv** (table containing samples and their respective groups)
- **contrasts_GLmicroarray.csv** (table containing all pairwise comparisons)

<br>

### 9b. Perform Individual Probeset Level DE

```R
# Calculate results
res <- lm_fit_pairwise(probeset_level_data, design)

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

- `probeset_level_data` (R object containing probeset level expression values after summarization of normalized probeset level data, output from [Step 7](#7-probeset-summarization))
- `design` (R object containing the limma study design matrix, indicating the group that each sample belongs to, output from [Step 9a](#9a-generate-design-matrix) above)

**Output Data:**

- INTERIM.csv (statistical values from individual probeset level DE analysis, including:
  - Log2fc between all pairwise comparisons
  - T statistic for all pairwise comparison tests
  - P value for all pairwise comparison tests
  - Adjusted P value for all pairwise comparison tests)

<br>

### 9c. Add Annotation and Stats Columns and Format DE Table

```R
## Reformat Table for consistency across DE analyses tables within GeneLab ##

# Read in DE table 
df_interim <- read.csv("INTERIM.csv")

# Bind columns from gene mapped expression table
df_interim <- df_interim %>% 
  dplyr::bind_cols(probeset_expression_matrix.gene_mapped)

df_interim <- df_interim %>% dplyr::rename_with(reformat_names, .cols = matches('\\.condition'), group_name_mapping = design_data$mapping)


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

print("Remove extra columns from final table")

# These columns are data mapped to column PROBEID as per the original Manufacturer and can be linked as needed
colnames_to_remove = c(
  "AveExpr" # Replaced by 'All.mean' column
)

df_interim <- df_interim %>% dplyr::select(-any_of(colnames_to_remove))

PROBE_INFO_COLUMN_ORDER = c(
  "ProbesetID",
  "count_gene_mappings",
  "gene_mapping_source"
)

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

- `design_data` (a list of R objects containing the sample information and metadata, output from [Step 9a](#9a-generate-design-matrix) above)
- INTERIM.csv (statistical values from individual probeset level DE analysis, output from [Step 9b](#9b-perform-individual-probeset-level-de) above)
- `probeset_expression_matrix.gene_mapped` (R object containing probeset level expression values after summarization of normalized probeset level data combined with gene annotations specified by [Biomart](https://bioconductor.org/packages/3.14/bioc/html/biomaRt.html) or custom annotations, output from [Step 8a](#8a-get-probeset-annotations) above)

**Output Data:**

- **differential_expression_GLmicroarray.csv** (table containing normalized probeset expression values for each sample, group statistics, Limma probeset DE results for each pairwise comparison, and gene annotations. The ProbesetID is the unique index column.)

> All steps of the Microarray pipeline are performed using R markdown and the completed R markdown is rendered (via Quarto) as an html file (**NF_MAAffymetrix_v\*_GLmicroarray.html**) and published in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/) for the respective dataset.
