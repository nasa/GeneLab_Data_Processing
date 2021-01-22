# GeneLab bioinformatics processing pipeline for two-channel microarray data

> **This page holds an overview and example code of how GeneLab processes two-channel Microarray datasets. Exact processing code and GL-DPPD-710X revision used for specific datasets are available in the [GLDS_Processing_Scripts](https://developer.nasa.gov/asaravia/GeneLab_Data_Processing/tree/master/Microarray/2-channel_arrays/GLDS_Processing_Scripts) sub-directory and are also provided with their processed data in the [GeneLab Data Systems (GLDS) repository](https://genelab-data.ndc.nasa.gov/genelab/projects).**  

---

**Date:** Month day, 2020  
**Revision:** -   
**Document Number:** GL-DPPD-710X  

**Submitted by:**  
Homer Fogle (GeneLab Data Processing Team)

**Approved by:**  
Sylvain Costes (GeneLab Project Manager)  
Samrawit Gebre (GeneLab Deputy Project Manager and Interim GeneLab Configuration Manager)  
Amanda Saravia-Butler (GeneLab Data Processing Representative)  
Jonathan Galazka (GeneLab Project Scientist)

---

# Table of contents  

- [**Software used**](#software-used)
- [**General processing overview with example code**](#general-processing-overview-with-example-code)
  - [**Input Data and Parameters**](#input-data-and-parameters)
  - [**Output Data**](#output-data)
  
---

# Software used  

### Required R Packages available from the CRAN repository:  

|Program   | Relevant Links|
|:---------|:--------------|
|shiny     | [https://cran.r-project.org/web/packages/shiny/index.html](https://cran.r-project.org/web/packages/shiny/index.html)|
|rmarkdown | [https://cran.r-project.org/web/packages/rmarkdown/index.html](https://cran.r-project.org/web/packages/rmarkdown/index.html)|
|knitr     | [https://cran.r-project.org/web/packages/knitr/index.html](https://cran.r-project.org/web/packages/knitr/index.html)|
|xfun      | [https://cran.r-project.org/web/packages/xfun/index.html](https://cran.r-project.org/web/packages/xfun/index.html)|
|DT        | [https://cran.r-project.org/web/packages/DT/index.html](https://cran.r-project.org/web/packages/DT/index.html)|
|R.utils   | [https://cran.r-project.org/web/packages/R.utils/index.html](https://cran.r-project.org/web/packages/R.utils/index.html)|
|dplyr     | [https://cran.r-project.org/web/packages/dplyr/index.html](https://cran.r-project.org/web/packages/dplyr/index.html)|
|tidyr     | [https://cran.r-project.org/web/packages/tidyr/index.html](https://cran.r-project.org/web/packages/tidyr/index.html)|
|statmod   | [https://cran.r-project.org/web/packages/statmod/index.html](https://cran.r-project.org/web/packages/statmod/index.html)|


### Required R Packages available from the Bioconductor repository:

|Program       | Relevant Links|
|:-------------|:--------------|
|Risa          | [https://www.bioconductor.org/packages/release/bioc/html/Risa.html](https://www.bioconductor.org/packages/release/bioc/html/Risa.html)|
|limma         | [https://bioconductor.org/packages/release/bioc/html/limma.html](https://bioconductor.org/packages/release/bioc/html/limma.html)|
|GEOquery      | [https://bioconductor.org/packages/release/bioc/html/GEOquery.html](https://bioconductor.org/packages/release/bioc/html/GEOquery.html)|
|ArrayExpress  | [https://www.bioconductor.org/packages/release/bioc/html/ArrayExpress.html](https://www.bioconductor.org/packages/release/bioc/html/ArrayExpress.html)|
|STRINGdb      | [https://www.bioconductor.org/packages/release/bioc/html/STRINGdb.html](https://www.bioconductor.org/packages/release/bioc/html/STRINGdb.html)|
|AnnotationDbi | [https://www.bioconductor.org/packages/release/bioc/html/AnnotationDbi.html](https://www.bioconductor.org/packages/release/bioc/html/AnnotationDbi.html)|


> Exact versions are available along with the processing code for each specific dataset in the [GLDS_Processing_Scripts](https://developer.nasa.gov/asaravia/GeneLab_Data_Processing/tree/master/Microarray/2-channel_arrays/GLDS_Processing_Scripts) sub-directory. 

---

# General processing overview with example code  

The [glds_two_channel_arrays.Rmd](https://developer.nasa.gov/asaravia/GeneLab_Data_Processing/blob/master/Microarray/2-channel_arrays/glds_two_channel_arrays.Rmd) file is an R Markdown Interactive Document for processing GeneLab curated datasets that contain two-color gene expression assays. It can be run within the RStudio IDE with the Knit HTML button or on an R console with the command:  
rmarkdown::render("glds_two_channel_arrays.Rmd")
- The [map_annotation.csv](https://developer.nasa.gov/asaravia/GeneLab_Data_Processing/blob/master/Microarray/2-channel_arrays/map_annotation.csv), [organisms.csv](https://developer.nasa.gov/asaravia/GeneLab_Data_Processing/blob/master/Microarray/2-channel_arrays/organisms.csv) and [microarray_functions.R](https://developer.nasa.gov/asaravia/GeneLab_Data_Processing/blob/master/Microarray/2-channel_arrays/microarray_functions.R) files must be in the working directory to run the R Markdown Interactive Document. 
- Note that the [.Rprofile](https://developer.nasa.gov/asaravia/GeneLab_Data_Processing/blob/master/Microarray/2-channel_arrays/.Rprofile) file must also be present in the working directory in order for the interactive session to accept input files larger than 5 MB.

An html report file and dynamic site will be generated from parameters defined in the manual_entry code block. Interactive parameter selection is available using the Knit with Parameters option. Processed data files will also be exported to the working directory.  

Details about design model and normalization function parameters can be found in the [Limma package documentation](https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf).  

***
## Input Data and Parameters

* **GLDS Accession #** --- GeneLab repository accession number of the dataset or any string to create a seperate output folder (e.g. GLDS-8, GLDS8, Test1)
* **Organism** --- Supported species annotations are listed in the YAML header
* **Microarray Data Format** --- The Platform or Scanner file formatting of the raw data. This info can be found in the [GeneLab Data Repository](https://genelab-data.ndc.nasa.gov/genelab/projects) under 'PROTOCOLS', in the 'DESCRIPTION' of the 'image_aquisition' 'TYPE', and the imaging software used can be found in the 'DESCRIPTION' of the 'image_aquisition? array scanning protocol?' 'TYPE'. Gzip compressed files will automatically be extracted.
  + Agilent GenePix --- Agilent microarray platform imaged with GenePix software
  + Agilent --- Agilent microarray platform imaged with Agilent software
  + GenePix --- GenePix software formatting and intensity estimation
  + Spot --- SPOT software formatting with mean foreground estimation and morph background estimation
  + Imagene --- Imagene software formatting and intensity estimation
  + Bluefuse --- Bluefuse software formatting and intensity estimation
  + Generic --- Generic spot array
* **Background Correction** --- Options supported by limma. This info can be found in the [GeneLab Data Repository](https://genelab-data.ndc.nasa.gov/genelab/projects) under 'PROTOCOLS', in the 'DESCRIPTION' of the 'feature_extraction - not always there' 'TYPE'
  + normexp --- Adaptive background correction recommended as the default method.
  + subtract --- Subtract the background intensity from the foreground intensity for each spot.
  + morph --- Stabilize the variability of the M-values as a function of intensity.
* **Within Array Normalization** --- Options supported by limma. This info can be found in the [GeneLab Data Repository](https://genelab-data.ndc.nasa.gov/genelab/projects) under 'PROTOCOLS', in the 'DESCRIPTION' of the 'feature_extraction - is this true for every dataset?' 'TYPE'
  + loess --- Loess within array normalization recommended for Agilent arrays
  + printtiploess --- Print-tip recommended in general for others
  + robustspline ---  An empirical Bayes compromise between print-tip and global loess normalization
* **Between Array Normalization** --- Options supported by limma. This info can be found in the [GeneLab Data Repository](https://genelab-data.ndc.nasa.gov/genelab/projects) under 'PROTOCOLS', in the 'DESCRIPTION' of the 'feature_extraction - is this true for every dataset?' 'TYPE'
  + Aquantile --- Recommended as default. Applies quantile normalization to A-values between arrays
  + quantile --- normalizes directly to the individual red and green intensities
* **Primary Annotation Keytype** --- Gene identifier type for probe annotation from [Ann.dbi](https://www.bioconductor.org/packages/release/data/annotation/) sources. Mapped keytypes can be inferred from column headers in the probe annotation file or raw data files.
* **Raw Data Files** --- Select the raw data files to be analyzed. Raw data files for each dataset can be found in the [GeneLab Data Repository](https://genelab-data.ndc.nasa.gov/genelab/projects) under 'STUDY FILES' in the /GeneLab Processed Microarray Data Files/Raw Data/ directory.
* **Probe Annotation File** --- Select the Annotation file provided by the dataset submitter. This file can be downloaded from the [GeneLab Data Repository](https://genelab-data.ndc.nasa.gov/genelab/projects) for each microarray dataset under 'STUDY FILES' -> 'Microarray Data Files' -> &ast;adf.txt (from [ArrayExpress](https://www.ebi.ac.uk/arrayexpress/)) or GPL*.soft (from [GEO](https://www.ncbi.nlm.nih.gov/geo/)), or a platform specific file supplied by the dataset submitter. If multiple annotation formats are available, GPL is generally preferred.
* **ISA Metadata File** --- GeneLab ISAtab study metadata file, which can be downloaded from the [GeneLab Data Repository](https://genelab-data.ndc.nasa.gov/genelab/projects) under 'STUDY FILES' -> 'Study Metadata Files' -> *ISA.zip


***
## Output Data 

* list and define all output files generated from the script (similar to how all input parameters are listed and defined above)

