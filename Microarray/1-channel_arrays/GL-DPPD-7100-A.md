# GeneLab bioinformatics processing pipeline for one-channel microarray data

> **This page holds an overview and example code of how GeneLab processes one-channel Microarray datasets. Exact processing code and GL-DPPD-7100 revision used for specific datasets are available in the [GLDS_Processing_Scripts](https://developer.nasa.gov/asaravia/GeneLab_Data_Processing/tree/master/Microarray/1-channel_arrays/GLDS_Processing_Scripts) sub-directory and are also provided with their processed data in the [GeneLab Data Systems (GLDS) repository](https://genelab-data.ndc.nasa.gov/genelab/projects).**  

---

**Date:** Month day, 2020  
**Revision:** A   
**Document Number:** GL-DPPD-7100-A  

**Submitted by:**  
Homer Fogle (GeneLab Data Processing Team)

**Approved by:**  
Sylvain Costes (GeneLab Project Manager)  
Samrawit Gebre (GeneLab Deputy Project Manager and Interim GeneLab Configuration Manager)  
Amanda Saravia-Butler (GeneLab Data Processing Representative)  
Jonathan Galazka (GeneLab Project Scientist)

---

## Updates from previous revision

This version adds processing capability for additional Affymetrix ST type platforms as well as Nimblegen, Illumina Beadchip, and Genepix formatted Agilent expression arrays. Automated parsing of ISAtab metadata identifies paired sample data, technical replicates, timecourse factors, split-channel datafiles, and sample factor groups. Differential gene expression analysis for all pairwise group comparisons and gene level probe annotation have been added to the pipeline. 

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
|oligo         | [https://www.bioconductor.org/packages/release/bioc/html/oligo.html](https://www.bioconductor.org/packages/release/bioc/html/oligo.html)|
|affy          | [http://bioconductor.org/packages/release/bioc/html/affy.html](http://bioconductor.org/packages/release/bioc/html/affy.html)|
|illuminaio    | [http://bioconductor.org/packages/release/bioc/html/illuminaio.html](http://bioconductor.org/packages/release/bioc/html/illuminaio.html)|
|HELP          | [https://www.bioconductor.org/packages/release/bioc/html/HELP.html](https://www.bioconductor.org/packages/release/bioc/html/HELP.html)|


> Exact versions are available along with the processing code for each specific dataset in the [GLDS_Processing_Scripts](https://developer.nasa.gov/asaravia/GeneLab_Data_Processing/tree/master/Microarray/1-channel_arrays/GLDS_Processing_Scripts) sub-directory. 

---

# General processing overview with example code  

The [glds_one_channel_arrays.Rmd](https://developer.nasa.gov/asaravia/GeneLab_Data_Processing/blob/master/Microarray/1-channel_arrays/glds_one_channel_arrays.Rmd) file is an R Markdown Interactive Document for processing GeneLab curated datasets that contain one-color gene expression assays. It can be run within the RStudio IDE with the Knit HTML button or on an R console with the command:  
rmarkdown::render("glds_one_channel_arrays.Rmd")
- The [map_annotation.csv](https://developer.nasa.gov/asaravia/GeneLab_Data_Processing/blob/master/Microarray/1-channel_arrays/map_annotation.csv),[organisms.csv](https://developer.nasa.gov/asaravia/GeneLab_Data_Processing/blob/master/Microarray/1-channel_arrays/organisms.csv) and [microarray_functions.R](https://developer.nasa.gov/asaravia/GeneLab_Data_Processing/blob/master/Microarray/1-channel_arrays/microarray_functions.R) files must be in the working directory to run the R Markdown Interactive Document. 
- Note that the .Rprofile file must also be present in the working directory in order for the interactive session to accept input files larger than 5 MB.

An html report file and dynamic site will be generated from parameters defined in the manual_entry code block. Interactive parameter selection is available using the Knit with Parameters option. Processed data files will also be exported to the working directory.  

Details about design model and normalization function parameters can be found in the [Limma package documentation](https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf).  

***
## Input Data and Parameters

* **GLDS Accession #** --- GeneLab repository accession number of the dataset or any string to create a seperate output folder (e.g. GLDS-8, GLDS8, Test1)
* **Organism** --- Supported species annotations are listed in the YAML header
* **Microarray Data Format** --- The Platform or Scanner file formatting of the raw data. This info can be found in the [GeneLab Data Repository](https://genelab-data.ndc.nasa.gov/genelab/projects) under 'PROTOCOLS', in the 'DESCRIPTION' of the 'nucleic acid hybridization' 'TYPE'. Gzip compressed files will automatically be extracted.
  + Affymetrix Expression --- 3'-IVT type arrays; have probesets targeting the 3' end of transcripts
  + Affymetrix ST --- Gene or Exon type arrays; have probesets with greater coverage of transcripts
  + NimbleGen --- NimbleGen platform data in XYS or PAIR format
  + Illumina BeadChip --- Illumina platform expression data in IDAT or TXT formats
  + Agilent GenePix --- Agilent microarray platform in GPR format
  + Agilent --- Agilent microarray platform in TXT format
* **Study Design** --- Model to use when analyzing differential gene expression
  + Group Contrast --- Groups are defined by concatenation of study factor variables for each sample. Differential gene expression analysis is performed between all group pairs. This is suggested as default.
  + Repeated Measures --- Group contrast analysis with a paired sample parameter
* **Primary Annotation Keytype** --- Gene identifier type for probe annotation from [Ann.dbi](https://www.bioconductor.org/packages/release/data/annotation/) sources. Mapped keytypes can be inferred from column headers in the probe annotation file or raw data files. 
* **Raw Data Files** --- Select the raw data files to be analyzed. Raw data files for each dataset can be found in the [GeneLab Data Repository](https://genelab-data.ndc.nasa.gov/genelab/projects) under 'STUDY FILES' in the /GeneLab Processed Microarray Data Files/Raw Data/ directory
* **Probe Annotation File** --- Select the Annotation file provided by the dataset submitter. This file can be downloaded from the [GeneLab Data Repository](https://genelab-data.ndc.nasa.gov/genelab/projects) for each microarray dataset under 'STUDY FILES' -> 'Microarray Data Files' -> &ast;adf.txt (from [ArrayExpress](https://www.ebi.ac.uk/arrayexpress/)) or GPL*.soft (from [GEO](https://www.ncbi.nlm.nih.gov/geo/)), or a platform specific file supplied by the dataset submitter.
* **ISA Metadata File** --- GeneLab ISAtab study metadata file, which can be downloaded from the [GeneLab Data Repository](https://genelab-data.ndc.nasa.gov/genelab/projects) under 'STUDY FILES' -> 'Study Metadata Files' -> *ISA.zip

***
## Output Data 

* list and define all output files generated from the script (similar to how all input parameters are listed and defined above)

