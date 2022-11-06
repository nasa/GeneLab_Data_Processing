# GL_RefAnnotTable Workflow Information and Usage Instructions

## General workflow info
The current GeneLab Reference Annotation Table (GL_RefAnnotTable) pipeline is implemented as an R workflow that can be run from a command line interface (CLI) using bash. The workflow can be used even if you are unfamiliar with R, but if you want to learn more about R, visit the [R-project about page here](https://www.r-project.org/about.html). Additionally, an introduction to R along with installation help and information about using R for bioinformatics can be found [here at Happy Belly Bioinformatics](https://astrobiomike.github.io/R/basics).  

## Utilizing the workflow

1. [Install R and R packages](#1-install-r-and-r-packages)  
2. [Download the workflow files](#2-download-the-workflow-files)  
3. [Setup Execution Permission for Workflow Scripts](#3-setup-execution-permission-for-workflow-scripts)
4. [Run the workflow](#4-run-the-workflow)  

<br>

### 1. Install R and R packages

We recommend installing R via the [Comprehensive R Archive Network (CRAN)](https://cran.r-project.org/) as follows: 

1. Select the [CRAN Mirror](https://cran.r-project.org/mirrors.html) closest to your location.
2. Click the link under the "Download and Install R" section that's consistent with your machine.
3. Clink on the R-4.2.1 package consistent with your machine to download.
4. Double click on the R-4.2.1.pkg downloaded in step 3 and follow the installation instructions.

Once R is installed, open a CLI terminal and run the following command to activate R:

```bash
R
```

Within an active R environment, run the following commands to install the required R packages:

```R
install.packages("tidyverse", version = 1.3.2, repos = "http://cran.us.r-project.org")

install.packages("BiocManager", version = 3.15, repos = "http://cran.us.r-project.org")

BiocManager::install("STRINGdb", version = 3.15)
BiocManager::install("PANTHER.db", version = 3.15)
BiocManager::install("rtracklayer", version = 3.15)
```

<br>

### 2. Download the Workflow Files

All files required for utilizing the GL_RefAnnotTable workflow for generating reference annotation tables are in the [workflow_code](workflow_code) directory. To get a copy of latest GL_RefAnnotTable version on to your system, run the following command:

```bash
curl -LO https://github.com/nasa/GeneLab_Data_Processing/releases/download/GL_RefAnnotTable_1.0.0/GL_RefAnnotTable_1.0.0.zip
``` 

<br>

### 3. Setup Execution Permission for Workflow Scripts

Once you've downloaded the GL_RefAnnotTable workflow directory as a zip file, unzip the workflow then `cd` into the GL_RefAnnotTable_1.0.0 directory on the CLI. Next, run the following command to set the execution permissions for the R script:

```bash
chmod -R u+x *R
```

<br>

### 4. Run the Workflow

While in the GL_RefAnnotTable workflow directory, you are now able to run the workflow. Below is an example of how to run the workflow to build an annotation table for Mus musculus (mouse):

```bash
Rscript GL-DPPD-7110_build-genome-annots-tab.R MOUSE
```

**Input data:**

- No input files required, but a target organism must be specified as a positional command line argument, `MOUSE` is used in the example above. Run `Rscript GL-DPPD-7110_build-genome-annots-tab.R` with no positional arugments to see the list of currently available organisms. 

**Output data:**

- *-GL-annotations.tsv (Tab delineated table of gene annotations)
- *-GL-build-info.txt (Text file containing information used to create the annotation table, including tool and tool versions and date of creation)

