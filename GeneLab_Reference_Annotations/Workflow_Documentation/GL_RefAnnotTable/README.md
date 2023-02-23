# GL_RefAnnotTable Workflow Information and Usage Instructions

## General workflow info
The current GeneLab Reference Annotation Table (GL_RefAnnotTable) pipeline is implemented as an R workflow that can be run from a command line interface (CLI) using bash. The workflow can be used even if you are unfamiliar with R, but if you want to learn more about R, visit the [R-project about page here](https://www.r-project.org/about.html). Additionally, an introduction to R along with installation help and information about using R for bioinformatics can be found [here at Happy Belly Bioinformatics](https://astrobiomike.github.io/R/basics).  

## Utilizing the workflow

1. [Install conda, R, and R packages](#1-install-conda-r-and-r-packages)
2. [Download the workflow files](#2-download-the-workflow-files)  
3. [Run the workflow](#3-run-the-workflow)  

<br>

### 1. Install conda, R, and R packages

This details installing R in a [conda](https://docs.conda.io/en/latest/) environment. If you don't have conda yet, we recommend installing a Miniconda, Python3 version appropriate for your system – as exemplified [here](https://astrobiomike.github.io/unix/conda-intro#getting-and-installing-conda) along with an introduction to conda overall.

Once conda is installed on your system, we recommend installing [mamba](https://github.com/mamba-org/mamba#mamba), as it generally allows for much faster conda installations:

```bash
conda install -n base -c conda-forge mamba
```

Then the starting conda environment can be created as follows: 

```bash
mamba create -n GL-gen-ref-annots -c conda-forge -c bioconda -c defaults \
             r-base=4.2.0 r-biocmanager=1.30.18 bioconductor-panther.db=1.0.11 \
             bioconductor-delayedarray=0.24.0 r-tidyverse=1.3.2 r-remotes=2.4.2 \
             r-systemfonts=1.0.4 -y
```

When that finishes, activate the environment (using `conda`):

```bash
conda activate GL-gen-ref-annots
```

Once in the environment, enter R by running the following command: 

```bash
R
```

Within the active R environment, run the following commands to install the required R packages (this may take a few minutes):

```R
BiocManager::install(version = 3.15, ask = FALSE, update = FALSE, force = TRUE)

BiocManager::install(c("STRINGdb", "rtracklayer"), version = 3.15, update = FALSE)

# when finished, we can exit R
quit("no")
```

<br>

### 2. Download the Workflow Files

All files required for utilizing the GL_RefAnnotTable workflow for generating reference annotation tables are in the [workflow_code](workflow_code) directory. To get a copy of latest GL_RefAnnotTable version on to your system, run the following command:

```bash
# this link won't work until i update the release to 1.0.1
# curl -L https://github.com/nasa/GeneLab_Data_Processing/releases/download/GL_RefAnnotTable_1.0.1/GL-DPPD-7110_build-genome-annots-tab-1.0.1.R.txt > GL-DPPD-7110_build-genome-annots-tab.R

## until i pack up 1.0.1, this will work
curl -LO https://raw.githubusercontent.com/nasa/GeneLab_Data_Processing/gl-reference-annotations/GeneLab_Reference_Annotations/Workflow_Documentation/GL_RefAnnotTable/workflow_code/GL-DPPD-7110_build-genome-annots-tab.R
``` 

<br>

### 3. Run the Workflow

While in the GL_RefAnnotTable workflow directory, you are now able to run the workflow. Below is an example of how to run the workflow to build an annotation table for *C. elegans*:

```bash
Rscript GL-DPPD-7110_build-genome-annots-tab.R WORM
```

**Input data:**

- No input files required, but a target organism must be specified as a positional command line argument, `MOUSE` is used in the example above. Run `Rscript GL-DPPD-7110_build-genome-annots-tab.R` with no positional arugments to see the list of currently available organisms. 

**Output data:**

- *-GL-annotations.tsv (Tab delineated table of gene annotations)
- *-GL-build-info.txt (Text file containing information used to create the annotation table, including tool and tool versions and date of creation)

