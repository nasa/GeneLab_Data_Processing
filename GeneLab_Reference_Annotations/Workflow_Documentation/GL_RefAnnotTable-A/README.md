# GL_RefAnnotTable-A Workflow Information and Usage Instructions <!-- omit in toc -->

## Table of Contents <!-- omit in toc -->
- [General Workflow Info](#general-workflow-info)
- [Utilizing the Workflow](#utilizing-the-workflow)
  - [Approach 1: Using Apptainer](#approach-1-using-apptainer)
    - [1. Install Apptainer](#1-install-apptainer)
    - [2. Download the Workflow Files](#2-download-the-workflow-files)
    - [3. Fetch Apptainer Image](#3-fetch-apptainer-image)
    - [4. Run the Workflow](#4-run-the-workflow)
    - [5. Run the Annotations Database Creation Function as a Stand-Alone Script](#5-run-the-annotations-database-creation-function-as-a-stand-alone-script)
  - [Approach 2: Using a Local R Environment](#approach-2-using-a-local-r-environment)
    - [1. Install R and Required R Packages](#1-install-r-and-required-r-packages)
    - [2. Download the Workflow Files](#2-download-the-workflow-files-1)
    - [3. Set Execution Permissions for Workflow Scripts](#3-set-execution-permissions-for-workflow-scripts)
    - [4. Run the Workflow](#4-run-the-workflow-1)
    - [5. Run the Annotations Database Creation Function as a Stand-Alone Script](#5-run-the-annotations-database-creation-function-as-a-stand-alone-script-1)

<br>

---

## General Workflow Info

The current GeneLab Reference Annotation Table (GL_RefAnnotTable-A) pipeline is implemented as an R workflow that can be run from a command line interface (CLI) using bash. The workflow can be executed using either a Apptainer (formerly Singularity) container or a local R environment. The workflow can be used even if you are unfamiliar with R, but if you want to learn more about R, visit the [R-project about page here](https://www.r-project.org/about.html). Additionally, an introduction to R along with installation help and information about using R for bioinformatics can be found [here at Happy Belly Bioinformatics](https://astrobiomike.github.io/R/basics).

<br>

---

## Utilizing the Workflow

The GL_RefAnnotTable-A workflow can be run using two approaches:

1. **[Using Apptainer](#approach-1-using-apptainer)**.  

2. **[Using a local R environment](#approach-2-using-a-local-r-environment)**.  

Please follow the instructions for the approach that best matches your setup and preferences. Each method is explained in the sections below.

<br>

---

### Approach 1: Using Apptainer

This approach allows you to run the workflow within a containerized environment, ensuring consistency and reproducibility.

<br>

---

#### 1. Install Apptainer

Apptainer can be installed either through [Anaconda](https://anaconda.org/conda-forge/singularity) or as documented on the [Apptainer documentation page](https://apptainer.org/docs/admin/main/installation.html).

> **Note**: If you prefer to use Anaconda, we recommend installing Miniconda for your system, as instructed by [Happy Belly Bioinformatics](https://astrobiomike.github.io/unix/conda-intro#getting-and-installing-conda).
>
> Once conda is installed on your system, you can install Apptainer by running:
>
> ```bash
> conda install -c conda-forge apptainer
> ```

<br>

---

#### 2. Download the Workflow Files

Download the latest version of the GL_RefAnnotTable-A workflow:

```bash
curl -LO https://github.com/nasa/GeneLab_Data_Processing/releases/download/GL_RefAnnotTable-A_1.1.0/GL_RefAnnotTable-A_1.1.0.zip
unzip GL_RefAnnotTable-A_1.1.0.zip
cd GL_RefAnnotTable-A_1.1.0
```

<br>

---

#### 3. Fetch Apptainer Image

To fetch the Apptainer images needed for the workflow, run:

```bash
bash bin/prepull_apptainer.sh config/software/by_docker_image.config
```
> Note: This command should be run in the directory containing the GL_RefAnnotTable-A_1.1.0 folder downloaded in [step 2](#2-download-the-workflow-files). Depending on your network speed, this may take approximately 20 minutes.

Once complete, an apptainer folder containing the Apptainer images will be created. Export this folder as an Apptainer configuration environment variable:

```bash
export APPTAINER_CACHEDIR=$(pwd)/apptainer
```

<br>

---

#### 4. Run the Workflow

While in the `GL_RefAnnotTable-A_1.1.0` directory, you can now run the workflow. Below is an example for generating an annotation table for Mus musculus (mouse):

```bash
apptainer exec -B $(pwd):/work \
$APPTAINER_CACHEDIR/quay.io-nasa_genelab-gl-refannottable-a-1.1.0.img \
bash -c "cd /work && Rscript GL-DPPD-7110-A_build-genome-annots-tab.R 'Mus musculus'"
```

**Input data:**

- No input files are required.
- Specify the target organism using a positional command line argument. `Mus musculus` is used in the example above.
- To see a list of all available organisms, run `Rscript GL-DPPD-7110-A_build-genome-annots-tab.R` without positional arguments. The correct argument for each organism can also be found in the 'species' column of the [GL-DPPD-7110-A_annotations.csv](../../Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv)
- Optional: a reference table CSV can be supplied as a second positional argument instead of using the default [GL-DPPD-7110-A_annotations.csv](../../Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv)

**Output data:**

- *-GL-annotations.tsv (Tab delineated table of gene annotations)
- *-GL-build-info.txt (Text file containing information used to create the annotation table, including tool and tool versions and date of creation)

<br>

---

#### 5. Run the Annotations Database Creation Function as a Stand-Alone Script

If the reference table does not specify an annotations database for the target organism in the annotations column, the `install_annotations` function (defined in `install-org-db.R`) will be executed. This function can also be run as a stand-alone script:

```bash
apptainer exec -B $(pwd):/work \
  $APPTAINER_CACHEDIR/quay.io-nasa_genelab-gl-refannottable-a-1.1.0.img \
  bash -c "cd /work && Rscript install-org-db.R 'Bacillus subtilis'"
```

**Input data:**

- The target organism must be specified as the first positional command line argument. `Bacillus subtilis` is used in the example above. The correct argument for each organism can be found in the 'species' column of [GL-DPPD-7110-A_annotations.csv](https://raw.githubusercontent.com/nasa/GeneLab_Data_Processing/master/GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv)
- Optional: A local reference table can be supplied as a second positional argument. If not provided, the script will download the current version of GL-DPPD-7110-A_annotations.csv from Github by default.

**Output data:**

- org.*.eg.db/ (Species-specific annotation database, as a local R package)

<br>

---

### Approach 2: Using a Local R Environment  

This approach allows you to run the workflow directly in your local R environment without using Apptainer containers.

<br>

---

#### 1. Install R and Required R Packages

We recommend installing R via the [Comprehensive R Archive Network (CRAN)](https://cran.r-project.org/):

1. Select the [CRAN Mirror](https://cran.r-project.org/mirrors.html) closest to your location.  
2. Navigate to the download page for your operating system.  
3. Download and install R (e.g., R-4.4.0).  

Once R is installed, open a terminal and start R:

```bash
R
```

Within an active R environment, run the following commands to install the required R packages:

```R
install.packages("tidyverse")

install.packages("BiocManager")

BiocManager::install("STRINGdb")
BiocManager::install("PANTHER.db")
BiocManager::install("rtracklayer")
BiocManager::install("AnnotationForge")
BiocManager::install("biomaRt")
BiocManager::install("GO.db")
```

<br>

---

#### 2. Download the Workflow Files

All files required for utilizing the GL_RefAnnotTable-A workflow for generating reference annotation tables are in the [workflow_code](workflow_code) directory. To get a copy of latest GL_RefAnnotTable version on to your system, run the following command:

```bash
curl -LO https://github.com/nasa/GeneLab_Data_Processing/releases/download/GL_RefAnnotTable-A_1.1.0/GL_RefAnnotTable-A_1.1.0.zip
``` 

<br>

---

#### 3. Set Execution Permissions for Workflow Scripts

Once you've downloaded the GL_RefAnnotTable-A workflow directory as a zip file, unzip the workflow then `cd` into the GL_RefAnnotTable-A_1.1.0 directory on the CLI. Next, run the following command to set the execution permissions for the R script:

```bash
unzip GL_RefAnnotTable-A_1.1.0.zip
cd GL_RefAnnotTable-A_1.1.0
chmod -R u+x *R
```

<br>

---

#### 4. Run the Workflow 

While in the GL_RefAnnotTable workflow directory, you are now able to run the workflow. Below is an example of how to run the workflow to build an annotation table for Mus musculus (mouse):

```bash
Rscript GL-DPPD-7110-A_build-genome-annots-tab.R 'Mus musculus'
```

**Input data:**

- No input files are required. Specify the target organism using a positional command line argument. `Mus musculus` is used in the example above. To see a list of all available organisms, run `Rscript GL-DPPD-7110-A_build-genome-annots-tab.R` without positional arguments. The correct argument for each organism can also be found in the 'species' column of [GL-DPPD-7110-A_annotations.csv](../../Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv)

- Optional: a reference table CSV can be supplied as a second positional argument instead of using the default [GL-DPPD-7110-A_annotations.csv](../../Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv)

**Output data:**

- *-GL-annotations.tsv (Tab delineated table of gene annotations)
- *-GL-build-info.txt (Text file containing information used to create the annotation table, including tool and tool versions and date of creation)

<br>

---

#### 5. Run the Annotations Database Creation Function as a Stand-Alone Script

If the reference table does not specify an annotations database for the target organism in the 'annotations' column, the `install_annotations` function (defined in `install-org-db.R`) will be executed. This function can also be run as a stand-alone script:

```bash
Rscript install-org-db.R 'Bacillus subtilis'
```

**Input data:**

- The target organism must be specified as the first positional command line argument. `Bacillus subtilis` is used in the example above. The correct argument for each organism can be found in the 'species' column of [GL-DPPD-7110-A_annotations.csv](https://raw.githubusercontent.com/nasa/GeneLab_Data_Processing/master/GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv)
- Optional: A local reference table can be supplied as a second positional argument. If not provided, the script will download the current version of GL-DPPD-7110-A_annotations.csv from Github by default.

**Output data:**

- org.*.eg.db/ (species-specific annotation database, as a local R package)

<br>
