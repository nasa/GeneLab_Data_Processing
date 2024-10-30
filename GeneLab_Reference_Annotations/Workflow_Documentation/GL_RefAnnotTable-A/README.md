# GL_RefAnnotTable-A Workflow Information and Usage Instructions <!-- omit in toc -->

## Table of Contents <!-- omit in toc -->

- [General Workflow Information](#general-workflow-information)
- [Utilizing the Workflow](#utilizing-the-workflow)
  - [1. Download the Workflow Files](#1-download-the-workflow-files)
  - [2. Run the Workflow](#2-run-the-workflow)
    - [Approach 1: Using Singularity](#approach-1-using-singularity)
      - [Step 1: Install Singularity](#step-1-install-singularity)
      - [Step 2: Fetch the Singularity Image](#step-2-fetch-the-singularity-image)
      - [Step 3: Run the Workflow](#step-3-run-the-workflow)
    - [Approach 2: Using a Local R Environment](#approach-2-using-a-local-r-environment)
      - [Step 1: Install R and Required R Packages](#step-1-install-r-and-required-r-packages)
      - [Step 2: Run the Workflow](#step-2-run-the-workflow)
    - [Workflow Input/Output Data](#workflow-input-output-data)
  - [3. Run the Annotations Database Creation Function as a Stand-Alone Script](#3-run-the-annotations-database-creation-function-as-a-stand-alone-script)
    - [Using Singularity](#using-singularity)
    - [Using a Local R Environment](#using-a-local-r-environment)

<br> 

---

## General Workflow Information

The current GeneLab Reference Annotation Table (GL_RefAnnotTable-A) pipeline is implemented as an R workflow that can be run from a command line interface (CLI) using bash. The workflow can be executed using either a Singularity container or a local R environment. The workflow can be used even if you are unfamiliar with R, but if you want to learn more about R, visit the [R-project about page here](https://www.r-project.org/about.html). Additionally, an introduction to R along with installation help and information about using R for bioinformatics can be found [here at Happy Belly Bioinformatics](https://astrobiomike.github.io/R/basics).

<br> 

---

## Utilizing the Workflow

### 1. Download the Workflow Files

Download the latest version of the GL_RefAnnotTable-A workflow:

```bash
curl -LO https://github.com/nasa/GeneLab_Data_Processing/releases/download/GL_RefAnnotTable-A_1.1.0/GL_RefAnnotTable-A_1.1.0.zip
unzip GL_RefAnnotTable-A_1.1.0.zip
```

<br> 

---

### 2. Run the Workflow

The GL_RefAnnotTable-A workflow can be run using one of two approaches:

- **[Approach 1: Using Singularity](#approach-1-using-singularity)**
- **[Approach 2: Using a Local R Environment](#approach-2-using-a-local-r-environment)**

Please follow the instructions for the approach that best matches your setup and preferences. Each method is explained in detail below. 

> **Note**: If you encounter timeout errors, you can increase the default timeout (3600 seconds) by modifying the `options(timeout=3600)` line at the top of the `GL-DPPD-7110-A_build-genome-annots-tab.R` script.

<br> 

---

### Approach 1: Using Singularity

This approach allows you to run the workflow within a containerized environment, ensuring consistency and reproducibility.

<br>

#### Step 1: Install Singularity

Singularity is a containerization platform for running applications portably and reproducibly. We use container images hosted on Quay.io to encapsulate all the necessary software and dependencies required by the GL_RefAnnotTable-A workflow. This setup allows you to run the workflow without installing any software directly on your system. 
  
> ***Note**: Other containerization tools like Docker or Apptainer can also be used to pull and run these images.*
   

We recommend installing Singularity system-wide as per the official [Singularity installation documentation](https://docs.sylabs.io/guides/3.10/admin-guide/admin_quickstart.html).  
 

> ***Note**: While Singularity is also available through [Anaconda](https://anaconda.org/conda-forge/singularity), we recommend installing Singularity system-wide following the official installation documentation.*

<br> 

#### Step 2: Fetch the Singularity Image

To pull the Singularity image needed for the workflow, you can use the provided script as directed below or pull the image directly.  

> ***Note**: This command should be run in the location containing the `GL_RefAnnotTable-A_1.1.0` directory that was downloaded in [step 1](#1-download-the-workflow-files). Depending on your network speed, fetching the images will take approximately 20 minutes.*  
 

```bash
bash GL_RefAnnotTable-A_1.1.0/bin/prepull_singularity.sh GL_RefAnnotTable-A_1.1.0/config/software/by_docker_image.config
```
 
Once complete, a `singularity` folder containing the Singularity images will be created. Run the following command to export this folder as an environment variable:  
 

```bash
export SINGULARITY_CACHEDIR=$(pwd)/singularity
```

<br> 

#### Step 3: Run the Workflow

While in the directory containing the `GL_RefAnnotTable-A_1.1.0` folder that was downloaded in [step 1](#1-download-the-workflow-files), you can now run the workflow. Below is an example for generating the annotation table for *Mus musculus* (mouse): 
 

```bash
singularity exec -B $(pwd)/GL_RefAnnotTable-A_1.1.0:/work \
$SINGULARITY_CACHEDIR/quay.io-nasa_genelab-gl-refannottable-a-1.1.0.img \
Rscript /work/GL-DPPD-7110-A_build-genome-annots-tab.R 'Mus musculus'
```

<br> 

---

### Approach 2: Using a Local R Environment

This approach allows you to run the workflow directly in your local R environment without using containers.

<br>

#### Step 1: Install R and Required R Packages

We recommend installing R via the [Comprehensive R Archive Network (CRAN)](https://cran.r-project.org/):

1. Select the [CRAN Mirror](https://cran.r-project.org/mirrors.html) closest to your location.
   
2. Navigate to the download page for your operating system.
   
3. Download and install R (e.g., R-4.4.0). 

Once R is installed, install the required R packages as follows:

Open a terminal and start R: 
 

```bash
R
``` 

 
Within the R environment, run the following commands to install the required packages: 
 

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

#### Step 2: Run the Workflow

While in the directory containing the `GL_RefAnnotTable-A_1.1.0` folder that was downloaded in [step 1](#1-download-the-workflow-files), you can now run the workflow. Below is an example of how to run the workflow to build an annotation table for *Mus musculus* (mouse):  
 

```bash
Rscript GL_RefAnnotTable-A_1.1.0/GL-DPPD-7110-A_build-genome-annots-tab.R 'Mus musculus'
```

 <br> 

 --- 

 ### Workflow Input/Output Data 

The input and output data are the same for both [Approach 1: Using Singularity](#approach-1-using-singularity) and [Approach 2: Using a Local R Environment](#approach-2-using-a-local-r-environment).

 <br>

**Input data:**

- No input files are required. Specify the species name of the target organism using a positional command line argument. `Mus musculus` is used in both the Singularity and the local R environment examples above.
  > **Notes**:  
  > - To see a list of all available organisms, run `Rscript GL-DPPD-7110-A_build-genome-annots-tab.R` without positional arguments.  
  > - The correct argument for each organism can also be found in the 'species' column of the [GL-DPPD-7110-A_annotations.csv](../../Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv)
  
- *Optional*: A local reference table CSV file can be supplied as a second positional argument. If not provided, the script will download the current version of the [GL-DPPD-7110-A_annotations.csv](../../Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv) table by default.
   

**Output data:**

- *-GL-annotations.tsv (Tab delineated table of gene annotations)
  
- *-GL-build-info.txt (Text file containing information used to create the annotation table, including tool and tool versions and date of creation)

<br>

---

### 3. Run the Annotations Database Creation Function as a Stand-Alone Script 

If the reference table does not specify an annotations database for the target organism in the 'annotations' column of the [GL-DPPD-7110-A_annotations.csv](../../Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv) file, the `install_annotations` function (defined in `install-org-db.R`) will be executed by default. This function can also be run as a stand-alone script:

> **Note**: If you encounter timeout errors, you can increase the default timeout (3600 seconds) by modifying the `options(timeout=3600)` line at the top of the `install-org-db.R` script.

<br>

#### Using Singularity 

```bash
singularity exec -B $(pwd)/GL_RefAnnotTable-A_1.1.0:/work \
$SINGULARITY_CACHEDIR/quay.io-nasa_genelab-gl-refannottable-a-1.1.0.img \
Rscript /work/install-org-db.R 'Bacillus subtilis'
```
 
<br>

#### Using a Local R Environment

```bash
Rscript GL_RefAnnotTable-A_1.1.0/install-org-db.R 'Bacillus subtilis'
```

<br>

**Input data:**

- The species name of the target organism must be specified as the first positional command line argument. `Bacillus subtilis` is used in both the Singularity and local R examples above.
  > **Note**: The correct argument for each organism can also be found in the 'species' column of the [GL-DPPD-7110-A_annotations.csv](../../Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv)
  
- *Optional*: A local reference table CSV file can be supplied as a second positional argument. If not provided, the script will download the current version of the [GL-DPPD-7110-A_annotations.csv](../../Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv) table by default.  


**Output data:**

- org.*.eg.db/ (Species-specific annotation database, as a local R package)

<br>

---
