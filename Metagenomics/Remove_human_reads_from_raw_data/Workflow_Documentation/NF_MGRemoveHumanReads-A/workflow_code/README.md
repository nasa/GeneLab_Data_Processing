# NF_MGRemoveHumanReads-A Workflow Information and Usage Instructions


## General workflow info
The current pipeline for how GeneLab identifies and removes human DNA in Illumina metagenomics sequencing data (MGRemoveHumanReads), [GL-DPPD-7105-A.md](../../Pipeline_GL-DPPD-7105_Versions/GL-DPPD-7105-A.md), is implemented as a [NextFlow](https://www.nextflow.io/docs/stable/index.html) DSL2 workflow and utilizes [Singularity](https://docs.sylabs.io/guides/3.10/user-guide/introduction.html) run all tools in containers. This workflow (NF_MGRemoveHumanReads-A) is run using the command line interface (CLI) of any unix-based system. The workflow can be used even if you are unfamiliar with NextFlow and Singularity, but if you want to learn more about those, [this NextFlow tutorial](https://training.nextflow.io/basic_training/) within [NextFlow's documentation](https://www.nextflow.io/docs/stable/index.html) is a good place to start for that.

## Utilizing the workflow

1. [Install conda, mamba, and `genelab-utils` package](#1-install-conda-mamba-and-genelab-utils-package) 
2. [Install NextFlow and Singularity](#2-install-NextFlow-Singularity)
   2a. [Install Nextflow](#2a-install-nextflow)  
   2b. [Install Singularity](#2b-install-singularity)
3. [Download the workflow template files](#3-download-the-workflow-template-files)  
4. [Modify the variables in the Remove_Human_Reads.config file](#4-modify-the-variables-in-the-config-file)  
5. [Run the workflow](#5-run-the-workflow)  


<br>

---

### 1. Install conda, mamba, and `genelab-utils` package
We recommend installing a Miniconda, Python3 version appropriate for your system, as exemplified in [the Happy Belly Bioinformatics conda tutorial](https://astrobiomike.github.io/unix/conda-intro#getting-and-installing-conda).  

Once conda is installed on your system, we recommend installing [mamba](https://github.com/mamba-org/mamba#mamba), as it generally allows for much faster conda installations:

```bash
conda install -n base -c conda-forge mamba
```

> You can read a quick intro to mamba [here](https://astrobiomike.github.io/unix/conda-intro#bonus-mamba-no-5) if wanted.

Once mamba is installed, you can install the genelab-utils conda package in a new environment with the following command:

```bash
mamba create -n genelab-utils -c conda-forge -c bioconda -c defaults -c astrobiomike 'genelab-utils>=1.1.02'
```

The environment then needs to be activated:

```bash
conda activate genelab-utils
```

<br>

---

### 2. Install Nextflow and Singularity 

#### 2a. Install Nextflow

Nextflow can be installed either through [Anaconda](https://anaconda.org/bioconda/nextflow) or as documented on the [Nextflow documentation page](https://www.nextflow.io/docs/latest/getstarted.html).

> Note: If you want to install Anaconda, we recommend installing a Miniconda, Python3 version appropriate for your system, as instructed by [Happy Belly Bioinformatics](https://astrobiomike.github.io/unix/conda-intro#getting-and-installing-conda).  
> 
> Once conda is installed on your system, you can install the latest version of Nextflow by running the following commands:
> 
> ```bash
> conda install -c bioconda nextflow
> nextflow self-update
> ```

<br>

#### 2b. Install Singularity

Singularity is a container platform that allows usage of containerized software. This enables the GeneLab RCP workflow to retrieve and use all software required for processing without the need to install the software directly on the user's system.

We recommend installing Singularity on a system wide level as per the associated [documentation](https://docs.sylabs.io/guides/3.10/admin-guide/admin_quickstart.html).

> Note: Singularity is also available through [Anaconda](https://anaconda.org/conda-forge/singularity).

<br>

---


### 3. Download the workflow template files
All workflow files for removing human reads from metagenomics data are in the [workflow_code](workflow_code) directory. To get a copy of the latest NF_MGRemoveHumanReads-A version on to your system, run the following command:

```bash
GL-get-workflow NF_MGRemoveHumanReads-A
```

This downloaded the workflow into a directory called `NF_MGRemoveHumanReads-*/`, with the workflow version number at the end.

> Note: If wanting an earlier version, the wanted version can be provided as an optional argument like so:
> ```bash
> GL-get-workflow NF_MGRemoveHumanReads-A --wanted-version 1.0.0
> ```

<br>

---

### 4. Modify the variables in the Remove_Human_Reads.config file
Once you've downloaded the workflow template, you can modify the variables in the [Remove_Human_Reads.config](workflow_code/Remove_Human_Reads.config) file as needed. For example, you will have to provide a text file containing a single-column list of unique sample identifiers (see an example of how to set this up below - if you are running the example dataset, this file is provided in the [workflow_code](workflow_code) directory [here](workflow_code/unique-sample-IDs.txt)). You will also need to indicate the path to your input data (raw reads) and the root directory for where the kraken2 reference database should be stored (it will be setup automatically). Additionally, if necessary, you'll need to modify each variable in the [Remove_Human_Reads.config](workflow_code/Remove_Human_Reads.config) file to be consistent with the study you want to process and the machine you're using. 

> Note: If you are unfamiliar with how to specify paths, one place you can learn more is [here](https://astrobiomike.github.io/unix/getting-started#the-unix-file-system-structure).  

**Example for how to create a single-column list of unique sample identifiers from your raw data file names**

For example, if you only want to process a subset of the read files within the reads directory and have paired-end read data for 2 samples located in `../Raw_Sequence_Data/` relative to your workflow directory, that would look like this:

```bash
ls ../Raw_Sequence_Data/
```

```
Sample-1_R1.fastq.gz
Sample-1_R2.fastq.gz
Sample-2_R1.fastq.gz
Sample-2_R2.fastq.gz
```

You would set up your `unique-sample-IDs.txt` file as follows:

```bash
cat unique-sample-IDs.txt
```

```
Sample-1
Sample-2
```

<br>

---

### 5. Run the workflow

While in the directory holding the NextFlow file, .config file, and other workflow files that you downloaded in [step 3](#3-download-the-workflow-template-files), here is one example command of how to run the workflow:

```bash
nextflow run *path/to/Remove_Human_Reads.nf* -ansi-log false -specify_reads false
```

* `-ansi-log false` – specifies to print out each command being run to the screen
* `-resume` – continues to run the workflow using cached data from the previous run
* `-specify_reads false` - processes all reads in the working directory, without requiring a sample ID list


See `nextflow -h` and [NextFlow's documentation](https://www.nextflow.io/docs/master/index.html) for more options and details.

A quick example can be run with the files included in the [workflow_code](workflow_code) directory after specifying a location for the reference database in the [Remove_Human_Reads.config](workflow_code/Remove_Human_Reads.config) file.

---

## Reference database info
The database we use was built with kraken2 v2.1.1 as detailed below, and can be downloaded to run with this workflow (it's ~4.3 GB uncompressed).

---
