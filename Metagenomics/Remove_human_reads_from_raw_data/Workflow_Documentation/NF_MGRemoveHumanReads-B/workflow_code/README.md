# NF_MGRemoveHumanReads-B Workflow Information and Usage Instructions


## General workflow info
The current pipeline for how GeneLab identifies and removes human DNA in Illumina metagenomics sequencing data (MGRemoveHumanReads), [GL-DPPD-7105-A.md](../../Pipeline_GL-DPPD-7105_Versions/GL-DPPD-7105-A.md), is implemented as a [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow and utilizes [conda](https://docs.conda.io/en/latest/) environments to install/run all tools. This workflow (NF_MGRemoveHumanReads-A) is run using the command line interface (CLI) of any unix-based system. The workflow can be used even if you are unfamiliar with Snakemake and conda, but if you want to learn more about those, [this Snakemake tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html) within [Snakemake's documentation](https://snakemake.readthedocs.io/en/stable/) is a good place to start for that, and an introduction to conda with installation help and links to other resources can be found [here at Happy Belly Bioinformatics](https://astrobiomike.github.io/unix/conda-intro). 

## Utilizing the workflow

1. [Install conda and NextFlow](#1-install-conda-NextFlow)  
2. [Download the workflow template files](#2-download-the-workflow-template-files)  
3. [Modify the variables in the Remove_Human_Reads.config file](#3-modify-the-variables-in-the-config-file)  
4. [Run the workflow](#4-run-the-workflow)  

### 1. Install conda, mamba, and `genelab-utils` package
We recommend installing a Miniconda, Python3 version appropriate for your system, as exemplified in [the above link](https://astrobiomike.github.io/unix/conda-intro#getting-and-installing-conda).  

Once Conda is installed, you can install NextFlow into your specified directory using the following code:

```bash
curl -s https://get.nextflow.io | bash

sudo mv nextflow /usr/local/bin
```

### 2. Download the workflow template files
All workflow files for removing human reads from metagenomics data are in the [workflow_code](workflow_code) directory. To get a copy of the latest SW_MGRemoveHumanReads-A version on to your system, run the following command:

```bash
GL-get-workflow MG-remove-human-reads-B
```

This downloaded the workflow into a directory called `NF_MGRemoveHumanReads-*/`, with the workflow version number at the end.

> Note: If wanting an earlier version, the wanted version can be provided as an optional argument like so:
> ```bash
> GL-get-workflow MG-remove-human-reads --wanted-version 1.0.0
> ```

### 3. Modify the variables in the Remove_Human_Reads.config file
Once you've downloaded the workflow template, you can modify the variables in the [Remove_Human_Reads.config](workflow_code/Remove_Human_Reads.config) file as needed. For example, you will have to provide a text file containing a single-column list of unique sample identifiers (see an example of how to set this up below - if you are running the example dataset, this file is provided in the [workflow_code](workflow_code) directory [here](workflow_code/unique-sample-IDs.txt)). You will also need to indicate the path to your input data (raw reads) and the root directory for where the kraken2 reference database should be stored (it will be setup automatically). Additionally, if necessary, you'll need to modify each variable in the [config.yaml](workflow_code/config.yaml) file to be consistent with the study you want to process and the machine you're using. 

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

### 4. Run the workflow

While in the directory holding the NextFlow file, .config file, and other workflow files that you downloaded in [step 2](#2-download-the-workflow-template-files), here is one example command of how to run the workflow:

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
