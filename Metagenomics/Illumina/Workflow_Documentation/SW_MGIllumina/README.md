# SW_MGIllumina Workflow Information and Usage Instructions


## General workflow info
The current GeneLab Illumina metagenomics sequencing data processing pipeline (MGIllumina), [GL-DPPD-7107.md](../../Pipeline_GL-DPPD-7107_Versions/GL-DPPD-7107.md), is implemented as a [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow and utilizes [conda](https://docs.conda.io/en/latest/) environments to install/run all tools. This workflow (SW_MGIllumina) is run using the command line interface (CLI) of any unix-based system. The workflow can be used even if you are unfamiliar with Snakemake and conda, but if you want to learn more about those, [this Snakemake tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html) within [Snakemake's documentation](https://snakemake.readthedocs.io/en/stable/) is a good place to start for that, and an introduction to conda with installation help and links to other resources can be found [here at Happy Belly Bioinformatics](https://astrobiomike.github.io/unix/conda-intro).  

> **Note on reference databases**  
> Many reference databases are relied upon throughout this workflow. They will be installed and setup automatically the first time the workflow is run. All together, after installed and unpacked, they will take up about 240 GB of storage, but they may also require up to 500GB during installation and initial un-packing, so be sure there is enough room on your system before running the workflow.

## Utilizing the workflow

1. [Install conda and `genelab-utils` package](#1-install-conda-and-genelab-utils-package)  
2. [Download the workflow template files](#2-download-the-workflow-template-files)  
3. [Modify the variables in the config.yaml file](#3-modify-the-variables-in-the-configyaml-file)  
4. [Run the workflow](#4-run-the-workflow)  

### 1. Install conda and `genelab-utils` package
We recommend installing a Miniconda, Python3 version appropriate for your system, as exemplified in [the above link](https://astrobiomike.github.io/unix/conda-intro#getting-and-installing-conda).  

Once conda is installed on your system, you can install the genelab-utils conda package in a new environment with the following command:

```bash
conda create -n genelab-utils -c conda-forge -c bioconda -c defaults -c astrobiomike 'genelab-utils>=1.1.02'
```

The environment then needs to be activated:

```bash
conda activate genelab-utils
```

### 2. Download the workflow template files
All files required for utilizing the GeneLab workflow for processing Illumina metagenomics sequencing data are in the [workflow_code](workflow_code) directory. To get a copy of the latest SW_MGIllumina version on to your system, run the following command:

```bash
GL-get-workflow MG-Illumina
```

This downloaded the workflow into a directory called `SW_MGIllumina_*/`, with the workflow version number at the end.

### 3. Modify the variables in the config.yaml file
Once you've downlonaded the workflow template, you can modify the variables in your downloaded version of the [config.yaml](workflow_code/SW_MGIllumina_1.0.0/config.yaml) file as needed in order to match your dataset and system setup. For example, you will have to provide a text file containing a single-column list of unique sample identifiers (see an example of how to set this up below). You will also need to indicate the paths to your input data (raw reads) and the root directory for where the reference databases should be stored (they will be setup automatically). Additionally, if necessary, you'll need to modify each variable in the config.yaml file to be consistent with the study you want to process and the machine you're using. 

> Note: If you are unfamiliar with how to specify paths, one place you can learn more is [here](https://astrobiomike.github.io/unix/getting-started#the-unix-file-system-structure).  

**Example for how to create a single-column list of unique sample identifiers from your raw data file names**

For example, if you have paired-end read data for 2 samples located in `../Raw_Data/` relative to your workflow directory, that would look like this:

```bash
ls ../Raw_Data/
```

```
Sample-1_R1_raw.fastq.gz
Sample-1_R2_raw.fastq.gz
Sample-2_R1_raw.fastq.gz
Sample-2_R2_raw.fastq.gz
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

While in the directory holding the Snakefile, config.yaml, and other workflow files that you downloaded in [step 2](#2-download-the-workflow-template-files), here is one example command of how to run the workflow:

```bash
snakemake --use-conda --conda-prefix ${CONDA_PREFIX}/envs -j 2 -p
```

* `--use-conda` – specifies to use the conda environments included in the workflow (these are specified in the files in the workflow `envs/` directory)
* `--conda-prefix` – indicates where the needed conda environments will be stored. Adding this option will also allow the same conda environments to be re-used when processing additional datasets, rather than making new environments each time you run the workflow. The value listed for this option, `${CONDA_PREFIX}/envs`, points to the default location for conda environments (note: the variable `${CONDA_PREFIX}` will be expanded to the appropriate location on whichever system it is run on).
* `-j` – assigns the number of jobs Snakemake should run concurrently (keep in mind that many of the thread and cpu parameters set in the config.yaml file will be multiplied by this)
* `-p` – specifies to print out each command being run to the screen

See `snakemake -h` and [Snakemake's documentation](https://snakemake.readthedocs.io/en/stable/) for more options and details.

---
