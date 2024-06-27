# Workflow Information and Usage Instructions

## General Workflow Info

### Implementation Tools

The current GeneLab Illumina amplicon sequencing data processing pipeline (AmpIllumina), [GL-DPPD-7104-B.md](../../Pipeline_GL-DPPD-7104_Versions/GL-DPPD-7104-B.md), is implemented as a [Nextflow](https://nextflow.io/) DSL2 workflow and utilizes [Singularity](https://docs.sylabs.io/guides/3.10/user-guide/introduction.html) containers or [conda](https://docs.conda.io/en/latest/) environments to install/run all tools. This workflow is run using the command line interface (CLI) of any unix-based system.  While knowledge of creating workflows in nextflow is not required to run the workflow as is, [the Nextflow documentation](https://nextflow.io/docs/latest/index.html) is a useful resource for users who want to modify and/or extend this workflow.   

## Utilizing the Workflow

1. [Install nextflow, conda and singularity](#1-install-nextflow-conda-and-singularity)  
   1a. [Install nextflow and conda](#1a-install-nextflow-and-conda)  
   1b. [Install singularity](#1b-install-singularity)  

2. [Download the workflow files](#2-download-the-workflow-files)  

3. [Run the workflow](#3-run-the-workflow)  
   3a. [Approach 1: Run slurm jobs in singularity containers with OSD accession as input](#3a-approach-1-run-slurm-jobs-in-singularity-containers-with-osd-accession-as-input)   
   3b. [Approach 2: Run slurm jobs in singularity containers with a csv file as input](#3b-approach-2-run-slurm-jobs-in-singularity-containers-with-a-csv-file-as-input)  
   3c. [Approach 3: Run jobs locally in conda environments and specify the path to one or more existing conda environments](#3c-approach-run-jobs-locally-in-conda-environments-and-specify-the-path-to-one-or-more-existing-conda-environments)  
   3d. [Modify parameters and cpu resources in the nextflow config file](#3d-modify-parameters-and-cpu-resources-in-the-nextflow-config-file)  

4. [Workflow outputs](#4-workflow-outputs)  
   4a. [Main outputs](#4a-main-outputs)  
   4b. [Resource logs](#4b-resource-logs)  

<br>

### 1. Install nextflow, conda and singularity



####  1a. Install nextflow and conda

Nextflow can be installed either through [Anaconda](https://anaconda.org/bioconda/nextflow) or as documented on the [Nextflow documentation page](https://www.nextflow.io/docs/latest/getstarted.html).

> Note: If you want to install anaconda, we recommend installing a miniconda, python3 version appropriate for your system, as instructed by [Happy Belly Bioinformatics](https://astrobiomike.github.io/unix/conda-intro#getting-and-installing-conda).  

We recommend installing a miniconda, python3 version appropriate for your system, as exemplified in [the above link](https://astrobiomike.github.io/unix/conda-intro#getting-and-installing-conda).

Once conda is installed on your system, we recommend installing [mamba](https://github.com/mamba-org/mamba#mamba), as it generally allows for much faster conda installations.

```bash
conda install -n base -c conda-forge mamba
```

> You can read a quick intro to mamba [here](https://astrobiomike.github.io/unix/conda-intro#bonus-mamba-no-5).

Once mamba is installed, you can install the genelab-utils conda package which contains nextflow with the following command:

```bash
mamba create -n genelab-utils -c conda-forge -c bioconda -c defaults -c astrobiomike genelab-utils
```

The environment then needs to be activated:

```bash
conda activate genelab-utils

# Test that nextflow is installed
nextflow -h

# Update nextflow
nextflow self-update
```

<br>

#### 1b. Install singularity

Singularity is a container platform that allows usage of containerized software. This enables the GeneLab workflow to retrieve and use all software required for processing without the need to install the software directly on the user's system.

We recommend installing singularity on a system wide level as per the associated [documentation](https://docs.sylabs.io/guides/3.10/admin-guide/admin_quickstart.html).

<br>

### 2. Download the workflow files

All files required for utilizing the NF_XXX GeneLab workflow for processing amplicon illumina data are in the [workflow_code](workflow_code) directory. To get a copy of latest *NF_XXX* version on to your system, the code can be downloaded as a zip file from the release page then unzipped after downloading by running the following commands: 

```bash
wget https://github.com/nasa/GeneLab_Data_Processing/releases/download/NF_AmpIllumina/NF_AmpIllumina.zip
unzip NF_AmpIllumina.zip &&  cd NF_XXX-X_X.X.X
```

OR by using the genelab-utils conda package

```bash
GL-get-workflow  Amplicon-Illumina
```

<br>

### 3. Run the Workflow

For options and detailed help on how to run the workflow, run the following command:

```bash
nextflow run main.nf --help
```

> Note: Nextflow commands use both single hyphen arguments (e.g. -help) that denote general nextflow arguments and double hyphen arguments (e.g. --csv_file) that denote workflow specific parameters.  Take care to use the proper number of hyphens for each argument.

<br>

#### 3a. Approach 1: Run slurm jobs in singularity containers with OSD accession as input

```bash
nextflow run main.nf -resume -profile slurm,singularity --GLDS_accession GLDS-487 --target_region 16S
```

<br>

#### 3b. Approach 2: Run slurm jobs in singularity containers with a csv file as input

```bash
nextflow run main.nf -resume -profile slurm,singularity  --csv_file PE_file.csv --target_region 16S --F_primer AGAGTTTGATCCTGGCTCAG --R_primer CTGCCTCCCGTAGGAGT 
```

<br>

#### 3c. Approach 3: Run jobs locally in conda environments and specify the path to one or more existing conda environment(s)

```bash
nextflow run main.nf -resume -profile conda --csv_file SE_file.csv --target_region 16S --F_primer AGAGTTTGATCCTGGCTCAG --R_primer CTGCCTCCCGTAGGAGT --conda.qc <path/to/existing/conda/environment>
```

<br>

**Required Parameters For All Approaches:**

* `-run main.nf` - Instructs nextflow to run the NF_XXX workflow 
* `-resume` - Resumes  workflow execution using previously cached results
* `-profile` – Specifies the configuration profile(s) to load, `singularity` instructs nextflow to setup and use singularity for all software called in the workflow
* `--target_region` – Specifies the amplicon target region to be analyzed, 16S, 18S or ITS.
  

  *Required only if you would like to pull and process data directly from OSDR*

* `--GLDS_accession` – A Genelab / OSD accession number e.g. GLDS-487.

*Required only if --GLDS_accession is not passed as an argument*

* `--csv_file` –  A 2-column input file with these headers [sample_id, read]. Please see the sample `file.csv` in this repository for an example on how to format this file.

* `--F_primer` – Forward primer sequence.

* `--R_primer` – Reverse primer sequence.

> See `nextflow run -h` and [Nextflow's CLI run command documentation](https://nextflow.io/docs/latest/cli.html#run) for more options and details on how to run nextflow.

<br>

#### 3d. Modify parameters and cpu resources in the nextflow config file

Additionally, the parameters and workflow resources can be directly specified in the nextflow.config file. For detailed instructions on how to modify and set parameters in the nextflow.config file, please see the [documentation here](https://www.nextflow.io/docs/latest/config.html).

Once you've downloaded the workflow template, you can modify the parameters in the `params` scope and cpus/memory requirements in the `process` scope in your downloaded version of the [nextflow.config](workflow_code/nextflow.config) file as needed in order to match your dataset and system setup. For example, you can directly set the the full paths to available conda environments in the `conda` scope within the `params`  scope. Additionally, if necessary, you'll need to modify each variable in the nexflow.config file to be consistent with the study you want to process and the machine you're using.

### 4.  Workflow outputs

#### 4a. Main outputs

The outputs from this pipeline are documented in the [GL-DPPD-7106](../../Pipeline_GL-DPPD-7104_Versions/GL-DPPD-7104-B.md) processing protocol.

#### 4b. Resource logs

Standard nextflow resource usage logs are also produced as follows:

- Output:
  - Resource_Usage/execution_report_{timestamp}.html (an html report that includes metrics about the workflow execution including computational resources and exact workflow process commands)
  - Resource_Usage/execution_timeline_{timestamp}.html (an html timeline for all processes executed in the workflow)
  - Resource_Usage/execution_trace_{timestamp}.txt (an execution tracing file that contains information about each process executed in the workflow, including: submission time, start time, completion time, cpu and memory used, machine-readable output)

> Further details about these logs can also found within [this Nextflow documentation page](https://www.nextflow.io/docs/latest/tracing.html#execution-report).

