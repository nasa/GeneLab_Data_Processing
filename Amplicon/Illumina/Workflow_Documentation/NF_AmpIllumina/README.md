# Workflow Information and Usage Instructions

## General Workflow Info

### Implementation Tools

The current GeneLab Illumina amplicon sequencing data processing pipeline (AmpIllumina), [GL-DPPD-7104-C.md](../../Pipeline_GL-DPPD-7104_Versions/GL-DPPD-7104-C.md), is implemented as a [Nextflow](https://nextflow.io/) DSL2 workflow and utilizes [Singularity](https://docs.sylabs.io/guides/3.10/user-guide/introduction.html) containers, [Docker](https://docs.docker.com/get-started/) containers, or [conda](https://docs.conda.io/en/latest/) environments to install/run all tools. This workflow is run using the command line interface (CLI) of any unix-based system.  While knowledge of creating workflows in Nextflow is not required to run the workflow as is, [the Nextflow documentation](https://nextflow.io/docs/latest/index.html) is a useful resource for users who want to modify and/or extend this workflow.   

### Resource Requirements <!-- omit in toc -->

The table below details the default maximum resource allocations for individual Nextflow processes.

| CPU Cores | Memory |
|-------------------|----------------|
| 2                 | 5 GB           |

> **Note:** These per-process resource allocations are defaults. They can be adjusted by modifying `cpus` and `memory` directives in the configuration file: [`nextflow.config`](workflow_code/nextflow.config).

## Utilizing the Workflow

1. [Install Nextflow, Singularity, and Conda](#1-install-nextflow-singularity-and-conda)  
   1a. [Install Nextflow and Conda](#1a-install-nextflow-and-conda)  
   1b. [Install Singularity](#1b-install-singularity)  

2. [Download the Workflow Files](#2-download-the-workflow-files)  

3. [Fetch Singularity Images](#3-fetch-singularity-images)  

4. [Run the Workflow](#4-run-the-workflow)  
   4a. [Approach 1: Start with OSD or GLDS accession as input](#4a-approach-1-start-with-an-osd-or-glds-accession-as-input)  
   4b. [Approach 2: Start with a runsheet csv file as input](#4b-approach-2-start-with-a-runsheet-csv-file-as-input)  
   4c. [Modify parameters and compute resources in the Nextflow config file](#4c-modify-parameters-and-compute-resources-in-the-nextflow-config-file)  

5. [Workflow Outputs](#5-workflow-outputs)  
   5a. [Main outputs](#5a-main-outputs)  
   5b. [Resource logs](#5b-resource-logs)  

6. [Post-processing](#7-post-processing)  

<br>

---

### 1. Install Nextflow, Singularity, and Conda

#### 1a. Install Nextflow and Conda

Nextflow can be installed either through [Anaconda](https://anaconda.org/bioconda/nextflow) or as documented on the [Nextflow documentation page](https://www.nextflow.io/docs/latest/getstarted.html).

> Note: If you want to install Anaconda, we recommend installing a Miniconda version appropriate for your system, as instructed by [Happy Belly Bioinformatics](https://astrobiomike.github.io/unix/conda-intro#getting-and-installing-conda).  
> 
> Once conda is installed on your system, you can install the latest version of Nextflow by running the following commands:
> 
> ```bash
> conda install -c bioconda nextflow
> nextflow self-update
> ```
> You may also install [mamba](https://mamba.readthedocs.io/en/latest/index.html) first which is a faster implementation of conda and can be used as a drop-in replacement:
> ```bash
> conda install -c conda-forge mamba
> ```

<br>

#### 1b. Install Singularity

Singularity is a container platform that allows usage of containerized software. This enables the GeneLab workflow to retrieve and use all software required for processing without the need to install the software directly on the user's system.

We recommend installing Singularity on a system wide level as per the associated [documentation](https://docs.sylabs.io/guides/3.10/admin-guide/admin_quickstart.html).

> Note: Singularity is also available through [Anaconda](https://anaconda.org/conda-forge/singularity).

> Note: Alternatively, Docker can be used in place of Singularity. To get started with Docker, see the [Docker CE installation documentation](https://docs.docker.com/engine/install/).

<br>

---

### 2. Download the Workflow Files

All files required for utilizing the NF_AmpIllumina GeneLab workflow for processing amplicon Illumina data are in the [workflow_code](workflow_code) directory. To get a copy of the latest *NF_AmpIllumina* version on to your system, the code can be downloaded as a zip file from the release page then unzipped after downloading by running the following commands: 

```bash
wget https://github.com/nasa/GeneLab_Data_Processing/releases/download/NF_AmpIllumina_1.0.0/NF_AmpIllumina_1.0.0.zip
unzip NF_AmpIllumina_1.0.0.zip && cd NF_AmpIllumina_1.0.0
```

<br>

---

### 3. Fetch Singularity Images

Although Nextflow can fetch Singularity images from a url, doing so may cause issues as detailed [here](https://github.com/nextflow-io/nextflow/issues/1210).

To avoid this issue, run the following command to fetch the Singularity images prior to running the NF_AmpIllumina workflow:

> Note: This command should be run in the location containing the `NF_AmpIllumina_1.0.0` directory that was downloaded in [step 2](#2-download-the-workflow-files) above. Depending on your network speed, fetching the images will take ~20 minutes. Approximately 4GB of RAM is needed to download and build the Singularity images.

```bash
bash ./bin/prepull_singularity.sh nextflow.config
```

Once complete, a `singularity` folder containing the Singularity images will be created. Run the following command to export this folder as a Nextflow configuration environment variable to ensure Nextflow can locate the fetched images:

```bash
export NXF_SINGULARITY_CACHEDIR=$(pwd)/singularity
```

<br>

---

### 4. Run the Workflow

> ***Note:** All the commands in this step must be run from within the `NF_AmpIllumina_1.0.0` directory that was downloaded in [step 2](#2-download-the-workflow-files) above.*

For options and detailed help on how to run the workflow, run the following command:

```bash
nextflow run main.nf --help
```

> Note: Nextflow commands use both single hyphen arguments (e.g. -help) that denote general Nextflow arguments and double hyphen arguments (e.g. --input_file) that denote workflow specific parameters.  Take care to use the proper number of hyphens for each argument.

> Note: This workflow assumes that all your raw reads end with the same suffix. If they don't, please modify your input read filenames to have the same suffix as shown in [SE_file.csv](workflow_code/SE_file.csv) and [PE_file.csv](workflow_code/PE_file.csv).

<br>

#### 4a. Approach 1: Start with OSD or GLDS accession as input

```bash
nextflow run main.nf \
   -resume \
   -profile singularity \
   --target_region 16S \
   --accession OSD-487 
```

<br>

#### 4b. Approach 2: Start with a runsheet csv file as input

```bash
nextflow run main.nf \
   -resume \
   -profile singularity \
   --target_region 16S \
   --input_file PE_file.csv \
   --F_primer AGAGTTTGATCCTGGCTCAG \
   --R_primer CTGCCTCCCGTAGGAGT 
```

<br>


**Required Parameters For All Approaches:**

* `main.nf` - Instructs Nextflow to run the NF_AmpIllumina workflow 
* `-resume` - Resumes  workflow execution using previously cached results
* `-profile` – Specifies the configuration profile(s) to load (multiple options can be provided as a comma-separated list)
   * Software environment profile options (choose one):
      * `singularity` - instructs Nextflow to use Singularity container environments
      * `docker` - instructs Nextflow to use Docker container environments
      * `conda` - instructs Nextflow to use conda environments via the conda package manager
        > *Note: By default, Nextflow will create environments at runtime using the yaml files in the [workflow_code/envs](workflow_code/envs/) folder. You can change this behavior by using the `--conda_*` workflow parameters or by editing the [nextflow.config](workflow_code/nextflow.config) file to specify a centralized conda environments directory via the `conda.cacheDir` parameter.*
      * `mamba` - instructs Nextflow to use conda environments via the mamba package manager 
   * Other option (can be combined with the software environment option above using a comma, e.g. `-profile slurm,singularity`):
      * `slurm` - instructs Nextflow to use the [Slurm cluster management and job scheduling system](https://slurm.schedmd.com/overview.html) to schedule and run the jobs on a Slurm HPC cluster
* `--target_region` - Specifies the amplicon target region to be analyzed, available option are: 16S, 18S, or ITS

  *Required only if you would like to pull and process data directly from OSDR*

* `--accession` - A GeneLab / OSD accession number e.g. GLDS-487.


**Additional Required Parameters For Approach 1:** 

* `--accession` - The OSD or GLDS accession number specifying the [OSDR](https://osdr.nasa.gov/bio/repo/) dataset to process, e.g. OSD-487 or GLDS-487
  > *Note: Not all datasets have the same OSD and GLDS number, so make sure the correct OSD or GLDS number is specified*


**Additional Required Parameters For Approach 2:** 

* `--input_file` –  A single-end or paired-end runsheet csv file containing assay metadata for each sample, including sample_id, forward (path to forward read), [reverse (path to reverse read, for paired-end only),] paired (boolean, TRUE | FALSE), groups (specifies sample treatment group name). Please see the [runsheet documentation](./examples/runsheet) in this repository for examples on how to format this file.

* `--F_primer` - Forward primer sequence

* `--R_primer` - Reverse primer sequence


**Additional [Optional] Parameters For ALL Approaches**
> *Note: See `nextflow run -h` and [Nextflow's CLI run command documentation](https://nextflow.io/docs/latest/cli.html#run) for more options and details on how to run Nextflow.*

* `--errorStrategy` - Error handling strategy for Nextflow processes . If processes fail, use "ignore" to allow the workflow to continue running (default: "terminate")
* `--primers_linked` - Whether forward and reverse primers are linked (default: TRUE)
* `--anchored_primers` - Whether primers are anchored at the start of reads (default: TRUE)
* `--min_cutadapt_len` - Minimum length of reads to keep after Cutadapt trimming (default: 130) 
* `--discard_untrimmed` - Whether to discard untrimmed reads (default: TRUE)
* `--left_trunc` - Truncate forward reads after this many bases. Reads shorter than this are discarded (default: 0; no truncation)
* `--right_trunc` - Truncate reverse reads after this many bases. Reads shorter than this are discarded (default: 0; no truncation)
* `--left_maxEE` - Maximum expected errors allowed in forward reads (default: 1)
* `--right_maxEE` - Maximum expected errors allowed in reverse reads (default: 1)
* `--concatenate_reads_only` - Whether to concatenate paired reads end-to-end instead of merging based on overlapping regions (default: FALSE)
* `--diff_abund_method` - Differential abundance testing method to use: "all", "ancombc1", "ancombc2", or "deseq2" (default: "all")
* `--group` - Column name in input CSV file containing groups to be compared (default: "groups")
* `--samples_column` - Column name in input CSV file containing sample names (default: "sample_id")
* `--remove_struc_zeros` - Whether to remove structural zeros when running ANCOMBC (default: false)
* `--remove_rare` - Whether to filter out rare features and samples with low library sizes (default: false)
* `--prevalence_cutoff` - Taxa with prevalence below this fraction will be excluded.(default: 0.15)
* `--library_cutoff` - Samples with library sizes below this threshold will be excluded. (default: 100)
* `--conda_cutadapt` - Path to existing Cutadapt conda environment (default: null; creates new environment)
* `--conda_diversity` - Path to existing R diversity analysis conda environment (default: null; creates new environment)
* `--conda_dp_tools` - Path to existing dp_tools conda environment (default: null; creates new environment)
* `--conda_fastqc` - Path to existing FastQC conda environment (default: null; creates new environment)
* `--conda_multiqc` - Path to existing MultiQC conda environment (default: null; creates new environment)
* `--conda_R` - Path to existing R conda environment (default: null; creates new environment)
* `--conda_zip` - Path to existing zip conda environment (default: null; creates new environment)

<br>

#### 4c. Modify parameters and compute resources in the Nextflow config file

Additionally, all parameters and workflow resources can be directly specified in the [nextflow.config](./workflow_code/nextflow.config) file. For detailed instructions on how to modify and set parameters in the config file, please see the [documentation here](https://www.nextflow.io/docs/latest/config.html).

Once you've downloaded the workflow template, you can modify the parameters in the `params` scope and cpus/memory requirements in the `process` scope in your downloaded version of the [nextflow.config](workflow_code/nextflow.config) file as needed in order to match your dataset and system setup. Additionally, if necessary, you can modify each variable in the [nextflow.config](workflow_code/nextflow.config) file to be consistent with the study you want to process and the computer you're using for processing.

<br>

---

### 5. Workflow Outputs

#### 5a. Main Outputs

The outputs from this pipeline are documented in the [GL-DPPD-7104-C](../../Pipeline_GL-DPPD-7104_Versions/GL-DPPD-7104-C.md) processing protocol.

#### 5b. Resource Logs

Standard Nextflow resource usage logs are also produced as follows:

**Nextflow Resource Usage Logs**
   - Resource_Usage/execution_report_{timestamp}.html (an html report that includes metrics about the workflow execution including computational resources and exact workflow process commands)
   - Resource_Usage/execution_timeline_{timestamp}.html (an html timeline for all processes executed in the workflow)
   - Resource_Usage/execution_trace_{timestamp}.txt (an execution tracing file that contains information about each process executed in the workflow, including: submission time, start time, completion time, cpu and memory used, machine-readable output)

> Further details about these logs can also found within [this Nextflow documentation page](https://www.nextflow.io/docs/latest/tracing.html#execution-report).

<br>

---

### 6. Post-processing

> Please note that to run the post-processing workflow successfully, you MUST run the processing workflow above via the [launch.sh](workflow_code/launch.sh) script first. Please see the [script](workflow_code/launch.sh) for how to run it and make sure to edit the place holders before running it.

The post-processing workflow generates a README file, a protocols file, an md5sums table, and a file association table suitable for uploading to OSDR.

For options and detailed help on how to run the post-processing workflow, run the following command:

```bash
nextflow run post_processing.nf --help
```

To generate the post-processing files after running the main processing workflow successfully, modify and set the parameters in [post_processing.config](workflow_code/post_processing.config), then run the following command:

```bash
nextflow run post_processing.nf \
   -c post_processing.config \ 
   -resume \
   -profile singularity
``` 

The outputs of the post-processing workflow are described below:
> *Note: The outputs will be in a directory called `Post_Processing` by default*

**Post processing workflow output files** 
 - Post_processing/FastQC_Outputs/filtered_multiqc_GLAmpSeq_report.zip (Filtered sequence MultiQC report with paths purged) 
 - Post_processing/FastQC_Outputs/raw_multiqc_GLAmpSeq_report.zip (Raw sequence MultiQC report with paths purged)
 - Post_processing/<GLDS_accession>_-associated-file-names.tsv (File association table for OSDR curation)
 - Post_processing/<GLDS_accession>_amplicon-validation.log (Automated verification and validation log file)
 - Post_processing/processed_md5sum_GLAmpSeq.tsv (md5sums for the files published on OSDR)
 - Post_processing/processing_info_GLAmpSeq.zip  (Zip file containing all files used to run the workflow and required logs with paths purged) 
 - Post_processing/protocol.txt  (File describing the methods used by the workflow)
 - Post_processing/README_GLAmpSeq.txt (README file listing and describing the outputs of the workflow)
