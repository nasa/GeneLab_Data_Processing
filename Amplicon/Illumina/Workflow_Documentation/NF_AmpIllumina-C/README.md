# Workflow Information and Usage Instructions

## General Workflow Info

### Implementation Tools

The current GeneLab Illumina amplicon sequencing data processing pipeline (AmpIllumina), [GL-DPPD-7104-C.md](../../Pipeline_GL-DPPD-7104_Versions/GL-DPPD-7104-C.md), is implemented as a [Nextflow](https://nextflow.io/) DSL2 workflow and utilizes [Singularity](https://docs.sylabs.io/guides/3.10/user-guide/introduction.html) containers or [conda](https://docs.conda.io/en/latest/) environments to install/run all tools. This workflow is run using the command line interface (CLI) of any unix-based system.  While knowledge of creating workflows in Nextflow is not required to run the workflow as is, [the Nextflow documentation](https://nextflow.io/docs/latest/index.html) is a useful resource for users who want to modify and/or extend this workflow.   

### Resource Requirements <!-- omit in toc -->

The table below details the default maximum resource allocations for individual Nextflow processes.

| Default CPU Cores | Default Memory |
|-------------------|----------------|
| 2                 | 5 GB           |

> **Note:** These per-process resource allocations are defaults. They can be adjusted by modifying `cpus` and `memory` directives in the configuration file: [`nextflow.config`](workflow_code/nextflow.config).

## Utilizing the Workflow

1. [Install Nextflow and Singularity](#1-install-nextflow-and-singularity)  
   1a. [Install Nextflow](#1a-install-nextflow)  
   1b. [Install Singularity](#1b-install-singularity)  

2. [Download the Workflow Files](#2-download-the-workflow-files)  

3. [Fetch Singularity Images](#3-fetch-singularity-images)  

4. [Run the Workflow Directly with Nextflow](#4-run-the-workflow-directly-with-nextflow)  
   4a. [Approach 1: Run Slurm jobs in Singularity containers with OSD or GLDS accession as input](#4a-approach-1-run-slurm-jobs-in-singularity-containers-with-osd-or-glds-accession-as-input)  
   4b. [Approach 2: Run Slurm jobs in Singularity containers with a csv file as input](#4b-approach-2-run-slurm-jobs-in-singularity-containers-with-a-csv-file-as-input)  
   4c. [Approach 3: Run jobs locally in conda or mamba environments and specify the path to one or more existing conda environments](#4c-approach-3-run-jobs-locally-in-conda-or-mamba-environments-and-specify-the-path-to-one-or-more-existing-conda-environments)  
   4d. [Modify parameters and cpu resources in the Nextflow config file](#4d-modify-parameters-and-cpu-resources-in-the-nextflow-config-file)  

5. [Run the Workflow Indirectly Using the Python Wrapper Script](#5-run-the-workflow-indirectly-using-the-python-wrapper-script)  
   5a. [Approach 1: Use an OSD or GeneLab acession as input](#5a-approach-1-use-an-osd-or-genelab-acession-as-input)  
   5b. [Approach 2: Use a csv file as input to the workflow](#5b-approach-2-use-a-csv-file-as-input-to-the-workflow)  
   5c. [Approach 3: Use a csv file as input to the workflow and supply extra arguments to Nextflow run](#5c-approach-3-use-a-csv-file-as-input-to-the-workflow-and-supply-extra-arguments-to-nextflow-run)  
   5d. [Approach 4: Create an edited Nextflow config file but do not run the workflow](#5d-approach-4-create-an-edited-nextflow-config-file-but-do-not-run-the-workflow)  

6. [Workflow Outputs](#6-workflow-outputs)  
   6a. [Main Outputs](#6a-main-outputs)  
   6b. [Resource Logs](#6b-resource-logs)  

7. [Post-processing](#7-post-processing)  

<br>

---

### 1. Install Nextflow and Singularity

#### 1a. Install Nextflow

Nextflow can be installed either through [Anaconda](https://anaconda.org/bioconda/nextflow) or as documented on the [Nextflow documentation page](https://www.nextflow.io/docs/latest/getstarted.html).

> Note: If you want to install Anaconda, we recommend installing a Miniconda, Python3 version appropriate for your system, as instructed by [Happy Belly Bioinformatics](https://astrobiomike.github.io/unix/conda-intro#getting-and-installing-conda).  
> 
> Once conda is installed on your system, you can install the latest version of Nextflow by running the following commands:
> 
> ```bash
> conda install -c bioconda nextflow
> nextflow self-update
> ```
> You can also install [mamba](https://mamba.readthedocs.io/en/latest/index.html) which is a faster implementation of conda like so:
> ```bash
> conda install -c conda-forge mamba
> ```

<br>

#### 1b. Install Singularity

Singularity is a container platform that allows usage of containerized software. This enables the GeneLab workflow to retrieve and use all software required for processing without the need to install the software directly on the user's system.

We recommend installing Singularity on a system wide level as per the associated [documentation](https://docs.sylabs.io/guides/3.10/admin-guide/admin_quickstart.html).

> Note: Singularity is also available through [Anaconda](https://anaconda.org/conda-forge/singularity).

> Note: Alternatively, Docker can be used in place of Singularity. See the [Docker CE installation documentation](https://docs.docker.com/engine/install/).

<br>

---

### 2. Download the Workflow Files

All files required for utilizing the NF_AmpIllumina GeneLab workflow for processing amplicon illumina data are in the [workflow_code](workflow_code) directory. To get a copy of the latest NF_AmpIllumina version on to your system, the code can be downloaded as a zip file from the release page then unzipped after downloading by running the following commands: 

```bash
wget https://github.com/nasa/GeneLab_Data_Processing/releases/download/NF_AmpIllumina_1.0.0/NF_AmpIllumina_1.0.0.zip
unzip NF_AmpIllumina_1.0.0.zip &&  cd NF_AmpIllumina_1.0.0
```

<br>

---

### 3. Fetch Singularity Images

Although Nextflow can fetch Singularity images from a url, doing so may cause issues as detailed [here](https://github.com/nextflow-io/nextflow/issues/1210).

To avoid this issue, run the following command to fetch the Singularity images prior to running the NF_AmpIllumina workflow:

> Note: This command should be run in the location containing the `NF_AMPIllumina` directory that was downloaded in [step 2](#2-download-the-workflow-files) above. Depending on your network speed, fetching the images will take ~20 minutes. Approximately 4GB of RAM is needed to download and build the Singularity images.

```bash
bash ./bin/prepull_singularity.sh nextflow.config
```

Once complete, a `singularity` folder containing the Singularity images will be created. Run the following command to export this folder as a Nextflow configuration environment variable to ensure Nextflow can locate the fetched images:

```bash
export NXF_SINGULARITY_CACHEDIR=$(pwd)/singularity
```

<br>

---

### 4. Run the Workflow Directly with Nextflow

While in the location containing the `NF_AmpIllumina_1.0.0` directory that was downloaded in [step 2](#2-download-the-workflow-files), you are now able to run the workflow.

For options and detailed help on how to run the workflow, run the following command:

```bash
nextflow run main.nf --help
```

> Note: Nextflow commands use both single hyphen arguments (e.g. -help) that denote general Nextflow arguments and double hyphen arguments (e.g. --input_file) that denote workflow specific parameters.  Take care to use the proper number of hyphens for each argument.

> Note: This workflow assumes that all your raw reads end with the same suffix. If they don't, please modify your input read filenames to have the same suffix as shown in [SE_file.csv](workflow_code/SE_file.csv) and [PE_file.csv](workflow_code/PE_file.csv).

> Note: To use Docker instead of Singularity, use `-profile docker` in the Nextflow run command. Nextflow will automatically pull images as needed.

<br>

#### 4a. Approach 1: Run Slurm jobs in Singularity containers with OSD or GLDS accession as input

```bash
nextflow run main.nf \
   -resume \
   -profile slurm,singularity \
   --accession GLDS-487 \
   --target_region 16S
```

<br>

#### 4b. Approach 2: Run Slurm jobs in Singularity containers with a csv file as input

```bash
nextflow run main.nf \
   -resume \
   -profile slurm,singularity \
   --input_file PE_file.csv \
   --target_region 16S \
   --F_primer AGAGTTTGATCCTGGCTCAG \
   --R_primer CTGCCTCCCGTAGGAGT 
```

<br>

#### 4c. Approach 3: Run jobs locally in conda or mamba environments and specify the path to one or more existing conda environment(s)

```bash
nextflow run main.nf \
   -resume \
   -profile mamba \
   --input_file SE_file.csv \
   --target_region 16S \
   --F_primer AGAGTTTGATCCTGGCTCAG \
   --R_primer CTGCCTCCCGTAGGAGT \
   --conda_R <path/to/existing/conda/environment>
```

<br>

**Required Parameters For All Approaches:**

* `main.nf` - Instructs Nextflow to run the NF_AmpIllumina workflow 
* `-resume` - Resumes  workflow execution using previously cached results
* `-profile` - Specifies the configuration profile(s) to load, `singularity` instructs Nextflow to setup and use Singularity for all software called in the workflow, `conda` or `mamba` instructs Nextflow to use conda environments
* `--target_region` - Specifies the amplicon target region to be analyzed, 16S, 18S or ITS.

  *Required only if you would like to pull and process data directly from OSDR*

* `--accession` - A GeneLab / OSD accession number e.g. GLDS-487.

**Required only if --accession is not passed as an argument**

* `--input_file` -  A 4-column (single-end) or 5-column (paired-end) input csv file with the following headers (sample_id, forward, [reverse,] paired, groups). Please see the sample [SE_file.csv](workflow_code/SE_file.csv) and [PE_file.csv](workflow_code/PE_file.csv) in this repository for examples on how to format this file.

* `--F_primer` - Forward primer sequence.

* `--R_primer` - Reverse primer sequence.

> See `nextflow run -h` and [Nextflow's CLI run command documentation](https://nextflow.io/docs/latest/cli.html#run) for more options and details on how to run Nextflow.

<br>

#### 4d. Modify parameters and cpu resources in the Nextflow config file

Additionally, the parameters and workflow resources can be directly specified in the nextflow.config file. For detailed instructions on how to modify and set parameters in the nextflow.config file, please see the [documentation here](https://www.nextflow.io/docs/latest/config.html).

Once you've downloaded the workflow template, you can modify the parameters in the `params` scope and cpus/memory requirements in the `process` scope in your downloaded version of the [nextflow.config](workflow_code/nextflow.config) file as needed in order to match your dataset and system setup. Additionally, if necessary, you'll need to modify each variable in the  [nextflow.config](workflow_code/nextflow.config) file to be consistent with the study you want to process and the machine you're using.

**Additional Parameters**

* `--errorStrategy` - Error handling strategy for Nextflow processes . If processes fail, use "ignore" to allow the workflow to continue running (default: "terminate")
* `--primers_linked` - Whether forward and reverse primers are linked (default: TRUE)
* `--anchored_primers` - Whether primers are anchored at the start of reads (default: TRUE)
* `--min_cutadapt_len` - Minimum length of reads to keep after Cutadapt primer trimming (default: 130) 
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
* `--conda_diversity` - Path to existing R diversity conda environment (default: null; creates new environment)
* `--conda_dp_tools` - Path to existing dp_tools conda environment (default: null; creates new environment)
* `--conda_fastqc` - Path to existing FastQC conda environment (default: null; creates new environment)
* `--conda_multiqc` - Path to existing MultiQC conda environment (default: null; creates new environment)
* `--conda_R` - Path to existing R conda environment (default: null; creates new environment)
* `--conda_zip` - Path to existing zip conda environment (default: null; creates new environment)

<br>

---

### 5. Run the Workflow Indirectly Using the Python Wrapper Script

For options and detailed help on how to run the workflow using the script, run the following command:

```bash
python run_workflow.py
```

#### 5a. Approach 1: Use an OSD or GeneLab acession as input

```bash
python run_workflow.py \
   --run \
   --target-region 16S \
   --accession GLDS-487 \
   --profile slurm,singularity 
```

#### 5b. Approach 2: Use a csv file as input to the workflow

```bash
python run_workflow.py \
   --run \
   --target-region 16S \
   --input-file PE_file.csv \
   --F-primer AGAGTTTGATCCTGGCTCAG \
   --R-primer CTGCCTCCCGTAGGAGT \
   --profile singularity 
```

#### 5c. Approach 3: Use a csv file as input to the workflow and supply extra arguments to Nextflow run 

Here were want to monitor our jobs with Nextflow Tower.

```bash
export TOWER_ACCESS_TOKEN=<ACCESS TOKEN>
export TOWER_WORKSPACE_ID=<WORKSPACE ID>
python run_workflow.py \
   --run \
   --target-region 16S \
   --input-file PE_file.csv \
   --F-primer AGAGTTTGATCCTGGCTCAG \
   --R-primer CTGCCTCCCGTAGGAGT \
   --profile slurm,conda \
   --extra 'with-tower'
```

#### 5d. Approach 4: Create an edited Nextflow config file but do not run the workflow

```bash
python run_workflow.py \
   --target-region 16S \
   --accession GLDS-487 \
   --profile slurm,singularity
```

> Note: When using the wrapper script, all outputs generated by the workflow will be in a directory specified by the `--output-dir` parameter. This will be the parent directory `..` by default.

<br>

---

### 6. Workflow Outputs

#### 6a. Main Outputs

The outputs from this pipeline are documented in the [GL-DPPD-7104-C](../../Pipeline_GL-DPPD-7104_Versions/GL-DPPD-7104-C.md) processing protocol.

#### 6b. Resource Logs

Standard Nextflow resource usage logs are also produced as follows:

- Output:
  - Resource_Usage/execution_report_{timestamp}.html (an html report that includes metrics about the workflow execution including computational resources and exact workflow process commands)
  - Resource_Usage/execution_timeline_{timestamp}.html (an html timeline for all processes executed in the workflow)
  - Resource_Usage/execution_trace_{timestamp}.txt (an execution tracing file that contains information about each process executed in the workflow, including: submission time, start time, completion time, cpu and memory used, machine-readable output)

> Further details about these logs can also found within [this Nextflow documentation page](https://www.nextflow.io/docs/latest/tracing.html#execution-report).

<br>

---

### 7. Post-processing

> Please note that to run the post-processing workflow successfully, you MUST run the processing workflow above via the [launch.sh](workflow_code/launch.sh) script first. Please see the [script](workflow_code/launch.sh) for how to run it and make sure to edit the place holders before running it.

For options and detailed help on how to run the post-processing workflow, run the following command:

```bash
nextflow run post_processing.nf --help
```

To generate a README file, a protocols file, a md5sums table and a file association table after running the processing workflow sucessfully, modify and set the parameters in [post_processing.config](workflow_code/post_processing.config) then run the following command:

```bash
nextflow run post_processing.nf \
   -c post_processing.config 
   -resume \
   -profile slurm,singularity
``` 

The outputs of the run will be in a directory called `Post_Processing` by default and they are as follows:
 - Post_processing/FastQC_Outputs/filtered_multiqc_GLAmpSeq_report.zip (Filtered sequence MultiQC report with paths purged) 
 - Post_processing/FastQC_Outputs/raw_multiqc_GLAmpSeq_report.zip (Raw sequence MultiQC report with paths purged)
 - Post_processing/<GLDS_accession>_-associated-file-names.tsv (File association table for curation)
 - Post_processing/<GLDS_accession>_amplicon-validation.log (Automatic verification and validation log file)
 - Post_processing/processed_md5sum_GLAmpSeq.tsv (md5sums for the files to be released on OSDR)
 - Post_processing/processing_info_GLAmpSeq.zip  (Zip file containing all files used to run the workflow and required logs with paths purged) 
 - Post_processing/protocol.txt  (File describing the methods used by the workflow)
 - Post_processing/README_GLAmpSeq.txt (README file listing and describing the outputs of the workflow)
