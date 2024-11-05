# Workflow Information and Usage Instructions

## General Workflow Info

### Implementation Tools

The current GeneLab Illumina metagenomics sequencing data processing pipeline (MGIllumina-A), [GL-DPPD-7107-A.md](../../Pipeline_GL-DPPD-7107_Versions/GL-DPPD-7107-A.md), is implemented as a [Nextflow](https://nextflow.io/) DSL2 workflow and utilizes [Singularity](https://docs.sylabs.io/guides/3.10/user-guide/introduction.html) containers or [conda](https://docs.conda.io/en/latest/) environments to install/run all tools. This workflow is run using the command line interface (CLI) of any unix-based system.  While knowledge of creating workflows in nextflow is not required to run the workflow as is, [the Nextflow documentation](https://nextflow.io/docs/latest/index.html) is a useful resource for users who want to modify and/or extend this workflow. 

> **Note on reference databases**  
> Many reference databases are relied upon throughout this workflow. They will be installed and setup automatically the first time the workflow is run. All together, after installed and unpacked, they will take up about about 340 GB of storage, but they may also require up to 500GB during installation and initial un-packing, so be sure there is enough room on your system before running the workflow.

<br>

## Utilizing the Workflow

1. [Install Nextflow and Singularity](#1-install-nextflow-and-singularity)  
   1a. [Install Nextflow](#1a-install-nextflow)  
   1b. [Install Singularity](#1b-install-singularity)  

2. [Download the workflow files](#2-download-the-workflow-files)  

3. [Fetch Singularity Images](#3-fetch-singularity-images)  

4. [Run the workflow](#4-run-the-workflow)  
   4a. [Approach 1: Run slurm jobs in singularity containers with OSD accession as input](#4a-approach-1-run-slurm-jobs-in-singularity-containers-with-osd-accession-as-input)   
   4b. [Approach 2: Run slurm jobs in singularity containers with a csv file as input](#4b-approach-2-run-slurm-jobs-in-singularity-containers-with-a-csv-file-as-input)  
   4c. [Approach 3: Run jobs locally in conda environments and specify the path to one or more existing conda environments](#4c-approach-run-jobs-locally-in-conda-environments-and-specify-the-path-to-one-or-more-existing-conda-environments)  
   4d. [Modify parameters and cpu resources in the nextflow config file](#4d-modify-parameters-and-cpu-resources-in-the-nextflow-config-file)  

5. [Workflow outputs](#5-workflow-outputs)  
   5a. [Main outputs](#5a-main-outputs)  
   5b. [Resource logs](#5b-resource-logs)  

6. [Post Processing](#6-post-processing)  

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

<br>

#### 1b. Install Singularity

Singularity is a container platform that allows usage of containerized software. This enables the GeneLab workflow to retrieve and use all software required for processing without the need to install the software directly on the user's system.

We recommend installing Singularity on a system wide level as per the associated [documentation](https://docs.sylabs.io/guides/3.10/admin-guide/admin_quickstart.html).

> Note: Singularity is also available through [Anaconda](https://anaconda.org/conda-forge/singularity).

<br>

---

### 2. Download the workflow files

All files required for utilizing the NF_MGIllumina-A GeneLab workflow for processing metagenomics Illumina data are in the [workflow_code](workflow_code) directory. To get a copy of latest *NF_MGIllumina-A* version on to your system, the code can be downloaded as a zip file from the release page then unzipped after downloading by running the following commands: 

```bash
wget https://github.com/nasa/GeneLab_Data_Processing/releases/download/NF_MGIllumina-A_1.0.0/NF_MGIllumina-A_1.0.0.zip
unzip NF_MGIllumina-A_1.0.0.zip &&  cd NF_MGIllumina-A_1.0.0
```

<br>

---

### 3. Fetch Singularity Images

Although Nextflow can fetch Singularity images from a url, doing so may cause issues as detailed [here](https://github.com/nextflow-io/nextflow/issues/1210).

To avoid this issue, run the following command to fetch the Singularity images prior to running the NF_MGIllumina-A workflow:

> Note: This command should be run from within the `NF_MGIllumina-A_1.0.0` directory that was downloaded in [step 2](#2-download-the-workflow-files) above.  

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

> ***Note:** All the commands in this step must be run from within the `NF_MGIllumina-A_1.0.0` directory that was downloaded in [step 2](#2-download-the-workflow-files) above.*

For options and detailed help on how to run the workflow, run the following command:

```bash
nextflow run main.nf --help
```

> Note: Nextflow commands use both single hyphen arguments (e.g. -help) that denote general nextflow arguments and double hyphen arguments (e.g. --csv_file) that denote workflow specific parameters.  Take care to use the proper number of hyphens for each argument.

<br>

#### 4a. Approach 1: Run slurm jobs in singularity containers with OSD or GLDS accession as input

```bash
nextflow run main.nf -resume -profile slurm,singularity --accession OSD-574
```

<br>

#### 4b. Approach 2: Run slurm jobs in singularity containers with a csv file as input

```bash
nextflow run main.nf -resume -profile slurm,singularity  --csv_file PE_file.csv
```

<br>

#### 4c. Approach 3: Run jobs locally in conda environments and specify the path to one or more existing conda environment(s)

```bash
nextflow run main.nf -resume -profile conda --csv_file SE_file.csv --conda.qc <path/to/existing/conda/environment>
```

<br>

**Required Parameters For All Approaches:**

* `-run main.nf` - Instructs nextflow to run the NF_MGIllumina-A workflow 

* `-resume` - Resumes workflow execution using previously cached results

* `-profile` – Specifies the configuration profile(s) to load, `singularity` instructs nextflow to setup and use singularity for all software called in the workflow

*Required only if you would like to pull and process data directly from OSDR*

* `--accession` – A Genelab / OSD accession number e.g. OSD-574.

*Required only if --accession is not passed as an argument*

* `--csv_file` –  A single-end or paired-end input csv file containing assay metadata for each sample, including sample_id, forward, reverse, and/or paired. Please see the sample [SE_file.csv](workflow_code/SE_file.csv) and [PE_file.csv](workflow_code/PE_file.csv) in this repository for examples on how to format this file.

> See `nextflow run -h` and [Nextflow's CLI run command documentation](https://nextflow.io/docs/latest/cli.html#run) for more options and details on how to run nextflow.

<br>

#### 4d. Modify parameters and cpu resources in the nextflow config file

Additionally, the parameters and workflow resources can be directly specified in the nextflow.config file. For detailed instructions on how to modify and set parameters in the nextflow.config file, please see the [documentation here](https://www.nextflow.io/docs/latest/config.html).

Once you've downloaded the workflow template, you can modify the parameters in the `params` scope and cpus/memory requirements in the `process` scope in your downloaded version of the [nextflow.config](workflow_code/nextflow.config) file as needed in order to match your dataset and system setup. For example, you can directly set the the full paths to available conda environments in the `conda` scope within the `params` scope. Additionally, if necessary, you'll need to modify each variable in the [nextflow.config](workflow_code/nextflow.config) file to be consistent with the study you want to process and the machine you're using.

<br>

---

### 5. Workflow outputs

#### 5a. Main outputs

The outputs from this pipeline are documented in the [GL-DPPD-7107-A](../../Pipeline_GL-DPPD-7107_Versions/GL-DPPD-7107-A.md) processing protocol.

#### 5b. Resource logs

Standard nextflow resource usage logs are also produced as follows:

- **Output:**
  - Resource_Usage/execution_report_{timestamp}.html (an html report that includes metrics about the workflow execution including computational resources and exact workflow process commands)
  - Resource_Usage/execution_timeline_{timestamp}.html (an html timeline for all processes executed in the workflow)
  - Resource_Usage/execution_trace_{timestamp}.txt (an execution tracing file that contains information about each process executed in the workflow, including: submission time, start time, completion time, cpu and memory used, machine-readable output)

> Further details about these logs can also found within [this Nextflow documentation page](https://www.nextflow.io/docs/latest/tracing.html#execution-report).

<br>

---

### 6. Post Processing

For options and detailed help on how to run the post-processing workflow, run the following command:

```bash
nextflow run post_processing.nf --help
```

To generate a README file, a protocols file, a md5sums table and a file association table after running the processing workflow sucessfully, modify and set the parameters in [post_processing.config](workflow_code/post_processing.config) then run the following command:

```bash
nextflow -C post_processing.config run post_processing.nf -resume -profile slurm,singularity
``` 

The outputs of the run will be in a directory called `Post_Processing` by default and they are as follows:

 - Post_processing/FastQC_Outputs/filtered_multiqc_GLmetagenomics_report.zip (Filtered sequence multiqc report with paths purged)  

 - Post_processing/FastQC_Outputs/raw_multiqc_GLmetagenomics_report.zip (Raw sequence multiqc report with paths purged)  

 - Post_processing/<GLDS_accession>_-associated-file-names.tsv (File association table for curation)  

 - Post_processing/<GLDS_accession>_metagenomics-validation.log (Automatic verification and validation log file)  

 - Post_processing/processed_md5sum_GLmetagenomics.tsv (md5sums for the files to be released on OSDR)  

 - Post_processing/processing_info_GLmetagenomics.zip  (Zip file containing all files used to run the workflow and required logs with paths purged)  

 - Post_processing/protocol.txt  (File describing the methods used by the workflow)  

 - Post_processing/README_GLmetagenomics.txt (README file listing and describing the outputs of the workflow)

