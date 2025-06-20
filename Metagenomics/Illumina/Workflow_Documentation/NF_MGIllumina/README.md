# Workflow Information and Usage Instructions

## General Workflow Info

### Implementation Tools

The current GeneLab Illumina metagenomics sequencing data processing pipeline (MGIllumina-A), [GL-DPPD-7107-A.md](../../Pipeline_GL-DPPD-7107_Versions/GL-DPPD-7107-A.md), is implemented as a [Nextflow](https://nextflow.io/) DSL2 workflow and utilizes [Singularity](https://docs.sylabs.io/guides/3.10/user-guide/introduction.html) containers, [Docker](https://docs.docker.com/get-started/) containers, or [conda](https://docs.conda.io/en/latest/) environments to install/run all tools. This workflow is run using the command line interface (CLI) of any unix-based system.  While knowledge of creating workflows in Nextflow is not required to run the workflow as is, [the Nextflow documentation](https://nextflow.io/docs/latest/index.html) is a useful resource for users who want to modify and/or extend this workflow. 

> **Note on reference databases**  
> Many reference databases are relied upon throughout this workflow. They will be installed and setup automatically the first time the workflow is run. All together, after installed and unpacked, they will take up about about 340 GB of storage, but they may also require up to 500GB during installation and initial un-packing, so be sure there is enough room on your system before running the workflow.

<br>

## Utilizing the Workflow

1. [Installing Nextflow, Singularity, and conda](#1-installing-nextflow-singularity-and-conda)  
   1a. [Install Nextflow and conda](#1a-install-nextflow-and-conda)  
   1b. [Install Singularity](#1b-install-singularity)  
2. [Download the workflow files](#2-download-the-workflow-files)  
3. [Fetch Singularity Images](#3-fetch-singularity-images)  
4. [Run the workflow](#4-run-the-workflow)  
   4a. [Approach 1: Start with OSD or GLDS accession as input](#4a-approach-1-start-with-an-osd-or-glds-accession-as-input)    
   4b. [Approach 2: Start with a runsheet csv file as input](#4b-approach-2-start-with-a-runsheet-csv-file-as-input)  
   4c. [Modify parameters and compute resources in the Nextflow config file](#4c-modify-parameters-and-compute-resources-in-the-nextflow-config-file)
5. [Workflow outputs](#5-workflow-outputs)  
   5a. [Main outputs](#5a-main-outputs)  
   5b. [Resource logs](#5b-resource-logs)  
6. [Post Processing](#6-post-processing)  

<br>

---

### 1. Installing Nextflow, Singularity, and conda

#### 1a. Install Nextflow and conda

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
> conda install -c bioconda nextflow
> nextflow self-update
> ```

<br>

#### 1b. Install Singularity

Singularity is a container platform that allows usage of containerized software. This enables the GeneLab workflow to retrieve and use all software required for processing without the need to install the software directly on the user's system.

We recommend installing Singularity on a system wide level as per the associated [documentation](https://docs.sylabs.io/guides/3.10/admin-guide/admin_quickstart.html).

> Note: Singularity is also available through [Anaconda](https://anaconda.org/conda-forge/singularity).

> Note: Alternatively, Docker can be used in place of Singularity. To get started with Docker, see the [Docker CE installation documentation](https://docs.docker.com/engine/install/).

<br>

---

### 2. Download the workflow files

All files required for utilizing the NF_MGIllumina GeneLab workflow for processing metagenomics Illumina data are in the [workflow_code](workflow_code) directory. To get a copy of latest *NF_MGIllumina* version on to your system, the code can be downloaded as a zip file from the release page then unzipped after downloading by running the following commands: 

```bash
wget https://github.com/nasa/GeneLab_Data_Processing/releases/download/NF_MGIllumina_1.0.0/NF_MGIllumina_1.0.0.zip
unzip NF_MGIllumina_1.0.0.zip &&  cd NF_MGIllumina_1.0.0
```

<br>

---

### 3. Fetch Singularity Images

Although Nextflow can fetch Singularity images from a url, doing so may cause issues as detailed [here](https://github.com/nextflow-io/nextflow/issues/1210).

To avoid this issue, run the following command to fetch the Singularity images prior to running the NF_MGIllumina workflow:

> Note: This command should be run from within the `NF_MGIllumina_1.0.0` directory that was downloaded in [step 2](#2-download-the-workflow-files) above.  

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

> ***Note:** All the commands in this step must be run from within the `NF_MGIllumina_1.0.0` directory that was downloaded in [step 2](#2-download-the-workflow-files) above.*

For options and detailed help on how to run the workflow, run the following command:

```bash
nextflow run main.nf --help
```

> Note: Nextflow commands use both single hyphen arguments (e.g. -help) that denote general Nextflow 
arguments and double hyphen arguments (e.g. --input_file) that denote workflow specific parameters.  
Take care to use the proper number of hyphens for each argument.

<br>

#### 4a. Approach 1: Start with an OSD or GLDS accession as input

```bash
nextflow run main.nf -resume -profile singularity --accession OSD-574
```

<br>

#### 4b. Approach 2: Start with a runsheet csv file as input

```bash
nextflow run main.nf -resume -profile singularity  --input_file PE_file.csv
```

<br>

**Required Parameters For All Approaches:**

* `-run main.nf` - Instructs Nextflow to run the NF_MGIllumina workflow 

* `-resume` - Resumes workflow execution using previously cached results

* `-profile` – Specifies the configuration profile(s) to load (multiple options can be provided as a comma-separated list)
   * Software environment profile options (choose one):
      * `singularity` - instructs Nextflow to use Singularity container environments
      * `docker` - instructs Nextflow to use Docker container environments
      * `conda` - instructs Nextflow to use conda environments via the conda package manager. By default, Nextflow will create environments at runtime using the yaml files in the [workflow_code/envs](workflow_code/envs/) folder. You can change this behavior by using the `--conda_*` workflow parameters or by editing the [nextflow.config](workflow_code/nextflow.config) file to specify a centralized conda environments directory via the `conda.cacheDir` parameter
      * `mamba` - instructs Nextflow to use conda environments via the mamba package manager. 
   * Other option (can be combined with the software environment option above):
      * `slurm` - instructs Nextflow to use the [Slurm cluster management and job scheduling system](https://slurm.schedmd.com/overview.html) to schedule and run the jobs on a Slurm HPC cluster.

* `--accession` – A Genelab / OSD accession number e.g. OSD-574.
   > *Required only if you would like to download and process data directly from OSDR*

* `--input_file` –  A single-end or paired-end runsheet csv file containing assay metadata for each sample, including sample_id, forward, reverse, and/or paired. Please see the [runsheet documentation](./examples/runsheet) in this repository for examples on how to format this file.
   > *Required only if `--accession` is not passed as an argument*

<br> 

> See `nextflow run -h` and [Nextflow's CLI run command documentation](https://nextflow.io/docs/latest/cli.html#run) for more options and details on how to run Nextflow.  
> For additional information on editing the `nextflow.config` file, see [Step 4c](#4c-modify-parameters-and-cpu-resources-in-the-nextflow-config-file) below.   


<br>

#### 4c. Modify parameters and compute resources in the Nextflow config file

Additionally, all parameters and workflow resources can be directly specified in the [nextflow.config](./workflow_code/nextflow.config) file. For detailed instructions on how to modify and set parameters in the config file, please see the [documentation here](https://www.nextflow.io/docs/latest/config.html).

Once you've downloaded the workflow template, you can modify the parameters in the `params` scope and cpus/memory requirements in the `process` scope in your downloaded version of the [nextflow.config](workflow_code/nextflow.config) file as needed in order to match your dataset and system setup. Additionally, if necessary, you can modify each variable in the [nextflow.config](workflow_code/nextflow.config) file to be consistent with the study you want to process and the computer you're using for processing.

<br>

---

### 5. Workflow outputs

#### 5a. Main outputs

> Note: The outputs from the GeneLab Illumina metagenomics sequencing data processing pipeline workflow are documented in the [GL-DPPD-7107-A.md](../../Pipeline_GL-DPPD-7107_Versions/GL-DPPD-7107-A.md) processing protocol.

#### 5b. Resource logs

Standard Nextflow resource usage logs are also produced as follows:

**Nextflow Resource Usage Logs**
   - Output:
      - Resource_Usage/execution_report_{timestamp}.html (an html report that includes metrics about the workflow execution including computational resources and exact workflow process commands)
      - Resource_Usage/execution_timeline_{timestamp}.html (an html timeline for all processes executed in the workflow)
      - Resource_Usage/execution_trace_{timestamp}.txt (an execution tracing file that contains information about each process executed in the workflow, including: submission time, start time, completion time, cpu and memory used, machine-readable output)

> Further details about these logs can also found within [this Nextflow documentation page](https://www.nextflow.io/docs/latest/tracing.html#execution-report).

<br>

---

### 6. Post Processing

The post-processing workflow generates a README file, a protocols file, an md5sums
table, and a file association table suitable for uploading to OSDR.

For options and detailed help on how to run the post-processing workflow, run the following command:

```bash
nextflow run post_processing.nf --help
```

To generate the post-processing files after running the main processing workflow successfully, modify and set the parameters in [post_processing.config](workflow_code/post_processing.config), then run the following command:

```bash
nextflow run post_processing.nf -C post_processing.config -resume -profile singularity
``` 

The outputs of the post-processing workflow are described below:

**Post processing workflow**
   - Output:
      - Post_processing/FastQC_Outputs/filtered_multiqc_GLmetagenomics_report.zip (Filtered sequence multiqc report with paths purged)
      - Post_processing/FastQC_Outputs/raw_multiqc_GLmetagenomics_report.zip (Raw sequence multiqc report with paths purged)
      - Post_processing/<GLDS_accession>_-associated-file-names.tsv (File association table for curation)
      - Post_processing/<GLDS_accession>_metagenomics-validation.log (Automated verification and validation log file)
      - Post_processing/processed_md5sum_GLmetagenomics.tsv (md5sums for the files to be released on OSDR)
      - Post_processing/processing_info_GLmetagenomics.zip  (Zip file containing all files used to run the workflow and required logs with paths purged)
      - Post_processing/protocol.txt  (File describing the methods used by the workflow)
      - Post_processing/README_GLmetagenomics.txt (README file listing and describing the outputs of the workflow)

