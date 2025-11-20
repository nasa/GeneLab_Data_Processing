# NF_MGRemoveHostReads Workflow Information and Usage Instructions


## General workflow info
The current GeneLab Host Identification and Removal pipeline for metagenomics sequencing (MGRemoveHostReads), [GL-DPPD-7105-B.md](../../Pipeline_GL-DPPD-7105_Versions/GL-DPPD-7105-B.md), is implemented as a [Nextflow](https://www.nextflow.io/docs/stable/index.html) DSL2 workflow and utilizes [Singularity](https://docs.sylabs.io/guides/3.10/user-guide/introduction.html) containers or [conda](https://docs.conda.io/en/latest/) environments to install/run all tools. This workflow (NF_MGRemoveHostReads) is run using the command line interface (CLI) of any unix-based system. While knowledge of creating or modifying Nextflow workflows is not required to run the workflow as-is, the [Nextflow documentation](https://www.nextflow.io/docs/stable/index.html) is a useful resource for users who wish to modify and/or extend the workflow.

<br>

## Utilizing the Workflow

1. [Install Nextflow, Singularity, and Conda](#1-install-nextflow-singularity-and-conda)  
   1a. [Install Nextflow and Conda](#1a-install-nextflow-and-conda)  
   1b. [Install Singularity](#1b-install-singularity)

2. [Download the Workflow Files](#2-download-the-workflow-files)  

3. [Fetch Singularity Images](#3-fetch-singularity-images)  

4. [Run the Workflow](#4-run-the-workflow)  
   4a. [Start with a sample ID list as input](#4a-start-with-a-sample-ID-list-as-input)  
   4b. [Modify parameters and compute resources in the Nextflow config file](#4b-modify-parameters-and-compute-resources-in-the-nextflow-config-file)  

5. [Workflow Outputs](#5-workflow-outputs)  
   5a. [Main outputs](#5a-main-outputs)  
   5b. [Resource logs](#5b-resource-logs)

<br>

---

### 1. Install Nextflow, Singularity, and Conda

#### 1a. Install Nextflow and Conda

Nextflow can be installed either through the [Anaconda bioconda channel](https://anaconda.org/bioconda/nextflow) or as documented on the [Nextflow documentation page](https://www.nextflow.io/docs/latest/getstarted.html).

> Note: If you wish to install conda, we recommend installing a Miniforge version appropriate for your system, as documented on the [conda-forge website](https://conda-forge.org/download/), where you can find basic binaries for most systems. More detailed miniforge documentation is available in the [miniforge github repository](https://github.com/conda-forge/miniforge).
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

> Note: Singularity is also available through the [Anaconda conda-forge channel](https://anaconda.org/conda-forge/singularity).

> Note: Alternatively, Docker can be used in place of Singularity. To get started with Docker, see the [Docker CE installation documentation](https://docs.docker.com/engine/install/).

<br>

---

### 2. Download the Workflow Files

All files required for utilizing the NF_MGRemoveHostReads GeneLab workflow for removing host reads from metagenomics sequencing data are in the [workflow_code](workflow_code) directory. To get a copy of the latest *NF_MGRemoveHostReads* version on to your system, the code can be downloaded as a zip file from the release page then unzipped after downloading by running the following commands: 

```bash
wget <insert URL here>
unzip NF_MGRemoveHostReads_1.0.0.zip && cd NF_MGRemoveHostReads_1.0.0
```

<br>

---

### 3. Fetch Singularity Images

Although Nextflow can fetch Singularity images from a url, doing so may cause issues as detailed [here](https://github.com/nextflow-io/nextflow/issues/1210).

To avoid such issues, the required Singularity images can be manually fetched as follows before running the workflow:

```bash
mkdir -p singularity
cd singularity

# Pull required containers
singularity pull kraken2_2.1.6.img docker://quay.io/biocontainers/kraken2:2.1.6--pl5321h077b44d_0

cd ..
```

Once complete, a `singularity` folder containing the Singularity images will be created. Run the following command to export this folder as a Nextflow configuration environment variable to ensure Nextflow can locate the fetched images:

```bash
export NXF_SINGULARITY_CACHEDIR=$(pwd)/singularity
```

<br>

---

### 4. Run the Workflow

> ***Note:** All the commands in this step assume that the workflow will be run from within the `NF_MGRemoveHostReads_1.0.0` directory that was downloaded in [step 2](#2-download-the-workflow-files) above. They may also be run from a different location by providing the full path to the main.nf workflow file in the `NF_MGRemoveHostReads_1.0.0` directory.*


This workflow can be run by providing the path to a text file containing a single-column list of unique sample identifiers, an example of which is shown [here](workflow_code/unique-sample-IDs.txt), along with the path to input data (raw reads of samples). 

It also requires setting the root directory for where kraken2 reference databases are (or will be) stored. The workflow assumes databases follow the naming convention `kraken2-<host>-db`. If a database for a specified host is not found in the provided root directory, the workflow automatically builds one from scratch and saves it in the same directory using that name convention. 

In cases where the workflow is to build kraken2 database from scratch, it is important to ensure that the host organism's details are present in hosts.csv table [here](workflow_code/assets/hosts.csv). If not, they should be appended to the table in the following format: `name,species,refseq_ID,genome_build,FASTA_URL`. 

Alternatively, a pre-built database can be manually downloaded and unpacked into the root directory, provided it follows the same naming convention. An example of which is available in the [reference database info page](https://github.com/nasa/GeneLab_Data_Processing/blob/master/Metagenomics/Remove_host_reads/Workflow_Documentation/SW_MGRemoveHumanReads-A/reference-database-info.md), which describes how the human database was generated for a previous version of this workflow and how to obtain it for reuse.  

> Note: Nextflow commands use both single hyphen arguments (e.g. -profile) that denote general Nextflow arguments and double hyphen arguments (e.g. --osd) that denote workflow specific parameters.  Take care to use the proper number of hyphens for each argument.

<br>

#### 4a. Start with a sample ID list as input

```bash
nextflow run main.nf \
   -resume \
   -profile singularity \
   --ref_dbs_Dir <Path to kraken2 reference database> \
   --sample_id_list unique_sample_ids.txt \
   --reads_dir <Path to reads directory>
```

<br>


**Required Parameters:**

* `main.nf` - Instructs Nextflow to run the NF_MGRemoveHostReads workflow. If running in a directory other than `NF_MGRemoveHostReads_1.0.0`, replace with the full path to the NF_MGRemoveHostReads main.nf workflow file.
* `-resume` - Resumes  workflow execution using previously cached results
* `-profile` – Specifies the configuration profile(s) to load (multiple options can be provided as a comma-separated list)
   * Software environment profile options (choose one):
      * `singularity` - instructs Nextflow to use Singularity container environments
      * `docker` - instructs Nextflow to use Docker container environments
      * `conda` - instructs Nextflow to use conda environments via the conda package manager
      * `mamba` - instructs Nextflow to use conda environments via the mamba package manager 
   * Other option (can be combined with the software environment option above using a comma, e.g. `-profile slurm,singularity`):
      * `slurm` - instructs Nextflow to use the [Slurm cluster management and job scheduling system](https://slurm.schedmd.com/overview.html) to schedule and run the jobs on a Slurm HPC cluster
* `--ref_dbs_Dir` - Specifies the path to the directory where kraken2 databases are or will be stored
* `--sample_id_list` -  path to a single-column file with unique sample identifiers (type: string, default: null)
  > *Note: An example of this file is provided in the [workflow_code](workflow_code) directory [here](workflow_code/unique-sample-IDs.txt).*
* `--reads_dir` -  path to raw reads directory (type: string, default: null)


**Optional Parameters:**
> *Note: See `nextflow run -h` and [Nextflow's CLI run command documentation](https://nextflow.io/docs/latest/cli.html#run) for more options and details on how to run Nextflow.*


* `--is_single` - whether data is single-end  (type: boolean, default: false)  	
* `--single_suffix` - raw reads suffix that follows the unique part of sample names  (type: string, default: "_raw.fastq.gz")
* `--R1_suffix` - raw forward reads suffix that follows the unique part of sample names  (type: string, default: "_R1_raw.fastq.gz")
* `--R2_suffix` - raw reverse reads suffix that follows the unique part of sample names  (type: string, default: "_R2_raw.fastq.gz")
* `--single_out_suffix` - host-removed reads suffix that follows the unique part of sample names  (type: string, default: "_HRremoved_raw.fastq.gz")
* `--R1_out_suffix` - host-removed forward reads suffix that follows the unique part of sample names  (type: string, default: "_R1_HRremoved_raw.fastq.gz")
* `--R2_out_suffix` - host-removed reverse reads suffix that follows the unique part of sample names  (type: string, default: "_R2_HRremoved_raw.fastq.gz")
* `--host` - host organism, should match the entry provided under `name` column in [hosts.csv](workflow_code/assets/hosts.csv) (type: string, default: "human")

<br>

#### 4b. Modify parameters and compute resources in the Nextflow config file

Additionally, all parameters and workflow resources can be directly specified in the [nextflow.config](./workflow_code/nextflow.config) file. For detailed instructions on how to modify and set parameters in the config file, please see the [documentation here](https://www.nextflow.io/docs/latest/config.html).

Once you've downloaded the workflow template, you can modify the parameters in the `params` scope and cpus/memory requirements in the `process` scope in your downloaded version of the [nextflow.config](workflow_code/nextflow.config) file as needed in order to match your dataset and system setup. Additionally, if necessary, you can modify each variable in the [nextflow.config](workflow_code/nextflow.config) file to be consistent with the study you want to process and the computer you're using for processing.

<br>

---

### 5. Workflow Outputs

#### 5a. Main Outputs

The outputs from this pipeline are documented in the [GL-DPPD-7105-B](https://github.com/nasa/GeneLab_Data_Processing/blob/master/Metagenomics/Remove_host_reads/Pipeline_GL-DPPD-7105_Versions/GL-DPPD-7105-B.md) processing protocol.


The workflow also outputs the following:
  - processing_info/protocol.txt (a text file describing the methods used by the workflow)

#### 5b. Resource Logs

Standard Nextflow resource usage logs are also produced as follows:

**Nextflow Resource Usage Logs**
   - processing_info/execution_report_{timestamp}.html (an html report that includes metrics about the workflow execution including computational resources and exact workflow process commands)
   - processing_info/execution_timeline_{timestamp}.html (an html timeline for all processes executed in the workflow)
   - processing_info/execution_trace_{timestamp}.txt (an execution tracing file that contains information about each process executed in the workflow, including: submission time, start time, completion time, cpu and memory used, machine-readable output)

> Further details about these logs can also found within [this Nextflow documentation page](https://www.nextflow.io/docs/latest/tracing.html#execution-report).

<br>

---