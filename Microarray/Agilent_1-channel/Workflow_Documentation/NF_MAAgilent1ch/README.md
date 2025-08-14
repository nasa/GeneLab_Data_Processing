# NF_MAAgilent1ch Workflow Information and Usage Instructions <!-- omit in toc -->

## General Workflow Info <!-- omit in toc -->

### Implementation Tools <!-- omit in toc -->

The current GeneLab Agilent 1 Channel Microarray consensus processing pipeline (NF_MAAgilent1ch), [GL-DPPD-7112](../../Pipeline_GL-DPPD-7112_Versions/GL-DPPD-7112.md), is implemented as a [Nextflow](https://nextflow.io/) DSL2 workflow and utilizes [Singularity](https://docs.sylabs.io/guides/3.10/user-guide/introduction.html) to run all tools in containers. This workflow (NF_MAAgilent1ch) is run using the command line interface (CLI) of any unix-based system.  While knowledge of creating workflows in Nextflow is not required to run the workflow as is, [the Nextflow documentation](https://nextflow.io/docs/latest/index.html) is a useful resource for users who want to modify and/or extend this workflow.   

### Workflow & Subworkflows <!-- omit in toc -->

---

<!-- TODO Add workflow image -->

---
The NF_MAAgilent1ch workflow is composed of three subworkflows as shown in the image above.
Below is a description of each subworkflow and the additional output files generated that are not already indicated in the [GL-DPPD-7112 pipeline document](../../Pipeline_GL-DPPD-7112_Versions/GL-DPPD-7112.md):

1. **Analysis Staging Subworkflow**

   - Description:
     - This subworkflow extracts the metadata parameters (e.g. organism, Array Design REF) needed for processing from the OSD/GLDS ISA archive and retrieves the raw reads files hosted on the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/).
       > *OSD/GLDS ISA archive*: ISA directory containing Investigation, Study, and Assay (ISA) metadata files for a respective GLDS dataset - the *ISA.zip file is located in the [OSDR](https://osdr.nasa.gov/bio/repo/) under 'Files' -> 'Study Metadata Files' for any GeneLab Data Set (GLDS) in the OSDR.

2. **Agilent 1 Channel Microarray Processing Subworkflow**

   - Description:
     - This subworkflow uses the staged raw data and metadata parameters from the Analysis Staging Subworkflow to generate processed data using the [GL-DPPD-7112 pipeline](../../Pipeline_GL-DPPD-7112_Versions/GL-DPPD-7112.md).

1. **V&V Pipeline Subworkflow**

   - Description:
     - This subworkflow performs validation and verification (V&V) on the raw and processed data files.  It performs a series of checks on the output files generated and flags the results, using the flag codes indicated in the table below, which are outputted into a log file.  
       **V&V Flags**:

       |Flag Codes|Flag Name|Interpretation|
       |:---------|:--------|:-------------|
       | 2    | MANUAL   | Special flag that indicates a manual check that is advised. Often used to advise what assess in QA plots. |
       | 20    | GREEN   | Indicates the check passed all validation conditions |
       | 30    | YELLOW  | Indicates the check was flagged for minor issues (e.g. slight outliers) |
       | 50    | RED     | Indicates the check was flagged for moderate issues (e.g. major outliers) |
       | 80    | HALT    | Indicates the check was flagged for severe issues that trigger a processing halt (e.g. missing data) |

<br>

---
## Utilizing the Workflow <!-- omit in toc -->

- [1. Install Nextflow and Singularity](#1-install-nextflow-and-singularity)
  - [1a. Install Nextflow](#1a-install-nextflow)
  - [1b. Install Singularity](#1b-install-singularity)
- [2. Download the Workflow Files](#2-download-the-workflow-files)
- [3. Run the Workflow](#3-run-the-workflow)
  - [3a. Approach 1: Run the workflow on a GeneLab Agilent 1 Channel Microarray dataset](#3a-approach-1-run-the-workflow-on-a-genelab-agilent-1-channel-microarray-dataset)
  - [3b. Approach 2: Run the workflow on a non-GLDS dataset using a user-created runsheet](#3b-approach-2-run-the-workflow-on-a-non-glds-dataset-using-a-user-created-runsheet)
- [4. Additional Output Files](#4-additional-output-files)

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

Singularity is a container platform that allows usage of containerized software. This enables the GeneLab NF_MAAgilent1ch workflow to retrieve and use all software required for processing without the need to install the software directly on the user's system.

We recommend installing Singularity on a system wide level as per the associated [documentation](https://docs.sylabs.io/guides/3.10/admin-guide/admin_quickstart.html).

> Note: Singularity is also available through [Anaconda](https://anaconda.org/conda-forge/singularity).

<br>

---

### 2. Download the Workflow Files

All files required for utilizing the NF_MAAgilent1ch GeneLab workflow for processing Agilent 1 Channel Microarray data are in the [workflow_code](workflow_code) directory. To get a copy of latest NF_MAAgilent1ch version on to your system, the code can be downloaded as a zip file from the release page then unzipped after downloading by running the following commands: 

```bash
wget https://github.com/nasa/GeneLab_Data_Processing/releases/download/NF_MAAgilent1ch_1.0.4/NF_MAAgilent1ch_1.0.4.zip

unzip NF_MAAgilent1ch_1.0.4.zip
```

<br>

---

### 3. Run the Workflow

While in the location containing the `NF_MAAgilent1ch_1.0.4` directory that was downloaded in [step 2](#2-download-the-workflow-files), you are now able to run the workflow. Below are three examples of how to run the NF_MAAgilent1ch workflow:
> Note: Nextflow commands use both single hyphen arguments (e.g. -help) that denote general nextflow arguments and double hyphen arguments (e.g. --ensemblVersion) that denote workflow specific parameters.  Take care to use the proper number of hyphens for each argument.

<br>

#### 3a. Approach 1: Run the workflow on a GeneLab Agilent 1 Channel Microarray dataset

```bash
nextflow run NF_MAAgilent1ch_1.0.4/main.nf \ 
   -profile singularity \
   --osdAccession OSD-548 \
   --gldsAccession GLDS-548 
```

<br>

#### 3b. Approach 2: Run the workflow on a non-GLDS dataset using a user-created runsheet

> Note: Specifications for creating a runsheet manually are described [here](examples/runsheet/README.md).

```bash
nextflow run NF_MAAgilent1ch_1.0.4/main.nf \ 
   -profile singularity \
   --runsheetPath </path/to/runsheet> 
```

<br>

**Required Parameters For All Approaches:**

* `NF_MAAgilent1ch_1.0.4/main.nf` - Instructs Nextflow to run the NF_MAAgilent1ch workflow 

* `-profile` - Specifies the configuration profile(s) to load, `singularity` instructs Nextflow to setup and use singularity for all software called in the workflow


<br>

**Additional Required Parameters For [Approach 1](#3a-approach-1-run-the-workflow-on-a-genelab-agilent-1-channel-microarray-dataset):**

* `--osdAccession OSD-###` – specifies the OSD ID to process through the NF_MAAgilent1ch workflow (replace ### with the OSD number)

* `--gldsAccession GLDS-###` – specifies the GLDS ID to process through the NF_MAAgilent1ch workflow (replace ### with the GLDS number)  

<br>

**Additional Required Parameters For [Approach 2](#3b-approach-2-run-the-workflow-on-a-non-glds-dataset-using-a-user-created-runsheet):**

* `--runsheetPath` - specifies the path to a local runsheet (Default: a runsheet is automatically generated using the metadata on the GeneLab Repository for the GLDS dataset being processed) 

<br>

**Optional Parameters:**

* `--skipVV` - skip the automated V&V processes (Default: the automated V&V processes are active) 

* `--outputDir` - specifies the directory to save the raw and processed data files (Default: files are saved in the launch directory)  

<br>

All parameters listed above and additional optional arguments for the NF_MAAgilent1ch workflow, including debug related options that may not be immediately useful for most users, can be viewed by running the following command:

```bash
nextflow run NF_MAAgilent1ch_1.0.4/main.nf --help
```

See `nextflow run -h` and [Nextflow's CLI run command documentation](https://nextflow.io/docs/latest/cli.html#run) for more options and details common to all nextflow workflows.

<br>

---

### 4. Additional Output Files

All R code steps and output are rendered within a Quarto document yielding the following:

   - Output:
     - NF_MAAgilent1ch_1.0.4.html (html report containing executed code and output including QA plots)
  

The outputs from the Analysis Staging and V&V Pipeline Subworkflows are described below:
> Note: The outputs from the Agilent 1 Channel Microarray Processing Subworkflow are documented in the [GL-DPPD-7112.md](../../../Pipeline_GL-DPPD-7112_Versions/GL-DPPD-7112.md) processing protocol.

**Analysis Staging Subworkflow**

   - Output:
     - \*_microarray_v1_runsheet.csv (table containing metadata required for processing, including the raw reads files location)
     - \*-ISA.zip (the ISA archive of the GLDS datasets to be processed, downloaded from the GeneLab Data Repository)
   
   
**V&V Pipeline Subworkflow**

   - Output:
     - VV_log_VV_AGILE1CH.tsv.MANUAL_CHECKS_PENDING (table containing V&V flags for all checks performed. Also contains rows indicating suggested manual checks focusing on QA plots embedded in the html report)

<br>

Standard Nextflow resource usage logs are also produced as follows:
> Further details about these logs can also found within [this Nextflow documentation page](https://www.nextflow.io/docs/latest/tracing.html#execution-report).

**Nextflow Resource Usage Logs**
   - Output:
     - Resource_Usage/execution_report_{timestamp}.html (an html report that includes metrics about the workflow execution including computational resources and exact workflow process commands)
     - Resource_Usage/execution_timeline_{timestamp}.html (an html timeline for all processes executed in the workflow)
     - Resource_Usage/execution_trace_{timestamp}.txt (an execution tracing file that contains information about each process executed in the workflow, including: submission time, start time, completion time, cpu and memory used, machine-readable output)

<br>
