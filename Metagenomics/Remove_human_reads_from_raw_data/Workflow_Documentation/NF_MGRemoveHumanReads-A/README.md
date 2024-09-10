# NF_MGRemoveHumanReads-A Workflow Information and Usage Instructions


## General workflow info
The current pipeline for how GeneLab identifies and removes human DNA in Illumina metagenomics sequencing data (MGRemoveHumanReads), [GL-DPPD-7105-A.md](../../Pipeline_GL-DPPD-7105_Versions/GL-DPPD-7105-A.md), is implemented as a [NextFlow](https://www.nextflow.io/docs/stable/index.html) DSL2 workflow and utilizes [Singularity](https://docs.sylabs.io/guides/3.10/user-guide/introduction.html) to run all tools in containers. This workflow (NF_MGRemoveHumanReads-A) is run using the command line interface (CLI) of any unix-based system.  While knowledge of creating or modifying Nextflow workflows is not required to run the workflow as-is, the [Nextflow documentation](https://www.nextflow.io/docs/stable/index.html) is a useful resource for users who wish to modify and/or extend the workflow.

## Utilizing the workflow

1. [Install NextFlow and Singularity](#1-install-nextflow-and-singularity)  
   1a. [Install Nextflow](#1a-install-nextflow)  
   1b. [Install Singularity](#1b-install-singularity)
2. [Download the workflow template files](#2-download-the-workflow-template-files)  
3. [Modify the variables in the nextflow.config file](#3-modify-the-variables-in-the-nextflowconfig-file)  
4. [Run the workflow](#4-run-the-workflow)


<br>

---

### 1. Install Nextflow and Singularity

#### 1a. Install Nextflow

To install NextFlow, follow the [Nextflow documentation page](https://www.nextflow.io/docs/latest/getstarted.html).
> 
> Download and install NextFlow directly:
>
> ```bash
> curl -s https://get.nextflow.io | bash
> sudo mv nextflow /usr/local/bin
> ```

Or, if Conda is installed on your system (and you prefer to use it) install and activate the “genelab-utils” Conda package including NextFlow by running the following commands:
> 
> ```bash
> conda install -n base -c conda-forge mamba
> ```
> 
> You can read a quick intro to mamba [here](https://astrobiomike.github.io/unix/conda-intro#bonus-mamba-no-5) if desired.
> 
> Once mamba is installed, you can install the genelab-utils conda package in a new environment > with the following command:
> 
> ```bash
> mamba create -n genelab-utils -c conda-forge -c bioconda -c defaults -c astrobiomike 'genelab-utils>=1.1.02'
> 
> conda activate genelab-utils
> ```
>






<br>

#### 1b. Install Singularity

Singularity is a container platform that allows usage of containerized software. This enables the workflow to retrieve and use all software required for processing without the need to install the software directly on the user's system.

We recommend installing Singularity on a system wide level as per the associated [documentation](https://docs.sylabs.io/guides/3.10/admin-guide/admin_quickstart.html).

> Note: Singularity is also available through [Anaconda](https://anaconda.org/conda-forge/singularity) if you are using Conda.

<br>

---


### 2. Download the workflow template files
All workflow files for removing human reads from metagenomics data are in the [workflow_code](workflow_code) directory. You can do this by either downloading the files for this workflow from GitHub or by [cloning](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository) the repository.


If you are using Conda with “genelab-utils”, you can copy the workflow files to your system using the “GL-get-workflow” command:
> ```bash
> GL-get-workflow NF_MGRemoveHumanReads-A
> ```
> 
> This downloaded the workflow into a directory called `NF_MGRemoveHumanReads-*/`, with the workflow version number at the end.
> 
> Note: If wanting an earlier version, the wanted version can be provided as an optional argument like so:
> ```bash
> GL-get-workflow NF_MGRemoveHumanReads-A --wanted-version 1.0.0
> ```

<br>

---

### 3. Modify the variables in the nextflow.config file
Once you've downloaded the workflow template, you can modify the variables in the [nextflow.config](workflow_code/nextflow.config) file as needed. You will also need to indicate the path to your input data (raw reads) and the root directory for where the kraken2 reference database should be stored (it will be set up automatically). Additionally, if necessary, you'll need to modify each variable in the [nextflow.config](workflow_code/nextflow.config) file to be consistent with the study you want to process and the machine you're using. 
Confirm the following variables are appropriate for your data:
- DL_kraken
- single_end
- specify_reads
- Sample_id_list
- reads_dir
- PE_reads_suffix or SE_reads_suffix
- PE_reads_out_suffix or SE_reads_out_suffix
- kraken_output_dir
- human_db_path


If you want to specify certain read files you will have to provide a text file containing a single-column list of unique sample identifiers (see an example of how to set this up below - if you are running the example dataset, this file is provided in the [workflow_code](workflow_code) directory [here](workflow_code/unique-sample-IDs.txt)).
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

### 4. Run the workflow

While in the directory holding the NextFlow file, nextflow.config file, and other workflow files that you downloaded in [step 2](#2-download-the-workflow-template-files), here is one example command of how to run the workflow:

```bash
nextflow run *path/to/Remove_Human_Reads.nf* -resume -ansi-log false --DL_kraken false
```


* `-resume` – continues to run the workflow using cached data from the previous run
* `-ansi-log false` – specifies to print out each command being run to the screen instead of dynamically updating the log
* `--specify_reads false` - processes all reads in the working directory, without requiring a sample ID list
* `--single_end true` – indicates reads are single-ended
* `--DL_kraken true` – runs a process before the rest of the workflow to download and install the necessary database.
>
> Note - Nextflow options use a single dash prefix, e.g. -resume, whereas pipeline parameters use double dash notation, e.g. --specify_reads. All of the pipeline parameters can and should be set from the nextflow.config file to avoid typos or errors.
>

See `nextflow -h` and [NextFlow's documentation](https://www.nextflow.io/docs/master/index.html) for more options and details.


---

## Reference database info
The database we use was built with kraken2 v2.1.1 as detailed below, and can be downloaded to run with this workflow (it's ~4.3 GB uncompressed). The steps for building it are described on the [reference database info page](https://github.com/nasa/GeneLab_Data_Processing/blob/master/Metagenomics/Remove_human_reads_from_raw_data/Workflow_Documentation/SW_MGRemoveHumanReads-A/reference-database-info.md).

---

