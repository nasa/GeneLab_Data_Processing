# SW_AmpIllumina-A Workflow Information and Usage Instructions <!-- omit in toc -->


## General workflow info <!-- omit in toc -->
The current GeneLab Illumina amplicon sequencing data processing pipeline (AmpIllumina), [GL-DPPD-7104-A.md](../../Pipeline_GL-DPPD-7104_Versions/GL-DPPD-7104-A.md), is implemented as a [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow and utilizes [conda](https://docs.conda.io/en/latest/) environments to install/run all tools. This workflow (SW_AmpIllumina-A) is run using the command line interface (CLI) of any unix-based system. The workflow can be used even if you are unfamiliar with Snakemake and conda, but if you want to learn more about those, [this Snakemake tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html) within [Snakemake's documentation](https://snakemake.readthedocs.io/en/stable/) is a good place to start for that, and an introduction to conda with installation help and links to other resources can be found [here at Happy Belly Bioinformatics](https://astrobiomike.github.io/unix/conda-intro).  

<br>

---

## Utilizing the workflow <!-- omit in toc -->

1. [Install conda, mamba, and `genelab-utils` package](#1-install-conda-mamba-and-genelab-utils-package)
2. [Download the workflow template files](#2-download-the-workflow-template-files)
3. [Run the workflow using `run_workflow.py`](#3-run-the-workflow-using-run_workflowpy)   
   3a. [Approach 1: Run the workflow on a GeneLab Amplicon (Illumina) sequencing dataset with automatic retrieval of raw read files and metadata](#3a-approach-1-run-the-workflow-on-a-genelab-amplicon-illumina-sequencing-dataset-with-automatic-retrieval-of-raw-read-files-and-metadata)   
   3b. [Approach 2: Run the workflow on a non-OSD dataset using a user-created runsheet](#3b-approach-2-run-the-workflow-on-a-non-osd-dataset-using-a-user-created-runsheet)   
4. [Additional output files](#additional-output-files)

<br>

___

### 1. Install conda, mamba, `genelab-utils`, and `dp-tools` package
We recommend installing a Miniconda, Python3 version appropriate for your system, as exemplified in [the above link](https://astrobiomike.github.io/unix/conda-intro#getting-and-installing-conda).  

Once conda is installed on your system, we recommend installing [mamba](https://github.com/mamba-org/mamba#mamba), as it generally allows for much faster conda installations:

```bash
conda install -n base -c conda-forge mamba
```

> You can read a quick intro to mamba [here](https://astrobiomike.github.io/unix/conda-intro#bonus-mamba-no-5) if wanted.

Once mamba is installed, you can install the genelab-utils conda package in a new environment with the following command:

```bash
mamba create -n genelab-utils -c conda-forge -c bioconda -c defaults -c astrobiomike 'genelab-utils>=1.1.02' git pip
```

The environment then needs to be activated by running the following command:

```bash
conda activate genelab-utils
```

Once the `genelab-utils` environment is active, use `pip` to install the `dp-tools` package by running the following command:

```bash
pip install git+https://github.com/torres-alexis/dp_tools.git@amplicon_updates
```

<br>

___

### 2. Download the workflow template files
All files required for utilizing the GeneLab workflow for processing Illumina amplicon sequencing data are in the [workflow_code](workflow_code) directory. To get a copy of latest SW_AmpIllumina-A version on to your system, the code can be downloaded as a zip file from the release page then unzipped after downloading by running the following commands:

```bash
## Change to NASA GitHub when ready 
wget https://github.com/asaravia-butler/GeneLab_Data_Processing_TESTING/releases/download/SW_AmpIllumina-A_1.2.0/SW_AmpIllumina-A_1.2.0.zip
#wget https://github.com/nasa/GeneLab_Data_Processing/releases/download/SW_AmpIllumina-A_1.2.0/SW_AmpIllumina-A_1.2.0.zip

unzip SW_AmpIllumina-A_1.2.0.zip
```

This downloaded the workflow into a directory called `SW_AmpIllumina-A_1.2.0`. To run the workflow, you will need to move into that directory by running the following command:

```bash
cd SW_AmpIllumina-A_1.2.0
```

<br>

___

### 3. Run the workflow using `run_workflow.py`

While in the `SW_AmpIllumina-A_1.2.0` directory that was downloaded in [step 2](#2-download-the-workflow-template-files), you are now able to run the workflow using the `run_workflow.py` script in the [scripts/](workflow_code/scripts) sub-directory to set up the configuration files needed to execute the workflow.

> Note: The commands to run the workflow in each approach listed below allows for two sets of options. The options specified outside of the quotation marks are specific to the `run_workflow.py` script, and the options specified within the quotation marks are specific to `snakemake`.

<br>

___

#### 3a. Approach 1: Run the workflow on a GeneLab Amplicon (Illumina) sequencing dataset with automatic retrieval of raw read files and metadata

> This approach processes data hosted on the [NASA Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/). Upon execution, the command downloads then parses the OSD ISA.zip file to create a runsheet containing link(s) to the raw reads and the metadata required for processing. The runsheet is then used to prepare the necessary configuration files before executing the workflow using the specified Snakemake run command.

```bash
python ./scripts/run_workflow.py --OSD OSD-487 --run "snakemake --use-conda --conda-prefix ${CONDA_PREFIX}/envs -j 2 -p"
```

<br>

___

#### 3b. Approach 2: Run the workflow on a non-OSD dataset using a user-created runsheet

> If processing a non-OSD dataset, you must manually create the runsheet for your dataset to run the workflow. Specifications for creating a runsheet manually are described [here](examples/runsheet/README.md).

```bash
python ./scripts/run_workflow.py --runsheetPath </path/to/runsheet> --run "snakemake --use-conda --conda-prefix ${CONDA_PREFIX}/envs -j 2 -p"
```

**Parameter Definitions for `run_workflow.py`:** 

* `--OSD OSD-###` - specifies the OSD dataset to process through the SW_AmpIllumina workflow (replace ### with the OSD number)
   > *Used for Approach 1 only.*

* `--target` - specifies the genomic target for the assay. Options: 16S, 18S, ITS. This determines which reference database is used for taxonomic classification, and it is used to select the appropriate dataset from an OSD study when multiple options are available.

* `--runsheetPath` - specifies the path to a local runsheet containing the metadata and raw reads location (as a link or local path), used for processing a non-OSD dataset through the SW_AmpIllumina workflow 
   > *Optionally used for Approach 2 only, the form can be used instead of providing a runsheet on the NASA EDGE platform.*

* `--run` - specifies the command used to execute the snakemake workflow; snakemake-specific parameters are defined below

* `--outputDir` - specifies the output directory for the output files generated by the workflow 
   > *This is an optional command that can be added outside the quotation marks in either approach to specify the output directory. If this option is not used, the output files will be printed to the current working directory, i.e. in the `SW_AmpIllumina-A_1.2.0` directory that was downloaded in [step 2](#2-download-the-workflow-template-files).*

* `--trim-primers TRUE/FALSE` - specifies to trim primers (TRUE) or not (FALSE). Default: TRUE
   > *Note: Primers should virtually always be trimmed from amplicon datasets. This option is here for cases where they have already been removed.*

* `--min_trimmed_length` - specifies the minimum length of trimmed reads during cutadapt filtering. For paired-end data: if one read gets filtered, both reads are discarded. Default: 130
   > *Note: For paired-end data, all filtering and trimming should leave a minimum of an 8-base overlap of forward and reverse reads.*

* `--primers-linked TRUE/FALSE` - if set to TRUE, instructs cutadapt to treat the primers as linked. Default: TRUE
   > *Note: See [cutadapt documentation here](https://cutadapt.readthedocs.io/en/stable/recipes.html#trimming-amplicon-primers-from-paired-end-reads) for more info.*

* `--anchor-primers TRUE/FALSE` - indicates if primers should be anchored (TRUE) or not (FALSE) when provided to cutadapt. Default: TRUE
   > *Note: See [cutadapt documention here](https://cutadapt.readthedocs.io/en/stable/guide.html#anchored-5adapters) for more info.*

* `--discard-untrimmed TRUE/FALSE` - if set to TRUE, instructs cutadapt to remove reads if the primers were not found in the expected location; if set to FALSE, these reads are kept. Default: TRUE

* `--left-trunc` - dada2 parameter that specifies to truncate the forward reads to this length, bases beyond this length will be removed and reads shorter than this length are discarded. Default: 0 (no truncation)
   > *Note: See dada2 [filterAndTrim documentation](https://rdrr.io/bioc/dada2/man/filterAndTrim.html) for more info.*

* `--right-trunc` - dada2 parameter that specifies to truncate the reverse reads, bases beyond this length will be truncated and reads shorter than this length are discarded. Default: 0 (no truncation)
   > *Note: See dada2 [filterAndTrim documentation](https://rdrr.io/bioc/dada2/man/filterAndTrim.html) for more info.*

* `--left-maxEE` - dada2 parameter that specifies the maximum expected error (maxEE) allowed for each forward read, reads with a higher maxEE than provided will be discarded. Default: 1 
   > *Note: See dada2 [filterAndTrim documentation](https://rdrr.io/bioc/dada2/man/filterAndTrim.html) for more info.*

* `--right-maxEE` - dada2 parameter that specifies the maximum expected error (maxEE) allowed for each forward read, reads with a higher maxEE than provided will be discarded. Default: 1 
   > *Note: See dada2 [filterAndTrim documentation](https://rdrr.io/bioc/dada2/man/filterAndTrim.html) for more info.*

* `--concatenate_reads_only TRUE/FALSE` - if set to TRUE, specifies to concatenate forward and reverse reads only with dada2 instead of merging paired reads. Default: FALSE

* `--output-prefix ""` - specifies the prefix to use on all output files to distinguish multiple primer sets, leave as an empty string if only one primer set is being processed (if used, be sure to include a connecting symbol, e.g. "ITS-"). Default: ""

* `--specify-runsheet` - specifies the runsheet to use when multiple runsheets are generated.
   > *Optional argument used in Approach 1 for datasets that have multiple assays for the same amplicon target (e.g. [OSD-249](https://osdr.nasa.gov/bio/repo/data/studies/OSD-249)).*


**Parameter Definitions for `snakemake`**

* `--use-conda` – specifies to use the conda environments included in the workflow (these are specified in the [envs](workflow_code/envs) directory)
* `--conda-prefix` – indicates where the needed conda environments will be stored. Adding this option will also allow the same conda environments to be re-used when processing additional datasets, rather than making new environments each time you run the workflow. The value listed for this option, `${CONDA_PREFIX}/envs`, points to the default location for conda environments (note: the variable `${CONDA_PREFIX}` will be expanded to the appropriate location on whichever system it is run on).
* `-j` – assigns the number of jobs Snakemake should run concurrently
* `-p` – specifies to print out each command being run to the screen

See `snakemake -h` and [Snakemake's documentation](https://snakemake.readthedocs.io/en/stable/) for more options and details.

<br>

___

### Additional Output Files

The outputs from the `run_workflow.py` and differential abundance analysis (DAA) / visualizations scripts are described below:
> Note: Outputs from the Amplicon Seq - Illumina pipeline are documented in the [GL-DPPD-7104-A.md](../../Pipeline_GL-DPPD-7104_Versions/GL-DPPD-7104-A.md) processing protocol.

- **Metadata Outputs:**
  - \*_AmpSeq_v1_runsheet.csv (table containing metadata required for processing, including the raw reads files location)
  - \*-ISA.zip (the ISA archive of the OSD datasets to be processed, downloaded from the OSDR)
  - config.yaml (configuration file containing the metadata from the runsheet (\*_AmpSeq_v1_runsheet.csv), required for running the SW_AmpIllumina workflow) 
  - unique-sample-IDs.txt (text file containing the IDs of each sample used, required for running the SW_AmpIllumina workflow)
- **DAA and Visualization Outputs:**  
  - dendrogram_by_group.png (Dendrogram of euclidean distance - based hierarchical clustering of the samples, colored by experimental groups) 
  - PCoA_w_labels.png (Principle Coordinates Analysis plot of VST transformed ASV counts, with sample labels)
  - PCoA_without_labels.png (Principle Coordinates Analysis plot of VST transformed ASV counts, without labels)
  - Rarefaction.png (Rarefaction plot visualizing species richness for each sample)
  - richness_by_sample.png (Chao1 richness estimates and Shannon diversity estimates for each sample)
  - richness_by_group.png (Chao1 richness estimates and Shannon diversity estimates for each group)
  - relative_classes.png (Bar plot taxonomic summary of proportions of phyla identified in each group, by class)
  - relative_phyla.png (Bar plot taxonomic summary of proportions of phyla identified in each group, by phyla)
  - {group1}\_vs_{group2}.csv (Differential abundance tables for all pairwise contrasts of groups)
  - volcano\_{group1}\_vs_{group2}.png (Volcano plots for all pairwise contrasts of groups)
