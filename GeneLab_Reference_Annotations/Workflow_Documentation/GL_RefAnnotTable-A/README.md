# GL_RefAnnotTable-A Workflow Information and Usage Instructions <!-- omit in toc -->

## General Workflow Info <!-- omit in toc -->

### Implementation Tools <!-- omit in toc -->

The current GeneLab Reference Annotation Table (GL_RefAnnotTable-A) pipeline is implemented as an R workflow and utilizes [Singularity](https://docs.sylabs.io/guides/3.10/user-guide/introduction.html) to run all tools in a containerized environment. This workflow is run using the command line interface (CLI) of any unix-based system.

<br>

---
## Utilizing the Workflow

1. [Install Singularity](#1-install-singularity)
2. [Download the Workflow Files](#2-download-the-workflow-files)  
3. [Fetch Singularity Images](#3-fetch-singularity-images)  
4. [Run the Workflow](#4-run-the-workflow)  
5. [Run the annotations database creation function as a stand-alone script](#5-run-the-annotations-database-creation-function-as-a-stand-alone-script)

<br>

---

### 1. Install Singularity

Singularity is a container platform that allows usage of containerized software. This enables the GL_RefAnnotTable-A workflow to retrieve and use all software required for processing without the need to install the software directly on the user's system.

We recommend installing Singularity on a system wide level as per the associated [documentation](https://docs.sylabs.io/guides/3.10/admin-guide/admin_quickstart.html).

> Note: Singularity is also available through [Anaconda](https://anaconda.org/conda-forge/singularity).

<br>

---

### 2. Download the Workflow Files

All files required for utilizing the GL_RefAnnotTable-A workflow for generating reference annotation tables are in the [workflow_code](workflow_code) directory. To get a copy of latest GL_RefAnnotTable-A version on to your system, run the following commands:

```bash
curl -LO https://github.com/nasa/GeneLab_Data_Processing/releases/download/GL_RefAnnotTable-A_1.1.0/GL_RefAnnotTable-A_1.1.0.zip
unzip GL_RefAnnotTable-A_1.1.0.zip
```

<br>

---

### 3. Fetch Singularity Images

Although Singularity can fetch images from a url, doing so may cause issues as detailed [here](https://github.com/nextflow-io/nextflow/issues/1210).

To avoid this issue, run the following command to fetch the Singularity images prior to running the GL_RefAnnotTable-A workflow:
> Note: This command should be run in the location containing the `GL_RefAnnotTable-A_1.1.0` directory that was downloaded in [step 2](#2-download-the-workflow-files) above. Depending on your network speed, fetching the images will take ~20 minutes.  

```bash
bash GL_RefAnnotTable-A_1.1.0/bin/prepull_singularity.sh GL_RefAnnotTable-A_1.1.0/config/software/by_docker_image.config
```

Once complete, a `singularity` folder containing the Singularity images will be created. Run the following command to export this folder as a Singularity configuration environment variable:

```bash
export SINGULARITY_CACHEDIR=$(pwd)/singularity
```

<br>

---

### 4. Run the Workflow

While in the location containing the `GL_RefAnnotTable-A_1.1.0` directory that was downloaded in [step 2](#2-download-the-workflow-files), you are now able to run the workflow. Below is an example of how to run the workflow to build an annotation table for Mus musculus (mouse):

```bash
singularity exec -B $(pwd)/GL_RefAnnotTable-A_1.1.0:/work \
  $SINGULARITY_CACHEDIR/gl-refannottable_v1.0.0.sif \
  bash -c "cd /work && Rscript GL-DPPD-7110-A_build-genome-annots-tab.R 'Mus musculus'"
```

**Input data:**

- No input files are required. Specify the target organism using a positional command line argument. `Mus musculus` is used in the example above. To see a list of all available organisms, run the command without positional arguments. The correct argument for each organism can also be found in the 'species' column of the [GL-DPPD-7110-A_annotations.csv](../../Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv)

- Optional: a reference table CSV can be supplied as a second positional argument instead of using the default [GL-DPPD-7110-A_annotations.csv](../../Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv)

**Output data:**

- *-GL-annotations.tsv (Tab delineated table of gene annotations)
- *-GL-build-info.txt (Text file containing information used to create the annotation table, including tool and tool versions and date of creation)

<br>

---

### 5. Run the annotations database creation function as a stand-alone script

When the workflow is run, if the reference table does not specify an annotations database for the target_organism in the `annotations` column, the `install_annotations` function, defined in the `install-org-db.R` script, will be executed. This script will locally create and install an annotations database R package using AnnotationForge. This function can also be run as a stand-alone script from the command line:

```bash
singularity exec -B $(pwd)/GL_RefAnnotTable-A_1.1.0:/work \
  $SINGULARITY_CACHEDIR/gl-refannottable_v1.0.0.sif \
  bash -c "cd /work && Rscript install-org-db.R 'Bacillus subtilis' /path/to/GL-DPPD-7110-A_annotations.csv"
```

**Input data:**

- The target organism must be specified as the first positional command line argument, `Bacillus subtilis` is used in the example above. The correct argument for each organism can be found in the 'species' column of the [GL-DPPD-7110-A_annotations.csv](../../Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv)

- The path to a local reference table must also be supplied as the second positional argument

**Output data:**

- org.*.eg.db/ (species-specific annotation database, as a local R package)

<br>