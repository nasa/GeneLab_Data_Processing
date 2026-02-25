# Instructions for installing the tools required to process small RNAseq data using the GeneLab pipeline

> **This directory holds instructions for installing conda environments for the mirdeep2 and miRDP2 tools.**  

---

### Install prerequisites

  * **Anaconda**  
    To install conda environments, you'll first have to install [Anaconda](https://www.anaconda.com/). Click [here](https://docs.anaconda.com/anaconda/install/) for installation instructions.

<br>

### Install the **mirdeep2** conda environment by running the following command:

  ```
  conda create -n mirdeep2-2.0.1.3 mirdeep2=2.0.1.3 -c bioconda
  ```

  Activate the mirdeep2 conda environment with the following command:
  > This environment needs to be activated to run steps 1-6 of the [small RNAseq processing pipeline](../Pipeline_GL-DPPD-7101_Versions/GL-DPPD-7101-E.md)
  
  ```
  conda activate mirdeep2-2.0.1.3
  ``` 

<br>

### Install the **miRDP2** conda environment by running the following command:

  ```
  conda env create -n mirdp2-1.1.4 mirdeep-p2=1.1.4 -c bioconda
  ```

  Activate the RNAseq_R_tools conda environment with the following command:
  > This environment needs to be activated to run step 7 of the [RNAseq processing pipeline](../Pipeline_GL-DPPD-7101_Versions/GL-DPPD-7101-E.md)
  
  ```
  conda activate mirdp2-1.1.4
  ``` 
  Note: The tools in the RNAseq_R_tools conda environment can be run on a standard laptop.

  Once the environment is setup, install the rfam reference files. miRDP2 expects these files to exist within the install location.

  1. Download the ncRNA_rfam.tar.gz reference file [https://sourceforge.net/projects/mirdp2/files/latest_version/](https://sourceforge.net/projects/mirdp2/files/latest_version/)
  2. Decompress the downloaded file and build the bowtie index (note that the bowtie indexes MUST be placed in the miRDP2 install location and named as indicated)

  ```bash
  conda activate mirdp2-1.1.4
  tar -zxf ncRNA_rfam.tar.gz
  bowtie-build -f ${CONDA_PREFIX}/bin/scripts/index/rfam_index
  ```

<br>

### Troubleshooting conda environment installation

  If you run into any issues while installing the RNAseq conda environments, update conda using the following command then try re-creating the RNAseq conda environments using the installation instructions above.
  ```
  conda update conda
  ```

<br>