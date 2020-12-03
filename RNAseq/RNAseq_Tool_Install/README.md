# Instructions for installing the tools required to process RNAseq data using the GeneLab pipeline

> **This directory holds yaml files and instructions for how to install two conda environments, RNAseq_fq_to_counts_tools and RNAseq_R_tools, that contain all the software programs needed to process RNAseq data using the GeneLab pipeline.**  

---

### Install prerequisites
  * Anaconda  
    To install conda environments, you'll first have to install [Anaconda](https://www.anaconda.com/). CLick [here](https://docs.anaconda.com/anaconda/install/) for installation instructions.

### Install the **RNAseq_fq_to_counts_tools** conda environment by running the following command:

  ```
  conda env create -f RNAseq_fq_to_counts_tools.yml
  ```

  Activate the RNAseq_fq_to_counts_tools conda environment with the following command:
  
  ```
  conda activate RNAseq_fq_to_counts
  ``` 
