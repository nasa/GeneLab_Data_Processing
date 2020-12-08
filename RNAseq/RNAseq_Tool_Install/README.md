# Instructions for installing the tools required to process RNAseq data using the GeneLab pipeline

> **This directory holds yaml files and instructions for how to install two conda environments, RNAseq_fq_to_counts_tools and RNAseq_R_tools, that contain all the software programs needed to process RNAseq data using the GeneLab pipeline.**  

---

### Install prerequisites

  * **Anaconda**  
    To install conda environments, you'll first have to install [Anaconda](https://www.anaconda.com/). Click [here](https://docs.anaconda.com/anaconda/install/) for installation instructions.


### Install the **RNAseq_fq_to_counts_tools** conda environment by running the following command:

  ```
  conda env create -f RNAseq_fq_to_counts_tools.yml
  ```

  Activate the RNAseq_fq_to_counts_tools conda environment with the following command:
  > This environment needs to be activated to run steps 1-6 of the [RNAseq processing pipeline](https://github.com/nasa/GeneLab_Data_Processing/blob/master/RNAseq/GL-DPPD-7101-C.md)
  
  ```
  conda activate RNAseq_fq_to_counts
  ``` 


### Install the **RNAseq_R_tools** conda environment by running the following command:

  ```
  conda env create -f RNAseq_R_tools.yml
  ```

  Activate the RNAseq_R_tools conda environment with the following command:
  > This environment needs to be activated to run step 7 of the [RNAseq processing pipeline](https://github.com/nasa/GeneLab_Data_Processing/blob/master/RNAseq/GL-DPPD-7101-C.md)
  
  ```
  conda activate RNAseq_R_tools
  ``` 


### Create RNAseq output directory structure

  To generate the same directory structure that GeneLab uses to organize the RNAseq processing pipeline outputs follow the instructions below:
  1. Use the `cd` command to navigate to the location on your device where you want to output the RNAseq processing pipeline.
  
