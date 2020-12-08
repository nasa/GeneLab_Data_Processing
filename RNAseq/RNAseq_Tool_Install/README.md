# Instructions for installing the tools required to process RNAseq data using the GeneLab pipeline

> **This directory holds yaml files and instructions for how to install two conda environments, RNAseq_fq_to_counts_tools and RNAseq_R_tools, that contain all the software programs needed to process RNAseq data using the GeneLab pipeline. Additionally, instructions and the shell script for creating the directory structure that GeneLab uses to organize the RNAseq processing pipeline output files are provided.**  

---

### Install prerequisites

  * **Anaconda**  
    To install conda environments, you'll first have to install [Anaconda](https://www.anaconda.com/). Click [here](https://docs.anaconda.com/anaconda/install/) for installation instructions.

<br>

### Install the **RNAseq_fq_to_counts_tools** conda environment by running the following command:

  ```
  conda env create -f RNAseq_fq_to_counts_tools.yml
  ```

  Activate the RNAseq_fq_to_counts_tools conda environment with the following command:
  > This environment needs to be activated to run steps 1-6 of the [RNAseq processing pipeline](https://github.com/nasa/GeneLab_Data_Processing/blob/master/RNAseq/GL-DPPD-7101-C.md)
  
  ```
  conda activate RNAseq_fq_to_counts
  ``` 
  Note: At least 45GB of RAM is required to run the tools in the RNAseq_fq_to_counts_tools conda environment.

<br>

### Install the **RNAseq_R_tools** conda environment by running the following command:

  ```
  conda env create -f RNAseq_R_tools.yml
  ```

  Activate the RNAseq_R_tools conda environment with the following command:
  > This environment needs to be activated to run step 7 of the [RNAseq processing pipeline](https://github.com/nasa/GeneLab_Data_Processing/blob/master/RNAseq/GL-DPPD-7101-C.md)
  
  ```
  conda activate RNAseq_R_tools
  ``` 
  Note: The tools in the RNAseq_R_tools conda environment can be run on a standard laptop.

<br>

### Troubleshooting conda environment installation

  If you run into any issues while installing the RNAseq conda environments, update conda using the following command then try re-creating the RNAseq conda environments using the installation instructions above.
  ```
  conda update conda
  ```

<br>
  

### Create RNAseq output directory structure

  To generate the same directory structure that GeneLab uses to organize the RNAseq processing pipeline output files, follow the instructions below:
  1. Use the `cd` command to navigate to the location on your device where you want to create the RNAseq directory structure.
  2. Set up the top level directory for the GLDS dataset you want to process using the `mkdir` command as follows (the example below is for GLDS-245):
  ```
  mkdir GLDS-245
  ```  
  3. Go into the top level directory that you created in step 2 using the `cd` command as shown in the example below:
  ```
  cd GLDS-245
  ``` 
  4. Create all sub-directories within the top level directory using the GL_RNAseq_mkdir.sh script by executing the following command (make sure you first download the [GL_RNAseq_mkdir.sh]() file and save it in the topy level direcroty you made in step 2):
  ```
  ./GL_RNAseq_mkdir.sh
  ``` 
  5. Your directory structure should now be set up. 
  
