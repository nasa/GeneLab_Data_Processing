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
  
  2. Set up the top level (aka parent) directory for the GLDS dataset you want to process using the `mkdir` command as follows (the example below is for GLDS-245):
  ```
  mkdir GLDS-245
  ```  
  
  3. Go into the top level directory that you created in step 2 using the `cd` command as shown in the example below:
  ```
  cd GLDS-245
  ``` 
  
  4. Create all sub-directories (aka child and subsequent grandchild directories) within the top level directory using the GL_RNAseq_mkdir.sh script by executing the following command (make sure you first download the [GL_RNAseq_mkdir.sh](https://github.com/nasa/GeneLab_Data_Processing/blob/master/RNAseq/RNAseq_Tool_Install/GL_RNAseq_mkdir.sh) file and save it in the topy level direcroty you made in step 2):
  ```
  bash GL_RNAseq_mkdir.sh
  ``` 

5. Your directory structure should now be set up. If you wish to process the dataset using the same scripts that were used to create the processed data for your select dataset in the [GeneLab Repository](https://genelab-data.ndc.nasa.gov/genelab/projects), you can do so by following the steps below:  

   1. Download the scripts and metadata files that were used to processes your select dataset into the respective subdirectories (aka grandchild directories) of the `processing_scripts` child directory you made in step 4. Note: The processing scripts for all GeneLab RNAseq processed datasets are provided in the [RNAseq/GLDS_Processing_Scripts](https://github.com/nasa/GeneLab_Data_Processing/tree/master/RNAseq/GLDS_Processing_Scripts) directory of this repository.  

   2. Use a text editor such as [nano](https://www.nano-editor.org/) to change the paths indicated in each of the processing scripts you downloaded in step 5i to match the path to the top level (aka parent) directory you created in step 2.

      > **Note1:** If your system uses the [slurm](https://slurm.schedmd.com/overview.html) job scheduler, you will also have to customize the #SBATCH options to be consistent with your system's slurm settings (consult your system administrator regarding the settings needed for your system).  
      >
      > **Note2:** If you are not using the slurm job scheduler, you will have to remove the #SBATCH options from each slurm script and replace them with the equivalant options for the job scheduler your system uses.  
   
   3. If you wish to start processing from the beginning of the pipeline using raw sequence data, download the raw fastq files from the [GeneLab Repository](https://genelab-data.ndc.nasa.gov/genelab/projects) for your select dataset by clicking on 'STUDY FILES' in the left panel then click the arrowhead next to 'GeneLab Processed RNA-Seq Files' then click on 'Merged sequence data' then select and download all the *raw.fastq.gz files into the `/GLDS-#/00-RawData/Fastq` directory you made in steps 3 and 4.  
   
      > **Note1:** The *raw.fastq.gz files for datasets sequenced by GeneLab are stored under 'STUDY FILES' -> 'RNA-Seq' -> 'Raw sequence data'.  
      >
      > **Note2:** GeneLab provides processed data files from each step of the RNAseq processing pipeline in the [GeneLab Repository](https://genelab-data.ndc.nasa.gov/genelab/projects), so you may start processing from any step of the pipeline by downloading the appropriate input data for that step.  
   
   4. Once you've downloaded the scripts and input data you need to start processing (steps 5i and 5iii, respectively) and have revised all of your processing scripts to indicate the correct paths and settings for your machine (step 5ii), you may use the `cd` command to navigate to the appropriare subdirectory within the `processing_scripts` child directory and begin executing the scripts in the same order as indicated in the [RNAseq processing pipeline](https://github.com/nasa/GeneLab_Data_Processing/blob/master/RNAseq/GL-DPPD-7101-C.md).  
