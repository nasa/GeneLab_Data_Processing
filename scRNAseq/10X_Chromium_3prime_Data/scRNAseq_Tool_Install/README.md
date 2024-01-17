# Instructions for installing the tools required to process single cell RNAseq data using the GeneLab pipeline

> **This directory holds yaml files and instructions for how to install the conda environments that contain all the software programs needed to process single cell RNAseq data using the GeneLab pipeline. 
Additionally, instructions and the shell scripts for creating the directory structure that GeneLab uses to organize the single cell RNAseq processing pipeline output files are provided.**  

---

### Install Conda

To install conda environments, you'll first have to install [Anaconda](https://www.anaconda.com/). We recommend installing a Miniconda, Python3 version appropriate for your system, as instructed by [Happy Belly Bioinformatics](https://astrobiomike.github.io/unix/conda-intro#getting-and-installing-conda). 

<br>

---

### Download the yaml files

All the yaml files required to install the conda environments needed to run the scRNAseq pipeline are in the [scRNAseq_yaml_files](scRNAseq_yaml_files) sub-directory. To get a 
copy of the yaml files on to your system, copy the github web address of the [scRNAseq_yaml_files directory](scRNAseq_yaml_files), then paste it into [GitZip here](http://kinolien.github.io/gitzip/), and click download.


After download, unzip the scRNAseq_yaml_files.zip file and cd into the directory, where you will find 3 yaml files:
* [fastqc.yml](scRNAseq_yaml_files/fastqc.yml)
* [multiqc.yml](scRNAseq_yaml_files/multiqc.yml)
* [star.yml](scRNAseq_yaml_files/star.yml)

<br>

---

### Install the GeneLab scRNAseq tools via conda environments

#### Install the **fastqc** conda environment by running the following command:

  ```
  conda env create -f fastqc.yml
  ```

  Activate the fastqc conda environment with the following command:
  > This environment needs to be activated to run steps [1a](../Pipeline_GL-DPPD-7111_Versions/GL-DPPD-7111.md#1a-raw-data-qc) of the [scRNAseq processing pipeline](../Pipeline_GL-DPPD-7111_Versions/GL-DPPD-7111.md).
  
  ```
  conda activate fastqc
  ```   
  
<br>

#### Install the **multiqc** conda environment by running the following command:

  ```
  conda env create -f multiqc.yml
  ```

  Activate the multiqc conda environment with the following command:
  > This environment needs to be activated to run steps [1b](../Pipeline_GL-DPPD-7111_Versions/GL-DPPD-7111.md#1b-compile-raw-data-qc) and [3b](../Pipeline_GL-DPPD-7111_Versions/GL-DPPD-7111.md#3b-compile-alignment-logs) of the [scRNAseq processing pipeline](../Pipeline_GL-DPPD-7111_Versions/GL-DPPD-7111.md).
  
  ```
  conda activate multiqc
  ```   
 
  
<br>

#### Install the **star** conda environment by running the following command:

  ```
  conda env create -f star.yml
  ```

  Activate the star conda environment with the following command:
  > This environment needs to be activated to run steps [2](../Pipeline_GL-DPPD-7111_Versions/GL-DPPD-7111.md#2-build-star-reference) and [3a](../Pipeline_GL-DPPD-7111_Versions/GL-DPPD-7111.md#3a-align-reads-to-reference-genome-with-starsolo) of the [scRNAseq processing pipeline](../Pipeline_GL-DPPD-7111_Versions/GL-DPPD-7111.md).
  
  ```
  conda activate star
  ```   
  
<br>


---

<br>

### Troubleshooting conda environment installation

  If you run into any issues while installing the RNAseq conda environments, update conda using the following command then try re-creating the RNAseq conda environments using the installation instructions above.
  
  ```
  conda update conda
  ```

<br>
  
---

### Create RNAseq output directory structure

To generate the same directory structure that GeneLab uses to organize the scRNAseq processing pipeline output files, follow the instructions below:

  1. Use the `cd` command to navigate to the location on your device where you want to create the RNAseq directory structure.
  
  2. Set up the top level (aka parent) directory for the GLDS dataset you want to process using the `mkdir` command as follows (the example below is for GLDS-405):
  ```
  mkdir GLDS-405
  ```  
  
  3. Go into the top level directory that you created in step 2 using the `cd` command as shown in the example below:
  ```
  cd GLDS-405
  ``` 
  
  4. Download the [GL_scRNAseq_mkdir.sh](GL_scRNAseq_mkdir.sh) bash script into the directory you created in step 2. To download the bash script, copy the github web address of the [GL_scRNAseq_mkdir.sh bash script](GL_scRNAseq_mkdir.sh), then paste it into [GitZip here](http://kinolien.github.io/gitzip/), and click download.
     > Make the downloaded GL_scRNAseq_mkdir.sh bash script executable by running the following command:  
     >   
     > `chmod u+x GL_scRNAseq_mkdir.sh`
  
  5. Create all sub-directories (aka child and subsequent grandchild directories) within the top level directory (created in step 2) by executing the GL_scRNAseq_mkdir.sh script as follows:
  ```
  bash GL_scRNAseq_mkdir.sh
  ``` 

5. Your directory structure should now be set up. If you wish to process the dataset using the same scripts that were used to create the processed data for your select dataset in the [GeneLab 
Repository](https://genelab-data.ndc.nasa.gov/genelab/projects), you can do so by following the steps below:  

   1. Download the scripts and metadata files that were used to processes your select dataset into the respective subdirectories (aka grandchild directories) of the `processing_scripts` child directory you made in step 5. 
      > Note: The processing scripts 
for all GeneLab scRNAseq processed datasets are provided in the [GLDS_Processing_Scripts](../GLDS_Processing_Scripts) directory of this repository.  

   2. Use a text editor such as [nano](https://www.nano-editor.org/) to change the paths indicated in each of the processing scripts you downloaded in step 5i to match the path to the top level (aka parent) directory you created in step 2.

      > **Note1:** If your system uses the [slurm](https://slurm.schedmd.com/overview.html) job scheduler, you will also have to customize the #SBATCH options to be consistent with your system's slurm settings (consult your system administrator 
regarding the settings needed for your system).  
      >
      > **Note2:** If you are not using the slurm job scheduler, you will have to remove the #SBATCH options from each slurm script and replace them with the equivalent options for the job scheduler your system uses.  
   
   3. If you wish to start processing from the beginning of the pipeline using raw sequence data, download the raw fastq files from the [GeneLab Repository](https://genelab-data.ndc.nasa.gov/genelab/projects) for your select dataset by clicking on 
'STUDY FILES' in the left panel then click the arrowhead next to 'GeneLab Processed RNA-Seq Files' then click on 'Merged sequence data' then select and download all the *raw.fastq.gz files into the `/GLDS-#/00-RawData/Fastq` directory you made in 
steps 2-4.  
   
      > **Note1:** The *raw.fastq.gz files for datasets sequenced by GeneLab are stored under 'STUDY FILES' -> 'RNA-Seq' -> 'Raw sequence data'.  
      >
      > **Note2:** GeneLab provides processed data files from each step of the scRNAseq processing pipeline in the [GeneLab Repository](https://genelab-data.ndc.nasa.gov/genelab/projects), so you may start processing from any step of the pipeline by 
downloading the appropriate input data for that step.  
   
   4. Once you've downloaded the scripts and input data you need to start processing (steps 5i and 5iii, respectively) and have revised all of your processing scripts to indicate the correct paths and settings for your machine (step 5ii), you may begin executing the scripts in the same order as indicated in the [scRNAseq processing 
pipeline](../Pipeline_GL-DPPD-7111_Versions/GL-DPPD-7111.md).  

