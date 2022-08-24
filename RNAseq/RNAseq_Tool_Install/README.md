# Instructions for installing the tools required to process RNAseq data using the GeneLab pipeline

> **This directory holds yaml files and instructions for how to install two conda environments, RNAseq_fq_to_counts_tools and RNAseq_R_tools, that contain all the software programs needed to process RNAseq data using the GeneLab pipeline. Additionally, instructions and the shell scripts for creating the directory structure that GeneLab uses to organize the RNAseq processing pipeline output files are provided.**  

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
  > This environment needs to be activated to run steps 1-6 of the [RNAseq processing pipeline](../Pipeline_GL-DPPD-7101_Versions/GL-DPPD-7101-E.md)
  
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
  > This environment needs to be activated to run step 7 of the [RNAseq processing pipeline](../Pipeline_GL-DPPD-7101_Versions/GL-DPPD-7101-E.md)
  
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
  

### Create directory structure for Genome and STAR and RSEM Indices 

  To generate the same directory structure that GeneLab uses to organize genome reference files and STAR and RSEM indices needed for processing RNAseq data, follow the instructions below:
  1. Use the `cd` command to navigate to the location on your device where you want to create the genome and index directory structure.
  
  2. Download the [GL_genome_index_mkdir.sh](https://github.com/nasa/GeneLab_Data_Processing/blob/master/RNAseq/RNAseq_Tool_Install/GL_genome_index_mkdir.sh) file and save it to the location on your device you selected in step 1, then create the genome and index directories and sub-directories by executing the following command:
  ```
  bash GL_genome_index_mkdir.sh
  ```  
  
3. Your directory structure should now be set up. You'll next have to download the fasta and gtf files for your organism(s) of interest from [Ensembl](https://www.ensembl.org/) and save the fasta and gtf file in each respective organism's subdirectory in the 'Genomes' directory you made in step 2.  

   1. For animals you can find these Ensembl files [here](https://uswest.ensembl.org/index.html), plant Ensembl references can be found [here](https://plants.ensembl.org/index.html), and Ensembl reference files for microbes can be found [here](https://bacteria.ensembl.org/index.html). The reference files for the ERCC genes GeneLab uses can be downloaded [here](https://assets.thermofisher.com/TFS-Assets/LSG/manuals/ERCC92.zip). 

      > **Note:** Ensembl displays the current release of each organism's genome and annotation files; to use the exact release that was used to process data hosted in the [GeneLab Repository](https://genelab-data.ndc.nasa.gov/genelab/projects), download the files using the links provided in the [GeneLab_Reference_and_Annotation_Files](https://github.com/nasa/GeneLab_Data_Processing/tree/master/RNAseq/GeneLab_Reference_and_Annotation_Files) tables.
  
   2. If samples from the dataset you want to process were spiked with ERCC genes, you'll have to create a reference fasta and gtf file containing both the organism of interest and ERCC genes by concatenating your organism of interest's fasta file with the ERCC fasta file, and your organism of interest's gtf file with the ERCC gtf file. Below are example commands for how to do this with *Mus musculus* (mouse) Ensembl release 101:
      ```
      gunzip ./Genomes/Mus_musculus/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
      ```
      ```
      cat ./Genomes/Mus_musculus/Mus_musculus.GRCm38.dna.primary_assembly.fa ./Genomes/ERCC/ERCC92.fa > ./Genomes/Mus_musculus/Mus_musculus.GRCm38.dna.primary_assembly_and_ERCC92.fa
      ```
      ```
      gunzip ./Genomes/Mus_musculus/Mus_musculus.GRCm38.101.gtf.gz
      ```
      ```
      cat ./Genomes/Mus_musculus/Mus_musculus.GRCm38.101.gtf ./Genomes/ERCC/ERCC92.gtf > ./Genomes/Mus_musculus/Mus_musculus.GRCm38.101_and_ERCC92.gtf
      ```
  
  4. You can now set up your scripts for creating STAR indices for each organism of interest in the `./STAR_Indices/STAR_index_scripts` subdirectory you created in step 2 by following the instructions in [step 3](../Pipeline_GL-DPPD-7101_Versions/GL-DPPD-7101-E.md#3-build-star-reference) of the GL RNAseq processing pipeline. 
  
     > **Note:** Prior to creating a STAR index for your organism of interest, you will need to create a subdirectory in the `./STAR_Indices` directory you created in step 2, where the STAR index will be outputted. At GeneLab, we title these directories with the organism name and the read length of the raw sequence data. Below is an example for making a STAR index output directory for Mus musculus samples with ERCC spike-in having a read length of 149: 
     ```
     mkdir ./STAR_Indices/Mus_musculus_w_ERCC_RL-149
     ```
  
  5. You can now set up your scripts for creating RSEM indices for each organism of interest in the `./RSEM_Indices/RSEM_index_scripts` subdirectory you created in step 2 by following the instructions in [step 5](../Pipeline_GL-DPPD-7101_Versions/GL-DPPD-7101-E.md#5-build-rsem-reference) of the GL RNAseq processing pipeline. 
  
     > **Note:** Prior to creating a RSEM index for your organism of interest, you will need to create a subdirectory in the `./RSEM_Indices` directory you created in step 2, where the RSEM index will be outputted. At GeneLab, we title these directories with the organism name. Below is an example for making a RSEM index output directory for Mus musculus samples with ERCC spike-in: 
     ```
     mkdir ./RSEM_Indices/Mus_musculus_w_ERCC
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
  
  4. Create all sub-directories (aka child and subsequent grandchild directories) within the top level directory using the GL_RNAseq_mkdir.sh script by executing the following command (make sure you first download the [GL_RNAseq_mkdir.sh](https://github.com/nasa/GeneLab_Data_Processing/blob/master/RNAseq/RNAseq_Tool_Install/GL_RNAseq_mkdir.sh) file and save it in the top level directory you made in step 2):
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
   
   4. Once you've downloaded the scripts and input data you need to start processing (steps 5i and 5iii, respectively) and have revised all of your processing scripts to indicate the correct paths and settings for your machine (step 5ii), you may use the `cd` command to navigate to the appropriate subdirectory within the `processing_scripts` child directory and begin executing the scripts in the same order as indicated in the [RNAseq processing pipeline](../Pipeline_GL-DPPD-7101_Versions/GL-DPPD-7101-E.md).  
