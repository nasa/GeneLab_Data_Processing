# Estimate Host Reads Workflow Information and Usage Instructions

## General Workflow Information
The workflow for estimating host DNA in Illumina metagenomics sequencing data, as specified in the "Estimate Host Reads" workflow, is implemented using Nextflow. This workflow is intended to be run on any Unix-based system using Docker containers to ensure consistency and reproducibility of the computational environment.

## Utilizing the Workflow

1. **Install Docker**
2. **Download the workflow template files**
3. **Modify the variables in the Nextflow config file**
4. **Run the workflow**

### 1. Install Docker
We recommend installing Docker to handle all dependencies within containers. This simplifies the setup on any system and avoids compatibility issues:
```bash
# Install Docker following the official guidelines:
https://docs.docker.com/get-docker/
```

### 2. Download the Workflow Template Files
Clone the repository or download the workflow files from the designated repository or link. Ensure you have all required files, including the Nextflow script (.nf) and the associated configuration files.

### 3. Modify the Variables in the Nextflow Config File
Adjust the variables in the `nextflow.config` file to match your specific needs. This includes paths to input data, Docker containers, and output directories.

Once you've downloaded the workflow template, you can modify the variables in the [Estimate_Host_Reads.config](workflow_code/Estimate_Host_Reads.config) file as needed. For example, you will have to provide a text file containing a single-column list of unique sample identifiers (see an example of how to set this up below - if you are running the example dataset, this file is provided in the [workflow_code](workflow_code) directory [here](workflow_code/unique-sample-IDs.txt)). You will also need to indicate the path to your input data (raw reads) and the root directory for where the kraken2 reference database should be stored (it will be setup automatically). Additionally, if necessary, you'll need to modify each variable in the [config.yaml](workflow_code/config.yaml) file to be consistent with the study you want to process and the machine you're using. 

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

### 4. Run the Workflow
Navigate to the directory containing the Nextflow script and config file. Here is an example command to run the workflow:
```bash
nextflow run Estimate_Host_Reads.nf -profile docker
```
- `-profile docker` specifies that Docker containers should be used to run the tools.

See `nextflow run -help` and [Nextflow's documentation](https://www.nextflow.io/docs/latest/index.html) for more options and detailed information.

## Reference Database Information
The database used for host estimation is maintained and updated periodically. Links to download the latest version of the database can be found in the workflow documentation.

---

This workflow is designed for flexibility and can be adapted for various datasets and research needs. For any additional information or support, refer to the official Nextflow documentation or contact the support team.
