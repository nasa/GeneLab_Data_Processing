# Snakemake RNASeq pipeline
RNASeq Snakemake workflow. 

-- Incomplete --  
Pre-processing to RSEM counts.

### Install
Clone repo, cd into dir and create environment
```bash
conda env create -f rna_pipeline.yaml
conda activate rna_pipeline
```
### Prepare analysis dir
Create an analysis directory and copy `Snakefile, config.json, envs/` files/dir from cloned dir.  
Edit parameters as required in `config.json`.   
### Usage
For cluster environment supply `--cluster` param and `job_script.sh`
Start pipeline locally with :  
`` -j {int} `` Number of jobs to start
```bash
snakemake --use-conda --conda-frontend conda --conda-prefix ${CONDA_PREFIX}/envs -j 1 -p
```