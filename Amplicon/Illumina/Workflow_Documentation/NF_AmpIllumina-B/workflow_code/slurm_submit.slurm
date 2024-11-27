#!/bin/bash

#SBATCH --job-name="nf_master" ## Replace job_name with the name of the job you are running ##
#SBATCH --output=nf_master.o.%j ## Replace job_name with the name of the job you are running ##
#SBATCH --error=nf_master.e.%j ## Replace job_name with the name of the job you are running ##
#SBATCH --partition=normal ## Specifies the job queue to use, for urgent jobs change normal to priority ##
#SBATCH --mem=20G ## Memory required to run the job in MB, this example is showing 10,000 MB or 10GB, change this number based on how much RAM you need ##
#SBATCH --cpus-per-task=1 ## Number of CPUs to run the job, this example is showing 5 CPUs, change this number based on how many CPUs you need ##
#SBATCH --mail-user=name@domain ## Specifies the e-mail address to e-mail when the job is complete, replace this e-mail address with your NASA e-mail address ##
#SBATCH --mail-type=END ## Tells slurm to e-mail the address above when the job has completed ##

. ~/.profile


echo "nf_master" ## Replace job_name with the name of the job you are running ##
echo ""


## Add a time-stamp at the start of the job ##
start=$(date +%s)
echo "start time: $start"

## Print the name of the compute node executing the job ##
echo $HOSTNAME


## Activate the conda environemnt containing the tools you need to run your job ##
## You can see a list of all available environments by running the command: conda env list ##

source activate /path/to/envs/nextflow/  ## Replace conda_env_name with the name of the environment ##


## Print the version of the tool you are using to ensure the tool version is recorded ##
echo ""
echo "Nextflow version: " ## Replace Tool with the name of the tool you are using ##
nextflow -v ## Replace this command with the command the tool uses to print its version ##
echo ""


## The command(s) that you want to run in this slurm job ##
export NXF_SINGULARITY_CACHEDIR=singularity/
export TOWER_ACCESS_TOKEN=<ACCESS TOKEN>
export TOWER_WORKSPACE_ID=<WORKSPACE ID>
#nextflow run main.nf -resume -profile slurm,singularity --input_file PE_file.csv --target_region 16S --F_primer AGAGTTTGATCCTGGCTCAG --R_primer CTGCCTCCCGTAGGAGT  -with-tower ## Replace command with the command(s) you want to run ##
#nextflow run main.nf -resume -profile slurm,conda --accession GLDS-487 --target_region 16S -with-tower 
#python3 run_workflow.py --run --target-region 16S --input-file PE_file.csv --F-primer GTGCCAGCMGCCGCGGTAA --R-primer GGACTACHVGGGTWTCTAAT --profile slurm,singularity --singularity-cacheDir /path/to/singularity_images/ --extra 'with-tower' --R-memory '20 GB'
python3 run_workflow.py --run --target-region 16S --accession GLDS-487 --profile slurm,singularity --singularity-cacheDir /path/to/singularity_images/ --extra 'with-tower' --R-memory '20 GB'
## Add a time-stamp at the end of the job then calculate how long the job took to run in seconds, minutes, and hours ##
echo ""
end=$(date +%s)
echo "end time: $end"
runtime_s=$(echo $(( end - start )))
echo "total run time(s): $runtime_s"
sec_per_min=60
sec_per_hr=3600
runtime_m=$(echo "scale=2; $runtime_s / $sec_per_min;" | bc)
echo "total run time(m): $runtime_m"
runtime_h=$(echo "scale=2; $runtime_s / $sec_per_hr;" | bc)
echo "total run time(h): $runtime_h"
echo ""


## Print the slurm job ID so you have it recorded and can view slurm job statistics if needed ##
echo "slurm job ID: ${SLURM_JOB_ID}"
