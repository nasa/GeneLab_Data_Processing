#!/bin/bash

#SBATCH --job-name="nf_master" ## Replace job_name with the name of the job you are running ##
#SBATCH --output=nf_master.o.%j ## Replace job_name with the name of the job you are running ##
#SBATCH --error=nf_master.e.%j ## Replace job_name with the name of the job you are running ##
#SBATCH --partition=<patition_name> ## Specifies the job queue to use, for urgent jobs change normal to priority ##
#SBATCH --mem=20G ## Memory required to run the job in MB, this example is showing 10,000 MB or 10GB, change this number based on how much RAM you need ##
#SBATCH --cpus-per-task=1 ## Number of CPUs to run the job, this example is showing 5 CPUs, change this number based on how many CPUs you need ##
#SBATCH --mail-user=user@domain.com ## Specifies the e-mail address to e-mail when the job is complete, replace this e-mail address with your NASA e-mail address ##
#SBATCH --mail-type=END ## Tells slurm to e-mail the address above when the job has completed ##

# Load your user profile
. ~/.profile


echo "nf_master" ## Replace job_name with the name of the job you are running ##
echo ""


## Add a time-stamp at the start of the job ##
start=$(date +%s)
echo "start time: $start"

## Print the name of the compute node executing the job ##
echo $HOSTNAME

# Launch the workflow via the launch
# The workflow must be run this way if you'd like to run the post-processing workflow after it completes 
WORKFLOW_DIR='/path/to/the/workflow_code'

# Make sure you edit the following lines of the 'SETUP START' section 
# in the ./launch.sh script before submitting this job:
#conda activate /path/to/conda/envs/nextflow
#export NXF_SINGULARITY_CACHEDIR=/path/to/singularity_images/
#export TOWER_ACCESS_TOKEN=<YOUR ACCESS TOKEN>
#export TOWER_WORKSPACE_ID=<YOUR WORKSPACE ID>

# Processing
bash ./launch.sh processing ${WORKFLOW_DIR}/main.nf ${WORKFLOW_DIR}/nextflow.config '--accession GLDS-487'

# Post Processing
#bash ./launch.sh post_processing  ${WORKFLOW_DIR}/post_processing.nf  ${WORKFLOW_DIR}/post_processing.config \
#          '--name First M. Last --email name@domain.com --GLDS_accession GLDS-487 --OSD_accession ODS-487 --isa_zip  ../GeneLab/OSD-487_metadata_GLDS-487-ISA.zip --runsheet ../GeneLab/GLfile.csv'




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
