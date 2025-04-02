#!/usr/bin/env bash
set  -euo pipefail

# Script to launch the nextflow workflow on a slurm compute cluster

# Usage: bash ./launch.sh '[extra arguments]' [main.nf]
# Example: bash ./launch.sh '--accession GLDS-487 --target_region 16S' path/to/main.nf 


EXTRA=${1:-""}  # extra arguments to the nextflow run command
MAIN=${2:-"main.nf"} # path to the main.nf nextflow script


#==============================================================================
# SETUP START
#==============================================================================
eval "$(conda shell.bash hook)"
conda activate /path/to/envs/nextflow
export NXF_SINGULARITY_CACHEDIR=/path/to/singularity_images/
export TOWER_ACCESS_TOKEN=<YOUR ACCESS TOKEN>
export TOWER_WORKSPACE_ID=<YOUR WORKSPACE ID>

#==============================================================================
# UMASK CONFIGURATION 
#==============================================================================
echo "Setting umask to enable group read-access by default"
umask u=rwx,g=rx
echo "Umask settings for this launch: $(umask -S)"


#==============================================================================
# NEXTFLOW COMMAND START
#==============================================================================

RUN_NAME=MAIN_$(date +%Y%m%d%H%M%S)

RUN_COMMAND="nextflow run \
                    -name ${RUN_NAME} \
                    ${MAIN} \
                    -resume \
                    -profile slurm,singularity \
                    -with-tower \
                    -process.queue 'normal' \
                    -ansi-log false \
                    ${EXTRA}"

echo "Running command: ${RUN_COMMAND}"
echo ""
[ -d processing_scripts ] || mkdir processing_scripts
eval ${RUN_COMMAND} && echo ${RUN_COMMAND} > processing_scripts/command.txt

# Save the nextflow log to a file
echo "Creating Nextflow processing info file..."
nextflow log ${RUN_NAME} -f name,script > processing_scripts/nextflow_processing_info_GLAmpliseq.txt
echo nextflow log ${RUN_NAME} -f name,script >> processing_scripts/nextflow_processing_info_GLAmpliseq.txt
echo "Nextflow processing info written to processing_scripts/nextflow_processing_info_GLAmpliseq.txt"
 
# Set permissions on launch directory
echo ""
echo "Setting permissions on launch directory..."
chmod -R 755 .
echo "Permissions set to 755 recursively on launch directory"
