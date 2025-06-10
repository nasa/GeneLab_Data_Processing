!/usr/bin/env bash
set  -euo pipefail

# Script to launch a nextflow workflow on slurm cluster using singularity

# Usage: bash ./launch.sh [mode] [main.nf]  [config] '[extra arguments]'
# Examples

# Processing:
#      bash ./launch.sh processing path/to/main.nf path/to/nextflow.config '--accession GLDS-487 --target_region 16S'

# Postprocessing:
#      bash ./launch.sh post_processing  path/to/post_processing.nf  path/to/post_processing.config \
#           '--name First M. Last --email name@domain.com --GLDS_accession GLDS-487 --OSD_accession ODS-487 --isa_zip  ../GeneLab/OSD-487_metadata_GLDS-487-ISA.zip --runsheet ../GeneLab/GLfile.csv'



MODE=${1:-''} # Script run mode i.e. processing or post_processing.
MAIN=${2:-''}  # Path to main.nf or post_processing.nf nextflow script for processing and post_processing, respectively.
CONFIG=${3:-''} # nextflow config file i.e. nextflow.config or post_processing.config.
EXTRA=${4:-''}  # extra arguments to the nextflow run command particularly those in the params sections of the .config files.


#==============================================================================
# SETUP START
#==============================================================================
eval "$(conda shell.bash hook)"
conda activate /path/to/conda/envs/nextflow
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
if [ ${MODE} == "processing" ]; then

    RUN_NAME=MAIN_$(date +%Y%m%d%H%M%S)

    RUN_COMMAND="nextflow -C ${CONFIG}
                    run \
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


elif [ ${MODE} == "post_processing" ];then


    RUN_NAME=POST_$(date +%Y%m%d%H%M%S)

    RUN_COMMAND="nextflow -C ${CONFIG}
                    run \
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
    eval ${RUN_COMMAND}

else
    echo 'Please provide a valid mode to run the workflow.'
    echo 'Either processing or post_processing for running the processing or post_processing workflows, respectively.'
    exit 1
fi


# Set permissions on launch directory
echo ""
echo "Setting permissions on launch directory..."
chmod -R 755 .
echo "Permissions set to 755 recursively on launch directory"
                                                                                            
