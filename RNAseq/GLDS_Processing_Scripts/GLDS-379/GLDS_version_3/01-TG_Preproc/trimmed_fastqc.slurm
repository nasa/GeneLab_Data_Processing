#!/bin/bash
#
#SBATCH --job-name="GLDS-379_trimmed_fastqc"
#SBATCH --output=trimmed_fastqc_out_logs/trimmed_fastqc_%A_%a.out
#SBATCH --partition=normal
#SBATCH --mem=6000
#SBATCH --array=1-141
#
#SBATCH --mail-user=user@nasa.gov
#SBATCH --mail-type=END

. ~/.profile

echo "GLDS-379_trimmed_fastqc"
echo ""

start=$(date +%s)
echo "start time: $start"
echo $HOSTNAME
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

source activate /miniconda38_admin/envs/fastqc_03-2022

in_dir=/GLDS-379/01-TG_Preproc/Fastq
out_dir=/GLDS-379/01-TG_Preproc/FastQC_Reports

sample=$(cat samples.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)

echo ""
echo "SAMPLE: ${sample}"

echo ""
echo "FastQC version: "
fastqc --version
echo ""

call="fastqc -o $out_dir $in_dir/${sample}_R*_trimmed.fastq.gz"

echo $call
echo ""
eval $call


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
