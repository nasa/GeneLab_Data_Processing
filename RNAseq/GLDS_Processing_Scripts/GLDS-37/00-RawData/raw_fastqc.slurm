#!/bin/bash
#
#SBATCH --job-name="GLDS-37_raw_fastqc"
#SBATCH --output=raw_fastqc_out_logs/GLDS-37_raw_fastqc_%A_%a.out
#SBATCH --mem=6000
#SBATCH --array=1-56
#

start=$(date +%s)
echo "start time: $start"
echo $HOSTNAME
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

in_dir=/path/to/GLDS-37/00-RawData/Fastq
out_dir=/path/to/GLDS-37/00-RawData/FastQC_Reports

sample=$(cat samples.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)

echo "SAMPLE: ${sample}"

call="fastqc -o $out_dir $in_dir/${sample}*fastq.gz"

echo $call
eval $call

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
