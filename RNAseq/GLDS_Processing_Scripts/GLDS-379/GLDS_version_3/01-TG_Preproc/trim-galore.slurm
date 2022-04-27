#!/bin/bash
#
#SBATCH --job-name="trimgalore"
#SBATCH --output=TG_out_logs/trimgalore_%A_%a.out
#SBATCH --partition=normal
#SBATCH --mem=24000
#SBATCH --cpus-per-task=6
#SBATCH --array=1-141
#
#SBATCH --mail-user=user@nasa.gov
#SBATCH --mail-type=END

. ~/.profile


echo "trimgalore"
echo ""

start=$(date +%s)
echo "start time: $start"
echo $HOSTNAME

source activate /miniconda38_admin/envs/trim_galore_03-2022

out_dir=/GLDS-379/01-TG_Preproc/Fastq
in_dir=/GLDS-379/00-RawData/Fastq
trimRepdir=/GLDS-379/01-TG_Preproc/Trimming_Reports

sample=$(cat samples.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)

echo ""
echo "SAMPLE: ${sample}"

echo ""
echo "TrimGalore! version: "
trim_galore --version
echo ""
echo "cutadapt version: "
cutadapt --version
echo ""


call1="trim_galore --gzip \
	--cores 6 \
	--phred33 \
	--illumina \
	--paired \
	--output_dir $out_dir \
	$in_dir/${sample}_R1_raw.fastq.gz $in_dir/${sample}_R2_raw.fastq.gz"

echo $call1
echo ""
eval $call1


call2="mv $out_dir/${sample}_*trimming_report.txt $trimRepdir/"

echo $call2
echo ""
eval $call2


call3="mv $out_dir/${sample}_R1_raw_val_1.fq.gz $out_dir/${sample}_R1_trimmed.fastq.gz"

echo $call3
echo ""
eval $call3


call4="mv $out_dir/${sample}_R2_raw_val_2.fq.gz $out_dir/${sample}_R2_trimmed.fastq.gz"

echo $call4
echo ""
eval $call4


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
