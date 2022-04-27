#!/bin/bash
#
#SBATCH --job-name="GLDS-314_RSEM_count"
#SBATCH --output=rsem_count_out_logs/RSEM_counts_%A_%a.out
#SBATCH --partition=normal
#SBATCH --mem=40000
#SBATCH --cpus-per-task=14
#SBATCH --array=1-17
#
#SBATCH --mail-user=user@nasa.gov
#SBATCH --mail-type=END

. ~/.profile


echo "GLDS-314_RSEM_count"
echo ""

start=$(date +%s)
echo "start time: $start"
echo $HOSTNAME
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo ""

source activate /miniconda38_admin/envs/rsem_03-2022

input_bam=/GLDS-314/02-STAR_Alignment
counts_output=/GLDS-314/03-RSEM_Counts
rsem_ref=/RSEM_Indices/ensembl_plants_release48/Arabidopsis_thaliana/Atha

sample=$(cat samples.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)

echo ""
echo "SAMPLE: ${sample}"

echo ""
echo "RSEM version: "
rsem-calculate-expression --version
echo ""

call="rsem-calculate-expression --num-threads 14 \
	--alignments \
	--bam \
	--seed 12345 \
	--seed-length 20 \
	--estimate-rspd \
	--no-bam-output \
	--strandedness reverse \
	$input_bam/${sample}/${sample}_Aligned.toTranscriptome.out.bam \
	$rsem_ref \
	$counts_output/${sample}"

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
