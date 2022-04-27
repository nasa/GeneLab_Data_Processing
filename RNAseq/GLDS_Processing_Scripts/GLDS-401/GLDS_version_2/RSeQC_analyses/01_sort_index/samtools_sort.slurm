#!/bin/bash
#
#SBATCH --job-name="GLDS-401_sort"
#SBATCH --output=samtools_sort_out_logs/sort_%A_%a.out
#SBATCH --mem=30000
#SBATCH --cpus-per-task=10
#SBATCH --partition=normal
#SBATCH --array=1-12
#
#SBATCH --mail-user=user@nasa.gov
#SBATCH --mail-type=END

. ~/.profile

echo "samtools_sort"
echo ""

start=$(date +%s)
echo "start time: $start"
echo $HOSTNAME
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo ""

source activate /miniconda38_admin/envs/samtools_03-2022

outpath="/GLDS-401/02-STAR_Alignment"

sample=$(cat samples.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)

echo ""
echo "SAMPLE: ${sample}"
echo ""

echo "samtools version: "
samtools --version
echo ""

call="samtools sort -m 3G \
	--threads 10 \
	-o $outpath/${sample}/${sample}_Aligned.sortedByCoord_sorted.out.bam \
	$outpath/${sample}/${sample}_Aligned.sortedByCoord.out.bam"

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
