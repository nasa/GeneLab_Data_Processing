#!/bin/bash
#
#SBATCH --job-name="GLDS-314_align"
#SBATCH --output=STAR_align_out_logs/align_%A_%a.out
#SBATCH --partition=normal
#SBATCH --mem=60000
#SBATCH --cpus-per-task=18
#SBATCH --array=1-17
#
#SBATCH --mail-user=user@nasa.gov
#SBATCH --mail-type=END

. ~/.profile


echo "STAR_align"
echo ""

start=$(date +%s)
echo "start time: $start"
echo $HOSTNAME
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo ""

source activate /miniconda38_admin/envs/star_03-2022

in_dir=/GLDS-314/01-TG_Preproc/Fastq
out_dir=/GLDS-314/02-STAR_Alignment
REF=/STAR_Indices/ensembl_plants_release48/Arabidopsis_thaliana_RL-50

sample=$(cat samples.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)

echo ""
echo "SAMPLE: ${sample}"

echo ""
echo "STAR version: "
STAR --version
echo ""

mkdir $out_dir/${sample}

call="STAR --twopassMode Basic \
	--limitBAMsortRAM 55000000000 \
	--genomeDir $REF \
	--outSAMunmapped Within \
	--outFilterType BySJout \
	--outSAMattributes NH HI AS NM MD MC \
	--outFilterMultimapNmax 20 \
	--outFilterMismatchNmax 999 \
	--outFilterMismatchNoverReadLmax 0.04 \
	--alignIntronMin 20 \
	--alignIntronMax 1000000 \
	--alignSJoverhangMin 8 \
	--alignSJDBoverhangMin 1 \
	--sjdbScore 1 \
	--readFilesCommand zcat \
	--runThreadN 18 \
	--outSAMtype BAM SortedByCoordinate \
	--quantMode TranscriptomeSAM GeneCounts \
	--outSAMheaderHD @HD VN:1.4 SO:coordinate \
	--outFileNamePrefix $out_dir/${sample}/${sample}_ \
	--readFilesIn $in_dir/${sample}_trimmed.fastq.gz"

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

