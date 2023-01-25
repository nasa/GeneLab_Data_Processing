#!/bin/bash
#
#SBATCH --job-name="GLDS-104_STAR_align"
#SBATCH --output=STAR_align_out_logs/GLDS-104_STAR_align_%A_%a.out
#SBATCH --mem=45000
#SBATCH --array=1-12
#
#SBATCH --mail-user=user@nasa.gov
#SBATCH --mail-type=END

start=$(date +%s)
echo "start time: $start"
echo $HOSTNAME
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

out_dir="/GLDS-104/02-STAR_Alignment"
in_dir="/GLDS-104/01-TG_Preproc/Fastq"

REF="/STAR_Indices/Mus_musculus_RL-100"
GTF="/Genomes/Mus_musculus/Mus_musculus.GRCm38.96.gtf"

sample=$(cat samples.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)

mkdir $out_dir/${sample}
echo "SAMPLE: ${sample}"

call="STAR --twopassMode Basic \
	--limitBAMsortRAM 30000000000 \
	--genomeDir $REF \
	--outSAMunmapped Within \
	--outFilterType BySJout \
	--outSAMattributes NH HI AS NM MD MC \
	--outFilterMultimapNmax 20 \
	--outFilterMismatchNmax 999 \
	--outFilterMismatchNoverReadLmax 0.04 \
	--alignIntronMin 20 \
	--alignIntronMax 1000000 \
	--alignMatesGapMax 1000000 \
	--alignSJoverhangMin 8 \
	--alignSJDBoverhangMin 1 \
	--sjdbScore 1 \
	--readFilesCommand zcat \
	--runThreadN 15 \
	--outSAMtype BAM SortedByCoordinate \
	--quantMode TranscriptomeSAM \
	--outSAMheaderHD @HD VN:1.4 SO:coordinate \
	--outFileNamePrefix $out_dir/${sample}/${sample}_ \
	--readFilesIn $in_dir/${sample}_R1_trimmed.fastq.gz $in_dir/${sample}_R2_trimmed.fastq.gz"

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

