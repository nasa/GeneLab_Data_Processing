#!/bin/bash
#
#SBATCH --job-name="GLDS-245_rsem_count"
#SBATCH --output=rsem_count_out_logs/GLDS-245_rsem_count_%A_%a.out
#SBATCH --mem=40000
#SBATCH --array=1-59
#
#SBATCH --mail-type=END

start=$(date +%s)
echo "start time: $start"
echo $HOSTNAME
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

input_bam="/path/to/GLDS-245/02-STAR_Alignment"
rsem_ref="/path/to/RSEM_Indices/Mus_musculus_w_ERCC/Mmus_ERCC"
counts_output="/path/to/GLDS-245/03-RSEM_Counts"

sample=$(cat samples.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)

echo "SAMPLE: ${sample}"

call="rsem-calculate-expression --num-threads 6 \
	--alignments \
	--bam \
	--paired-end \
	--seed 12345 \
	--estimate-rspd \
	--no-bam-output \
	--strandedness reverse \
	$input_bam/${sample}/${sample}_Aligned.toTranscriptome.out.bam \
	$rsem_ref \
	$counts_output/${sample}"

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
