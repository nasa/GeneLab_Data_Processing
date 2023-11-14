#!/bin/bash
#
#SBATCH --job-name="GLDS-120_rsem_count"
#SBATCH --output=rsem_count_out_logs/GLDS-120_rsem_count_%A_%a.out
#SBATCH --mem=40000
#SBATCH --array=1-36
#
#SBATCH --mail-user=user@nasa.gov
#SBATCH --mail-type=END

. ~/.profile

echo "GLDS-120_RSEM_count"
echo ""


start=$(date +%s)
echo "start time: $start"
echo ""
echo $HOSTNAME
echo ""
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo ""

source activate /conda/envs/rnaseq_v1.0

input_bam="/GLDS-120/02-STAR_Alignment"
rsem_ref="/RSEM_Indices/ensembl_release44/Arabidopsis_thaliana/Atha"
counts_output="/GLDS-120/03-RSEM_Counts_fwd"

sample=$(cat samples.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)

echo "SAMPLE: ${sample}"
echo ""
echo "Pull Aligned.toTranscriptome.out.bam file from s3"
echo ""

aws s3 cp s3://genelab-repo-prod/genelab-data/GLDS-120/rna_seq/GLDS-120_rna_seq_${sample}_Aligned.toTranscriptome.out.bam /GLDS-120/02-STAR_Alignment/${sample}_Aligned.toTranscriptome.out.bam
echo ""

echo ""
echo "RSEM version: "
rsem-calculate-expression --version
echo ""


call="rsem-calculate-expression --num-threads 8 \
	--alignments \
	--bam \
	--paired-end \
	--seed 12345 \
	--estimate-rspd \
	--no-bam-output \
	--strandedness none \
	$input_bam/${sample}_Aligned.toTranscriptome.out.bam \
	$rsem_ref \
	$counts_output/${sample}"

echo $call
echo ""
eval $call
echo ""

echo "Remove bam file"
rm /GLDS-120/02-STAR_Alignment/${sample}_Aligned.toTranscriptome.out.bam

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

