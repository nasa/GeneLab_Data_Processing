#!/bin/bash
#
#SBATCH --job-name="GLDS-120_infer_expt"
#SBATCH --output=03_out_logs/GLDS-120_infer_expt_%A_%a.out
#SBATCH --mem=20000
#SBATCH --array=1-3
#
#SBATCH --mail-user=user@nasa.gov
#SBATCH --mail-type=END

. ~/.profile

echo "infer_experiment"
echo ""

start=$(date +%s)
echo "start time: $start"
echo $HOSTNAME
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo ""

source activate /conda/envs/rnaseq_v1.0

in_dir="/GLDS-120/02-STAR_Alignment"
out_dir="/GLDS-120/RSeQC_analyses/03_infer_experiment"
genome_files="/Genome_GTF_BED_Files/ensembl_release44"

sample=$(cat samples.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)

echo ""
echo "SAMPLE: ${sample}"
echo ""

echo "Pull Aligned.sortedByCoord.out.bam file from s3"
echo ""

aws s3 cp s3://genelab-repo-prod/genelab-data/GLDS-120/rna_seq/GLDS-120_rna_seq_${sample}_Aligned.sortedByCoord.out.bam /GLDS-120/02-STAR_Alignment/${sample}_Aligned.sortedByCoord.out.bam
echo ""

echo "Index Aligned.sortedByCoord.out.bam file"
echo ""

echo "samtools version: "
samtools --version
echo ""

call1="samtools index -@ 4 \
	$in_dir/${sample}_Aligned.sortedByCoord.out.bam"

echo $call1
echo ""
eval $call1
echo ""

echo "infer_experiment.py version: "
infer_experiment.py --version
echo ""

call="infer_experiment.py -r $genome_files/Arabidopsis_thaliana.TAIR10.44.gtf.bed -i $in_dir/${sample}_Aligned.sortedByCoord.out.bam -s 15000000 > $out_dir/${sample}_infer_expt.out"

echo $call
echo ""
eval $call
echo ""

echo "Remove bam file"
rm /GLDS-120/02-STAR_Alignment/${sample}_Aligned.sortedByCoord.out.bam

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

