#!/bin/bash
#
#SBATCH --job-name="GLDS-379_read_dist"
#SBATCH --output=04_out_logs/read_dist_%A_%a.out
#SBATCH --partition=normal
#SBATCH --mem=30000
#SBATCH --array=1-141
#
#SBATCH --mail-user=user@nasa.gov
#SBATCH --mail-type=END

. ~/.profile

echo "inner_distance"
echo ""

start=$(date +%s)
echo "start time: $start"
echo $HOSTNAME
echo ""

source activate /miniconda38_admin/envs/rseqc_03-2022

in_dir=/GLDS-379/02-STAR_Alignment
out_dir=/GLDS-379/RSeQC_analyses/04_inner_distance
genome_files=/Genome_GTF_BED_Files/ensembl_release101/Mus_musculus

sample=$(cat samples.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)

echo ""
echo "SAMPLE: ${sample}"
echo ""

echo "inner_distance.py version: "
inner_distance.py --version
echo ""

mkdir $out_dir/${sample}

call="inner_distance.py -r $genome_files/Mus_musculus.GRCm38.101_and_ERCC92.gtf.bed -i $in_dir/${sample}/${sample}_Aligned.sortedByCoord_sorted.out.bam -k 15000000 -l -150 -u 350 -o $out_dir/${sample}/${sample}"

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
