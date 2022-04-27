#!/bin/bash
#
#SBATCH --job-name="GLDS-314_infer_expt"
#SBATCH --output=03_out_logs/infer_expt_%A_%a.out
#SBATCH --partition=normal
#SBATCH --mem=20000
#SBATCH --array=1-17
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

source activate /miniconda38_admin/envs/rseqc_03-2022

in_dir="/GLDS-314/02-STAR_Alignment"
out_dir="/GLDS-314/RSeQC_analyses/03_infer_experiment"
genome_files="/Genome_GTF_BED_Files/ensembl_plants_release48/Arabidopsis_thaliana"

sample=$(cat samples.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)

echo ""
echo "SAMPLE: ${sample}"
echo ""

echo "infer_experiment.py version: "
infer_experiment.py --version
echo ""

call="infer_experiment.py -r $genome_files/Arabidopsis_thaliana.TAIR10.48.gtf.bed -i $in_dir/${sample}/${sample}_Aligned.sortedByCoord_sorted.out.bam -s 15000000 > $out_dir/${sample}_infer_expt.out"

echo $call
echo ""
eval $call
echo ""

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
