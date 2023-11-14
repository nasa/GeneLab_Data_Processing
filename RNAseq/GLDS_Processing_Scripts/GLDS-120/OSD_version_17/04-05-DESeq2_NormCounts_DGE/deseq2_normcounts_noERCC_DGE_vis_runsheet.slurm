#!/bin/bash
#
#SBATCH --job-name="deseq2_normcounts_noERCC_DGE_vis_runsheet"
#SBATCH --output=deseq2_normcounts_noERCC_DGE_vis_runsheet.out
#SBATCH --partition=priority
#SBATCH --mem=20000
#
#SBATCH --mail-user=user@nasa.gov
#SBATCH --mail-type=END

. ~/.profile


echo "deseq2_normcounts_noERCC_DGE_vis_runsheet"
echo ""

start=$(date +%s)
echo "start time: $start"
echo $HOSTNAME

source activate /miniconda38_admin/envs/RNAseq_Rtools_03-2022

echo ""
echo "R version: "
R --version
echo ""


Rscript --save deseq2_normcounts_noERCC_DGE_vis_runsheet.R


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
