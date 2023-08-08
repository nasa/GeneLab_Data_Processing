#!/bin/bash
#
#SBATCH --job-name="deseq2_normcounts_noERCC_DGE_vis"
#SBATCH --output=deseq2_normcounts_noERCC_DGE_vis.out
#SBATCH --mem=10000
#
#SBATCH --mail-user=user@nasa.gov
#SBATCH --mail-type=END

. ~/.profile


echo "deseq2_normcounts_noERCC_DGE_vis"
echo ""

start=$(date +%s)
echo "start time: $start"
echo $HOSTNAME

source activate /conda/envs/RNAseq_Rtools

echo ""
echo "R version: "
R --version
echo ""


Rscript --save deseq2_normcounts_noERCC_DGE_vis.R


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

