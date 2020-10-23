#!/bin/bash
#
#SBATCH --job-name="GLDS-235_deseq2_normcounts_wERCC_DGE_vis"
#SBATCH --output=deseq2_normcounts_wERCC_DGE_vis.out
#SBATCH --mem=40000
#
#SBATCH --mail-type=END

start=$(date +%s)
echo "start time: $start"
echo $HOSTNAME

Rscript --save deseq2_normcounts_wERCC_DGE_vis.R

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
