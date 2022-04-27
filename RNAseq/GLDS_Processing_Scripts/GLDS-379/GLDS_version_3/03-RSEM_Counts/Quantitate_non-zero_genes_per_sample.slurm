#!/bin/bash
#
#SBATCH --job-name="Quantitate_non-zero_genes_per_sample"
#SBATCH --output=Quantitate_non-zero_genes_per_sample.out
#SBATCH --partition=normal
#SBATCH --mem=10000
#
#SBATCH --mail-user=user@nasa.gov
#SBATCH --mail-type=END

. ~/.profile

echo "Quantitate_non-zero_genes_per_sample"
echo ""

start=$(date +%s)
echo "start time: $start"
echo $HOSTNAME
echo ""

source activate /miniconda38_admin/envs/RNAseq_Rtools_03-2022

echo "Location of R: "
which R
echo ""
echo "R version: "
R --version
echo ""

Rscript --save Quantitate_non-zero_genes_per_sample.R

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

