#!/bin/bash
#
#SBATCH --job-name="GLDS-401_raw_multiqc_report"
#SBATCH --output=GLDS-401_raw_multiqc_report.out
#SBATCH --partition=normal
#SBATCH --mem=10000
#
#SBATCH --mail-user=user@nasa.gov
#SBATCH --mail-type=END

. ~/.profile


echo "GLDS-401_raw_multiqc_report"
echo ""

start=$(date +%s)
echo "start time: $start"
echo $HOSTNAME

source activate /miniconda38_admin/envs/multiqc_03-2022

in_dir=/GLDS-401/00-RawData/FastQC_Reports
out_dir=/GLDS-401/00-RawData/FastQC_Reports/raw_multiqc_report

echo ""
echo "MultiQC version: "
multiqc --version
echo ""

call="multiqc --interactive -n raw_multiqc -o $out_dir $in_dir/"

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
