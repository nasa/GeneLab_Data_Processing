#!/bin/bash

. ~/.profile


echo "GLDS-379_raw_multiqc_report"
echo ""

start=$(date +%s)
echo "start time: $start"
echo $HOSTNAME

source activate /conda/envs/rnaseq_v1.0

in_dir=/GLDS-379/00-RawData/FastQC_Reports
out_dir=/GLDS-379/00-RawData/FastQC_Reports/raw_multiqc_report

echo ""
echo "MultiQC version: "
multiqc --version
echo ""


call="multiqc -n raw_multiqc -o $out_dir $in_dir/"

echo $call
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

