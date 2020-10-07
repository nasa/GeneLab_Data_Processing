#!/bin/bash
#
#SBATCH --job-name="GLDS-38_trimmed_multiqc_report"
#SBATCH --output=GLDS-38_trimmed_multiqc_report.out
#SBATCH --mem=16000
#

start=$(date +%s)
echo "start time: $start"
echo $HOSTNAME

in_dir=/path/to/GLDS-38/01-TG_Preproc/FastQC_Reports
out_dir=/path/to/GLDS-38/01-TG_Preproc/FastQC_Reports/trimmed_multiqc_report

call="multiqc -o $out_dir $in_dir/"

echo $call
eval $call

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
