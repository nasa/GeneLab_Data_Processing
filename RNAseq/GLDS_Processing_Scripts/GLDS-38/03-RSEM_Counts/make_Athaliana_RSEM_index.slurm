#!/bin/bash
#
#SBATCH --job-name="make_Athaliana_RSEM_index"
#SBATCH --output=make_Athaliana_RSEM_index.out
#SBATCH --mem=60000
#

start=$(date +%s)
echo "start time: $start"
echo $HOSTNAME

rsem-prepare-reference --gtf /path/to/Genomes/Arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.44.gtf \
/path/to/Genomes/Arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa \
/path/to/RSEM_Indices/Arabidopsis_thaliana/Athaliana

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
