#!/bin/bash
#
#SBATCH --job-name="make_Bsubtilis-168_RSEM_index"
#SBATCH --output=make_Bsubtilis-168_RSEM_index.out
#SBATCH --mem=30000
#
#SBATCH --mail-type=END

start=$(date +%s)
echo "start time: $start"
echo $HOSTNAME

rsem-prepare-reference --gtf /path/to/Genomes/Bacillus_subtilis/Bacillus_subtilis_subsp_subtilis_str_168.ASM904v1.46.gtf \
/path/to/Genomes/Bacillus_subtilis/Bacillus_subtilis_subsp_subtilis_str_168.ASM904v1.dna.toplevel.fa \
/path/to/RSEM_Indices/Bacillus_subtilis_168/Bsubtilis-168

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
