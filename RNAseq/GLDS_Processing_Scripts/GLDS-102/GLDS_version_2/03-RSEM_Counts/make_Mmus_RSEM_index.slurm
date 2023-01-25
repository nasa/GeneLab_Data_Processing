#!/bin/bash
#
#SBATCH --job-name="make_Mmus_RSEM_index"
#SBATCH --output=make_Mmus_RSEM_index.out
#SBATCH --mem=60000
#
#SBATCH --mail-type=END

start=$(date +%s)
echo "start time: $start"
echo $HOSTNAME

rsem-prepare-reference --gtf /path/to/Genomes/Mus_musculus/Mus_musculus.GRCm38.96.gtf \
/path/to/Genomes/Mus_musculus/Mus_musculus.GRCm38.dna.toplevel.fa \
/path/to/RSEM_Indices/Mus_musculus/Mmus

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
