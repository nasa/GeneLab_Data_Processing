#!/bin/bash
#
#SBATCH --job-name="make_Athaliana_RL-75_STAR_index"
#SBATCH --output=make_Athaliana_RL-75_STAR_index.out
#SBATCH --mem=60000
#

start=$(date +%s)
echo "start time: $start"
echo $HOSTNAME

STAR --runThreadN 15 \
--runMode genomeGenerate \
--limitGenomeGenerateRAM 40000000000 \
--genomeDir /path/to/STAR_Indices/Arabidopsis_thaliana_RL-75 \
--genomeFastaFiles /path/to/Genomes/Arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa \
--sjdbGTFfile /path/to/Genomes/Arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.44.gtf \
--sjdbOverhang 74

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
