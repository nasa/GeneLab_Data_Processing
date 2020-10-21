#!/bin/bash
#
#SBATCH --job-name="make_Mmus_RL-50_STAR_index"
#SBATCH --output=make_Mmus_RL-50_STAR_index.out
#SBATCH --mem=60000
#
#SBATCH --mail-type=END

start=$(date +%s)
echo "start time: $start"
echo $HOSTNAME

STAR --runThreadN 15 \
--runMode genomeGenerate \
--limitGenomeGenerateRAM 40000000000 \
--genomeDir /path/to/STAR_Indices/Mus_musculus_RL-50 \
--genomeFastaFiles /path/to/Genomes/Mus_musculus/Mus_musculus.GRCm38.dna.toplevel.fa \
--sjdbGTFfile /path/to/Genomes/Mus_musculus/Mus_musculus.GRCm38.96.gtf \
--sjdbOverhang 49

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
