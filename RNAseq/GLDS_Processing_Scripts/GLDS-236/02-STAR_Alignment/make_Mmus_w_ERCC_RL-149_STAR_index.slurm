#!/bin/bash
#
#SBATCH --job-name="make_Mmus_w_ERCC_RL-149_STAR_index"
#SBATCH --output=make_Mmus_w_ERCC_RL-149_STAR_index.out
#SBATCH --mem=60000
#
#SBATCH --mail-type=END

start=$(date +%s)
echo "start time: $start"
echo $HOSTNAME

STAR --runThreadN 15 \
--runMode genomeGenerate \
--limitGenomeGenerateRAM 40000000000 \
--genomeDir /path/to/STAR_Indices/Mus_musculus_w_ERCC_RL-149 \
--genomeFastaFiles /path/to/Genomes/Mus_musculus/Mus_musculus.GRCm38.dna.toplevel_and_ERCC92.fa \
--sjdbGTFfile /path/to/Genomes/Mus_musculus/Mus_musculus.GRCm38.96_and_ERCC92.gtf \
--sjdbOverhang 148

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
