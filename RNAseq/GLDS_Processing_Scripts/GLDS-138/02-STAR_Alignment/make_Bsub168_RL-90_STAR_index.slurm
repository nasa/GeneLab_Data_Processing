#!/bin/bash
#
#SBATCH --job-name="make_Bsub168_RL-90_STAR_index"
#SBATCH --output=make_Bsub168_RL-90_STAR_index.out
#SBATCH --mem=74000
#
#SBATCH --mail-type=END

start=$(date +%s)
echo "start time: $start"
echo $HOSTNAME

STAR --runThreadN 12 \
--runMode genomeGenerate \
--limitGenomeGenerateRAM 72000000000 \
--genomeSAindexNbases 10 \
--genomeDir /path/to/STAR_Indices/Bacillus_subtilis_168_RL-90/ \
--genomeFastaFiles /path/to/Genomes/Bacillus_subtilis/Bacillus_subtilis_subsp_subtilis_str_168.ASM904v1.dna.toplevel.fa \
--sjdbGTFfile /path/to/Genomes/Bacillus_subtilis/Bacillus_subtilis_subsp_subtilis_str_168.ASM904v1.46.gtf \
--sjdbOverhang 89

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
