#!/bin/bash
#
#SBATCH --job-name="generate_runsheet"
#SBATCH --output=generate_runsheet.out
#SBATCH --partition=priority
#SBATCH --mem=20000
#
#SBATCH --mail-user=user@nasa.gov
#SBATCH --mail-type=END

. ~/.profile


echo "generate_runsheet"
echo ""

start=$(date +%s)
echo "start time: $start"
echo $HOSTNAME


echo ""
echo "dp_tools version: "
echo ""
singularity run --bind $PWD:$PWD docker://quay.io/j_81/dp_tools:1.1.6 python -c "import dp_tools; print(dp_tools.__version__)"
echo ""


echo ""

call1="singularity run --bind $PWD:$PWD docker://quay.io/j_81/dp_tools:1.1.6 dpt-get-isa-archive \
	--accession GLDS-120"


echo $call1
echo ""
eval $call1


call2="singularity run --bind $PWD:$PWD docker://quay.io/j_81/dp_tools:1.1.6 dpt-isa-to-runsheet --accession GLDS-120 \
	--config-type bulkRNASeq \
	--config-version Latest \
	--isa-archive GLDS-120_metadata_GLDS-120-ISA.zip"


echo $call2
echo ""
eval $call2

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
