#!/usr/bin/env bash

sample_IDs_file=${1}
host=${2}
input_dir=${3}
output_file=${4}

# starting output file
printf "Sample_ID\tTotal_fragments_before\tTotal_fragments_after\tPercent_${host}_reads_removed\n" > ${output_file}

# looping through all input files and generating columns for final table
for sample in $(cat ${sample_IDs_file})
do
    cat ${input_dir}/${sample}-removal-info.tmp >> ${output_file}
done