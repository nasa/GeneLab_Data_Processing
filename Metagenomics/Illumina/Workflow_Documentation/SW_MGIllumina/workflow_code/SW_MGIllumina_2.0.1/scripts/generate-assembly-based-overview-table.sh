#!/usr/bin/env bash

sample_IDs_file=${1}
assemblies_dir=${2}
genes_dir=${3}
mapping_dir=${4}
bins_dir=${5}
MAGs_dir=${6}
output=${7}

# starting output file
printf "Sample_ID\tassembly_produced\tgene_calls_identified\tread_mapping_successful\tbins_recovered\tMAGs_recovered\n" > ${output}

# looping through all input files and generating columns for final table
for sample in $(cat ${sample_IDs_file})
do

    # checking assembly
    if [ ! -s ${assemblies_dir}/${sample}-assembly.fasta ]; then
        printf "No\n" >> assembly-status.tmp

        # removing empty output fasta
        rm -rf ${assemblies_dir}/${sample}-assembly.fasta

    else
        printf "Yes\n" >> assembly-status.tmp
    fi

    # checking gene calls
    if [ ! -s ${genes_dir}/${sample}-genes.faa ]; then
        printf "No\n" >> genes-status.tmp

        # removing empty output files
        rm -rf ${genes_dir}/${sample}-genes.faa ${genes_dir}/${sample}-genes.fasta ${genes_dir}/${sample}-genes.gff

    else
        printf "Yes\n" >> genes-status.tmp
    fi

    # checking read-mapping outputs
    if [ ! -s ${mapping_dir}/${sample}.bam ]; then
        printf "No\n" >> mapping-status.tmp

        # removing empty output files
        rm -rf ${mapping_dir}/${sample}.bam ${mapping_dir}/${sample}-metabat-assembly-depth.tsv

    else
        printf "Yes\n" >> mapping-status.tmp
    fi

    # getting number of bins recovered if any produced
    if compgen -G "${bins_dir}*.fasta" > /dev/null; then
        num_bins=$(ls ${bins_dir}*.fasta | grep -c "${sample}-bin.[0-9]*.fasta")
        printf "${num_bins}\n" >> bins-status.tmp
    else
        printf "0\n" >> bins-status.tmp
    fi

    # getting number of MAGs recovered
    if compgen -G "${MAGs_dir}*.fasta" >/dev/null; then
        num_MAGs=$(ls ${MAGs_dir}*.fasta | grep -c "${sample}-MAG-[0-9]*.fasta")
        printf "${num_MAGs}\n" >> MAGs-status.tmp
    else
        printf "0\n" >> MAGs-status.tmp
    fi

done

# combining, adding to output file and removing intermediates
cat <( paste ${sample_IDs_file} assembly-status.tmp \
             genes-status.tmp mapping-status.tmp \
             bins-status.tmp MAGs-status.tmp ) >> ${output}

rm assembly-status.tmp genes-status.tmp mapping-status.tmp bins-status.tmp MAGs-status.tmp
