#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/**************************************************************************************** 
*********************  Read mapping to contig assembly using Bowtie2 ********************
****************************************************************************************/

// This process builds the bowtie2 index and runs the mapping for each sample
process MAPPING {

    tag "Mapping ${sample_id}-s reads to its assembly ${assembly}..."
    label "mapping"

    input:
        tuple val(sample_id), path(assembly), path(reads), val(isPaired)
    output:
        tuple val(sample_id), path("${sample_id}.sam"), path("${sample_id}-mapping-info.txt"), emit: sam
        path("versions.txt"), emit: version
    script:
        """
        if [ ${isPaired}  == 'true' ]; then
            # Only running if the assembly produced anything
            if [ -s ${assembly} ]; then

                bowtie2-build ${assembly} ${sample_id}-index 
                bowtie2 --mm -q --threads ${task.cpus} \\
                        -x ${sample_id}-index  -1 ${reads[0]} -2 ${reads[1]} \\
                        --no-unal > ${sample_id}.sam  2> ${sample_id}-mapping-info.txt 
            rm ${sample_id}-index*
            else

                touch ${sample_id}.sam
                printf "Mapping not performed for ${sample_id} because the assembly didn't produce anything.\\n"

            fi
        # Single-end
        else

            # Only running if the assembly produced anything
            if [ -s ${assembly} ]; then

                bowtie2-build ${assembly} ${sample_id}-index 
                bowtie2 --mm -q --threads ${task.cpus} \\
                        -x ${sample_id}-index -r ${reads[0]} \\
                        --no-unal > ${sample_id}.sam  2> ${sample_id}-mapping-info.txt
                        
                rm ${sample_id}-index*
            else

                touch ${sample_id}.sam
                printf "Mapping not performed for ${sample_id} because the assembly didn't produce anything.\\n"

            fi

        fi
        bowtie2 --version  | head -n 1 | sed -E 's/.*(bowtie2-align-s version.+)/\\1/' > versions.txt
        """
}



// This process builds the bowtie2 index and runs the mapping for each sample
process SAM_TO_BAM {

    tag "Sorting and converting ${sample_id}-s sam to bam files..."
    label "mapping"

    input:
        tuple val(sample_id), path(sam), path(mapping_info)
    output:
        tuple val(sample_id), path("${sample_id}.bam"), emit: bam
        path("versions.txt"), emit: version
    script:
        """
        # Only running if the assembly produced anything
        if [ -s ${sam} ]; then

            samtools sort -@ ${task.cpus} ${sam} > ${sample_id}.bam 2> /dev/null

        else

            touch ${sample_id}.bam
            printf "Sorting and converting not performed for ${sample_id} because read mapping didn't produce anything.\\n"

        fi
        samtools --version | head -n1 > versions.txt
        """
}
