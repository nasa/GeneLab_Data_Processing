#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
//params.paired = false
//params.max_mem = 100e9

/**************************************************************************************** 
**************************  Sequence assembly and summary *******************************
****************************************************************************************/

// This process handles running the assembly for each individual sample.
process ASSEMBLE {

    tag "Assembling ${sample_id}-s reads using megahit..."

    input:
        tuple val(sample_id), path(reads), val(isPaired)
    output:
        tuple val(sample_id), path("${sample_id}_final.contigs.fa"), emit: contigs
        path("${sample_id}-assembly.log"), emit: log
        path("versions.txt"), emit: version
    script:
        """
        # Removing output directory if exists already but process still needs to be 
        # run (because there is no --force option to megahit i dont't think):        
        [ -d ${sample_id}-megahit-out/ ] && rm -rf ${sample_id}-megahit-out/

        if [ ${isPaired} == true ]; then
       
            BASENAME_FORWARD=`basename -s '.gz' ${reads[0]}`
            BASENAME_REVERSE=`basename -s '.gz' ${reads[1]}` 

            zcat  ${reads[0]} > \${BASENAME_FORWARD}
            zcat  ${reads[1]} > \${BASENAME_REVERSE}

            megahit -1 \${BASENAME_FORWARD} -2 \${BASENAME_REVERSE} \\
               -m ${params.max_mem} -t ${task.cpus} \\
               --min-contig-len 500 -o ${sample_id}-megahit-out > ${sample_id}-assembly.log 2>&1
         
        else

            BASENAME=`basename -s '.gz' ${reads[0]}`
            zcat ${reads[0]} > \${BASENAME}
            megahit -r \${BASENAME} -m ${params.max_mem} -t  ${task.cpus} \\
                --min-contig-len 500 -o ${sample_id}-megahit-out > ${sample_id}-assembly.log 2>&1
        fi
        
        mv ${sample_id}-megahit-out/final.contigs.fa ${sample_id}_final.contigs.fa
        megahit -v > versions.txt
        """
}


process RENAME_HEADERS {

    tag "Renaming ${sample_id}-s assembly fasta file-s headers..."
    label "bit"
    label "assembly"

    input:
        tuple val(sample_id), path(assembly)
    output:
        tuple val(sample_id), path("${sample_id}-assembly.fasta"), emit: contigs
        path("versions.txt"), emit: version
    script:
        """
        bit-rename-fasta-headers -i ${assembly} \\
                                 -w c_${sample_id} \\
                                 -o ${sample_id}-assembly.fasta

        # Checking the assembly produced anything (megahit can run, produce 
        # the output fasta, but it will be empty if no contigs were assembled)
        if [ ! -s ${sample_id}-assembly.fasta ]; then
            printf "${sample_id}\\tNo contigs assembled\\n" > Failed-assemblies.tsv
        fi
        bit-version |grep "Bioinformatics Tools"|sed -E 's/^\\s+//' > versions.txt
        """
}


// This process summarizes and reports general stats for all individual sample assemblies in one table.
process SUMMARIZE_ASSEMBLIES {

    tag "Generating a summary of all the assemblies..."
    label "bit"
    label "assembly"

    input:
        path(assemblies)      
    output:
        path("${params.additional_filename_prefix}assembly-summaries${params.assay_suffix}.tsv"), emit: summary
        path("versions.txt"), emit: version
    script:
        """
        bit-summarize-assembly \\
                 -o ${params.additional_filename_prefix}assembly-summaries${params.assay_suffix}.tsv \\
                 ${assemblies}
        bit-version |grep "Bioinformatics Tools"|sed -E 's/^\\s+//' > versions.txt
        """
}
