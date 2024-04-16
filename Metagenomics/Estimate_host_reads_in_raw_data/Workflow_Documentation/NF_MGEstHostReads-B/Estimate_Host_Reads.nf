// Initial logging of the workflow's parameters for tracking and debug purposes
log.info """\
    ESTIMATE HOST READS
    ===================================
    Download DB:        ${params.DL_kraken}
    Single end reads:   ${params.single_end}
    projectDir:         ${projectDir}"""
    .stripIndent()


// Process for paired-end reads using Kraken2
process PE_kraken2 {
    container params.kraken2container 
    tag "$sample_id"
    publishDir "$params.kraken_output_dir", pattern: "*.{txt,tsv}"

    input:
    path database 
    tuple val(sample_id), path(reads_ch)
    

    output:
    path "${sample_id}-kraken2-output.txt"
    path "${sample_id}-kraken2-report.tsv"

    script:
    """
    kraken2 --db $database --gzip-compressed \
            --threads 2 --use-names --paired  \
            --output ${sample_id}-kraken2-output.txt \
            --report ${sample_id}-kraken2-report.tsv \
            ${reads_ch[0]} ${reads_ch[1]}

    """
}

// Process for single-end reads using Kraken2
process SE_kraken2 {
    
    container params.kraken2container 
    tag "$sample_id"
    publishDir "$params.kraken_output_dir", pattern: "*.{txt,tsv}"

    input:
    path database 
    tuple val(sample_id), path(reads_ch)

    output:
    path "${sample_id}-kraken2-output.txt"
    path "${sample_id}-kraken2-report.tsv"
    path "${sample_id}${params.SE_reads_out_suffix}.gz"

    script:
    """
    kraken2 --db $database --gzip-compressed --threads 2 --use-names  \
            --output ${sample_id}-kraken2-output.txt \
            --report ${sample_id}-kraken2-report.tsv \
            ${reads_ch[0]}

    """
}



workflow {


    // Log the database path being used
    log.info "\nAccessing previous host reads database"
    database_ch = Channel.value(params.host_db_path)
    database_ch.view{"database path: ${it}"}

    // Conditional execution for single-end or paired-end data
    if(params.single_end == true) {
        log.info "\nReading Single-end data from ${params.reads_dir}\n"

        if (params.specify_reads) {
            reads_ch = Channel
            .fromPath("${params.sample_id_list}")
            .splitText()
            .map { it.trim() }
            .map { sample_id ->
                def files = file("${params.reads_dir}${sample_id}${params.SE_reads_suffix}")
                return [sample_id, files]
            }
        }
        else {
        reads_ch = Channel
            .fromPath("${params.reads_dir}/*${params.SE_reads_suffix}", checkIfExists: true)
            .map { readfile ->
                def sampleId = readfile.name.replaceAll("${params.SE_reads_suffix}\$", "")
                return tuple(sampleId, readfile)
            }
        }
        reads_ch.view{"reads: ${it}"}
        output_ch = SE_kraken2(database_ch, reads_ch)
        
    }
    else {
        log.info "\nReading Paired-end data from ${params.reads_dir}\n"
        // Load specific reads if specified
        if (params.specify_reads) {
            reads_ch = Channel
            .fromPath("${params.sample_id_list}")
            .splitText()
            .map { it.trim() }
            .map { sample_id ->
                def files = file("${params.reads_dir}${sample_id}${params.PE_reads_suffix}").toList().sort()
                return [sample_id, files]
            }
        }
        else {
            reads_ch = Channel.fromFilePairs(params.reads_dir + "*" + params.PE_reads_suffix, checkIfExists: true)
        }
        reads_ch.view{"reads: ${it}"}
        output_ch = PE_kraken2(database_ch, reads_ch)
    }
    // Calculate and log the final percentage of unclassified reads
    final_percent = output_ch[1]
        .collect{(it.text[0..5]).toFloat()}
        .average().trunc(2)
        .view{"\nRESULT: ${it}% of input reads were unclassified, available in ${params.kraken_output_dir}/reads "}

}
