

log.info """\
    REMOVING HUMAN READS
    ===================================
    Download DB:        ${params.DL_kraken}
    Single end reads:   ${params.single_end}
    Use SampleID file:  ${params.specify_reads}
    Outputs:            ${params.human_db_path}
    """
    .stripIndent()

// Process to set up the human reads database
process set_up_human_db {
    tag "Downloading human reads database to ${params.human_db_path}\n"
    publishDir path: "$projectDir" 

    output:
    path "${params.human_db_name}/"

    script:
    """
    curl -L -o ${params.human_db_name}.tar.gz https://ndownloader.figshare.com/files/25627058

    tar -xzvf ${params.human_db_name}.tar.gz
    rm ${params.human_db_name}.tar.gz

    """

}

// Process for paired-end (PE) read analysis with Kraken2
process PE_kraken2 {

    container params.kraken2container
    tag "$sample_id"
    publishDir "$params.kraken_output_dir", pattern: "*.{txt,tsv}"
    publishDir "$params.kraken_output_dir/reads", pattern: "*.fastq.gz"

    input:
    path database 
    tuple val(sample_id), path(reads_ch)
    

    output:
    path "${sample_id}-kraken2-output.txt"
    path "${sample_id}-kraken2-report.tsv"
    path "${sample_id}*.gz"

    script:
    """
    kraken2 --db $database --gzip-compressed \
            --threads 2 --use-names --paired  \
            --output ${sample_id}-kraken2-output.txt \
            --report ${sample_id}-kraken2-report.tsv \
            --unclassified-out "${sample_id}${params.PE_reads_out_suffix}" \
            ${reads_ch[0]} ${reads_ch[1]}

    gzip ${sample_id}*.fastq
    """
}

// Process for single-end (SE) read analysis with Kraken2
process SE_kraken2 {
    
    container params.kraken2container 
    tag "$sample_id"
    publishDir "$params.kraken_output_dir", pattern: "*.{txt,tsv}"
    publishDir "$params.kraken_output_dir/reads", pattern: "*.fastq.gz"

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
            --unclassified-out "${sample_id}${params.SE_reads_out_suffix}" \
            ${reads_ch[0]}

    gzip ${sample_id}${params.SE_reads_out_suffix}
    """
}


// Main workflow logic
workflow {

// Conditionally download the human reads database if requested
if(params.DL_kraken == true){
    log.info "\nPreparing to download new human reads database"
    database_ch =  set_up_human_db()
    database_ch.view{"database path: ${it}"}
}

else {
    log.info "\nAccessing previous human reads database"
    database_ch = Channel.value(params.human_db_path)
    database_ch.view{"database path: ${it}"}
}

// Process reads based on whether they are single-end or paired-end
if(params.single_end == true) {
    log.info "\nReading Single-end data from ${params.reads_dir}\n"
    
    // Specified reads handling (mimics Channel.fromFilePairs() method's output, but with SE)
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

    // Specified reads handling (mimics Channel.fromFilePairs() method's output)
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

// Calculate and log the final result
final_percent = output_ch[1]
    .collect{(it.text[0..5]).toFloat()}
    .average().trunc(2)
    .view{"\nRESULT: ${it}% of input reads were unclassified, available in ${params.kraken_output_dir}/reads "}
}