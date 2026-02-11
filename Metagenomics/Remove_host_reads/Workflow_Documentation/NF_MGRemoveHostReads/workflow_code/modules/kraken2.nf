process KRAKEN_2 {

    tag "${meta.id}"
    publishDir "${params.reads_outdir}", pattern: '*.fastq.gz'
    publishDir "${params.outdir}/results/kraken2-output", pattern: '*.{txt,tsv}'

    input:
    path database
    tuple val(meta), path(reads)
    
    output:
    path "${meta.id}*_HRremoved_*.gz", emit: host_removed
    path "${meta.id}-kraken2-output.txt", emit: output
    path "${meta.id}-kraken2-report.tsv", emit: report
    path("versions.txt"), emit: version

    script:
    def pe_flag = meta.paired_end ? "--paired" : ""
    def input   = meta.paired_end ? "${reads[0]} ${reads[1]}" : "${reads}"
    def output  = meta.paired_end ? "${meta.id}_R#.fastq" : "${meta.id}.fastq"
    """
    kraken2 --db $database --gzip-compressed \
            --threads ${task.cpus} --use-names ${pe_flag} \
            --output ${meta.id}-kraken2-output.txt \
            --report ${meta.id}-kraken2-report.tsv \
            --unclassified-out ${output} \
            ${input}

    # Compress intermediate FASTQ files
    gzip ${input}

    # Rename compressed files to final output names
    mv ${meta.id}_R_1.fastq ${meta.id}${params.R1_out_suffix}
    mv ${meta.id}_R_2.fastq ${meta.id}${params.R2_out_suffix}

    echo "Kraken2 \$(kraken2 -version | head -n 1 | awk '{print \$3}')" >> versions.txt            
    """
}