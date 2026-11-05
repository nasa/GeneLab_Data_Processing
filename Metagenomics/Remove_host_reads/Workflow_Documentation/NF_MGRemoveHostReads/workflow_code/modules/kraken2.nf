process KRAKEN_2 {

    tag "${meta.id}"

    input:
    path database
    tuple val(meta), path(reads)
    val out_suffix
    
    output:
    path "${meta.id}*${out_suffix}.fastq.gz", emit: host_removed
    path("${meta.id}-kraken2-output.txt"), emit: output
    path("${meta.id}-kraken2-report.tsv"), emit: report
    path("versions.txt"), emit: version

    script:
    if (meta.paired_end)
        """
        kraken2 --db $database --gzip-compressed \
            --threads ${task.cpus} --use-names --paired \
            --output ${meta.id}-kraken2-output.txt \
            --report ${meta.id}-kraken2-report.tsv \
            --unclassified-out ${meta.id}_R#.fastq \
            ${reads[0]} ${reads[1]}

        # Rename and compress files to final output names
        mv "${meta.id}_R_1.fastq" "${meta.id}_R1${out_suffix}.fastq" && \
        gzip ${meta.id}_R1${out_suffix}.fastq
        
        mv "${meta.id}_R_2.fastq" "${meta.id}_R2${out_suffix}.fastq" && \
        gzip ${meta.id}_R2${out_suffix}.fastq

        echo "Kraken2 \$(kraken2 -version | head -n 1 | awk '{print \$3}')" >> versions.txt
        """
    else
        """
        kraken2 --db $database --gzip-compressed \
                --threads ${task.cpus} --use-names \
                --output ${meta.id}-kraken2-output.txt \
                --report ${meta.id}-kraken2-report.tsv \
                --unclassified-out ${meta.id}.fastq \
                ${reads}

        # Rename and compress files to final output names
        test -f "${meta.id}.fastq" && mv "${meta.id}.fastq" "${meta.id}${out_suffix}.fastq" && \
        gzip ${meta.id}${out_suffix}.fastq

        echo "Kraken2 \$(kraken2 -version | head -n 1 | awk '{print \$3}')" >> versions.txt            
        """
}