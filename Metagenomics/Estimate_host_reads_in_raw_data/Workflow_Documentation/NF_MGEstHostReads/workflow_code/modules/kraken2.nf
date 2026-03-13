process KRAKEN_2 {

    tag "${meta.id}"

    input:
    path database
    tuple val(meta), path(reads)
    
    output:
    path "${meta.id}-kraken2-output.txt", emit: output
    path "${meta.id}-kraken2-report.tsv", emit: report
    path("versions.txt"), emit: version

    script:
    def input = meta.paired_end ? "--paired ${reads[0]} ${reads[1]}" : "${reads}"
    """
    kraken2 --db $database --gzip-compressed \
            --threads ${task.cpus} --use-names \
            --output ${meta.id}-kraken2-output.txt \
            --report ${meta.id}-kraken2-report.tsv \
            ${input}

    echo "Kraken2 \$(kraken2 -version | head -n 1 | awk '{print \$3}')" >> versions.txt      
    """
}