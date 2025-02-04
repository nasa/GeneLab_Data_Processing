process SAMTOOLS_STATS {
    // Reports samtools stats for a sample
    tag "Sample: ${ meta.id }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple path(genome_fasta), path(genome_gtf)

    output:
    tuple val(meta), path("*.stats"), emit: stats
    path  "versions.yml"            , emit: versions

    script:
    """
    samtools \
        stats \
        --threads ${ task.cpus } \
        --reference ${ genome_fasta } \
        --most-inserts 1.0 \
        ${ bam } \
        > ${ meta.id }.stats

    # Capture versions
    echo '"${ task.process }":' > versions.yml
    echo "    samtools: \$(samtools --version | head -n1 | sed 's/^.*samtools //; s/Using.*\$//')" >> versions.yml
    """
}