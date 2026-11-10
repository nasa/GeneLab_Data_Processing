process COPY_READS {
    tag "${ meta.id }"

    input:
        tuple val(meta), path("?.gz")

    output:
        tuple val(meta), path("${meta.id}*.gz"), emit: raw_reads

    script:
        if ( meta.paired_end ) {
        """
        cp -P 1.gz ${meta.id}${params.R1_suffix}
        cp -P 2.gz ${meta.id}${params.R2_suffix}
        """
        } else {
        """
        cp -P 1.gz ${meta.id}${params.single_suffix}
        """
        }
}
