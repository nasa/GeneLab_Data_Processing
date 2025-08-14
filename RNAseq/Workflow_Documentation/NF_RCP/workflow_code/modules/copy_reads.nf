process COPY_READS {
    tag "Sample: ${ meta.id }"

    input:
        tuple val(meta), path("?.gz")

    output:
        tuple val(meta), path("${meta.id}*.gz"), emit: raw_reads

    script:
        if ( meta.paired_end ) {
        """
        cp -P 1.gz ${meta.id}_R1_raw.fastq.gz
        cp -P 2.gz ${meta.id}_R2_raw.fastq.gz
        """
        } else {
        """
        cp -P 1.gz  ${meta.id}_raw.fastq.gz
        """
        }
}