process ALIGN_BOWTIE2 {
  // Aligns reads against Bowtie2 index
  tag "Sample: ${ meta.id }"

  input:
    tuple val(meta), path(reads)
    path(bowtie2_index_dir)

  output:
    path("${ meta.id }/${ meta.id }*"), emit: publishables // used to ensure direct files are available for publishing directive
    path("${ meta.id }/${ meta.id }.bowtie2.log"), emit: alignment_logs
    tuple val(meta), path("${ meta.id }/${meta.id}.bam"), emit: bam
    path("${ meta.id }/${ meta.id }.unmapped.fastq"), emit: unmapped_reads
    path("versions.yml"), emit: versions

  script:
    def readArgs = meta.paired_end ? "-1 ${ reads[0] } -2 ${ reads[1] }" : "-U ${ reads }"

    """
    export BOWTIE2_INDEXES=${ bowtie2_index_dir }

    
    mkdir -p ${ meta.id }
    bowtie2 -x ${ BOWTIE2_INDEXES } \
    ${readArgs} \
    --threads ${ task.cpus } \
    --minins 0 \
    --maxins 500 \
    -k 1 \
    --un ${ meta.id }/${ meta.id }.unmapped.fastq \
    2> ${ meta.id }/${ meta.id }.bowtie2.log \
    | samtools view -bS --threads ${ task.cpus } -o ${ meta.id }/${ meta.id }.bam -

    echo '"${task.process}":' > versions.yml
    echo "    bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')" >> versions.yml
    """
}