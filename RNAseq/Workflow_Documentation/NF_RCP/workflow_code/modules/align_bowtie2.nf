/*
    Aligns reads against Bowtie2 index, uses SAMtools to convert SAM to BAM
    and filter out unmapped reads.
*/

process ALIGN_BOWTIE2 {
  tag "Sample: ${ meta.id }"

  input:
    tuple val(meta), path(reads)
    path(bowtie2_index_dir)

  output:
    tuple val(meta), path("${meta.id}.bam"), emit: bam
    path("${ meta.id}_*unmapped.fastq.gz"), emit: unmapped_reads
    path("${ meta.id }.bowtie2.log"), emit: alignment_logs
    path("versions.yml"), emit: versions

  script:
    def unaligned = meta.paired_end ? "--un-conc-gz ${meta.id}.unmapped.fastq.gz" : "--un-gz ${meta.id}.unmapped.fastq.gz"
    def readArgs = meta.paired_end ? "-1 ${reads[0]} -2 ${reads[1]}" : "-U ${reads}"

    """
    INDEX=\$(find -L ${bowtie2_index_dir} -name "*.rev.1.bt2" | sed "s/\\.rev.1.bt2\$//")

    bowtie2 -x "\$INDEX" \\
      -a \\
      ${readArgs} \\
      --threads ${task.cpus} \\
      --minins 10 \\
      --maxins 1000 \\
      ${unaligned} \\
      2> ${meta.id}.bowtie2.log | \\
      samtools view -bS -F 4 --threads ${task.cpus} - > ${meta.id}.bam

    # Rename unmapped reads
    if [ ${meta.paired_end} == true ]; then
      # For paired-end data
      mv "${meta.id}.unmapped.fastq.1.gz" "${meta.id}_R1_unmapped.fastq.gz"
      mv "${meta.id}.unmapped.fastq.2.gz" "${meta.id}_R2_unmapped.fastq.gz"
    else
      # For single-end data
      mv "${meta.id}.unmapped.fastq.gz" "${meta.id}_unmapped.fastq.gz"
    fi

    echo '"${task.process}":' > versions.yml
    echo "    bowtie2: \$(bowtie2 --version | head -n1 | awk '{print \$3}')" >> versions.yml
    echo "    samtools: \$(samtools --version | head -n1 | awk '{print \$2}')" >> versions.yml
    """
}