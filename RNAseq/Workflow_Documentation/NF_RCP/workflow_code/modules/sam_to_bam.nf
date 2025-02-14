process SAM_TO_BAM {
  // Convert SAM to BAM, keeping only mapped reads, without sorting
  tag "Sample: ${ meta.id }"

  input:
    tuple val(meta), path(sam)

  output:
    tuple val(meta), path("${meta.id}.bam"), emit: bam
    path("versions.yml"), emit: versions

  script:
    """
    # Convert SAM to BAM
    samtools view -bS -F 4 --threads ${task.cpus} ${sam} > ${meta.id}.bam

    # Capture versions
    echo '"${task.process}":' > versions.yml
    echo "    samtools: \$(samtools --version | head -n1 | awk '{print \$2}')" >> versions.yml
    """
}