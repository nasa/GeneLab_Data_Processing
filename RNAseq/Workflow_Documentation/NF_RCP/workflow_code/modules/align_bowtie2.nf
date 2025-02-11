process ALIGN_BOWTIE2 {
  // Aligns reads against Bowtie2 index
  tag "Sample: ${ meta.id }"

  input:
    tuple val(meta), path(reads)
    path(bowtie2_index_dir)

  output:
    tuple val(meta), path("${meta.id}.sam"), emit: sam
    path("${ meta.id }.Unmapped.fastq*"), emit: unmapped_reads, optional: true
    path("${ meta.id }.bowtie2.log"), emit: alignment_logs
    path("versions.yml"), emit: versions

  script:
    def unaligned = meta.paired_end ? "--un-conc-gz ${meta.id}.Unmapped.fastq.gz" : "--un-gz ${meta.id}.Unmapped.fastq.gz"
    def readArgs = meta.paired_end ? "-1 ${reads[0]} -2 ${reads[1]}" : "-U ${reads}"

    """
    INDEX=\$(find -L ${bowtie2_index_dir} -name "*.rev.1.bt2" | sed "s/\\.rev.1.bt2\$//")

    bowtie2 -x "\$INDEX" \\
      ${readArgs} \\
      --threads ${task.cpus} \\
      --minins 0 \\
      --maxins 500 \\
      --no-unal \\
      ${unaligned} \\ 
      -S ${meta.id}.sam \\
      2> ${meta.id}.bowtie2.log

    echo '"${task.process}":' > versions.yml
    echo "    bowtie2: \$(bowtie2 --version | head -n1 | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')" >> versions.yml
    """
}