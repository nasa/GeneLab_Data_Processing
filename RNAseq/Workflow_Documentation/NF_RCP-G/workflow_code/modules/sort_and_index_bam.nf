process SORT_AND_INDEX_BAM {
  tag "Sample: ${ meta.id }"

  input:
    tuple val(meta), path(bam_file)

  output:
    tuple val(meta), path(sorted_bam_fname), path("${ sorted_bam_fname }.bai"), emit: sorted_bam
    tuple path(sorted_bam_fname), path("${ sorted_bam_fname }.bai"), emit: bam_only_files
    //path("versions.yaml"), emit: versions

  script:
    sorted_bam_fname = bam_file.name.contains('.out.bam') ? 
                   bam_file.name.replaceAll('.out.bam', '_sorted.out.bam') : 
                   bam_file.name.replaceAll('.bam', '_sorted.bam')
    mem_MB_per_thread = (task.memory.toMega().intValue() * 0.8 / task.cpus).intValue()
    """    
    samtools sort -m ${ mem_MB_per_thread }M \
                  --threads ${ task.cpus } \
                  -o ${ sorted_bam_fname } \
                  ${ bam_file }

    samtools index -@ ${ task.cpus } ${ sorted_bam_fname }

    #echo '"${task.process}":' > versions.yml
    #echo "    samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')" >> versions.yml
    """
}