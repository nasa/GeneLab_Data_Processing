process COUNT_ALIGNED {
  // Generates gene and isoform counts from alignments
  tag "Sample: ${ meta.id }, Strandedness: ${ strandedness } "

  input:
    tuple val(meta), path("${meta.id}_Aligned.toTranscriptome.out.bam")
    path(RSEM_REF)
    val(strandedness)

  output:
    tuple val(meta), path("${ meta.id }*"), emit: counts
    path("${ meta.id }*"), emit: only_counts
    path("${ meta.id }.genes.results"), emit: genes_results
    path("versions.yml"), emit: versions

  script:
    strandedness_opt_map = ["sense":"forward","antisense":"reverse","unstranded":"none"]
    // e.g. 'mus_musculus' should become 'Mmus'
    organism_str = "${ meta.organism_sci.substring(0,1).toUpperCase() }${ meta.organism_sci.split('_')[1].substring(0,3).toLowerCase() }"
    """
    rsem-calculate-expression --num-threads $task.cpus \
      ${ meta.paired_end ? '--paired-end' : '' } \
      --bam \
      --alignments \
      --no-bam-output \
      --estimate-rspd \
	  --seed-length 20 \
      --seed 12345 \
      --strandedness ${ strandedness_opt_map.get(strandedness) } \
      ${meta.id}_Aligned.toTranscriptome.out.bam \
      ${ RSEM_REF }/${ organism_str } \
      ${ meta.id }

    echo '"${task.process}":' > versions.yml
    # RSEM version reporting is broken in RSEM 1.3.3...
    #echo "    rsem: \$(rsem-calculate-expression --version | sed 's/Current version: RSEM v//')" >> versions.yml
    echo "    rsem: 1.3.3" >> versions.yml
    """
}