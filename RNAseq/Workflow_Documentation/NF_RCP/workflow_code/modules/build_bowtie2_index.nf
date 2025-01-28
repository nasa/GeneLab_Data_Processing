process BUILD_BOWTIE2_INDEX {
  // Builds Bowtie 2 index, this is ercc-spike-in and organism specific
  tag "Refs: ${ genome_fasta.baseName }, ${ genome_gtf.baseName }, Source: ${reference_source}${reference_source.toLowerCase().contains('ensembl') ? ', Version: ' + reference_version : ''}${params.genome_subsample ? ', GenomeSubsample: ' + params.genome_subsample : ''}"
  storeDir "${ derived_store_path }/Bowtie2_Indices/${ reference_source }/${reference_source.toLowerCase().contains('ensembl') ? reference_version + '/' : ''}${ meta.organism_sci }"

  input:
    val(derived_store_path)
    val(organism_sci)
    val(reference_source)
    val(reference_version)
    tuple path(genome_fasta), path(genome_gtf)
    val(meta)


  output:
    path("${ genome_fasta.baseName }"), emit: index_dir
  script:
    """
    mkdir -p ${ genome_fasta.baseName }

    bowtie2-build --threads ${task.cpus} \
      -f ${ genome_fasta } \
      ${ genome_fasta.baseName }/${ genome_fasta.baseName }
    """
}