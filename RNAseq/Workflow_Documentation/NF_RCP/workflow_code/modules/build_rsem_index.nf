process BUILD_RSEM_INDEX {
  // Builds RSEM index, this is ercc-spike-in, organism, and ensembl version specific
  tag "Refs: ${ genome_fasta }, ${ genome_gtf }, Source: ${reference_source}${reference_source.toLowerCase().contains('ensembl') ? ', Version: ' + reference_version : ''}${ params.genome_subsample ? ', GenomeSubsample: ' + params.genome_subsample : ''}"
  storeDir "${ derived_store_path }/RSEM_Indices/${ reference_source }/${reference_source.toLowerCase().contains('ensembl') ? reference_version + '/' : ''}${ meta.organism_sci }"

  input:
    val(derived_store_path)
    val(organism_sci)
    val(reference_source)
    val(reference_version)
    tuple path(genome_fasta), path(genome_gtf)
    val(meta)

  output:
    path("${ genome_fasta.baseName }"), emit: index_dir
    path("${ genome_fasta.baseName }/${ organism_str }.grp") // to ensure check expected file contents exist
    //path("versions.yaml"), emit: versions


  script:
    // e.g. 'mus_musculus' should become 'Mmus'
    organism_str = "${ meta.organism_sci.substring(0,1).toUpperCase() }${ meta.organism_sci.split('_')[1].substring(0,3).toLowerCase() }"
    """
    mkdir  ${ genome_fasta.baseName }
    rsem-prepare-reference --gtf ${ genome_gtf } ${ genome_fasta } ${ genome_fasta.baseName }/${ organism_str }

    #echo '"${task.process}":' > versions.yml
    #echo "    rsem: \$(rsem-calculate-expression --version | sed -e 's/Current version: RSEM v//g')" >> versions.yml
    """

}