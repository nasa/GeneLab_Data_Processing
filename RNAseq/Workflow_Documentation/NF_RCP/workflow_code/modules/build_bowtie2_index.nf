process BUILD_BOWTIE2_INDEX {
  // Builds Bowtie 2 index, this is ERCC-spike-in and organism-specific
  tag "Refs: ${ genome_fasta.baseName }, ${ genome_gtf.baseName }, Source: ${reference_source}${reference_source.toLowerCase().contains('ensembl') ? ', Version: ' + reference_version : ''}${params.genome_subsample ? ', GenomeSubsample: ' + params.genome_subsample : ''}"
  //storeDir "${ derived_store_path }/Bowtie2_Indices/${ reference_source }/${reference_source.toLowerCase().contains('ensembl') ? reference_version + '/' : ''}${ meta.organism_sci }"

  input:
    val(derived_store_path)
    val(organism_sci)
    val(reference_source)
    val(reference_version)
    tuple path(genome_fasta), path(genome_gtf)
    val(meta)

  output:
    path("bowtie2/"), emit: index_dir

  script:
    """
    mkdir -p bowtie2

    bowtie2-build --threads ${task.cpus} \
      --seed 1 \
      -f ${ genome_fasta } \
      bowtie2/${ genome_fasta.baseName }

    echo '"${task.process}":' > versions.yml
    echo "    bowtie2: \$(bowtie2 --version | head -n1 | awk '{print \$3}')" >> versions.yml
    """
}