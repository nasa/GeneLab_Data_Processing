process DOWNLOAD_REFERENCES {
  // Download and decompress genome and gtf files
  tag "Organism: ${organism_sci}, Reference Source: ${reference_source}${reference_source.toLowerCase().contains('ensembl') ? ', Reference Version: ' + reference_version : ''}"
  storeDir "${reference_store_path}/${reference_source}/${reference_source.toLowerCase().contains('ensembl') ? reference_version + '/' : ''}${organism_sci}"

  input:
    val(reference_store_path)
    val(organism_sci)
    val(reference_source)
    val(reference_version)
    val(fasta_url)
    val(gtf_url)
  
  output:
    tuple path("{*.fa,*.fna}"), path("*.gtf")       , emit: reference_files

  script:
  """
  wget ${fasta_url}
  gunzip *.gz

  wget ${gtf_url}
  gunzip *.gz
  """
}