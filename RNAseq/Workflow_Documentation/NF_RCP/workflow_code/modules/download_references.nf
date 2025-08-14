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
    tuple path("{*.fa,*.fna}"), path("*.gtf"), emit: reference_files

  script:
  """
  # Create temp directories for processing
  mkdir -p temp_fasta temp_gtf

  # Handle fasta file
  if [[ "${fasta_url}" == http* ]]; then
    wget --directory-prefix temp_fasta "${fasta_url}"
  else
    # It's a file path - just copy it directly to main directory
    cp "${fasta_url}" temp_fasta
  fi
  
  if [[ "${gtf_url}" == http* ]]; then
    wget --directory-prefix temp_gtf "${gtf_url}"
  else
    # It's a file path - just copy it directly to main directory
    cp "${gtf_url}" temp_gtf
  fi
    
  # Handle decompression if needed
  if ls temp_fasta/*.gz &> /dev/null; then
    gunzip temp_fasta/*.gz
  fi
  # Move processed files to main directory
  mv temp_fasta/*.fa temp_fasta/*.fna ./ 2>/dev/null || true

  if ls temp_gtf/*.gz &> /dev/null; then
    gunzip temp_gtf/*.gz
  fi
  # Move processed files to main directory
  mv temp_gtf/*.gtf ./ 2>/dev/null || true
  
  # Clean up temp directories
  rm -rf temp_fasta temp_gtf
  """
}