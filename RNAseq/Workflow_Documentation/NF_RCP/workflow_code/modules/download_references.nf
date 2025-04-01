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
    cd temp_fasta
    wget "${fasta_url}"
    
    # Handle decompression if needed
    if ls *.gz &> /dev/null; then
      gunzip *.gz
    fi
    # Move processed files to main directory
    mv *.fa *.fna ../ 2>/dev/null || true
    cd ..
  else
    # It's a file path - just copy it directly to main directory
    cp "${fasta_url}" ./
  fi

  # Handle gtf file
  if [[ "${gtf_url}" == http* ]]; then
    cd temp_gtf
    wget "${gtf_url}"
    
    # Handle decompression if needed
    if ls *.gz &> /dev/null; then
      gunzip *.gz
    fi
    # Move processed files to main directory
    mv *.gtf ../ 2>/dev/null || true
    cd ..
  else
    # It's a file path - just copy it directly to main directory
    cp "${gtf_url}" ./
  fi

  # Clean up temp directories
  rm -rf temp_fasta temp_gtf
  """
}