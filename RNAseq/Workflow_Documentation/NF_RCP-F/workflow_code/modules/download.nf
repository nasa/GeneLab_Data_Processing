/*
 * Processes for downloading reference files
 */

process DOWNLOAD_GUNZIP_REFERENCES {
  // Download and decompress genome and annotation files
  tag "Organism: ${ organism_sci }  Ensembl Version: ${ensemblVersion}"
  label 'networkBound'
  storeDir "${params.referenceStorePath}/${ref_source}/${ensemblVersion}/${ organism_sci }"

  input:
    tuple val(organism_sci), val(fasta_url), val(gtf_url)
    tuple val(ensemblVersion), val(ref_source)
  
  output:
    tuple path("*.fa"), path("*.gtf")

  script:
  """
  wget ${fasta_url}
  gunzip *.gz

  wget ${gtf_url}
  gunzip *.gz
  """
}

process DOWNLOAD_ERCC {
  label 'networkBound'
  storeDir "${params.referenceStorePath}/ERCC_thermofisher"

  input:
    val(has_ercc)

  output:
    tuple path("ERCC92.fa"), path("ERCC92.gtf")

  when:
    has_ercc

  script:
    """
    wget --no-check-certificate --quiet \
    -O ERCC92.zip \
    https://assets.thermofisher.com/TFS-Assets/LSG/manuals/ERCC92.zip  \
    && \
    unzip ERCC92.zip
    """
}
