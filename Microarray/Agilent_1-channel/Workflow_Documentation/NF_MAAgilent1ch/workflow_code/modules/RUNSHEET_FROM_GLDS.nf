process RUNSHEET_FROM_GLDS {
  // Downloads isa Archive and creates runsheet using GeneLab API
  tag "${ gldsAccession }"
  publishDir "${ params.outputDir }/${ gldsAccession }/Metadata",
    mode: params.publish_dir_mode

  input:
    val(osdAccession)
    val(gldsAccession)

  output:
    path("${ gldsAccession }_microarray_v?_runsheet.csv"), emit: runsheet
    path("*.zip"), emit: isaArchive

  script:
    def injects = params.biomart_attribute ? "--inject biomart_attribute=${ params.biomart_attribute }" : ''
    """
    dpt-get-isa-archive --accession ${ osdAccession }

    dpt-isa-to-runsheet --accession ${ gldsAccession } \
      --config-type microarray --isa-archive *.zip ${ injects }
    """
}