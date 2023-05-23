process RUNSHEET_FROM_GLDS {
  // Downloads isa Archive and creates runsheet using GeneLab API
  tag "${ gldsAccession }"
  publishDir "${ params.resultsDir }/Metadata",
    pattern: "*.zip",
    mode: params.publish_dir_mode

  input:
    val(osdAccession)
    val(gldsAccession)
    path("dp_tools__affymetrix")

  output:
    path("${ osdAccession }_microarray_v?_runsheet.csv"), emit: runsheet
    path("*.zip"), emit: isaArchive

  script:
    def injects = params.biomart_attribute ? "--inject biomart_attribute='${ params.biomart_attribute }'" : ''
    """

    dpt-get-isa-archive --accession ${ osdAccession }
    ls dp_tools__affymetrix

    dpt-isa-to-runsheet --accession ${ osdAccession } \
      --plugin-dir dp_tools__affymetrix \
      --isa-archive *.zip ${ injects }
    """
}