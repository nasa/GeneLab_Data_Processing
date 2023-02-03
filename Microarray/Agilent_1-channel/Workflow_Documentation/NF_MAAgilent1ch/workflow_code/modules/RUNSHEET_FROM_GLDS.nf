process RUNSHEET_FROM_GLDS {
  // Downloads isa Archive and creates runsheet using GeneLab API
  tag "${ gldsAccession }"
  publishDir "${ params.outputDir }/${ gldsAccession }/Metadata",
    pattern: "*.zip",
    mode: params.publish_dir_mode

  input:
    val(osdAccession)
    val(gldsAccession)
    path("dp_tools__agilent_1_channel")

  output:
    path("${ gldsAccession }_microarray_v?_runsheet.csv"), emit: runsheet
    path("*.zip"), emit: isaArchive

  script:
    def injects = params.biomart_attribute ? "--inject biomart_attribute='${ params.biomart_attribute }'" : ''
    """

    dpt-get-isa-archive --accession ${ osdAccession }
    ls dp_tools__agilent_1_channel

    dpt-isa-to-runsheet --accession ${ gldsAccession } \
      --plugin-dir dp_tools__agilent_1_channel \
      --isa-archive *.zip ${ injects }
    """
}