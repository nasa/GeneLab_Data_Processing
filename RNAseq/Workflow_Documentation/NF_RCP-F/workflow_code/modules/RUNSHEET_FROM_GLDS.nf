process RUNSHEET_FROM_GLDS {
  // Downloads isa Archive and creates runsheet using GeneLab API
  tag "${ gldsAccession }"
  publishDir "${ params.outputDir }/${ gldsAccession }/Metadata",
    pattern: "*.{zip,csv}",
    mode: params.publish_dir_mode

  input:
    // TEMP: RESTORE ONCE OSD SUPPORT ADDED val(osdAccession)
    val(gldsAccession)
    path(dp_tools_plugin)

  output:
    path("${ gldsAccession }_*_v?_runsheet.csv"), emit: runsheet
    path("*.zip"), emit: isaArchive

  script:
    def injects = params.biomart_attribute ? "--inject biomart_attribute='${ params.biomart_attribute }'" : ''
    """

    dpt-get-isa-archive --accession ${ gldsAccession }
    ls ${dp_tools_plugin}

    dpt-isa-to-runsheet --accession ${ gldsAccession } \
      --plugin-dir ${dp_tools_plugin} \
      --isa-archive *.zip ${ injects }
    """
}