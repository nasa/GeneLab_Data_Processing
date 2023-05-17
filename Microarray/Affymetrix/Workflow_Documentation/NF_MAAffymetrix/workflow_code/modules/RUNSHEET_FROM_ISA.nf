process RUNSHEET_FROM_ISA {
  // Generates Runsheet using a path to an ISA archive
  tag "${ gldsAccession }"
  publishDir "${ params.outputDir }/${ gldsAccession }/Metadata",
    pattern: "*.zip",
    mode: params.publish_dir_mode

  input:
    val(osdAccession)
    val(gldsAccession)
    path(isaArchive)
    path("dp_tools__agilent_1_channel")

  output:
    path("${ osdAccession }_microarray_v?_runsheet.csv"), emit: runsheet

  script:
    """
    dpt-isa-to-runsheet --accession ${ osdAccession } \
      --plugin-dir dp_tools__agilent_1_channel \
      --isa-archive ${ isaArchive }
    """
}