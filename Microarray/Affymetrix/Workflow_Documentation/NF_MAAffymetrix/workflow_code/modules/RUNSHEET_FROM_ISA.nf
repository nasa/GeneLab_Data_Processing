process RUNSHEET_FROM_ISA {
  // Generates Runsheet using a path to an ISA archive
  tag "${ gldsAccession }"
  publishDir "${ params.resultsDir }/Metadata",
    pattern: "*.csv",
    mode: params.publish_dir_mode

  publishDir "${ params.resultsDir }/Metadata",
    pattern: "*.zip",
    mode: params.publish_dir_mode

  input:
    val(osdAccession)
    val(gldsAccession)
    path(isaArchive, stageAs: "input/*")
    path("dp_tools__affymetrix")

  output:
    path("${ osdAccession }_microarray_v?_runsheet.csv"), emit: runsheet
    path("*.zip"), emit: isaArchive, optional: true

  script:
    """
    dpt-isa-to-runsheet --accession ${ osdAccession } \
      --plugin-dir dp_tools__affymetrix \
      --isa-archive ${ isaArchive }

    ${params.isaArchivePath ? "cp ${ isaArchive } ." : ""} // Publish ISA.zip provided by --isaArchivePath 
    """
}