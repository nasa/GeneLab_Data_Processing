process VV_AGILE1CH {
  // Log publishing
  publishDir "${ params.outputDir }/${ params.gldsAccession }",
    pattern:  "VV_report.tsv.MANUAL_CHECKS_PENDING" ,
    mode: params.publish_dir_mode,
    saveAs: { "VV_Logs/VV_log_${ task.process.replace(":","-") }.tsv.MANUAL_CHECKS_PENDING" }
  // V&V'ed data publishing
  publishDir "${ params.outputDir }/${ params.gldsAccession }",
    pattern: '00-DE/**',
    mode: params.publish_dir_mode
  publishDir "${ params.outputDir }/${ params.gldsAccession }",
    pattern: 'Metadata/**',
    mode: params.publish_dir_mode

  label 'VV'

  input:
    path("VV_INPUT/Metadata/*") // While files from processing are staged, we instead want to use the files located in the publishDir for QC
    path("VV_INPUT/00-DE/*") // "While files from processing are staged, we instead want to use the files located in the publishDir for QC
    val(skipVV) // Skips running V&V but will still publish the files
    path("dp_tools__agilent_1_channel")
  
  output:
    path("Metadata/*_runsheet.csv"), emit: VVed_runsheet
    path("00-DE/*"), emit: VVed_raw_reads
    path("VV_report.tsv.MANUAL_CHECKS_PENDING"), optional: params.skipVV, emit: log

  script:
    """
    # move from VV_INPUT to task directory
    # This allows detection as output files for publishing
    mv VV_INPUT/* .

    # Run V&V unless user requests to skip V&V
    if ${ !skipVV} ; then
      dpt validation run dp_tools__agilent_1_channel . Metadata/*_runsheet.csv
    fi
    """
}