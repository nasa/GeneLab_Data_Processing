process VV_AGILE1CH {
  // Log publishing
  publishDir "${ params.outputDir }/${ params.gldsAccession }",
    pattern:  "VV_log.tsv" ,
    mode: params.publish_dir_mode,
    saveAs: { "VV_Logs/VV_log_${ task.process.replace(":","-") }.tsv" }
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
  
  output:
    path("Metadata/*_runsheet.csv"), emit: VVed_runsheet
    path("00-DE"), emit: VVed_raw_reads
    path("VV_log.tsv"), optional: params.skipVV, emit: log

  script:
    """
    # move from VV_INPUT to task directory
    # This allows detection as output files for publishing
    mv VV_INPUT/* .

    # Run V&V unless user requests to skip V&V
    if ${ !params.skipVV} ; then
      VV_data_assets.py   --root-path . \\
                          --accession ${ params.gldsAccession } \\
                          --runsheet-path Metadata/*_runsheet.csv \\
                          --data-asset-sets  \\
                            'Metadata' \\
                          --run-components \\
                            'Metadata' \\
                            'Raw Reads' \\
                            'Raw Reads By Sample' \\
                          --max-flag-code ${ params.max_flag_code }
    fi
    """
}