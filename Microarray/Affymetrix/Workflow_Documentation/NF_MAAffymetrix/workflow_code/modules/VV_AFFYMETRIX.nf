process VV_AFFYMETRIX {
  // Log publishing
  publishDir "${ params.resultsDir }",
    pattern:  "VV_report.tsv.MANUAL_CHECKS_PENDING" ,
    mode: params.publish_dir_mode,
    saveAs: { "VV_Logs/VV_log_${ task.process.replace(":","-") }.tsv.MANUAL_CHECKS_PENDING" }
  // V&V'ed data publishing
  publishDir "${ params.resultsDir }",
    pattern: '00-RawData/**',
    mode: params.publish_dir_mode
  publishDir "${ params.resultsDir }",
    pattern: '01-oligo_NormExp/**',
    mode: params.publish_dir_mode
  publishDir "${ params.resultsDir }",
    pattern: '02-limma_DGE/**',
    mode: params.publish_dir_mode
  publishDir "${ params.resultsDir }",
    pattern: 'Metadata/**',
    mode: params.publish_dir_mode

  label 'VV'

  input:
    path("VV_INPUT/Metadata/*") // While files from processing are staged, we instead want to use the files located in the publishDir for QC
    path("VV_INPUT/*") // "While files from processing are staged, we instead want to use the files located in the publishDir for QC
    val(skipVV) // Skips running V&V but will still publish the files
    path("dp_tools__affymetrix_channel")
  
  output:
    path("Metadata/*_runsheet.csv"), emit: VVed_runsheet
    path("00-RawData/*"), emit: VVed_rawData
    path("01-oligo_NormExp/*"), emit: VVed_NormExp
    path("02-limma_DGE/*"), emit: VVed_DGE
    path("VV_report.tsv.MANUAL_CHECKS_PENDING"), optional: params.skipVV, emit: log
    path("versions.yml"), emit: versions

  script:
    """
    # move from VV_INPUT to task directory
    # This allows detection as output files for publishing
    mv VV_INPUT/* .

    # Run V&V unless user requests to skip V&V
    if ${ !skipVV} ; then
      dpt validation run dp_tools__affymetrix_channel . Metadata/*_runsheet.csv --max-flag-code ${ params.max_flag_code }
    fi

    # Export versions
    cat >> versions.yml <<END_OF_VERSIONS
    - name: dp_tools
      version: \$(python -c "import dp_tools; print(dp_tools.__version__)")
      homepage: https://github.com/J-81/dp_tools
      workflow task: ${task.process}
    END_OF_VERSIONS
    """
}