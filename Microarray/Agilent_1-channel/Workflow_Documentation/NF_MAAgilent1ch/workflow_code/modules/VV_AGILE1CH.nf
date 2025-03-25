process VV_AGILE1CH {
  // Log publishing
  publishDir "${ params.outputDir }/${ params.gldsAccession }",
    pattern:  "VV_report_GLmicroarray.tsv.MANUAL_CHECKS_PENDING" ,
    mode: params.publish_dir_mode,
    saveAs: { "VV_Logs/VV_log_${ task.process.replace(":","-") }_GLmicroarray.tsv.MANUAL_CHECKS_PENDING" }
  // V&V'ed data publishing
  publishDir "${ params.outputDir }/${ params.gldsAccession }",
    pattern: '00-RawData/**',
    mode: params.publish_dir_mode
  publishDir "${ params.outputDir }/${ params.gldsAccession }",
    pattern: '01-limma_NormExp/**',
    mode: params.publish_dir_mode
  publishDir "${ params.outputDir }/${ params.gldsAccession }",
    pattern: '02-limma_DGE/**',
    mode: params.publish_dir_mode
  publishDir "${ params.outputDir }/${ params.gldsAccession }",
    pattern: 'Metadata/**',
    mode: params.publish_dir_mode

  label 'VV'

  input:
    path("VV_INPUT/Metadata/*") // While files from processing are staged, we instead want to use the files located in the publishDir for QC
    path("VV_INPUT/*") // "While files from processing are staged, we instead want to use the files located in the publishDir for QC
    val(skipVV) // Skips running V&V but will still publish the files
    path(dp_tools_path)
  
  output:
    path("Metadata/*_runsheet.csv"), emit: VVed_runsheet
    path("02-limma_DGE/*"), emit: VVed_DGE, optional: true
    path("01-limma_NormExp/*"), emit: VVed_NormExp
    path("00-RawData/*"), emit: VVed_rawData
    path("VV_report_GLmicroarray.tsv.MANUAL_CHECKS_PENDING"), optional: params.skipVV, emit: log
    path("versions.yml"), emit: versions

  script:
    """
    # move from VV_INPUT to task directory
    # This allows detection as output files for publishing
    mv VV_INPUT/* .

    # Run V&V unless user requests to skip V&V
    if ${ !skipVV} ; then
      dpt validation run ${ dp_tools_path } . Metadata/*_runsheet.csv
      mv VV_report.tsv.MANUAL_CHECKS_PENDING VV_report_GLmicroarray.tsv.MANUAL_CHECKS_PENDING
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