process VV_AGILE1CH {
  // Log publishing
  publishDir "${ params.outputDir }/${ params.gldsAccession }",
    pattern:  "VV_report.tsv.MANUAL_CHECKS_PENDING" ,
    mode: params.publish_dir_mode,
    saveAs: { "VV_Logs/VV_log_${ task.process.replace(":","-") }.tsv.MANUAL_CHECKS_PENDING" }
  // V&V'ed data publishing
  publishDir "${ params.outputDir }/${ params.gldsAccession }",
    pattern: '00-RawData',
    mode: params.publish_dir_mode
  publishDir "${ params.outputDir }/${ params.gldsAccession }",
    pattern: '01-limma_NormExp',
    mode: params.publish_dir_mode
  publishDir "${ params.outputDir }/${ params.gldsAccession }",
    pattern: '02-limma_DGE',
    mode: params.publish_dir_mode
  publishDir "${ params.outputDir }/${ params.gldsAccession }",
    pattern: 'Metadata/**',
    mode: params.publish_dir_mode

  label 'VV'

  input:
    path("VV_INPUT/Metadata/*") // While files from processing are staged, we instead want to use the files located in the publishDir for QC
    tuple path("VV_INPUT/02-limma_DGE"),
          path("VV_INPUT/01-limma_NormExp"),
          path("VV_INPUT/00-RawData")
     // "While files from processing are staged, we instead want to use the files located in the publishDir for QC
    val(skipVV) // Skips running V&V but will still publish the files
    path("dp_tools__agilent_1_channel")
  
  output:
    path("Metadata/*_runsheet.csv"), emit: VVed_runsheet
    tuple path("02-limma_DGE"),
          path("01-limma_NormExp"),
          path("00-RawData"), emit: VVed_de_data
    path("VV_report.tsv.MANUAL_CHECKS_PENDING"), optional: params.skipVV, emit: log
    path("versions.yml"), emit: versions

  script:
    """
    # move from VV_INPUT to task directory
    # This allows detection as output files for publishing
    mv VV_INPUT/* .

    # Run V&V unless user requests to skip V&V
    if ${ !skipVV} ; then
      dpt validation run dp_tools__agilent_1_channel . Metadata/*_runsheet.csv
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