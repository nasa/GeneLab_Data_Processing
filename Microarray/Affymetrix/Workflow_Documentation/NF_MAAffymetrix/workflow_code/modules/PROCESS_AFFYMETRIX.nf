process PROCESS_AFFYMETRIX {
  publishDir "${ params.resultsDir }/GeneLab",
    pattern: "NF_MAAffymetrix_v${workflow.manifest.version}_GLmicroarray.html",
    mode: params.publish_dir_mode
  stageInMode 'copy'

  input:
    path(qmd) // quarto qmd file to render
    path(runsheet_csv) // runsheet to supply as parameter
    path(annotation_file_path)
    tuple val(ensemblVersion), val(ensemblSource)
    val(limit_biomart_query) // DEBUG option, limits biomart queries to the number specified if not set to false
    val(skipDE) // whether to skip DE

  output:
    path("NF_MAAffymetrix_v${workflow.manifest.version}_GLmicroarray.html"), emit: report

    tuple path("02-limma_DGE"),
          path("01-oligo_NormExp"),
          path("00-RawData"), emit: de

    path("versions.yml"), emit: versions 

  script:
    def limit_biomart_query_parameter = limit_biomart_query ? "-P DEBUG_limit_biomart_query:${limit_biomart_query}" : ''
    def run_DE = skipDE ? "-P run_DE:'false'" : ''
    """
        export HOME=\$PWD;

        quarto render \$PWD/${qmd} \
            -P 'workflow_version:${workflow.manifest.version}' \
            -P 'runsheet:${runsheet_csv}' \
            -P 'annotation_file_path:${annotation_file_path}' \
            -P 'ensembl_version:${ensemblVersion}' \
            -P 'local_annotation_dir:${params.referenceStorePath}' \
            ${limit_biomart_query_parameter} \
            ${run_DE}

        # Rename report
        mv Affymetrix.html NF_MAAffymetrix_v${workflow.manifest.version}_GLmicroarray.html

        cat >> versions.yml <<END_OF_VERSIONS
        - name: quarto
          version: \$(quarto --version)
          homepage: https://quarto.org/
          workflow task: ${task.process}
        END_OF_VERSIONS
    """
}