process AGILE1CH {
  publishDir "${ params.outputDir }/${ params.gldsAccession }/GeneLab",
    pattern: "NF_MAAgilent1ch_v${workflow.manifest.version}_GLmicroarray.html",
    mode: params.publish_dir_mode
  stageInMode 'copy'

  input:
    path(qmd) // quarto qmd file to render
    path(runsheet_csv) // runsheet to supply as parameter
    path(annotation_file_path) // runsheet to supply as parameter
    tuple val(ensemblVersion), val(ensemblSource)
    val(limit_biomart_query) // DEBUG option, limits biomart queries to the number specified if not set to false

  output:
    path("NF_MAAgilent1ch_v${workflow.manifest.version}_GLmicroarray.html"), emit: report

    tuple path("02-limma_DGE/contrasts_GLmicroarray.csv"),
          path("02-limma_DGE/SampleTable_GLmicroarray.csv"),
          path("02-limma_DGE/differential_expression_GLmicroarray.csv"),
          path("01-limma_NormExp/normalized_expression_GLmicroarray.csv"),
          path("00-RawData/raw_intensities_GLmicroarray.csv"), emit: de_all_files

    tuple path("02-limma_DGE"),
          path("01-limma_NormExp"),
          path("00-RawData"), emit: de

    path("versions.yml"), emit: versions // Note: Quarto version captured in script body.  R versions captured during render (part of qmd code).

  script:
    def limit_biomart_query_parameter = limit_biomart_query ? "-P DEBUG_limit_biomart_query:${limit_biomart_query}" : ''
    """
        export HOME=\$PWD;
        
        quarto render \$PWD/${qmd} \
            -P 'workflow_version:${workflow.manifest.version}' \
            -P 'runsheet:${runsheet_csv}' \
            -P 'annotation_file_path:${annotation_file_path}' \
            -P 'ensembl_version:${ensemblVersion}' \
            ${limit_biomart_query_parameter}

        # Rename report
        mv Agile1CMP.html NF_MAAgilent1ch_v${workflow.manifest.version}_GLmicroarray.html

        cat >> versions.yml <<END_OF_VERSIONS
        - name: quarto
          version: \$(quarto --version)
          homepage: https://quarto.org/
          workflow task: ${task.process}
        END_OF_VERSIONS
    """
}