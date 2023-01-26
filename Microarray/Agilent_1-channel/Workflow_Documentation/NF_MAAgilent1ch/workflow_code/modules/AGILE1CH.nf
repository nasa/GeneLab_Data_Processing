process AGILE1CH {
  stageInMode 'copy'

  input:
    path(qmd) // quarto qmd file to render
    path(runsheet_csv) // runsheet to supply as parameter
    val(biomart_attribute) // biomart attribute to supply as parameter (overridden by runsheet 'Array_Design_REF' column)
    path(annotation_file_path) // runsheet to supply as parameter
    val(organism) // runsheet to supply as parameter
    val(limit_biomart_query) // DEBUG option, limits biomart queries to the number specified if not set to false

  output:
    path("Agile1CMP.html"), emit: report

    tuple path("contrasts.csv"),
          path("SampleTable.csv"),
          path("differential_expression.csv"),
          path("visualization_output_table.csv"),
          path("visualization_PCA_table.csv"), emit: de

    path("versions.txt"), emit: version // TODO: convert version reporting to json format

  script:
    def limit_biomart_query_parameter = limit_biomart_query ? "-P DEBUG_limit_biomart_query:${limit_biomart_query}" : ''
    """
        quarto render \$PWD/${qmd} \
            -P 'runsheet:${runsheet_csv}' \
            -P 'biomart_attribute:${biomart_attribute}' \
            -P 'annotation_file_path:${annotation_file_path}' \
            -P 'organism:${organism}' \
            ${limit_biomart_query_parameter}

        echo "quarto version: \$(quarto --version)" >> versions.txt
    """
}