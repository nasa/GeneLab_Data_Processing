process PROCESS_AFFYMETRIX {
  publishDir "${ publishdir }/GeneLab",
    pattern: "NF_MAAffymetrix_v${workflow.manifest.version}_GLmicroarray.html",
    mode: params.publish_dir_mode
  stageInMode 'copy'

  input:
    val(publishdir)
    path(qmd) // quarto qmd file to render
    path(runsheet_csv) // runsheet to supply as parameter
    path(array_data_files) // staged, locally-named, decompressed raw array data files
    path(annotation_file_path)
    tuple val(ensemblVersion), val(ensemblSource)
    path(referenceStorePath) // path to custom annotation references
    path(array_annot_path) // path to custom array design info file
    val(skipDE) // whether to skip DE

  output:
    path("NF_MAAffymetrix_v${workflow.manifest.version}_GLmicroarray.html"), emit: report

    tuple path("02-limma_DGE"),
          path("01-oligo_NormExp"),
          path("00-RawData"), emit: de

    path("versions.yml"), emit: versions 

  script:
    def run_DE = skipDE ? "-P run_DE:'false'" : ''
    """
        export HOME=\$PWD;

        quarto render \$PWD/${qmd} \
            -P 'workflow_version:${workflow.manifest.version}' \
            -P 'runsheet:${runsheet_csv}' \
            -P 'annotation_file_path:${annotation_file_path}' \
            -P 'ensembl_version:${ensemblVersion}' \
            -P 'local_annotation_dir:${referenceStorePath}' \
            -P 'array_annot_path:${array_annot_path}' \
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