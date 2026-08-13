process PURGE_PROCESSING_INFO {

  publishDir "${ data_dir }/GeneLab",
    mode: params.publish_dir_mode

  input:
    path(data_dir)
    path(processing_info, stageAs: "unpurged_processing_info_GLmicroarray.txt")

  output:
    path("nextflow_processing_info_GLmicroarray.txt"), emit: purged_processing_info

  script:
    """
    array_annot_basename=\$(basename '${params.array_annot_path}')
    sed "s|[^ ]*/\${array_annot_basename}|\${array_annot_basename}|g" ${processing_info} > nextflow_processing_info_GLmicroarray.txt
    """
}