process GENERATE_MD5SUMS {
  // Generates tabular data indicating genelab standard publishing files, md5sum generation, and tool version table formatting
  publishDir "${ data_dir }/GeneLab",
    mode: params.publish_dir_mode

  input:
    path(data_dir)
    val(done_token) // ensures process runs after purging processing_info.txt

  output:
    path("processed_md5sum_GLmicroarray.tsv"), emit: processed_md5sum

  script:
    """
    generate_md5sum_files.py  --outdir ${data_dir} --workflow_version ${workflow.manifest.version}
    """
}