process DUMP_META {
  publishDir "${ publishdir }/GeneLab",
    mode: params.publish_dir_mode
  
  input:
    val(publishdir)
    val(meta)
  
  output:
    path("meta.sh")
  
  script:
  """
  # Write the meta file
  reformat_meta.sh '${ meta }' > meta.sh
  """
}