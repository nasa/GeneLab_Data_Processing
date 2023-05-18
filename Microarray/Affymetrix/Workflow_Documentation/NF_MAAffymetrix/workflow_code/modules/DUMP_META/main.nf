process DUMP_META {
  publishDir "${ params.outputDir }/${ params.gldsAccession }/GeneLab",
    mode: params.publish_dir_mode
  
  input:
    val(meta)
  
  output:
    path("meta.sh")
  
  script:
  """
  # Write the meta file
  reformat_meta.sh '${ meta }' > meta.sh
  """
}