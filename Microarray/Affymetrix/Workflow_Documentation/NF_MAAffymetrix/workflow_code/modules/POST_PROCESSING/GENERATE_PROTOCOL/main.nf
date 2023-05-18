process GENERATE_PROTOCOL {
  tag "${ params.gldsAccession }"
  publishDir "${ params.outputDir }/${ params.gldsAccession }/GeneLab",
    mode: params.publish_dir_mode

  input:
    path("software_versions.md")
    path("meta.sh")
  
  output:
    path("PROTOCOL.txt")
  
  script:
  """
  generate_protocol.sh
  """
}