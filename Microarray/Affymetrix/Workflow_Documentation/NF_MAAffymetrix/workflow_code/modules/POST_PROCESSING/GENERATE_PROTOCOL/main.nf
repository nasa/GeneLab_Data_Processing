process GENERATE_PROTOCOL {
  tag "${ params.gldsAccession }"
  publishDir "${ params.resultsDir }/GeneLab",
    mode: params.publish_dir_mode

  input:
    path("software_versions_GLmicroarray.md")
    path("meta.sh")
  
  output:
    path("PROTOCOL_GLmicroarray.txt")
  
  script:
  """
  generate_protocol.sh
  """
}