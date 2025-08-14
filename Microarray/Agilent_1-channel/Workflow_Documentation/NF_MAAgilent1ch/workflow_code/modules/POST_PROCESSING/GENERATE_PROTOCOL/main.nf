process GENERATE_PROTOCOL {
  tag "${ params.gldsAccession }"
  publishDir "${ params.outputDir }/${ params.gldsAccession }/GeneLab",
    mode: params.publish_dir_mode

  input:
    path("software_versions_GLmicroarray.md")
    val(organism)
  
  output:
    path("PROTOCOL_GLmicroarray.txt")
  
  script:
  """
  generate_protocol.sh $workflow.manifest.version \"$organism\"
  """
}