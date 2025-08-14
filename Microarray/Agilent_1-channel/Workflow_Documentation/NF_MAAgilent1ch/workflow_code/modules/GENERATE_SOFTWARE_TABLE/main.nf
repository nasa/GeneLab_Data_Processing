process GENERATE_SOFTWARE_TABLE {
  publishDir "${ params.outputDir }/${ params.gldsAccession }/GeneLab",
    pattern: "software_versions_GLmicroarray.md",
    mode: params.publish_dir_mode

  input:
    path("software_versions.yaml")
    val(filename)
  
  output:
    path("software_versions_GLmicroarray.md")
  
  script:
    """
    SoftwareYamlToMarkdownTable.py software_versions.yaml \"$filename\"
    """
}