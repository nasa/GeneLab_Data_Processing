process GENERATE_SOFTWARE_TABLE {
  publishDir "${ params.resultsDir }/GeneLab",
    pattern: "software_versions_GLmicroarray.md",
    mode: params.publish_dir_mode

  input:
    path("software_versions.yaml")
  
  output:
    path("software_versions_GLmicroarray.md")
  
  script:
    """
    SoftwareYamlToMarkdownTable.py software_versions.yaml
    """
}