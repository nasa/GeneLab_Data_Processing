process GENERATE_SOFTWARE_TABLE {
  publishDir "${ params.resultsDir }/GeneLab",
    pattern: "software_versions.md",
    mode: params.publish_dir_mode

  input:
    path("software_versions.yaml")
  
  output:
    path("software_versions.md")
  
  script:
    """
    SoftwareYamlToMarkdownTable.py software_versions.yaml
    """
}