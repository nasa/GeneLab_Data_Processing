process GENERATE_MD5SUMS {
  // Generates tabular data indicating genelab standard publishing files, md5sum generation, and tool version table formatting
  tag "${ params.gldsAccession }"
  publishDir "${ params.outputDir }/${ params.gldsAccession }/GeneLab",
    mode: params.publish_dir_mode

  input:
    path(data_dir)
    path(runsheet)
    path("dp_tools__affymetrix")

  output:
    path("*md5sum*")
    path("Missing_md5sum_files.txt"), optional: true

  script:
    """
    generate_md5sum_files.py  --root-path ${ data_dir } \\
                              --runsheet-path ${ runsheet } \\
                              --plug-in-dir "dp_tools__affymetrix"
    """
}