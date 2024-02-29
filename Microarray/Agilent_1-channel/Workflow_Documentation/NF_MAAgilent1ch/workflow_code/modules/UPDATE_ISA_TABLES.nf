process UPDATE_ISA_TABLES {
  // Generates tabular data indicating genelab standard publishing files, md5sum generation, and tool version table formatting
  tag "${ params.gldsAccession }"
  publishDir "${ params.outputDir }/${ params.gldsAccession }/GeneLab",
    mode: params.publish_dir_mode

  input:
    path(data_dir)
    path(runsheet)
    path(dp_tools__agilent_1_channel)

  output:
    path("updated_curation_tables") // directory containing extended ISA tables

  script:
    """
    update_curation_table.py  --root-path ${ data_dir } \\
                              --runsheet-path ${ runsheet } \\
                              --plug-in-dir ${ dp_tools__agilent_1_channel } \\
                              --isa-path ${ data_dir }/Metadata/*ISA*.zip

    # Update assay table with gldsAccession
    sed -i 's/${ params.osdAccession }/${ params.gldsAccession }/g' updated_curation_tables/a*.txt
    """
}