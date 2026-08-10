process UPDATE_ASSAY_TABLE {
  // Generates tabular data indicating genelab standard publishing files, md5sum generation, and tool version table formatting
  publishDir "${ data_dir }/GeneLab/updated_curation_tables",
    mode: params.publish_dir_mode

  input:
    path(data_dir)
    path(runsheet)
    path(isa_archive)
    val(glds_accession)

  output:
    path("a_*.txt"), emit: updated_assay_table

  script:
    """
    update_assay_table.py  --runsheet ${ runsheet } \
                              --glds_accession ${ glds_accession } \
                              --isa_zip ${ isa_archive }
    """
}