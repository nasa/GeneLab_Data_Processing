process UPDATE_ASSAY_TABLE {
    publishDir "${ch_processed_directory}/GeneLab/updated_curation_tables",
    mode: params.publish_dir_mode,
    pattern: "a_*.txt"

    input:
        path(ch_processed_directory)
    output:
        path("a_*.txt"), emit: assay_table

    script:
    """
    update_assay_table.py --outdir ${ch_processed_directory} --mode ${params.mode} --assay_suffix ${params.assay_suffix} --glds_accession ${params.accession}
    """
}