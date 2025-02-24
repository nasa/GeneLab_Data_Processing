process PARSE_QC_METRICS {
    publishDir "${ publishdir }/GeneLab",
    mode: params.publish_dir_mode

    input:
        val(publishdir)
        val(osd_accession)
        val(meta)
        path(isa_zip)
        path(all_multiqc_output)
        path(rsem_counts)

    output:
        path("qc_metrics${params.assay_suffix}.csv"), emit: file

    script:
        """
        parse_multiqc.py --osd-num ${ osd_accession } ${ meta.paired_end ? '--paired' : '' }
        """
}
