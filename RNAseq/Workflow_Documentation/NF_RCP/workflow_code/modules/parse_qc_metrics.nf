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
        path(runsheet_path)

    output:
        path("qc_metrics${params.assay_suffix}.csv"), emit: file

    script:
        def assay_suffix = params.assay_suffix ? "--assay_suffix ${params.assay_suffix}" : ""
        def mode_param = params.mode ? "--mode ${params.mode}" : ""
        """
        parse_multiqc.py --osd-num ${ osd_accession } ${ assay_suffix } ${ meta.paired_end ? '--paired' : '' } ${ mode_param } --runsheet ${ runsheet_path }
        """
}
