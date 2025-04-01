process CLEAN_MULTIQC {

    input:
        path(multiqc_report_zip)
        val(mqc_label)

    output:
        path("output/${mqc_label}_multiqc${params.assay_suffix}_data.zip"), emit: zipped_data

    script:
        """
        mkdir -p output
        clean_multiqc_paths.py ${multiqc_report_zip} output
        """
}