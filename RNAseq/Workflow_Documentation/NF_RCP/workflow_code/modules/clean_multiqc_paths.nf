process CLEAN_MULTIQC_PATHS {
    tag "${mqc_label}"

    input:
        path(multiqc_report_zip)
        val(mqc_label)

    output:
        path("output/${mqc_label}_multiqc${params.assay_suffix}_report.zip"), emit: zipped_report

    script:
        """
        # Setup
        mkdir -p output
        unzip -q ${multiqc_report_zip} -d ./output
        
        cd ./output
        zip -r ${multiqc_report_zip} . -i "*_data" "*.html"
        """
}