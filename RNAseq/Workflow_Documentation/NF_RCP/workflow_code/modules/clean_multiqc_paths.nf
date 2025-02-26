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
        
        cd output/${multiqc_report_zip.simpleName}/*_data/
        
        # Clean paths in all files
        for f in *; do
            sed -E 's|.*/work/[^/]*/[^/]*/mqc_in/([^/]*)|\\1|g' \$f > tmp && mv tmp \$f
            sed -E 's|.*/GLDS_Datasets/(.+)|\\1|g' \$f > tmp && mv tmp \$f
            sed -E 's|.+/miniconda.+/envs/[^/]*/||g' \$f > tmp && mv tmp \$f
            sed -E 's|/[^ ]*/GLDS-|GLDS-|g' \$f > tmp && mv tmp \$f
        done
        
        cd ../..
        zip -r ${multiqc_report_zip} ./*
        """
}