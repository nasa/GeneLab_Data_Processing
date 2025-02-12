process MULTIQC {
    // tag("Dataset-wide")
    
    input:
    path(sample_names)
    path("mqc_in/*") // any number of multiqc compatible files
    path(multiqc_config)
    

    output:
    // path("${ params.MQCLabel }_multiqc_GLbulkRNAseq_report.zip"), emit: zipped_report, to do: reimplement zip output w/ cleaned paths
    path("${ params.MQCLabel }_multiqc_GLbulkRNAseq_report"), emit: unzipped_report
    path("${ params.MQCLabel }_multiqc_GLbulkRNAseq_report.zip"), emit: zipped_report
    path("${ params.MQCLabel }_multiqc_GLbulkRNAseq_report/${ params.MQCLabel }_multiqc_GLbulkRNAseq.html"), emit: html
    path("${ params.MQCLabel }_multiqc_GLbulkRNAseq_report/${ params.MQCLabel }_multiqc_GLbulkRNAseq_data"), emit: data
    path("versions.yml"), emit: versions

    script:
    def config_arg = multiqc_config.name != "NO_FILE" ? "--config ${ multiqc_config }" : ""
    """
    multiqc \\
        --force \\
        --interactive \\
        -o ${ params.MQCLabel }_multiqc_GLbulkRNAseq_report \\
        -n ${ params.MQCLabel }_multiqc_GLbulkRNAseq \\
        ${ config_arg } \\
        .
    zip -r '${ params.MQCLabel }_multiqc_GLbulkRNAseq_report.zip' '${ params.MQCLabel }_multiqc_GLbulkRNAseq_report'

    echo "${task.process}:" > versions.yml
    echo "    multiqc: \$(multiqc --version | sed -e "s/multiqc, version //g")" >> versions.yml
    """
}
