process MULTIQC {
    // tag("Dataset-wide")
    
    input:
    path(sample_names)
    path("mqc_in/*") // any number of multiqc compatible files
    path(multiqc_config)
    val(mqc_label)
    

    output:
    // path("${ mqc_label }multiqc_GLbulkRNAseq_report.zip"), emit: zipped_report, to do: reimplement zip output w/ cleaned paths
    path("${ mqc_label }multiqc${ params.assay_suffix }_report"), emit: unzipped_report
    path("${ mqc_label }multiqc${ params.assay_suffix }_report.zip"), emit: zipped_report
    path("${ mqc_label }multiqc${ params.assay_suffix }_report/${ mqc_label }multiqc${ params.assay_suffix }.html"), emit: html
    path("${ mqc_label }multiqc${ params.assay_suffix }_report/${ mqc_label }multiqc${ params.assay_suffix }_data"), emit: data
    path("versions.yml"), emit: versions

    script:
    def config_arg = multiqc_config.name != "NO_FILE" ? "--config ${ multiqc_config }" : ""
    """
    multiqc \\
        --force \\
        --interactive \\
        -o ${ mqc_label }multiqc${ params.assay_suffix  }_report \\
        -n ${ mqc_label }multiqc${ params.assay_suffix  } \\
        ${ config_arg } \\
        .
    zip -r '${ mqc_label }multiqc${ params.assay_suffix  }_report.zip' '${ mqc_label }multiqc${ params.assay_suffix  }_report'

    echo '"${task.process}":' > versions.yml
    echo "    multiqc: \$(multiqc --version | sed "s/multiqc, version //")" >> versions.yml
    """
}
