#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process SOFTWARE_VERSIONS {

    tag "Writing out software versions..."
    label "fastqc"   //unix environment

    input:
        path(software_versions)

    output:
        path("software_versions.txt")

    script:
        """
        # Delete white spaces and write out unique software versions 
        grep -v "^\$" ${software_versions} | sort -u > software_versions.txt 
        """
}
