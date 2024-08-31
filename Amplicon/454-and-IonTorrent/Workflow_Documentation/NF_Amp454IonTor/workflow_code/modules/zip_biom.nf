#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


/**************************************************************************************** 
*********************  Zip Biom ********************************************************
****************************************************************************************/

process ZIP_BIOM {

    tag "Zipping the taxonomy counts...."

    input:
       path(taxonomy_and_counts_biom)
    output:
        path("${params.output_prefix}taxonomy-and-counts${params.assay_suffix}.biom.zip"), emit: biom
        path("versions.txt"), emit: version
    script:
        """
        zip -q ${params.output_prefix}taxonomy-and-counts${params.assay_suffix}.biom.zip \\
                ${taxonomy_and_counts_biom}

        zip -h | grep "Zip" | sed -E 's/(Zip.+\\)).+/\\1/' > versions.txt
        """
}
