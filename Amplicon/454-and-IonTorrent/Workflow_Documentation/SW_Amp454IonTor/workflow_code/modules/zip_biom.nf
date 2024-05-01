#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


/**************************************************************************************** 
*********************  Zip Biom ********************************************************
****************************************************************************************/

process ZIP_BIOM {

    tag "Zipping the taxonomy counts...."

    input:
       path(taxonomy_and_counts_biom)  // path("taxonomy-and-counts${params.assay_suffix}.biom")
    output:
        path("${params.output_prefix}taxonomy-and-counts${params.assay_suffix}.biom.zip")
    script:
        """
        zip -q ${params.output_prefix}taxonomy-and-counts${params.assay_suffix}.biom.zip \\
                ${taxonomy_and_counts_biom}
        """
}
