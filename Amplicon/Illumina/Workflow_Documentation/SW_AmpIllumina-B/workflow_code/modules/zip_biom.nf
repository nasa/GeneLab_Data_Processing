#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


/**************************************************************************************** 
*********************  Zip Biom ********************************************************
****************************************************************************************/

process ZIP_BIOM {

    tag "Zipping the taxonomy counts...."

    input:
       path(taxonomy_and_counts_biom)  // path("taxonomy-and-counts_${params.assay_suffix}.biom")
    output:
        path("taxonomy-and-counts_${params.assay_suffix}.biom.zip")
    script:
        """
        zip -j -q taxonomy-and-counts_${params.assay_suffix}.biom.zip \\
                ${taxonomy_and_counts_biom}
        """
}
