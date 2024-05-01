#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


/**************************************************************************************** 
*********************  Plot generation *************************************************
****************************************************************************************/

// This process generates R visualizations using trimmed read data and grouping info from the runsheet
process R_VISUALIZATION {

    tag "Generating ASV plots...."
     
    input:
        path(runsheet)
        path(sample_IDs_file)
        path(counts)
        path(taxonomy)
    output:
        path("plots/*")
    script:
        """
        [ -d plots/ ] || mkdir plots/
        Illumina-R-visualizations.R \\
                      "${runsheet}" \\
                      "${sample_IDs_file}" \\
                      "${counts}" \\
                      "${taxonomy}" \\
                      "plots/" \\
                      "${params.output_prefix}" \\
                      "${params.assay_suffix}"
        """
}
