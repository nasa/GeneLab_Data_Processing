#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

//params.assay_suffix = "_GLAmpSeq"
//params.group  = "groups"
//params.samples_column = "Sample Name"
//params.output_prefix = ""
//params.input_file = "../GeneLab/GLDS-487_amplicon_v1_runsheet.csv"
//params.asv_table = "../workflow_output/Final_Outputs/counts_GLAmpSeq.tsv"
//params.taxonomy = "../workflow_output/Final_Outputs/taxonomy_GLAmpSeq.tsv"


process PLOT_TAXONOMY  {

    tag "Making taxonomy plots..."
    label "visualization"

    input:
        val(meta)
        path(feature_table)
        path(taxonomy_table)
        path(metadata)

    output:
        path("taxonomy_plots/"), emit: output_dir
        path("versions.txt"), emit: version

    script:
        """
        plot_taxonomy.R \\
                  --metadata-table '${metadata}' \\
                  --feature-table '${feature_table}' \\
                  --taxonomy-table '${taxonomy_table}' \\
                  --group '${meta.group}' \\
                  --samples-column '${meta.samples}' \\
                  --assay-suffix  '${meta.assay_suffix}' \\
                  --output-prefix  '${meta.output_prefix}'
                  
        Rscript -e "VERSIONS=sprintf('tidyverse %s\\nglue %s\\ntools %s\\nhere %s\\n',  \\
                                    packageVersion('tidyverse'), \\
                                    packageVersion('glue'), \\
                                    packageVersion('tools'), \\
                                    packageVersion('here')); \\
                    write(x=VERSIONS, file='versions.txt', append=TRUE)"
        """

}


workflow {


     meta  = Channel.of(["samples": params.samples_column,
                            "group" : params.group,
                            "assay_suffix" : params.assay_suffix,
                            "output_prefix" : params.output_prefix
                            ])
                            
     
     asv_table       =  Channel.fromPath(params.asv_table, checkIfExists: true)
     taxonomy_table  =  Channel.fromPath(params.taxonomy, checkIfExists: true)
     metadata        =  Channel.fromPath(params.input_file, checkIfExists: true) 
     
     
     PLOT_TAXONOMY(meta, asv_table, taxonomy_table, metadata)
     
     
     emit:
         version = PLOT_TAXONOMY.out.version
}