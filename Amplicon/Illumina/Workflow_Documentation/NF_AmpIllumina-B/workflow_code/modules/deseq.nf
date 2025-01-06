#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


process DESEQ  {

    tag "Runnning diffrential abundance using DESEQ2..."
    label "visualization"

    input:
        val(meta)
        path(feature_table)
        path(taxonomy_table)
        path(metadata)
        path(version) // dummy path to ensure dependency between this step and the step that generates this file
    output:
        path("differential_abundance/deseq2/"), emit: output_dir
        path("versions.txt"), emit: version

    script:
        """
        run_deseq2.R \\
                  --metadata-table '${metadata}' \\
                  --feature-table '${feature_table}' \\
                  --taxonomy-table '${taxonomy_table}' \\
                  --group '${meta.group}' \\
                  --samples-column '${meta.samples}' \\
                  --assay-suffix  '${meta.assay_suffix}' \\
                  --output-prefix  '${meta.output_prefix}' \\
                  --target-region  '${meta.target_region}'
        
        Rscript -e "VERSIONS=sprintf('tidyverse %s\\nglue %s\\nDESeq2 %s\\nRColorBrewer %s\\n',  \\
                                    packageVersion('tidyverse'), \\
                                    packageVersion('glue'), \\
                                    packageVersion('DESeq2'), \\
                                    packageVersion('RColorBrewer')); \\
                    write(x=VERSIONS, file='versions.txt', append=TRUE)"     
        """
}


workflow {

    
    meta  = Channel.of(["samples": params.samples_column,
                        "group" : params.group,
                        "assay_suffix" : params.assay_suffix,
                        "output_prefix" : params.output_prefix
                        ])
                            
                            
    metadata  = Channel.fromPath(params.metadata, checkIfExists: true)
    asv_table = Channel.fromPath(params.asv_table, checkIfExists: true) 
    taxonomy  =  Channel.fromPath(params.taxonomy, checkIfExists: true)
    // Dummy file
    version  =  Channel.fromPath(params.taxonomy, checkIfExists: true)
    
    DESEQ(meta, metadata, asv_table, taxonomy, version)

    emit:
        version = DESEQ.out.version

}
