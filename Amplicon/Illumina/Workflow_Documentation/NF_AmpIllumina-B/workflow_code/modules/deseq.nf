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

    output:
        path("differential_abundance/"), emit: output_dir
        path("versions.txt"), emit: version

    script:
        """
        run_deseq.R \\
                  --metadata-table '${metadata}' \\
                  --feature-table '${feature_table}' \\
                  --taxonomy-table '${taxonomy_table}' \\
                  --group '${meta.group}' \\
                  --samples-column '${meta.samples}' \\
                  --assay-suffix  '${meta.assay_suffix}' \\
                  --output-prefix  '${meta.output_prefix}'
        
        Rscript -e "VERSIONS=sprintf('tidyverse %s\\nglue %s\\nhere %s\\nDESeq2 %s\\nRColorBrewer %s\\n',  \\
                                    packageVersion('tidyverse'), \\
                                    packageVersion('glue'), \\
                                    packageVersion('here'), \\
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

    
    DESEQ(meta, metadata, asv_table, taxonomy)

    emit:
        version = DESEQ.out.version

}
