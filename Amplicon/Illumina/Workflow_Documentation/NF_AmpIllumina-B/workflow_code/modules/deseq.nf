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
        path(dummy) // dummy path to ensure dependency between this step and the step that generates this file

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
                  --target-region  '${meta.target_region}' \\
                  --prevalence-cutoff ${meta.prevalence_cutoff} \\
                  --library-cutoff  ${meta.library_cutoff} ${meta.extra}

        
        Rscript -e "VERSIONS=sprintf('DESeq2 %s\\n', packageVersion('DESeq2')); \\
                    write(x=VERSIONS, file='versions.txt', append=TRUE)"     
        """
}


workflow {

    
    meta  = Channel.of(["samples": params.samples_column,
                        "group" : params.group,
                        "depth" : params.rarefaction_depth,
                        "assay_suffix" : params.assay_suffix,
                        "output_prefix" : params.output_prefix,
                        "target_region" : params.target_region,
                        "library_cutoff" : params.library_cutoff,
                        "prevalence_cutoff" : params.prevalence_cutoff,
                        "extra" : params.remove_rare ? "--remove-rare" : ""
                        ])
                            
                            
    metadata  = Channel.fromPath(params.metadata, checkIfExists: true)
    asv_table = Channel.fromPath(params.asv_table, checkIfExists: true) 
    taxonomy  =  Channel.fromPath(params.taxonomy, checkIfExists: true)
    // Dummy file
    dummy  =  Channel.fromPath(params.taxonomy, checkIfExists: true)

    DESEQ(meta, metadata, asv_table, taxonomy, dummy)

    emit:
        version = DESEQ.out.version

}
