#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
//params.assay_suffix = "_GLAmpSeq"
//params.group  = "groups"
//params.samples_column = "Sample Name"
//params.rarefaction_depth = 500
//params.output_prefix = ""
//params.input_file = "../GeneLab/GLDS-487_amplicon_v1_runsheet.csv"
//params.asv_table = "../workflow_output/Final_Outputs/counts_GLAmpSeq.tsv"
//params.taxonomy = "../workflow_output/Final_Outputs/taxonomy_GLAmpSeq.tsv"


process ALPHA_DIVERSITY {

    tag "Running alpha diversity analysis..."
    label "visualization"

    input:
        val(meta)
        path(feature_table)
        path(taxonomy_table)
        path(metadata)

    output:
        path("alpha_diversity/"), emit: output_dir
        path("versions.txt"), emit: version

    script:
        """
        alpha_diversity.R \\
                  --metadata-table '${metadata}' \\
                  --feature-table '${feature_table}' \\
                  --taxonomy-table '${taxonomy_table}' \\
                  --group '${meta.group}' \\
                  --samples-column '${meta.samples}' \\
                  --rarefaction-depth ${meta.depth}   \\
                  --assay-suffix  '${meta.assay_suffix}' \\
                  --output-prefix  '${meta.output_prefix}'
                  
        Rscript -e "VERSIONS=sprintf('tidyverse %s\\nglue %s\\nvegan %s\\nhere %s\\nphyloseq %s\\nFSA %s\\nmultcompView %s\\nrstatix %s\\npatchwork %s\\nRColorBrewer %s\\n',  \\
                                    packageVersion('tidyverse'), \\
                                    packageVersion('glue'), \\
                                    packageVersion('vegan'), \\
                                    packageVersion('here'), \\
                                    packageVersion('phyloseq'), \\
                                    packageVersion('FSA'), \\
                                    packageVersion('multcompView'), \\
                                    packageVersion('rstatix'), \\
                                    packageVersion('patchwork'), \\
                                    packageVersion('RColorBrewer')); \\
                    write(x=VERSIONS, file='versions.txt', append=TRUE)"
        """

}





process BETA_DIVERSITY {

    tag "Running beta diversity analysis..."
    label "visualization"

    input:
        val(meta)
        path(feature_table)
        path(taxonomy_table)
        path(metadata)

    output:
        path("beta_diversity/"), emit: output_dir
        path("versions.txt"), emit: version

    script:
        """
        beta_diversity.R \\
                  --metadata-table '${metadata}' \\
                  --feature-table '${feature_table}' \\
                  --taxonomy-table '${taxonomy_table}' \\
                  --group '${meta.group}' \\
                  --samples-column '${meta.samples}' \\
                  --rarefaction-depth ${meta.depth}   \\
                  --assay-suffix  '${meta.assay_suffix}' \\
                  --output-prefix  '${meta.output_prefix}'
        
        Rscript -e "VERSIONS=sprintf('tidyverse %s\\nglue %s\\nvegan %s\\nhere %s\\nphyloseq %s\\nDESeq2 %s\\nggdendro %s\\nbroom %s\\nRColorBrewer %s\\n',  \\
                                    packageVersion('tidyverse'), \\
                                    packageVersion('glue'), \\
                                    packageVersion('vegan'), \\
                                    packageVersion('here'), \\
                                    packageVersion('phyloseq'), \\
                                    packageVersion('DESeq2'), \\
                                    packageVersion('ggdendro'), \\
                                    packageVersion('broom'), \\
                                    packageVersion('RColorBrewer')); \\
                    write(x=VERSIONS, file='versions.txt', append=TRUE)"         
        """

}


workflow{


     meta  = Channel.of(["samples": params.samples_column,
                            "group" : params.group,
                            "depth" : params.rarefaction_depth,
                            "assay_suffix" : params.assay_suffix,
                            "output_prefix" : params.output_prefix
                            ])
                            
     
     asv_table       =  Channel.fromPath(params.asv_table, checkIfExists: true)
     taxonomy_table  =  Channel.fromPath(params.taxonomy, checkIfExists: true)
     metadata        =  Channel.fromPath(params.input_file, checkIfExists: true) 
                            
    ALPHA_DIVERSITY(meta, asv_table, taxonomy_table, metadata)
    BETA_DIVERSITY(meta, asv_table, taxonomy_table, metadata)

    // Capture software versions
    software_versions_ch = Channel.empty()
    ALPHA_DIVERSITY.out.version | mix(software_versions_ch) | set{software_versions_ch}
    BETA_DIVERSITY.out.version | mix(software_versions_ch) | set{software_versions_ch}
    
    emit:
        version = software_versions_ch

}