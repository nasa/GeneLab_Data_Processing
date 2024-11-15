#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
//params.diff_abund_method = "ancombc2"
//params.assay_suffix = "_GLAmpSeq"
//params.output_prefix = ""
//params.group  = "groups"
//params.samples_column = "Sample Name"
//params.input_file = "../GeneLab/GLDS-487_amplicon_v1_runsheet.csv"
//params.asv_table = "../workflow_output/Final_Outputs/counts_GLAmpSeq.tsv"
//params.taxonomy = "../workflow_output/Final_Outputs/taxonomy_GLAmpSeq.tsv"

process ANCOMBC {

    tag "Running ${params.diff_abund_method} for differential abundance testing..."
    label "visualization"

    input: 
        val(meta)
        path(feature_table)
        path(taxonomy)
        path(metadata)

    output:
        path("differential_abundance/"), emit: output_dir
        path("versions.txt"), emit: version

    script:
        def script_name = params.diff_abund_method == "ancombc1" ? "pairwise_ancombc1.R" : "pairwise_ancombc2.R"
        """
        ${script_name} \\
                  --metadata-table '${metadata}' \\
                  --feature-table '${feature_table}' \\
                  --taxonomy-table '${taxonomy}' \\
                  --group '${meta.group}' \\
                  --samples-column '${meta.samples}' \\
                  --assay-suffix  '${meta.assay_suffix}' \\
                  --output-prefix  '${meta.output_prefix}' \\
                  --cpus ${task.cpus}

	      Rscript -e "VERSION=sprintf('ANCOMBC %s',  packageVersion('ANCOMBC')); \\
                    write(x=VERSION, file='versions.txt', append=TRUE)"
                    
        Rscript -e "VERSIONS=sprintf('tidyverse %s\\nglue %s\\nANCOMBC %s\\nhere %s\\nphyloseq %s\\nmia %s\\ntaxize %s\\nDescTools %s\\npatchwork %s\\nggrepel %s\\n',  \\
                                    packageVersion('tidyverse'), \\
                                    packageVersion('glue'), \\
                                    packageVersion('ANCOMBC'), \\
                                    packageVersion('here'), \\
                                    packageVersion('phyloseq'), \\
                                    packageVersion('mia'), \\
                                    packageVersion('taxize'), \\
                                    packageVersion('DescTools'), \\
                                    packageVersion('patchwork'), \\
                                    packageVersion('ggrepel')); \\
                    write(x=VERSIONS, file='versions.txt', append=TRUE)"
        """

}



workflow {

    
    meta  = Channel.of(["samples": params.samples_column,
                        "group" : params.group,
                        "assay_suffix" : params.assay_suffix,
                        "output_prefix" : params.output_prefix
                        ])
                            
                            
    metadata  = Channel.fromPath(params.input_file, checkIfExists: true)
    asv_table = Channel.fromPath(params.asv_table, checkIfExists: true) 
    taxonomy  =  Channel.fromPath(params.taxonomy, checkIfExists: true)

    
    ANCOMBC(meta, asv_table, taxonomy, metadata)

    emit:
        version = ANCOMBC.out.version

}


