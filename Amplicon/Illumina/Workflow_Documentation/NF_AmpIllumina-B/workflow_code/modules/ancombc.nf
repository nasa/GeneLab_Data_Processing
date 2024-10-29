#!/usr/bin/env nextflow

nextflow.enable.dsl = 2
//params.ancombc_version = 2
//params.assay_suffix = "_GLAmpSeq"
//params.group  = "groups"
//params.samples_column = "Sample Name"
//params.metadata = "../GeneLab/GLDS-487_amplicon_v1_runsheet.csv"
//params.asv_table = "../workflow_output/Final_Outputs/counts_GLAmpSeq.tsv"
//params.taxonomy = "../workflow_output/Final_Outputs/taxonomy_GLAmpSeq.tsv"

process ANCOMBC {

    tag "Running ANCOMBC ${params.ancombc_version} for differential abundance testing..."

    input: 
        tuple val(group), val(samples_column)
        path(metadata)
        path(asv_table)
        path(taxonomy)

    output:
        path("Plots/"), emit: plots
        path("differential_abundance${params.assay_suffix}.csv"), emit: table
        path("versions.txt"), emit: version

    script:
        def script_name = params.ancombc_version == 1 ? "pairwise_ancombc1.R" : "pairwise_ancombc2.R"
        """
        ${script_name} \\
                  --out-file  "differential_abundance${params.assay_suffix}.csv"   \\
                  --metadata-table '${metadata}' \\
                  --feature-table '${asv_table}' \\
                  --taxonomy-table '${taxonomy}' \\
                  --group '${group}' \\
                  --samples-column '${samples_column}' \\
                  --cpus ${task.cpus}

	Rscript -e "VERSION=sprintf('ANCOMBC %s',  packageVersion('ANCOMBC')); \\
                    write(x=VERSION, file='versions.txt', append=TRUE)"
        """

}



workflow {


    data_ch   = Channel.of([params.group, params.samples_column])
    metadata  = Channel.fromPath(params.metadata, checkIfExists: true)
    asv_table = Channel.fromPath(params.asv_table, checkIfExists: true) 
    taxonomy   =  Channel.fromPath(params.taxonomy, checkIfExists: true)

    
    ANCOMBC(data_ch, metadata, asv_table, taxonomy)

    emit:
        version = ANCOMBC.out.version

}


