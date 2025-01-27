/*
 * DESeq2 Differential Gene Expression Analysis 
 */
// ERCC counts are removed before normalization

process DGE_DESEQ2 {

    input:
        val(meta)
        path(runsheet_path)
        path(gene_counts)

    output:
        tuple path("Normalized_Counts${params.output_suffix}.csv"),
              path(params.mode == "microbes" ? "FeatureCounts_Unnormalized_Counts${params.output_suffix}.csv" : 
                   "RSEM_Unnormalized_Counts${params.output_suffix}.csv"),                emit: norm_counts
        path("contrasts${params.output_suffix}.csv"),                                     emit: contrasts
        path("SampleTable${params.output_suffix}.csv"),                                   emit: sample_table      
        path("differential_expression_no_annotations${params.output_suffix}.csv"),        emit: dge_table
        path("VST_Normalized_Counts${params.output_suffix}.csv"),                         emit: vst_norm_counts
        path("summary.txt"),                                                              emit: summary
        path("versions.txt"),                                                             emit: versions_txt

    script:
        def output_filename_suffix = params.output_suffix ?: ""
        def microbes = params.mode == 'microbes' ? 'TRUE' : 'FALSE'
        def dge_rmd_file = "${projectDir}/bin/deseq2_dge.Rmd"
        def debug_dummy_counts = params.use_dummy_gene_counts ? 'TRUE' : 'FALSE'

        """
        Rscript -e "rmarkdown::render('${dge_rmd_file}', 
        output_file = 'DESeq2_DGE.html',
        output_dir = '\${PWD}',
            params = list(
                cpus = ${task.cpus},
                work_dir = '\${PWD}',
                output_directory = '\${PWD}',
                output_filename_suffix = '${output_filename_suffix}',
                runsheet_path = '${runsheet_path}',
                microbes = ${microbes},
                gene_id_type = '${meta.gene_id_type}',
                input_gene_results_dir = '\${PWD}',
                DEBUG_MODE_LIMIT_GENES = FALSE,
                DEBUG_MODE_ADD_DUMMY_COUNTS = ${debug_dummy_counts}
            ))"
        """
}