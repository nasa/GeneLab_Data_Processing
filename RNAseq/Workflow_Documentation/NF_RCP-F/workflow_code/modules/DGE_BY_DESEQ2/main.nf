/*
 * Different Gene Expression Analysis Processes
 */
process DGE_BY_DESEQ2 {

  input:
    path("runsheet.csv")
    path("Rsem_gene_counts/*")
    val(meta)
    path(annotation_file)
    path("dge_annotation_R_scripts.zip")

  output:
    tuple path("norm_counts_output/Normalized_Counts.csv"),
          path("norm_counts_output/RSEM_Unnormalized_Counts.csv"), emit: norm_counts

    tuple path("dge_output/contrasts.csv"),
          path("dge_output/SampleTable.csv"),
          path("dge_output/differential_expression.csv"),
          path("dge_output/visualization_output_table.csv"),
          path("dge_output/visualization_PCA_table.csv"), emit: dge

    path("norm_counts_output/ERCC_Normalized_Counts.csv"), optional: true, emit: norm_counts_ercc

    tuple path("dge_output_ercc/ERCCnorm_contrasts.csv"),
          path("dge_output_ercc/ERCCnorm_SampleTable.csv"),
          path("dge_output_ercc/ERCCnorm_differential_expression.csv"),
          path("dge_output_ercc/visualization_output_table_ERCCnorm.csv"),
          path("dge_output_ercc/visualization_PCA_table_ERCCnorm.csv"), optional: true, emit: dge_ercc

    path("dge_output/summary.txt"), emit: summary
    path("dge_output_ercc/ERCCnorm_summary.txt"), optional: true, emit: summary_ercc

    path("versions.txt"), emit: version

  script:
    """
    # Unzip r scripts
    unzip dge_annotation_R_scripts.zip

    Rscript --vanilla dge_annotation_R_scripts/dge_annotation_workflow.R \\
        --runsheet_path runsheet.csv \\
        ${ params.use_dummy_gene_counts ? '--DEBUG_MODE_ADD_DUMMY_COUNTS' : ''} \\
        --input_gene_results_dir "Rsem_gene_counts" \\
        --primary_keytype ${ meta.primary_keytype } \\
        --normalization 'default' \\
        --normalized_counts_output_prefix "norm_counts_output/" \\
        --dge_output_prefix "dge_output/" \\
        --annotation_file_path ${annotation_file} \\
        --extended_table_output_prefix "dge_output/"\\
        --extended_table_output_suffix ".csv"

    if ${ meta.has_ercc ? 'true' : 'false'}
    then
        Rscript --vanilla dge_annotation_R_scripts/dge_annotation_workflow.R \\
            --runsheet_path runsheet.csv \\
            ${ params.use_dummy_gene_counts ? '--DEBUG_MODE_ADD_DUMMY_COUNTS' : ''} \\
            --input_gene_results_dir "Rsem_gene_counts" \\
            --primary_keytype ${ meta.primary_keytype } \\
            --normalization 'ERCC-groupB' \\
            --normalized_counts_output_prefix "norm_counts_output/ERCC_" \\
            --dge_output_prefix "dge_output_ercc/ERCCnorm_" \\
            --annotation_file_path ${annotation_file} \\
            --extended_table_output_prefix "dge_output_ercc/"\\
            --extended_table_output_suffix "_ERCCnorm.csv"
    fi
    # bump
    """
}