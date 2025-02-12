/*
 *  Generate PCA Table from Normalized Counts
 */

process GENERATE_PCA_TABLE {
    // tag "Dataset-wide"

    input:
        path(normalized_counts)

    output:
        path("visualization_PCA_table${params.output_suffix}.csv"), emit: pca_table
        path("versions.txt"),    emit: versions_txt

    script:
        def output_filename_suffix = params.output_suffix ?: ""
        def pca_rmd_file = "${projectDir}/bin/generate_pca_table.Rmd"

        """
        Rscript -e "rmarkdown::render('${pca_rmd_file}', 
        output_file = 'Generate_PCA_Table.html',
        output_dir = '\${PWD}',
            params = list(
                work_dir = '\${PWD}',
                output_directory = '\${PWD}',
                output_filename_suffix = '${output_filename_suffix}',
                input_table_path = '${normalized_counts}'
            ))"
        """
}