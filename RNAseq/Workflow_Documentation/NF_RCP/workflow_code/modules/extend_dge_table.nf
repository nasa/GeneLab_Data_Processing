/*
 *  Extend Differential Expresssion Table To Generate Visualization Table
 */

process EXTEND_DGE_TABLE {
    tag "Dataset-wide"

    input:
        path(dge_table)

    output:
        path("visualization_output_table${params.output_suffix}.csv"), emit: dge_visualization_table
        path("versions.txt"),    emit: versions_txt

    script:
        def output_filename_suffix = params.output_suffix ?: ""
        def extend_rmd_file = "${projectDir}/bin/extend_dge_table.Rmd"

        """
        Rscript -e "rmarkdown::render('${extend_rmd_file}', 
        output_file = 'Extend_DGE_Table.html',
        output_dir = '\${PWD}',
            params = list(
                work_dir = '\${PWD}',
                output_directory = '\${PWD}',
                output_filename_suffix = '${output_filename_suffix}',
                input_table_path = '${dge_table}'
            ))"
        """
}