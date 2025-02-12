/*
 * Add GeneLab_Reference_Annotations annotations to DGE table
 */

process ADD_GENE_ANNOTATIONS {
    // tag "Dataset-wide"

    input:
        val(meta)
        val(gene_annotations_url)
        path(dge_no_annotations)

    output:
        path("differential_expression${params.output_suffix}.csv"),  emit: annotated_dge_table
        path("versions.txt"), emit: versions_txt

    script:
        def output_filename_suffix = params.output_suffix ?: ""
        def annotations_rmd_file = "${projectDir}/bin/add_gene_annotations.Rmd"

        """
        Rscript -e "rmarkdown::render('${annotations_rmd_file}', 
        output_file = 'DGE_Annotations.html',
        output_dir = '\${PWD}',
            params = list(
                work_dir = '\${PWD}',
                output_directory = '\${PWD}',
                output_filename_suffix = '${output_filename_suffix}',
                annotation_file_path = '${gene_annotations_url}',
                gene_id_type = '${meta.gene_id_type}',
                input_table_path = '${dge_no_annotations}'
            ))"
        """
}