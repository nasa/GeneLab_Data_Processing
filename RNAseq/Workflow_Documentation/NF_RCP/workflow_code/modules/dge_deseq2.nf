/*
 * DESeq2 Differential Gene Expression Analysis 
 */
// ERCC counts are removed before normalization

process DGE_DESEQ2 {

    input:
        val(meta)
        path(runsheet_path)
        path(gene_counts)
        path("dge_deseq2.Rmd")
        val(output_label)

    output:
        tuple path("Normalized_Counts${output_label}${params.assay_suffix}.csv"),
              path(params.mode == "microbes" ? "FeatureCounts_Unnormalized_Counts${output_label}${params.assay_suffix}.csv" : 
                   "RSEM_Unnormalized_Counts${output_label}${params.assay_suffix}.csv"),                               emit: norm_counts
        path("contrasts${output_label}${params.assay_suffix}.csv"),                                                    emit: contrasts
        path("SampleTable${output_label}${params.assay_suffix}.csv"),                                                  emit: sample_table      
        path("differential_expression_no_annotations${output_label}${params.assay_suffix}.csv"),                       emit: dge_table
        path("VST_Counts${output_label}${params.assay_suffix}.csv"),                                        emit: vst_norm_counts
        path("summary.txt"),                                                                            emit: summary
        path("versions2.txt"),                                                                          emit: versions

    script:
        def output_filename_label = output_label ?: ""
        def output_filename_suffix = params.assay_suffix ?: ""
        def microbes = params.mode == 'microbes' ? 'TRUE' : 'FALSE'
        def debug_dummy_counts = params.use_dummy_gene_counts ? 'TRUE'  : 'FALSE'
        def input_counts_path = params.mode == 'microbes' ? gene_counts : "gene_counts"

        """
        if [[ "${params.mode}" != "microbes" ]]; then
            mkdir -p gene_counts
            mv ${gene_counts} gene_counts/
        fi
        Rscript -e "rmarkdown::render('dge_deseq2.Rmd', 
            output_file = 'DGE_DESeq2.html',
            output_dir = '\${PWD}',
            params = list(
                cpus = ${task.cpus},
                parallel_config = '${params.dge_parallel_config ?: ""}',
                work_dir = '\${PWD}',
                output_directory = '\${PWD}',
                output_filename_label = '${output_filename_label}',
                output_filename_suffix = '${output_filename_suffix}',
                runsheet_path = '${runsheet_path}',
                microbes = ${microbes},
                gene_id_type = '${meta.gene_id_type}',
                input_counts = '${input_counts_path}',
                DEBUG_MODE_LIMIT_GENES = FALSE,
                DEBUG_MODE_ADD_DUMMY_COUNTS = ${debug_dummy_counts}
            ))"

        Rscript -e "versions <- c(); 
                    versions['R'] <- gsub(' .*', '', gsub('R version ', '', R.version\\\$version.string));
                    versions['BioConductor'] <- as.character(BiocManager::version()); 
                    pkg_list <- c('BiocParallel', 'DESeq2', 'tidyverse', 'dplyr', 'knitr', 'stringr', 'yaml');
                    if (${microbes} != TRUE) {
                        pkg_list <- c(pkg_list, 'tximport');
                    }
                    for(pkg in pkg_list) {
                        versions[pkg] <- as.character(packageVersion(pkg))
                    };
                    cat('"RNASEQ_DGE_DESEQ2":\\n', 
                        paste0('    ', names(versions), ': ', versions, collapse='\\n'), 
                        '\\n', sep='', file='versions2.txt')"
        """
}