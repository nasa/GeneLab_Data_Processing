here::i_am("dge_annotation_R_scripts/dge_annotation_workflow.R")
library("optparse")
library("here")
library("cli")

parser <- OptionParser()
parser <- add_option(parser, c("--skip_perform_dge"),
    action = "store_true", default = FALSE,
    help = "Skips running the DGE, this can be used when the output from the DGE already exist",
)
parser <- add_option(parser, c("--skip_gene_annotation"),
    action = "store_true", default = FALSE,
    help = "Skips running the gene annotation",
)
parser <- add_option(parser, c("--work_dir"),
    default = getwd(),
    help = "Working directory for all analysis",
)
parser <- add_option(parser, c("--input_gene_results_dir"),
    help = "Directory containing .genes.results files for each sample",
)
parser <- add_option(parser, c("--primary_keytype"),
    help = "Gene ID used for annotation mapping. Often 'ENSEMBL' or 'TAIR' for most GLDS datasets",
)
parser <- add_option(parser, c("--normalization"),
    help = "Normalization approach. Allowed Values: ['default','ERCC-groupB']",
)
parser <- add_option(parser, c("--normalized_counts_output_prefix"),
    help = "output directory and prefix for normalized counts and related files",
)
parser <- add_option(parser, c("--dge_output_prefix"),
    help = "output directory and prefix for dge tables and related files",
)
parser <- add_option(parser, c("--DEBUG_MODE_LIMIT_GENES"),
    default = FALSE, action = "store_true",
    help = "Performs analysis on last 150 genes",
)
parser <- add_option(parser, c("--DEBUG_MODE_ADD_DUMMY_COUNTS"),
    default = FALSE, action = "store_true",
    help = "Replaces all gene counts with random values from 0 to 5000",
)
parser <- add_option(parser, c("--runsheet_path"),
    help = "runsheet csv path, one of two allowed metadata inputs, exactly one metadata input must be supplied",
)
parser <- add_option(parser, c("--annotation_file_path"),
    help = "Annotation database file to use for adding gene annotations",
)
parser <- add_option(parser, c("--extended_table_output_prefix"),
    help = "Visualization table output prefix",
)
parser <- add_option(parser, c("--extended_table_output_suffix"),
    help = "Visualization table output suffix",
)

args <- parse_args(parser)

# TODO: add mandatory argument checking

cat_bullet(ansi_align(names(args)), " : ", args)

if (!args$skip_perform_dge) {
    cli_alert_warning("Running Perform_DGE.Rmd")
    rmarkdown::render(here("dge_annotation_R_scripts", "Perform_DGE.Rmd"),
        output_dir = args$work_dir,
        params = list(
            work_dir = args$work_dir,
            input_gene_results_dir = args$input_gene_results_dir,
            primary_keytype = args$primary_keytype,
            runsheet_path = args$runsheet_path,
            normalization = args$normalization,
            dge_output_prefix = args$dge_output_prefix,
            normalized_counts_output_prefix = args$normalized_counts_output_prefix,
            DEBUG_MODE_LIMIT_GENES = args$DEBUG_MODE_LIMIT_GENES,
            DEBUG_MODE_ADD_DUMMY_COUNTS = args$DEBUG_MODE_ADD_DUMMY_COUNTS
        )
    )
    cli_alert_success("Done running Perform_DGE.Rmd")
} else {
    cli_alert_warning("Skipping Perform_DGE.Rmd")
}

if (!args$skip_gene_annotation) {
    cli_alert_warning("Running Add_Gene_Annotations.Rmd")
    rmarkdown::render(here("dge_annotation_R_scripts", "Add_Gene_Annotations.Rmd"),
        output_dir = args$work_dir,
        params = list(
            input_table_path = paste0(args$dge_output_prefix, "differential_expression_no_annotations.csv"),
            work_dir = args$work_dir,
            annotation_file_path = args$annotation_file_path,
            primary_keytype = args$primary_keytype,
            annotated_output_prefix = args$dge_output_prefix
        )
    )
    cli_alert_success("Done running Add_Gene_Annotations.Rmd")
} else {
    cli_alert_warning("Skipping Add_Gene_Annotations.Rmd")
}

if (!args$skip_gene_annotation) {
    cli_alert_warning("Running Extend_DGE_Table.Rmd")
    rmarkdown::render(here("dge_annotation_R_scripts", "Extend_DGE_Table.Rmd"),
        output_dir = args$work_dir,
        params = list(
            input_table_path = paste0(args$dge_output_prefix, "differential_expression.csv"),
            work_dir = args$work_dir,
            extended_table_output_prefix = args$extended_table_output_prefix,
            extended_table_output_suffix = args$extended_table_output_suffix
        )
    )
    cli_alert_success("Done running Extend_DGE_Table.Rmd")
} else {
    cli_alert_warning("Skipping Extend_DGE_Table.Rmd")
}

if (!args$skip_gene_annotation) {
    cli_alert_warning("Running Generate_PCA_Table.Rmd")
    rmarkdown::render(here("dge_annotation_R_scripts", "Generate_PCA_Table.Rmd"),
        output_dir = args$work_dir,
        params = list(
            input_table_path = paste0(args$normalized_counts_output_prefix, "Normalized_Counts.csv"),
            work_dir = args$work_dir,
            pca_table_output_prefix = args$extended_table_output_prefix,
            pca_table_output_suffix = args$extended_table_output_suffix
        )
    )
    cli_alert_success("Done running Generate_PCA_Table.Rmd")
} else {
    cli_alert_warning("Skipping Generate_PCA_Table.Rmd")
}

cli_alert_success("All done!")