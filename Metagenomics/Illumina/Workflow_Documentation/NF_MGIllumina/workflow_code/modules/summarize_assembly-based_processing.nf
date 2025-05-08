#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/**************************************************************************************** 
*********************  Summarize Assembly based metagenomics processing **************************
****************************************************************************************/

process GENERATE_ASSEMBLY_PROCESSING_OVERVIEW_TABLE {

    tag "Summarizing the results of assemnly processing...."
    label "bit"

    input:
        path(sample_IDs_file)
        path(MAGs_overview)
        path(MAGs_dir)
        path(assemblies)
        path(genes_aa)
        path(metabat_assembly_depth_files)
        path(bins)
        path(bam_files)
    output:
        path("${params.additional_filename_prefix}Assembly-based-processing-overview${params.assay_suffix}.tsv")
    script:
        """
        mkdir assemblies_dir/ && mv *-assembly.fasta assemblies_dir/
        mkdir genes_dir/ && mv *-genes.faa genes_dir/ 
        mkdir mapping_dir/ && mv *-metabat-assembly-depth.tsv *.bam  mapping_dir/
        mkdir bins_dir/ && mv  *-bin*.fasta  bins_dir/
        bash generate-assembly-based-overview-table.sh \\
                ${sample_IDs_file} \\
                assemblies_dir/ \\
                genes_dir/ \\
                mapping_dir/ \\
                bins_dir/ \\
                ${MAGs_dir}/ \\
                ${params.additional_filename_prefix}Assembly-based-processing-overview${params.assay_suffix}.tsv
        """
}

