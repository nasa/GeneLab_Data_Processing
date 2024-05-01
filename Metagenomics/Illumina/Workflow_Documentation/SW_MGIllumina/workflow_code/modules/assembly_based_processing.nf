#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/**************************************************************************************** 
***************************  Assembly-based processing workflow *************************
****************************************************************************************/

// Assembly-based workflow
include { ASSEMBLE; RENAME_HEADERS; SUMMARIZE_ASSEMBLIES } from "./assembly.nf"
include { MAPPING; SAM_TO_BAM } from "./read_mapping.nf"
include { CALL_GENES; REMOVE_LINEWRAPS } from "./assembly_annotation.nf"
include { KO_ANNOTATION; FILTER_KFAMSCAN } from "./assembly_annotation.nf"
include { TAX_CLASSIFICATION } from "./assembly_annotation.nf"
include { GET_COV_AND_DET } from "./coverage.nf"
include { COMBINE_GENE_ANNOTS_TAX_AND_COVERAGE; MAKE_COMBINED_GENE_LEVEL_TABLES } from "./combine_contig_annotation.nf"
include { COMBINE_CONTIG_TAX_AND_COVERAGE; MAKE_COMBINED_CONTIG_TAX_TABLES } from "./combine_contig_annotation.nf"
include { METABAT_BINNING } from "./binning.nf"
include { summarize_bins } from "./summarize_bins.nf"
include { summarize_mags } from "./summarize_MAG.nf"
include { GENERATE_ASSEMBLY_PROCESSING_OVERVIEW_TABLE } from "./summarize_assembly-based_processing.nf"

workflow assembly_based {

    take:
        file_ch
        filtered_ch
        ko_db_dir
        cat_db
        gtdbtk_db_dir
        use_gtdbtk_scratch_location

    main:
        /*****************************************************
        ************* Assembly-based analysis ****************
        *****************************************************/
        // Assemble reads to contigs
        assembly_ch = ASSEMBLE(filtered_ch) | RENAME_HEADERS
        assemblies_ch = assembly_ch.map{
                                      sample_id, assembly -> file("${assembly}")
                                      }.collect()
        SUMMARIZE_ASSEMBLIES(assemblies_ch)
        
        // Map reads to assembly
        read_mapping_ch = MAPPING(assembly_ch.join(filtered_ch)) | SAM_TO_BAM

        // Annotate assembly
        genes_ch = CALL_GENES(assembly_ch) | REMOVE_LINEWRAPS
        if (ko_db_dir){
            annotations_ch = KO_ANNOTATION(assembly_ch.join(genes_ch), ko_db_dir) | FILTER_KFAMSCAN
        }else{
            SETUP_KOFAMSCAN_DB()
            annotations_ch = KO_ANNOTATION(assembly_ch.join(genes_ch), 
                                           SETUP_KOFAMSCAN_DB.out.ko_db_dir) | FILTER_KFAMSCAN
        }

        if (cat_db){         
            taxonomy_ch = TAX_CLASSIFICATION(assembly_ch.join(genes_ch), cat_db)
        }else{
            SETUP_CAT_DB(params.database.CAT_DB_LINK)
            taxonomy_ch = TAX_CLASSIFICATION(assembly_ch.join(genes_ch), SETUP_CAT_DB.out.cat_db)
        }

        // Calculate gene coverage and depth
        coverage_ch = GET_COV_AND_DET(read_mapping_ch
                                          .join(assembly_ch)
                                          .join(genes_ch))

        // Combine contig annotation
        tax_and_cov_ch = COMBINE_GENE_ANNOTS_TAX_AND_COVERAGE(coverage_ch
                                                                .join(annotations_ch)
                                                                .join(taxonomy_ch)
                                                                .join(genes_ch)
                                                                .join(assembly_ch))
                                                                
        gene_coverage_annotation_and_tax_files_ch = tax_and_cov_ch.map{
                                                         sample_id, coverage -> file("${coverage}")
                                                         }.collect()

        MAKE_COMBINED_GENE_LEVEL_TABLES(gene_coverage_annotation_and_tax_files_ch)

        combined_cov_ch = COMBINE_CONTIG_TAX_AND_COVERAGE(coverage_ch
                                                            .join(taxonomy_ch)
                                                            .join(genes_ch)
                                                            .join(assembly_ch))

        MAKE_COMBINED_CONTIG_TAX_TABLES(combined_cov_ch.map{
                                               sample_id, coverage -> file("${coverage}")
                                               }.collect())

        // Assembly binning
        binning_ch = METABAT_BINNING(assembly_ch.join(read_mapping_ch)) 
        binning_ch | summarize_bins
        metabat_assembly_depth_files_ch = binning_ch.map{
                                           sample_id, depth, bins -> file("${depth}")
                                           }.collect()
        bins_ch = binning_ch.map{
                         sample_id, depth, bins -> bins instanceof List ? bins.each{it}: bins 
                         }.flatten().collect()
   
         
        // Check Bins and Summarize MAGs
        if(gtdbtk_db_dir){
            summarize_mags(summarize_bins.out.bins_checkm_results,
                      bins_ch,
                      gtdbtk_db_dir, use_gtdbtk_scratch_location,
                      gene_coverage_annotation_and_tax_files_ch)
        }else{
            SETUP_GTDBTK_DB()
            summarize_mags(summarize_bins.out.bins_checkm_results,
                      bins_ch, 
                      SETUP_GTDBTK_DB.out.gtdbtk_db_dir, use_gtdbtk_scratch_location,
                      gene_coverage_annotation_and_tax_files_ch)          
        }

        // Get the predicted amino acids for all the samples
        genes_aa_ch = genes_ch.map{sample_id, aa, nt, gff -> file("${aa}")}.collect() 
    
        // Generating a file with sample ids on a new line
        file_ch.map{row -> "${row.sample_id}"}
              .collectFile(name: "${baseDir}/unique-sample-IDs.txt", newLine: true)
              .set{sample_ids_ch}


        bam_files = read_mapping_ch.map{sample_id, bam -> file("${bam}")}.collect()
        // Summarize Assembly-based analysis
        GENERATE_ASSEMBLY_PROCESSING_OVERVIEW_TABLE(sample_ids_ch, summarize_mags.out.MAGs_overview,
                                                    summarize_mags.out.MAGs_dir, assemblies_ch,
                                                    genes_aa_ch,
                                                    metabat_assembly_depth_files_ch,
                                                    bins_ch,
                                                    bam_files)
}
