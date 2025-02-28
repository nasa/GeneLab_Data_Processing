#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/**************************************************************************************** 
***************************  Assembly-based processing workflow *************************
****************************************************************************************/

// Assembly-based workflow
// Processes to create the required database(s) if not provided
include { SETUP_CAT_DB; SETUP_KOFAMSCAN_DB; SETUP_GTDBTK_DB;
          SETUP_CHOCOPHLAN} from "./database_creation.nf"

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

        software_versions_ch = Channel.empty()

        // Assemble reads to contigs
        ASSEMBLE(filtered_ch)
        ASSEMBLE.out.contigs | RENAME_HEADERS
        assembly_ch = RENAME_HEADERS.out.contigs
        assemblies_ch = assembly_ch.map{
                                      sample_id, assembly -> file("${assembly}")
                                      }.collect()
        SUMMARIZE_ASSEMBLIES(assemblies_ch)

        // Write failed assemblies to a Failed assemblies file
        failed_assemblies = RENAME_HEADERS.out.failed_assembly
        failed_assemblies
              .map{ it.text }
              .collectFile(name: "${params.assemblies_dir}/Failed-assemblies.tsv", cache: false)
        
        // Map reads to assembly
        MAPPING(assembly_ch.join(filtered_ch))
        MAPPING.out.sam | SAM_TO_BAM
        read_mapping_ch = SAM_TO_BAM.out.bam

        // Annotate assembly
        CALL_GENES(assembly_ch)
        CALL_GENES.out.genes | REMOVE_LINEWRAPS
        genes_ch = REMOVE_LINEWRAPS.out.genes

        if (ko_db_dir != null){

            KO_ANNOTATION(assembly_ch.join(genes_ch), ko_db_dir)
            KO_ANNOTATION.out.temp_table | FILTER_KFAMSCAN
            annotations_ch = FILTER_KFAMSCAN.out.ko_annotation

        }else{

            SETUP_KOFAMSCAN_DB()
            SETUP_KOFAMSCAN_DB.out.version | mix(software_versions_ch) | set{software_versions_ch}
            KO_ANNOTATION(assembly_ch.join(genes_ch), SETUP_KOFAMSCAN_DB.out.ko_db_dir)
            KO_ANNOTATION.out.temp_table | FILTER_KFAMSCAN
            annotations_ch = FILTER_KFAMSCAN.out.ko_annotation

        }

        if (cat_db != null){         

            TAX_CLASSIFICATION(assembly_ch.join(genes_ch), cat_db)
            taxonomy_ch = TAX_CLASSIFICATION.out.taxonomy

        }else{

            SETUP_CAT_DB(params.database.CAT_DB_LINK)
            SETUP_CAT_DB.out.version | mix(software_versions_ch) | set{software_versions_ch}
            TAX_CLASSIFICATION(assembly_ch.join(genes_ch), SETUP_CAT_DB.out.cat_db)
            taxonomy_ch = TAX_CLASSIFICATION.out.taxonomy
        }

        // Calculate gene coverage and depth
        GET_COV_AND_DET(read_mapping_ch
                           .join(assembly_ch)
                           .join(genes_ch))
        coverage_ch = GET_COV_AND_DET.out.coverages

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
        METABAT_BINNING(assembly_ch.join(read_mapping_ch))
        binning_ch = METABAT_BINNING.out.bins
        binning_ch | summarize_bins
        depth_ch = METABAT_BINNING.out.depth
        metabat_assembly_depth_files_ch = depth_ch.map{
                                           sample_id, depth -> file("${depth}")
                                           }.collect()
        bins_ch = binning_ch.map{
                         sample_id, bins -> bins instanceof List ? bins.each{it}: bins 
                         }.flatten().collect()
   
         
        // Check Bins and Summarize MAGs
        if(gtdbtk_db_dir != null){
            summarize_mags(summarize_bins.out.bins_checkm_results,
                      bins_ch,
                      gtdbtk_db_dir, use_gtdbtk_scratch_location,
                      gene_coverage_annotation_and_tax_files_ch)
        }else{
            SETUP_GTDBTK_DB(params.database.GTDBTK_LINK)
            SETUP_GTDBTK_DB.out.version | mix(software_versions_ch) | set{software_versions_ch}
            summarize_mags(summarize_bins.out.bins_checkm_results,
                      bins_ch, 
                      SETUP_GTDBTK_DB.out.gtdbtk_db_dir, use_gtdbtk_scratch_location,
                      gene_coverage_annotation_and_tax_files_ch)          
        }

        // Get the predicted amino acids for all the samples
        genes_aa_ch = genes_ch.map{sample_id, aa, nt -> file("${aa}")}.collect() 
    
        // Generating a file with sample ids on a new line
        file_ch.map{row -> "${row.sample_id}"}
              .collectFile(name: "${launchDir}/unique-sample-IDs.txt", newLine: true)
              .set{sample_ids_ch}


        bam_files = read_mapping_ch.map{sample_id, bam -> file("${bam}")}.collect()
        // Summarize Assembly-based analysis
        GENERATE_ASSEMBLY_PROCESSING_OVERVIEW_TABLE(sample_ids_ch, summarize_mags.out.MAGs_overview,
                                                    summarize_mags.out.MAGs_dir, assemblies_ch,
                                                    genes_aa_ch,
                                                    metabat_assembly_depth_files_ch,
                                                    bins_ch,
                                                    bam_files)
   
        // Capture software versions 
        ASSEMBLE.out.version | mix(software_versions_ch) | set{software_versions_ch}
        RENAME_HEADERS.out.version | mix(software_versions_ch) | set{software_versions_ch}
        SUMMARIZE_ASSEMBLIES.out.version | mix(software_versions_ch) | set{software_versions_ch}
        MAPPING.out.version | mix(software_versions_ch) | set{software_versions_ch}
        SAM_TO_BAM.out.version | mix(software_versions_ch) | set{software_versions_ch}
        CALL_GENES.out.version | mix(software_versions_ch) | set{software_versions_ch}
        REMOVE_LINEWRAPS.out.version | mix(software_versions_ch) | set{software_versions_ch}
        KO_ANNOTATION.out.version | mix(software_versions_ch) | set{software_versions_ch}
        FILTER_KFAMSCAN.out.version | mix(software_versions_ch) | set{software_versions_ch}
        TAX_CLASSIFICATION.out.version | mix(software_versions_ch) | set{software_versions_ch}
        GET_COV_AND_DET.out.version | mix(software_versions_ch) | set{software_versions_ch}
        MAKE_COMBINED_GENE_LEVEL_TABLES.out.version | mix(software_versions_ch) | set{software_versions_ch}
        MAKE_COMBINED_CONTIG_TAX_TABLES.out.version | mix(software_versions_ch) | set{software_versions_ch}
        METABAT_BINNING.out.version | mix(software_versions_ch) | set{software_versions_ch}
        summarize_bins.out.versions | mix(software_versions_ch) | set{software_versions_ch}
        summarize_mags.out.versions | mix(software_versions_ch) | set{software_versions_ch}


    emit:
        versions = software_versions_ch

}
