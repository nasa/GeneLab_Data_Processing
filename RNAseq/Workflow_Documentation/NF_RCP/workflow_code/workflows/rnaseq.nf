include { PARSE_RUNSHEET } from './parse_runsheet.nf'
include { PARSE_ANNOTATIONS_TABLE } from '../modules/parse_annotations_table.nf'
include { FETCH_ISA } from '../modules/fetch_isa.nf'
include { ISA_TO_RUNSHEET } from '../modules/isa_to_runsheet.nf'
include { GET_ACCESSIONS } from '../modules/get_accessions.nf'
// include { PREPARE_REFERENCES } from './prepare_references.nf'
include { DOWNLOAD_REFERENCES } from '../modules/download_references.nf'
include { SUBSAMPLE_GENOME } from '../modules/subsample_genome.nf'
include { DOWNLOAD_ERCC } from '../modules/download_ercc.nf'
include { CONCAT_ERCC } from '../modules/concat_ercc.nf'
include { GTF_TO_PRED } from '../modules/gtf_to_pred.nf'
include { PRED_TO_BED } from '../modules/pred_to_bed.nf'
include { STAGE_RAW_READS } from './stage_raw_reads.nf'
include { FASTQC as RAW_FASTQC } from '../modules/fastqc.nf'
include { GET_MAX_READ_LENGTH } from '../modules/get_max_read_length.nf'
include { TRIMGALORE } from '../modules/trimgalore.nf'
include { FASTQC as TRIMMED_FASTQC } from '../modules/fastqc.nf'
include { BUILD_STAR_INDEX } from '../modules/build_star_index.nf'
include { ALIGN_STAR } from '../modules/align_star.nf'
include { SORT_AND_INDEX_BAM } from '../modules/sort_and_index_bam.nf'
include { INFER_EXPERIMENT } from '../modules/rseqc.nf'
include { GENEBODY_COVERAGE } from '../modules/rseqc.nf'
include { INNER_DISTANCE } from '../modules/rseqc.nf'
include { READ_DISTRIBUTION } from '../modules/rseqc.nf'
include { ASSESS_STRANDEDNESS } from '../modules/assess_strandedness.nf'
include { BUILD_RSEM_INDEX } from '../modules/build_rsem_index.nf'
include { QUANTIFY_STAR_GENES } from '../modules/quantify_star_genes.nf'
include { COUNT_ALIGNED } from '../modules/count_aligned.nf' 
include { EXTRACT_RRNA } from '../modules/extract_rrna.nf'
include { REMOVE_RRNA } from '../modules/remove_rrna.nf'
include { QUANTIFY_RSEM_GENES } from '../modules/quantify_rsem_genes.nf'
include { PARSE_QC_METRICS } from '../modules/parse_qc_metrics.nf'

include { MULTIQC as RAW_READS_MULTIQC } from '../modules/multiqc.nf' addParams(MQCLabel:"raw")
include { MULTIQC as TRIMMED_READS_MULTIQC } from '../modules/multiqc.nf' addParams(MQCLabel:"trimmed")
include { MULTIQC as TRIMMING_MULTIQC } from '../modules/multiqc.nf' addParams(MQCLabel:"trimming")
include { MULTIQC as ALIGN_MULTIQC } from '../modules/multiqc.nf' addParams(MQCLabel:"align")
include { MULTIQC as GENEBODY_COVERAGE_MULTIQC } from '../modules/multiqc.nf' addParams(MQCLabel:"geneBody_cov") //PublishTo: "RSeQC_Analyses/02_geneBody_coverage", 
include { MULTIQC as INFER_EXPERIMENT_MULTIQC } from '../modules/multiqc.nf' addParams(MQCLabel:"infer_exp") //PublishTo: "RSeQC_Analyses/03_infer_experiment", 
include { MULTIQC as INNER_DISTANCE_MULTIQC } from '../modules/multiqc.nf' addParams(MQCLabel:"inner_dist") //PublishTo: "RSeQC_Analyses/04_inner_distance", 
include { MULTIQC as READ_DISTRIBUTION_MULTIQC } from '../modules/multiqc.nf' addParams(MQCLabel:"read_dist") //PublishTo: "RSeQC_Analyses/05_read_distribution",
include { MULTIQC as COUNT_MULTIQC } from '../modules/multiqc.nf' addParams(MQCLabel:"RSEM_count")
include { MULTIQC as QUALIMAP_MULTIQC } from '../modules/multiqc.nf' addParams(MQCLabel:"qualimap") 
include { MULTIQC as ALL_MULTIQC } from '../modules/multiqc.nf' addParams(MQCLabel:"all")
//include { QUALIMAP_BAM_QC } from '../modules/qualimap.nf' not implemented
//include { QUALIMAP_RNASEQ_QC } from '../modules/qualimap.nf' not implemented

include { DGE_DESEQ2 } from '../modules/dge_deseq2.nf'
include { ADD_GENE_ANNOTATIONS } from '../modules/add_gene_annotations.nf'
// include { EXTEND_DGE_TABLE } from '../modules/extend_dge_table.nf'
// include { GENERATE_PCA_TABLE } from '../modules/generate_pca_table.nf'


include { VV_RAW_READS;
          VV_TRIMMED_READS;
          VV_STAR_ALIGNMENTS;
          VV_RSEQC;
          VV_RSEM_COUNTS;
          VV_DESEQ2_ANALYSIS;
          VV_CONCAT_FILTER } from '../modules/vv.nf'

def colorCodes = [
    c_line: "┅" * 70,
    c_back_bright_red: "\u001b[41;1m",
    c_bright_green: "\u001b[32;1m",
    c_blue: "\033[0;34m",
    c_yellow: "\u001b[33;1m",
    c_reset: "\033[0m"
]

workflow RNASEQ {
    take:
        dp_tools_plugin
        annotations_csv_url_string
        accession
        isa_archive_path
        runsheet_path
        api_url
        force_single_end
        truncate_to
        reference_source
        reference_version
        reference_fasta
        reference_gtf
        reference_store_path
        derived_store_path
    main:
        publishdir = "results" // default path passed to publishDir, updated below to "GLDS-#" if processing and OSDR dataset

        // Set up runsheet
        if (runsheet_path == null) {
            GET_ACCESSIONS( accession, api_url ) //Get both OSD and GLDS accessions based on the input accession
            accessions_txt = GET_ACCESSIONS.out.accessions_txt // returns accessions.txt with line1 = osd_accession, line2 = glds_accession. 
            osd_accession = accessions_txt.map { it.readLines()[0].trim() }
            glds_accession = accessions_txt.map { it.readLines()[1].trim() }
            publishdir = accessions_txt.map { it.readLines()[1].trim() }
            //Fetch ISA archive if not provided
            if (isa_archive_path == null) {
                FETCH_ISA(osd_accession, glds_accession)
                isa_archive = FETCH_ISA.out.isa_archive
            }
            //Convert ISA archive to runsheet
            ISA_TO_RUNSHEET( osd_accession, glds_accession, isa_archive, dp_tools_plugin )
            runsheet_path = ISA_TO_RUNSHEET.out.runsheet
        }

        PARSE_RUNSHEET(runsheet_path)

        samples = PARSE_RUNSHEET.out.samples
        //samples | view
        samples | first 
                | map { meta, reads -> meta }
                | set { ch_meta }
        
        has_ercc = ch_meta.map { it.has_ercc }

        ch_meta | map { meta -> meta.organism_sci }
        | set { organism_sci }

        PARSE_ANNOTATIONS_TABLE(annotations_csv_url_string, organism_sci)
        gene_annotations_url = PARSE_ANNOTATIONS_TABLE.out.gene_annotations_url

        // Use manually provided reference genome files if provided. Reference source and version are optional.
        if (params.reference_fasta && params.reference_gtf) {
            genome_references_pre_subsample = Channel.fromPath([params.reference_fasta, params.reference_gtf], checkIfExists: true).toList()
            genome_references_pre_subsample | view
            Channel.value(params.reference_source) | set { reference_source }
            Channel.value(params.reference_version) | set { reference_version }
        } else{
            // Use annotations table to get genome reference files
            DOWNLOAD_REFERENCES(reference_store_path, organism_sci, PARSE_ANNOTATIONS_TABLE.out.reference_source, PARSE_ANNOTATIONS_TABLE.out.reference_version, PARSE_ANNOTATIONS_TABLE.out.reference_fasta_url, PARSE_ANNOTATIONS_TABLE.out.reference_gtf_url)
            genome_references_pre_subsample = DOWNLOAD_REFERENCES.out.reference_files
            reference_source = PARSE_ANNOTATIONS_TABLE.out.reference_source
            reference_version = PARSE_ANNOTATIONS_TABLE.out.reference_version
        }

        // Genomic region subsampling step is used only for debugging / testing 
        if ( params.genome_subsample ) {
            SUBSAMPLE_GENOME( derived_store_path, organism_sci, genome_references_pre_subsample, reference_source, reference_version )
            SUBSAMPLE_GENOME.out.build | flatten | toList | set { genome_references_pre_ercc }
        } else {
            genome_references_pre_subsample | flatten | toList | set { genome_references_pre_ercc }
        }

        // Add ERCC Fasta and GTF to genome files
        DOWNLOAD_ERCC(has_ercc, reference_store_path).ifEmpty([file("ERCC92.fa"), file("ERCC92.gtf")]) | set { ch_maybe_ercc_refs }
        CONCAT_ERCC( reference_store_path, organism_sci, reference_source, reference_version, genome_references_pre_ercc, ch_maybe_ercc_refs, has_ercc )
        .ifEmpty { genome_references_pre_ercc.value }  | set { genome_references }
        
        // Convert GTF file to RSeQC-compatible BED file
        GTF_TO_PRED(
            derived_store_path,
            organism_sci,
            reference_source,
            reference_version,
            genome_references | map { it[1] }
        )
        PRED_TO_BED( 
            derived_store_path,
            organism_sci,
            reference_source,
            reference_version,
            GTF_TO_PRED.out.genome_pred
        )
        genome_bed = PRED_TO_BED.out.genome_bed

        // Metadata and reference files are ready. Stage the raw reads, find the max read length, and build the STAR index.

        // Stage the raw or truncated reads.
        STAGE_RAW_READS( samples )
        raw_reads = STAGE_RAW_READS.out.raw_reads
        samples_txt = STAGE_RAW_READS.out.samples_txt
        //samples_txt | view

        RAW_FASTQC( raw_reads )

        RAW_FASTQC.out.fastqc | map { it -> [ it[1], it[2] ] }
                          | flatten
                          | collect        // Collect all zip files into a single list
                          | set { raw_fastqc_zip }     // Create a channel with all zip files
        
        GET_MAX_READ_LENGTH( raw_fastqc_zip )
        GET_MAX_READ_LENGTH.out.length 
        | map { it.toString().toInteger() }  // ensure it's an integer
        | set { max_read_length }
        //max_read_length.view { "Max read length: $it" }

        // Trim raw reads
        TRIMGALORE( raw_reads )
        trimmed_reads = TRIMGALORE.out.reads
        trimgalore_reports = TRIMGALORE.out.reports | collect

        // Run FastQC on trimmed reads
        TRIMMED_FASTQC( trimmed_reads )
        TRIMMED_FASTQC.out.fastqc | map { it -> [ it[1], it[2] ] } \
                              | flatten \
                              | unique \
                              | collect \
                              | set { trimmed_fastqc_zip }


        // Build STAR genome index
        BUILD_STAR_INDEX(derived_store_path, organism_sci, reference_source, reference_version, genome_references, ch_meta, max_read_length)
        star_index_dir = BUILD_STAR_INDEX.out.index_dir

        // STAR two-pass alignment
        ALIGN_STAR( trimmed_reads, star_index_dir )
        star_alignment_logs = ALIGN_STAR.out.alignment_logs | collect
        transcriptome_aligned_bam = ALIGN_STAR.out.bam_to_transcriptome // STAR transcriptome-aligned BAM file used for transcript-level quantification & abundance estimation with RSEM
        
        // Sort and index genome coordinate-aligned bam files
        SORT_AND_INDEX_BAM( ALIGN_STAR.out.bam_by_coord)
        sorted_bam = SORT_AND_INDEX_BAM.out.sorted_bam

        // RSeQC modules
        INFER_EXPERIMENT( sorted_bam, genome_bed )
        GENEBODY_COVERAGE( sorted_bam, genome_bed )
        INNER_DISTANCE( sorted_bam, genome_bed, max_read_length )
        READ_DISTRIBUTION( sorted_bam, genome_bed )
        infer_expt_out = INFER_EXPERIMENT.out.log | map { it[1] }
                               | collect
        // Parse RSeQC infer_experiment.py results using thresholds set in bin/assess_strandedness.py to determine the strandedness of the dataset
        ASSESS_STRANDEDNESS( infer_expt_out )
        strandedness = ASSESS_STRANDEDNESS.out | map { it.text.split(":")[0] }

        // Combine RSeQC module logs
        ch_rseqc_logs = Channel.empty()
        ch_rseqc_logs | mix(INFER_EXPERIMENT.out.log_only,
                         GENEBODY_COVERAGE.out.all_output,
                         INNER_DISTANCE.out.all_output,
                         READ_DISTRIBUTION.out.log_only)
                   | collect
                   | set{ ch_rseqc_logs }

        // Run Qualimap BAM QC and rnaseq - not implemented
        // QUALIMAP_BAM_QC( sorted_bam, genome_bed, strandedness )
        // qualimap_bamqc_outputs = QUALIMAP_BAM_QC.out.results.toList()
        // QUALIMAP_RNASEQ_QC( sorted_bam, genome_references | map { it[1] }, strandedness )
        // qualimap_rnaseq_qc_outputs = QUALIMAP_RNASEQ_QC.out.results.toList()
        // qualimap_outputs = qualimap_bamqc_outputs
        //                     | concat (qualimap_rnaseq_qc_outputs)

        // Quantify STAR gene counts
        QUANTIFY_STAR_GENES( samples_txt, ALIGN_STAR.out.reads_per_gene | toSortedList, strandedness)

        // Build RSEM transcriptome index
        BUILD_RSEM_INDEX(derived_store_path, organism_sci, reference_source, reference_version, genome_references, ch_meta)
        rsem_index_dir = BUILD_RSEM_INDEX.out.index_dir

        // Run RSEM on the transcriptome-aligned BAM from STAR to calculate transcript expression estimates and quantify gene counts
        COUNT_ALIGNED( transcriptome_aligned_bam, rsem_index_dir, strandedness )
        EXTRACT_RRNA ( organism_sci, genome_references | map { it[1] })
        REMOVE_RRNA ( ch_meta, EXTRACT_RRNA.out.rrna_ids, COUNT_ALIGNED.out.genes_results )
        rsem_counts = COUNT_ALIGNED.out.counts | map { it[1] } | collect
        QUANTIFY_RSEM_GENES( samples_txt, rsem_counts )


        // MultiQC
        ch_multiqc_config = params.multiqc_config ? Channel.fromPath( params.multiqc_config ) : Channel.fromPath("NO_FILE")
        RAW_READS_MULTIQC( samples_txt, raw_fastqc_zip, ch_multiqc_config )
        TRIMMING_MULTIQC( samples_txt, trimgalore_reports, ch_multiqc_config )
        TRIMMED_READS_MULTIQC( samples_txt, trimmed_fastqc_zip, ch_multiqc_config )
        ALIGN_MULTIQC( samples_txt, star_alignment_logs, ch_multiqc_config )
        INFER_EXPERIMENT_MULTIQC( samples_txt, INFER_EXPERIMENT.out.log | map { it[1] } | collect, ch_multiqc_config )
        GENEBODY_COVERAGE_MULTIQC( samples_txt, GENEBODY_COVERAGE.out.log | map { it[1] } | collect, ch_multiqc_config )
        INNER_DISTANCE_MULTIQC( samples_txt, INNER_DISTANCE.out.log | map { it[1] } | collect, ch_multiqc_config )
        READ_DISTRIBUTION_MULTIQC( samples_txt, READ_DISTRIBUTION.out.log | map { it[1] } | collect, ch_multiqc_config )
        // QUALIMAP_MULTIQC ( samples_txt, qualimap_outputs, ch_multiqc_config)
        COUNT_MULTIQC( samples_txt, rsem_counts, ch_multiqc_config )
        all_multiqc_input = raw_fastqc_zip
                    | concat(trimgalore_reports)
                    | concat(trimmed_fastqc_zip)
                    | concat(star_alignment_logs)
                    | concat(INFER_EXPERIMENT.out.log | map { it[1] } | collect)
                    | concat(GENEBODY_COVERAGE.out.log | map { it[1] } | collect)
                    | concat(INNER_DISTANCE.out.log | map { it[1] } | collect)
                    | concat(READ_DISTRIBUTION.out.log | map { it[1] } | collect)
                    // | concat(qualimap_outputs)
                    | concat(rsem_counts)
                    | collect
        ALL_MULTIQC( samples_txt, all_multiqc_input, ch_multiqc_config )
        

        // Normalize counts, DGE 
        DGE_DESEQ2( ch_meta, runsheet_path, COUNT_ALIGNED.out.genes_results | toSortedList )
        dge_table = DGE_DESEQ2.out.dge_table
        // Add annotations to DGE table
        ADD_GENE_ANNOTATIONS( ch_meta, gene_annotations_url, dge_table )
        annotated_dge_table = ADD_GENE_ANNOTATIONS.out.annotated_dge_table
        // Extend DGE table to generate visualization table
        // Step being removed on update
        //EXTEND_DGE_TABLE( annotated_dge_table )
        // Generate PCA table from normalized counts 
        // Step being removed on update
        //GENERATE_PCA_TABLE ( DGE_DESEQ2.out.norm_counts | map { it[1] })

        // Parse QC metrics
        all_multiqc_output = RAW_READS_MULTIQC.out.data
            | concat(TRIMMED_READS_MULTIQC.out.data)
            | concat(ALIGN_MULTIQC.out.data)
            | concat(GENEBODY_COVERAGE_MULTIQC.out.data)
            | concat(INFER_EXPERIMENT_MULTIQC.out.data)
            | concat(INNER_DISTANCE_MULTIQC.out.data)
            | concat(READ_DISTRIBUTION_MULTIQC.out.data)
            | concat(COUNT_MULTIQC.out.data)
            | collect

        PARSE_QC_METRICS(
            publishdir,
            osd_accession,
            ch_meta,
            isa_archive,
            all_multiqc_output,
            QUANTIFY_RSEM_GENES.out.publishables
        )

        // Run dp_tools V&V modules on workflow outputs and publish them
        // VV_RAW_READS( 
        //     publishdir,
        //     ch_meta,
        //     runsheet_path,
        //     raw_reads | map{ it -> it[1] } | collect,
        //     raw_fastqc_zip,
        //     RAW_READS_MULTIQC.out.unzipped_report,
        //     RAW_READS_MULTIQC.out.zipped_report,
        //     "${ projectDir }/bin/dp_tools__NF_RCP" // dp_tools plugin
        // )
        // VV_TRIMMED_READS( 
        //     publishdir,
        //     ch_meta,
        //     runsheet_path,
        //     trimmed_reads | map{ it -> it[1] } | collect,
        //     trimmed_fastqc_zip,
        //     TRIMMED_READS_MULTIQC.out.unzipped_report,
        //     TRIMMED_READS_MULTIQC.out.zipped_report,
        //     "${ projectDir }/bin/dp_tools__NF_RCP" // dp_tools plugin
        // )
        // VV_STAR_ALIGNMENTS(
        //     publishdir,
        //     runsheet_path,
        //     ALIGN_STAR.out.publishables | collect,
        //     QUANTIFY_STAR_GENES.out.publishables | collect,
        //     ALIGN_MULTIQC.out.zipped_report,
        //     ALIGN_MULTIQC.out.unzipped_report,
        //     SORT_AND_INDEX_BAM.out.bam_bed | collect,
        //     "${ projectDir }/bin/dp_tools__NF_RCP" // dp_tools plugin
        // )
        // VV_RSEQC(
        //     publishdir,
        //     ch_meta,
        //     runsheet_path,
        //     ch_rseqc_logs,
        //     GENEBODY_COVERAGE_MULTIQC.out.zipped_report
        //         .mix( GENEBODY_COVERAGE_MULTIQC.out.unzipped_report )
        //         .collect(),
        //     INFER_EXPERIMENT_MULTIQC.out.zipped_report
        //         .mix( INFER_EXPERIMENT_MULTIQC.out.unzipped_report )
        //         .collect(),
        //     INNER_DISTANCE_MULTIQC.out.zipped_report
        //         .mix( INNER_DISTANCE_MULTIQC.out.unzipped_report )
        //         .collect(),
        //     READ_DISTRIBUTION_MULTIQC.out.zipped_report
        //         .mix( READ_DISTRIBUTION_MULTIQC.out.unzipped_report )
        //         .collect(),
        //     "${ projectDir }/bin/dp_tools__NF_RCP" // dp_tools plugin
        // )
        // VV_RSEM_COUNTS(
        //     publishdir,
        //     runsheet_path,
        //     COUNT_ALIGNED.out.only_counts | collect,
        //     QUANTIFY_RSEM_GENES.out.publishables,
        //     COUNT_MULTIQC.out.zipped_report,
        //     COUNT_MULTIQC.out.unzipped_report,
        //     "${ projectDir }/bin/dp_tools__NF_RCP" // dp_tools plugin
        // )
        // VV_DESEQ2_ANALYSIS(
        //     publishdir,
        //     ch_meta,
        //     runsheet_path,
        //     QUANTIFY_RSEM_GENES.out.publishables,
        //     COUNT_MULTIQC.out.zipped_report,
        //     COUNT_MULTIQC.out.unzipped_report,
        //     DGE_DESEQ2.out.norm_counts,
        //     DGE_DESEQ2.out.contrasts
        //         .mix( DGE_DESEQ2.out.sample_table )
        //         .mix( annotated_dge_table )
        //         .mix( GENERATE_PCA_TABLE.out.pca_table ),
        //     DGE_DESEQ2.out.norm_counts_ercc | ifEmpty( { file("NO_FILES.placeholder") }),
        //     DGE_BY_DESEQ2.out.dge_ercc | ifEmpty( { file("NO_FILES.placeholder") }),
        //     "${ projectDir }/bin/dp_tools__NF_RCP" // dp_tools plugin
        //             )
    emit:
        PARSE_QC_METRICS.out.file
}