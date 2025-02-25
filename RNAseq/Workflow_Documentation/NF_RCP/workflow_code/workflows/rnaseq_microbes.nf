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
include { GTF_TO_BED } from '../modules/gtf_to_bed.nf'
include { STAGE_RAW_READS } from './stage_raw_reads.nf'
include { FASTQC as RAW_FASTQC } from '../modules/fastqc.nf'
include { GET_MAX_READ_LENGTH } from '../modules/get_max_read_length.nf'
include { TRIMGALORE } from '../modules/trimgalore.nf'
include { FASTQC as TRIMMED_FASTQC } from '../modules/fastqc.nf'
include { BUILD_BOWTIE2_INDEX } from '../modules/build_bowtie2_index.nf'
include { ALIGN_BOWTIE2 } from '../modules/align_bowtie2.nf'
include { SAM_TO_BAM } from '../modules/sam_to_bam.nf'
include { SORT_AND_INDEX_BAM } from '../modules/sort_and_index_bam.nf'
include { INFER_EXPERIMENT } from '../modules/rseqc.nf'
include { GENEBODY_COVERAGE } from '../modules/rseqc.nf'
include { INNER_DISTANCE } from '../modules/rseqc.nf'
include { READ_DISTRIBUTION } from '../modules/rseqc.nf'
include { ASSESS_STRANDEDNESS } from '../modules/assess_strandedness.nf'

include { MULTIQC as RAW_READS_MULTIQC } from '../modules/multiqc.nf'
include { MULTIQC as TRIMMED_READS_MULTIQC } from '../modules/multiqc.nf'
include { MULTIQC as TRIMMING_MULTIQC } from '../modules/multiqc.nf'
include { MULTIQC as ALIGN_MULTIQC } from '../modules/multiqc.nf'
include { MULTIQC as GENEBODY_COVERAGE_MULTIQC } from '../modules/multiqc.nf'
include { MULTIQC as INFER_EXPERIMENT_MULTIQC } from '../modules/multiqc.nf'
include { MULTIQC as INNER_DISTANCE_MULTIQC } from '../modules/multiqc.nf'
include { MULTIQC as READ_DISTRIBUTION_MULTIQC } from '../modules/multiqc.nf'
include { MULTIQC as COUNT_MULTIQC } from '../modules/multiqc.nf'
include { MULTIQC as ALL_MULTIQC } from '../modules/multiqc.nf'

include { CLEAN_MULTIQC_PATHS as CLEAN_RAW_READS_MULTIQC_PATHS } from '../modules/clean_multiqc_paths.nf'
include { CLEAN_MULTIQC_PATHS as CLEAN_TRIMMED_READS_MULTIQC_PATHS } from '../modules/clean_multiqc_paths.nf'
include { CLEAN_MULTIQC_PATHS as CLEAN_TRIMMING_MULTIQC_PATHS } from '../modules/clean_multiqc_paths.nf'
include { CLEAN_MULTIQC_PATHS as CLEAN_ALIGN_MULTIQC_PATHS } from '../modules/clean_multiqc_paths.nf'
include { CLEAN_MULTIQC_PATHS as CLEAN_GENEBODY_COVERAGE_MULTIQC_PATHS } from '../modules/clean_multiqc_paths.nf'
include { CLEAN_MULTIQC_PATHS as CLEAN_INFER_EXPERIMENT_MULTIQC_PATHS } from '../modules/clean_multiqc_paths.nf'
include { CLEAN_MULTIQC_PATHS as CLEAN_INNER_DISTANCE_MULTIQC_PATHS } from '../modules/clean_multiqc_paths.nf'
include { CLEAN_MULTIQC_PATHS as CLEAN_READ_DISTRIBUTION_MULTIQC_PATHS } from '../modules/clean_multiqc_paths.nf'
include { CLEAN_MULTIQC_PATHS as CLEAN_COUNT_MULTIQC_PATHS } from '../modules/clean_multiqc_paths.nf'
include { CLEAN_MULTIQC_PATHS as CLEAN_ALL_MULTIQC_PATHS } from '../modules/clean_multiqc_paths.nf'

// include { QUALIMAP_BAM_QC } from '../modules/qualimap.nf'
// include { QUALIMAP_RNASEQ_QC } from '../modules/qualimap.nf'
include { GET_GTF_FEATURES } from '../modules/get_gtf_features.nf'
include { FEATURECOUNTS } from '../modules/featurecounts.nf'
include { EXTRACT_RRNA } from '../modules/extract_rrna.nf'
include { REMOVE_RRNA_FEATURECOUNTS } from '../modules/remove_rrna_featurecounts.nf'
include { DGE_DESEQ2 } from '../modules/dge_deseq2.nf'
include { ADD_GENE_ANNOTATIONS } from '../modules/add_gene_annotations.nf'
//include { EXTEND_DGE_TABLE } from '../modules/extend_dge_table.nf'
//include { GENERATE_PCA_TABLE } from '../modules/generate_pca_table.nf'
include { SOFTWARE_VERSIONS } from '../modules/software_versions.nf'
include { MD5SUM as RAW_MD5SUM } from '../modules/md5sum.nf'
include { MD5SUM as PROCESSED_MD5SUM } from '../modules/md5sum.nf'
include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'
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

workflow RNASEQ_MICROBES {
    take:
        ch_outdir
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
        // Parse accession, structure output directory as:
        // params.outdir/
        //   ├── GLDS-#/results/   # Main pipeline results
        //   └── nextflow_info/    # Pipeline execution metadata
        if ( accession ) {
            GET_ACCESSIONS( accession, api_url )
            osd_accession = GET_ACCESSIONS.out.accessions_txt.map { it.readLines()[0].trim() }
            glds_accession = GET_ACCESSIONS.out.accessions_txt.map { it.readLines()[1].trim() }
            ch_outdir = ch_outdir.combine(glds_accession).map { outdir, glds -> "$outdir/$glds" }
        }
        else {
            ch_outdir = ch_outdir.map { it + "/results" }
        }

        // if runsheet_path is not provided, set it up from ISA input
        // If ISA input is not provided, use the accession to get the ISA
        if ( runsheet_path == null ) {
            if ( isa_archive_path == null ) {
                FETCH_ISA( ch_outdir, osd_accession, glds_accession )
                isa_archive = FETCH_ISA.out.isa_archive
            }
            ISA_TO_RUNSHEET( ch_outdir, osd_accession, glds_accession, isa_archive, dp_tools_plugin )
            runsheet_path = ISA_TO_RUNSHEET.out.runsheet
        }

        // Validate input parameters and runsheet
        validateParameters()

        // Print summary of supplied parameters
        log.info paramsSummaryLog(workflow)

        // Parse the runsheet
        PARSE_RUNSHEET( runsheet_path )

        // Get samples from runsheet
        samples = PARSE_RUNSHEET.out.samples
        //samples | view

        // Get dataset-wide metadata
        samples | first 
                | map { meta, reads -> meta }
                | set { ch_meta }

        // Set metadata 
        ch_meta | map { meta -> meta.organism_sci }
        | set { organism_sci }

        PARSE_ANNOTATIONS_TABLE( annotations_csv_url_string, organism_sci )

        // Use manually provided reference genome files if provided. Reference source and version are optional.
        if ( params.reference_fasta && params.reference_gtf ) {
            genome_references_pre_subsample = Channel.fromPath([params.reference_fasta, params.reference_gtf], checkIfExists: true ).toList()
            genome_references_pre_subsample | view
            Channel.value( params.reference_source ) | set { reference_source }
            Channel.value( params.reference_version ) | set { reference_version }
        } else{
            // Use annotations table to get genome reference files
            DOWNLOAD_REFERENCES( reference_store_path, organism_sci, PARSE_ANNOTATIONS_TABLE.out.reference_source, PARSE_ANNOTATIONS_TABLE.out.reference_version, PARSE_ANNOTATIONS_TABLE.out.reference_fasta_url, PARSE_ANNOTATIONS_TABLE.out.reference_gtf_url )
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
        DOWNLOAD_ERCC( ch_meta.map { it.has_ercc }, reference_store_path ).ifEmpty([file("ERCC92.fa"), file("ERCC92.gtf")]) | set { ch_maybe_ercc_refs }
        CONCAT_ERCC( reference_store_path, organism_sci, reference_source, reference_version, genome_references_pre_ercc, ch_maybe_ercc_refs, ch_meta.map { it.has_ercc } )
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

        // Stage the raw or truncated reads.
        STAGE_RAW_READS( samples )
        raw_reads = STAGE_RAW_READS.out.raw_reads
        samples_txt = STAGE_RAW_READS.out.samples_txt
        //samples_txt | view

        // Run FastQC on raw reads
        RAW_FASTQC( raw_reads )
        // Collect the raw read fastqc zip files
        RAW_FASTQC.out.fastqc | map { it -> [ it[1], it[2] ] }
        | flatten
        | collect        // Collect all zip files into a single list
        | set { raw_fastqc_zip }     // Create a channel with all zip files
        
        // Get the max read length by parsing the raw read fastqc zip files
        GET_MAX_READ_LENGTH( raw_fastqc_zip )
        max_read_length = GET_MAX_READ_LENGTH.out.length | map { it.toString().toInteger() }
        //max_read_length.view { "Max read length: $it" }

        // Trim raw reads
        TRIMGALORE( raw_reads )
        trimmed_reads = TRIMGALORE.out.reads
        trimgalore_reports = TRIMGALORE.out.reports | collect

        // Run FastQC on trimmed reads
        TRIMMED_FASTQC( trimmed_reads )
        TRIMMED_FASTQC.out.fastqc | map { it -> [ it[1], it[2] ] } 
        | flatten 
        | collect 
        | set { trimmed_fastqc_zip }

        // Build Bowtie 2 genome index
        BUILD_BOWTIE2_INDEX( derived_store_path, organism_sci, reference_source, reference_version, genome_references, ch_meta )
        bowtie2_index_dir = BUILD_BOWTIE2_INDEX.out.index_dir

        // Align reads using Bowtie2
        ALIGN_BOWTIE2( trimmed_reads, bowtie2_index_dir )
        bowtie2_alignment_logs = ALIGN_BOWTIE2.out.alignment_logs | collect
        
        // Convert Bowtie2 SAM to BAM (query-name order, matching FASTQ input order )
        SAM_TO_BAM( ALIGN_BOWTIE2.out.sam ) 

        // Sort and index BAM files to convert from query-name order to genome coordinate order
        SORT_AND_INDEX_BAM( SAM_TO_BAM.out.bam )
        sorted_bam = SORT_AND_INDEX_BAM.out.sorted_bam
        bams = sorted_bam.map { it[1] } | toSortedList()

        // RSeQC modules
        INFER_EXPERIMENT( sorted_bam, genome_bed )
        GENEBODY_COVERAGE( sorted_bam, genome_bed )
        INNER_DISTANCE( sorted_bam, genome_bed, max_read_length )
        READ_DISTRIBUTION( sorted_bam, genome_bed )
        infer_expt_out = INFER_EXPERIMENT.out.log | map { it[1] }
        | collect
        
        // Combine RSeQC module logs
        ch_rseqc_logs = Channel.empty()
        ch_rseqc_logs 
        | mix( INFER_EXPERIMENT.out.log_only,
                GENEBODY_COVERAGE.out.all_output,
                INNER_DISTANCE.out.all_output,
                READ_DISTRIBUTION.out.log_only )
                | collect
                | set{ ch_rseqc_logs }

        // Parse RSeQC infer_experiment.py results using thresholds set in bin/assess_strandedness.py to determine the strandedness of the dataset
        ASSESS_STRANDEDNESS( infer_expt_out )
        strandedness = ASSESS_STRANDEDNESS.out | map { it.text.split(":")[0] }

        // Generate gene counts table from genome-coordinate sorted Bowtie2-aligned BAMs using featureCounts
        GET_GTF_FEATURES( genome_references )
        gtf_features = GET_GTF_FEATURES.out.gtf_features.map { it.text.trim() }
        FEATURECOUNTS( ch_meta, genome_references, gtf_features, strandedness, bams )
        counts = FEATURECOUNTS.out.counts

        // Find rRNA sequences in the genome and remove them from the counts table
        EXTRACT_RRNA( organism_sci, genome_references | map { it[1] })
        REMOVE_RRNA_FEATURECOUNTS ( counts, EXTRACT_RRNA.out.rrna_ids )

        // Normalize counts, DGE 
        DGE_DESEQ2( ch_meta, runsheet_path, counts )

        // Add annotations to DGE table
        ADD_GENE_ANNOTATIONS( ch_meta, PARSE_ANNOTATIONS_TABLE.out.gene_annotations_url, DGE_DESEQ2.out.dge_table )
        annotated_dge_table = ADD_GENE_ANNOTATIONS.out.annotated_dge_table

        // MultiQC
        ch_multiqc_config = params.multiqc_config ? Channel.fromPath( params.multiqc_config ) : Channel.fromPath("NO_FILE")
        RAW_READS_MULTIQC(samples_txt, raw_fastqc_zip, ch_multiqc_config, "raw")
        TRIMMING_MULTIQC(samples_txt, trimgalore_reports, ch_multiqc_config, "trimming")
        TRIMMED_READS_MULTIQC(samples_txt, trimmed_fastqc_zip, ch_multiqc_config, "trimmed")
        ALIGN_MULTIQC(samples_txt, bowtie2_alignment_logs, ch_multiqc_config, "align")
        INFER_EXPERIMENT_MULTIQC(samples_txt, INFER_EXPERIMENT.out.log | map { it[1] } | collect, ch_multiqc_config, "infer_exp")
        GENEBODY_COVERAGE_MULTIQC(samples_txt, GENEBODY_COVERAGE.out.log | map { it[1] } | collect, ch_multiqc_config, "geneBody_cov")
        INNER_DISTANCE_MULTIQC(samples_txt, INNER_DISTANCE.out.log | map { it[1] } | collect, ch_multiqc_config, "inner_dist")
        READ_DISTRIBUTION_MULTIQC(samples_txt, READ_DISTRIBUTION.out.log | map { it[1] } | collect, ch_multiqc_config, "read_dist")
        COUNT_MULTIQC(samples_txt, FEATURECOUNTS.out.summary, ch_multiqc_config, "featureCounts")
        

        // Software Version Capturing
        nf_version = '"NEXTFLOW":\n    nextflow: '.concat("${nextflow.version}\n")
        ch_nextflow_version = Channel.value(nf_version)
        ch_software_versions = Channel.empty()
        // Mix in versions from each process
        ch_software_versions = ch_software_versions
            | mix(ISA_TO_RUNSHEET.out.versions)  
            | mix(GTF_TO_PRED.out.versions)
            | mix(PRED_TO_BED.out.versions)
            | mix(RAW_FASTQC.out.versions)
            | mix(TRIMGALORE.out.versions)
            | mix(ALIGN_BOWTIE2.out.versions)
            | mix(SAM_TO_BAM.out.versions)
            | mix(INFER_EXPERIMENT.out.versions)
            | mix(GENEBODY_COVERAGE.out.versions)
            | mix(INNER_DISTANCE.out.versions)
            | mix(READ_DISTRIBUTION.out.versions)
            | mix(FEATURECOUNTS.out.versions)
            | mix(RAW_READS_MULTIQC.out.versions)
            | mix(DGE_DESEQ2.out.versions)
            | mix(ch_nextflow_version)
        // Process the versions:
        ch_software_versions 
            | unique  
            | collectFile(
                newLine: true, 
                cache: false
            )
            | set { ch_final_software_versions }
        // Convert software versions combined yaml to markdown table
        SOFTWARE_VERSIONS(ch_final_software_versions)

        // Generate md5sums for raw and processed data
        RAW_MD5SUM( 
            STAGE_RAW_READS.out.ch_all_raw_reads
            | concat (RAW_READS_MULTIQC.out.zipped_report)
            | collect,
            "raw_"
        )

        // PROCESSED_MD5SUM(
        //     TRIMMED_READS_MULTIQC.out.zipped_report,
        //     "processed_"
        // )

        // Run Qualimap BAM QC and rnaseq
        // QUALIMAP_BAM_QC( sorted_bam, genome_bed, strandedness )
        // QUALIMAP_RNASEQ_QC( sorted_bam, genome_references | map { it[1] }, strandedness )
        // qualimap_outputs = QUALIMAP_BAM_QC.out.results
        //             // | concat(QUALIMAP_RNASEQ_QC.out.results )
        //             | collect

        // Clean paths in outputs before VVing & publishing
        CLEAN_RAW_READS_MULTIQC_PATHS(RAW_READS_MULTIQC.out.zipped_report, "raw")
        CLEAN_TRIMMED_READS_MULTIQC_PATHS(TRIMMED_READS_MULTIQC.out.zipped_report, "trimmed")
        CLEAN_TRIMMING_MULTIQC_PATHS(TRIMMING_MULTIQC.out.zipped_report, "trimming")
        CLEAN_ALIGN_MULTIQC_PATHS(ALIGN_MULTIQC.out.zipped_report, "align")
        CLEAN_INFER_EXPERIMENT_MULTIQC_PATHS(INFER_EXPERIMENT_MULTIQC.out.zipped_report, "infer_exp")
        CLEAN_GENEBODY_COVERAGE_MULTIQC_PATHS(GENEBODY_COVERAGE_MULTIQC.out.zipped_report, "geneBody_cov")
        CLEAN_INNER_DISTANCE_MULTIQC_PATHS(INNER_DISTANCE_MULTIQC.out.zipped_report, "inner_dist")
        CLEAN_READ_DISTRIBUTION_MULTIQC_PATHS(READ_DISTRIBUTION_MULTIQC.out.zipped_report, "read_dist")

        VV_RAW_READS(
            dp_tools_plugin,
            ch_outdir,
            ch_meta,
            runsheet_path,
            raw_reads | map{ it -> it[1] } | collect,
            raw_fastqc_zip,
            CLEAN_RAW_READS_MULTIQC_PATHS.out.zipped_report
        )
    emit:
        CLEAN_READ_DISTRIBUTION_MULTIQC_PATHS.out.zipped_report
}