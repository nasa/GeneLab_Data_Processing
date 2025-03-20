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
include { PARSE_QC_METRICS } from '../modules/parse_qc_metrics.nf'

// include { CLEAN_MULTIQC as CLEAN_RAW_READS_MULTIQC} from '../modules/clean_multiqc.nf'
// include { CLEAN_MULTIQC as CLEAN_TRIMMED_READS_MULTIQC } from '../modules/clean_multiqc.nf'
// include { CLEAN_MULTIQC as CLEAN_TRIMMING_MULTIQC } from '../modules/clean_multiqc.nf'
// include { CLEAN_MULTIQC as CLEAN_ALIGN_MULTIQC } from '../modules/clean_multiqc.nf'
// include { CLEAN_MULTIQC as CLEAN_GENEBODY_COVERAGE_MULTIQC } from '../modules/clean_multiqc.nf'
// include { CLEAN_MULTIQC as CLEAN_INFER_EXPERIMENT_MULTIQC } from '../modules/clean_multiqc.nf'
// include { CLEAN_MULTIQC as CLEAN_INNER_DISTANCE_MULTIQC } from '../modules/clean_multiqc.nf'
// include { CLEAN_MULTIQC as CLEAN_READ_DISTRIBUTION_MULTIQC } from '../modules/clean_multiqc.nf'
// include { CLEAN_MULTIQC as CLEAN_COUNT_MULTIQC } from '../modules/clean_multiqc.nf'
// include { CLEAN_MULTIQC as CLEAN_ALL_MULTIQC } from '../modules/clean_multiqc.nf'

// include { QUALIMAP_BAM_QC } from '../modules/qualimap.nf'
// include { QUALIMAP_RNASEQ_QC } from '../modules/qualimap.nf'
include { GET_GTF_FEATURES } from '../modules/get_gtf_features.nf'
include { FEATURECOUNTS } from '../modules/featurecounts.nf'
include { QUANTIFY_FEATURECOUNTS_GENES } from '../modules/quantify_featurecounts_genes.nf'
include { EXTRACT_RRNA } from '../modules/extract_rrna.nf'
include { REMOVE_RRNA_FEATURECOUNTS } from '../modules/remove_rrna_featurecounts.nf'
include { DGE_DESEQ2 } from '../modules/dge_deseq2.nf'
include { DGE_DESEQ2 as DGE_DESEQ2_RRNA_RM } from '../modules/dge_deseq2.nf'
include { ADD_GENE_ANNOTATIONS } from '../modules/add_gene_annotations.nf'
include { ADD_GENE_ANNOTATIONS as ADD_GENE_ANNOTATIONS_RRNA_RM } from '../modules/add_gene_annotations.nf'
//include { EXTEND_DGE_TABLE } from '../modules/extend_dge_table.nf'
//include { GENERATE_PCA_TABLE } from '../modules/generate_pca_table.nf'
include { SOFTWARE_VERSIONS } from '../modules/software_versions.nf'
include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'
include { VV_RAW_READS;
    VV_TRIMMED_READS;
    VV_BOWTIE2_ALIGNMENT;
    VV_RSEQC;
    VV_FEATURECOUNTS;
    VV_DGE_DESEQ2 } from '../modules/vv.nf'
include { GENERATE_PROTOCOL } from '../modules/generate_protocol.nf'

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
        //   ├── [GLDS-#|results]/ # Main pipeline results
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

        // Ensure ch_outdir is a proper value channel that can be reused multiple times
        ch_outdir = ch_outdir.first()

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
        // log.info paramsSummaryLog(workflow)

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
        //STAGE_RAW_READS( ch_outdir.map { it + "/00-RawData/Fastq" }, samples )
        raw_reads = STAGE_RAW_READS.out.raw_reads
        samples_txt = STAGE_RAW_READS.out.samples_txt
        //samples_txt | view

        // Run FastQC on raw reads
        RAW_FASTQC( raw_reads )
        // Collect the raw read fastqc zip files
        RAW_FASTQC.out.fastqc | map { it -> [ it[1], it[2] ] }
        //RAW_FASTQC( ch_outdir.map { it + "/00-RawData/FastQC_Reports" }, raw_reads )
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
        //TRIMMED_FASTQC( ch_outdir.map { it + "/01-TG_Preproc/FastQC_Reports" }, trimmed_reads )
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
        
        QUANTIFY_FEATURECOUNTS_GENES( samples_txt, FEATURECOUNTS.out.counts )

        // Use the GTF to find rRNA genes, remove them from the counts table
        EXTRACT_RRNA( organism_sci, genome_references | map { it[1] })
        REMOVE_RRNA_FEATURECOUNTS ( counts, EXTRACT_RRNA.out.rrna_ids )

        dge_script = "${projectDir}/bin/dge_deseq2.Rmd"
        dge_annotations_script = "${projectDir}/bin/add_gene_annotations.Rmd"

        // Normalize counts, DGE, Add annotations to DGE table
        DGE_DESEQ2( ch_meta, runsheet_path, counts, dge_script, "" )
        ADD_GENE_ANNOTATIONS( ch_meta, PARSE_ANNOTATIONS_TABLE.out.gene_annotations_url, DGE_DESEQ2.out.dge_table, dge_annotations_script, "" )
        annotated_dge_table = ADD_GENE_ANNOTATIONS.out.annotated_dge_table
        
        // For rRNArm counts: Normalize counts, DGE, Add annotations to DGE table
        DGE_DESEQ2_RRNA_RM( ch_meta, runsheet_path, REMOVE_RRNA_FEATURECOUNTS.out.counts_rrnarm, dge_script, "_rRNArm" )
        ADD_GENE_ANNOTATIONS_RRNA_RM( ch_meta, PARSE_ANNOTATIONS_TABLE.out.gene_annotations_url, DGE_DESEQ2_RRNA_RM.out.dge_table, dge_annotations_script, "_rRNArm" )
        annotated_dge_table_rrna_rm = ADD_GENE_ANNOTATIONS_RRNA_RM.out.annotated_dge_table

        // MultiQC
        ch_multiqc_config = params.multiqc_config ? Channel.fromPath( params.multiqc_config ) : Channel.fromPath("NO_FILE")
        RAW_READS_MULTIQC(samples_txt, raw_fastqc_zip, ch_multiqc_config, "raw_")
        TRIMMING_MULTIQC(samples_txt, trimgalore_reports, ch_multiqc_config, "trimming_")
        TRIMMED_READS_MULTIQC(samples_txt, trimmed_fastqc_zip, ch_multiqc_config, "trimmed_")
        ALIGN_MULTIQC(samples_txt, bowtie2_alignment_logs, ch_multiqc_config, "align_")
        INFER_EXPERIMENT_MULTIQC(samples_txt, INFER_EXPERIMENT.out.log | map { it[1] } | collect, ch_multiqc_config, "infer_exp_")
        GENEBODY_COVERAGE_MULTIQC(samples_txt, GENEBODY_COVERAGE.out.log | map { it[1] } | collect, ch_multiqc_config, "geneBody_cov_")
        INNER_DISTANCE_MULTIQC(samples_txt, INNER_DISTANCE.out.log | map { it[1] } | collect, ch_multiqc_config, "inner_dist_")
        READ_DISTRIBUTION_MULTIQC(samples_txt, READ_DISTRIBUTION.out.log | map { it[1] } | collect, ch_multiqc_config, "read_dist_")
        COUNT_MULTIQC(samples_txt, FEATURECOUNTS.out.summary, ch_multiqc_config, "featureCounts_")
        

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
        SOFTWARE_VERSIONS(ch_outdir, ch_final_software_versions)

        // Parse QC metrics
        all_multiqc_output = RAW_READS_MULTIQC.out.data
            | concat( TRIMMED_READS_MULTIQC.out.data )
            | concat( ALIGN_MULTIQC.out.data )
            | concat( GENEBODY_COVERAGE_MULTIQC.out.data )
            | concat( INFER_EXPERIMENT_MULTIQC.out.data )
            | concat( INNER_DISTANCE_MULTIQC.out.data )
            | concat( READ_DISTRIBUTION_MULTIQC.out.data )
            | concat( COUNT_MULTIQC.out.data )
            | collect

        // Generate QC metrics table
        PARSE_QC_METRICS(
            ch_outdir,
            osd_accession,
            ch_meta,
            isa_archive,
            all_multiqc_output,
            DGE_DESEQ2.out.norm_counts | map{ it -> it[1] }
        )

        // Clean paths in outputs before VVing & publishing
        // CLEAN_RAW_READS_MULTIQC(RAW_READS_MULTIQC.out.zipped_data, "raw")
        // CLEAN_TRIMMED_READS_MULTIQC(TRIMMED_READS_MULTIQC.out.zipped_data, "trimmed")
        // CLEAN_TRIMMING_MULTIQC(TRIMMING_MULTIQC.out.zipped_data, "trimming")
        // CLEAN_ALIGN_MULTIQC(ALIGN_MULTIQC.out.zipped_data, "align")
        // CLEAN_INFER_EXPERIMENT_MULTIQC(INFER_EXPERIMENT_MULTIQC.out.zipped_data, "infer_exp")
        // CLEAN_GENEBODY_COVERAGE_MULTIQC(GENEBODY_COVERAGE_MULTIQC.out.zipped_data, "geneBody_cov")
        // CLEAN_INNER_DISTANCE_MULTIQC(INNER_DISTANCE_MULTIQC.out.zipped_data, "inner_dist")
        // CLEAN_READ_DISTRIBUTION_MULTIQC(READ_DISTRIBUTION_MULTIQC.out.zipped_data, "read_dist")
        // CLEAN_COUNT_MULTIQC(COUNT_MULTIQC.out.zipped_data, "featureCounts")

        VV_RAW_READS(
            dp_tools_plugin,
            ch_outdir,
            ch_meta,
            runsheet_path,
            raw_reads | map{ it -> it[1] } | collect,
            raw_fastqc_zip,
            RAW_READS_MULTIQC.out.zipped_data,
            RAW_READS_MULTIQC.out.html
        )

        VV_TRIMMED_READS(
            dp_tools_plugin,
            ch_outdir,
            ch_meta,
            runsheet_path,
            trimmed_reads | map{ it -> it[1] } | collect,
            trimmed_fastqc_zip,
            TRIMMED_READS_MULTIQC.out.zipped_data,
            TRIMMED_READS_MULTIQC.out.html,
            TRIMGALORE.out.reports | collect,
            TRIMMING_MULTIQC.out.zipped_data,
            TRIMMING_MULTIQC.out.html
        )

        VV_BOWTIE2_ALIGNMENT(
            dp_tools_plugin,
            ch_outdir,
            ch_meta,
            runsheet_path,
            ALIGN_BOWTIE2.out.alignment_logs | collect,                       // log files
            ALIGN_BOWTIE2.out.unmapped_reads | collect,                       // unmapped reads
            ALIGN_MULTIQC.out.zipped_data,                                  // MultiQC report
            ALIGN_MULTIQC.out.html,
            SORT_AND_INDEX_BAM.out.sorted_bam | map{ it -> it[1] } | collect, // sorted BAM files
            SORT_AND_INDEX_BAM.out.sorted_bam | map{ it -> it[2] } | collect  // BAM index files
        )

        VV_RSEQC(
            dp_tools_plugin,
            ch_outdir,
            ch_meta,
            runsheet_path,
            ch_rseqc_logs,
            GENEBODY_COVERAGE_MULTIQC.out.zipped_data,
            GENEBODY_COVERAGE_MULTIQC.out.html,
            INFER_EXPERIMENT_MULTIQC.out.zipped_data,
            INFER_EXPERIMENT_MULTIQC.out.html,
            Channel.empty() | mix(INNER_DISTANCE_MULTIQC.out.zipped_data) | collect | ifEmpty({ file("PLACEHOLDER1") }),
            Channel.empty() | mix(INNER_DISTANCE_MULTIQC.out.html) | collect | ifEmpty({ file("PLACEHOLDER2") }),
            READ_DISTRIBUTION_MULTIQC.out.zipped_data,
            READ_DISTRIBUTION_MULTIQC.out.html
        )

        VV_FEATURECOUNTS(
            dp_tools_plugin,
            ch_outdir,
            ch_meta,
            runsheet_path,
            FEATURECOUNTS.out.counts,
            FEATURECOUNTS.out.summary,
            QUANTIFY_FEATURECOUNTS_GENES.out.num_non_zero_genes,
            COUNT_MULTIQC.out.zipped_data,
            COUNT_MULTIQC.out.html
        )

        VV_DGE_DESEQ2(
            dp_tools_plugin,
            ch_outdir,
            ch_meta,
            runsheet_path,
            DGE_DESEQ2.out.norm_counts,
            DGE_DESEQ2.out.vst_norm_counts,
            DGE_DESEQ2.out.sample_table,
            DGE_DESEQ2.out.contrasts,
            ADD_GENE_ANNOTATIONS.out.annotated_dge_table,
            DGE_DESEQ2_RRNA_RM.out.norm_counts,
            DGE_DESEQ2_RRNA_RM.out.vst_norm_counts,
            DGE_DESEQ2_RRNA_RM.out.sample_table,
            DGE_DESEQ2_RRNA_RM.out.contrasts,
            ADD_GENE_ANNOTATIONS_RRNA_RM.out.annotated_dge_table
        )


        GENERATE_PROTOCOL(ch_outdir,
            ch_meta,
            strandedness,
            SOFTWARE_VERSIONS.out.software_versions_yaml,
            reference_source,
            reference_version,
            genome_references_pre_ercc
        )

    emit:
        outdir = ch_outdir
}