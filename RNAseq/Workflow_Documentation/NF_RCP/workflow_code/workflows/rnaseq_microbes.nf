include { PARSE_RUNSHEET } from './parse_runsheet.nf'
include { PARSE_ANNOTATIONS_TABLE } from '../modules/parse_annotations_table.nf'
include { FETCH_ISA } from '../modules/fetch_isa.nf'
include { ISA_TO_RUNSHEET } from '../modules/isa_to_runsheet.nf'
include { GET_ACCESSIONS } from '../modules/get_accessions.nf'
include { DOWNLOAD_REFERENCES } from '../modules/download_references.nf'
include { SUBSAMPLE_GENOME } from '../modules/subsample_genome.nf'
include { DOWNLOAD_ERCC } from '../modules/download_ercc.nf'
include { CONCAT_ERCC } from '../modules/concat_ercc.nf'
include { GTF_TO_PRED } from '../modules/gtf_to_pred.nf'
include { PRED_TO_BED } from '../modules/pred_to_bed.nf'
include { STAGE_RAW_READS } from './stage_raw_reads.nf'
include { STAGE_ANALYSIS } from './stage_analysis.nf'
include { FASTQC as RAW_FASTQC } from '../modules/fastqc.nf'
include { GET_MAX_READ_LENGTH } from '../modules/get_max_read_length.nf'
include { TRIMGALORE } from '../modules/trimgalore.nf'
include { FASTQC as TRIMMED_FASTQC } from '../modules/fastqc.nf'
include { BUILD_BOWTIE2_INDEX } from '../modules/build_bowtie2_index.nf'
include { ALIGN_BOWTIE2 } from '../modules/align_bowtie2.nf'
include { SORT_AND_INDEX_BAM } from '../modules/sort_and_index_bam.nf'
include { INFER_EXPERIMENT } from '../modules/rseqc.nf'
include { GENEBODY_COVERAGE } from '../modules/rseqc.nf'
include { INNER_DISTANCE } from '../modules/rseqc.nf'
include { READ_DISTRIBUTION } from '../modules/rseqc.nf'
include { ASSESS_STRANDEDNESS } from '../modules/assess_strandedness.nf'
include { GET_GTF_FEATURES } from '../modules/get_gtf_features.nf'
include { FEATURECOUNTS } from '../modules/featurecounts.nf'
include { QUANTIFY_FEATURECOUNTS_GENES } from '../modules/quantify_featurecounts_genes.nf'
include { EXTRACT_RRNA } from '../modules/extract_rrna.nf'
include { REMOVE_RRNA_FEATURECOUNTS } from '../modules/remove_rrna_featurecounts.nf'
include { DGE_DESEQ2 } from '../modules/dge_deseq2.nf'
include { DGE_DESEQ2 as DGE_DESEQ2_RRNA_RM } from '../modules/dge_deseq2.nf'
include { ADD_GENE_ANNOTATIONS } from '../modules/add_gene_annotations.nf'
include { ADD_GENE_ANNOTATIONS as ADD_GENE_ANNOTATIONS_RRNA_RM } from '../modules/add_gene_annotations.nf'
include { 
    MULTIQC as RAW_READS_MULTIQC;
    MULTIQC as TRIMMED_READS_MULTIQC;
    MULTIQC as ALIGN_MULTIQC;
    MULTIQC as GENEBODY_COVERAGE_MULTIQC;
    MULTIQC as INFER_EXPERIMENT_MULTIQC;
    MULTIQC as INNER_DISTANCE_MULTIQC;
    MULTIQC as READ_DISTRIBUTION_MULTIQC;
    MULTIQC as COUNT_MULTIQC
    MULTIQC as ALL_MULTIQC
} from '../modules/multiqc.nf'
include { PARSE_QC_METRICS } from '../modules/parse_qc_metrics.nf'
include { VV_RAW_READS;
    VV_TRIMMED_READS;
    VV_BOWTIE2_ALIGNMENT;
    VV_RSEQC;
    VV_FEATURECOUNTS;
    VV_DGE_DESEQ2;
    VV_CONCAT_FILTER } from '../modules/vv.nf'
include { SOFTWARE_VERSIONS } from '../modules/software_versions.nf'
include { GENERATE_PROTOCOL } from '../modules/generate_protocol.nf'
include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

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
        // Stage analysis setup (directory structure, inputs, and raw reads)
        STAGE_ANALYSIS(
            ch_outdir,
            dp_tools_plugin,
            accession,
            isa_archive_path,
            runsheet_path,
            api_url
        )
        ch_outdir = STAGE_ANALYSIS.out.ch_outdir
        samples = STAGE_ANALYSIS.out.samples
        raw_reads = STAGE_ANALYSIS.out.raw_reads
        samples_txt = STAGE_ANALYSIS.out.samples_txt
        runsheet_path = STAGE_ANALYSIS.out.runsheet_path
        isa_archive = STAGE_ANALYSIS.out.isa_archive
        osd_accession = STAGE_ANALYSIS.out.osd_accession
        glds_accession = STAGE_ANALYSIS.out.glds_accession

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

        // Run FastQC on raw reads  
        RAW_FASTQC( raw_reads ) 
        RAW_FASTQC.out.fastqc | map { it -> [ it[1], it[2] ] } // Collect the raw read fastqc zip files
        | flatten
        | collect // Collect all zip files into a single list
        | set { raw_fastqc_zip } // Create a channel with all zip files
        
        // Get the max read length by parsing the raw read fastqc zip files
        GET_MAX_READ_LENGTH( raw_fastqc_zip )
        max_read_length = GET_MAX_READ_LENGTH.out.length | map { it.toString().toInteger() }

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

        // Align reads using Bowtie2 (output as BAM)
        ALIGN_BOWTIE2( trimmed_reads, bowtie2_index_dir )
        bowtie2_alignment_logs = ALIGN_BOWTIE2.out.alignment_logs | collect

        // Sort and index BAM files 
        SORT_AND_INDEX_BAM( ALIGN_BOWTIE2.out.bam )
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

        // Create FeatureCounts table from BAMs
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
        // For rRNArm counts: Normalize counts, DGE, Add annotations to DGE table
        DGE_DESEQ2_RRNA_RM( ch_meta, runsheet_path, REMOVE_RRNA_FEATURECOUNTS.out.counts_rrnarm, dge_script, "_rRNArm" )
        ADD_GENE_ANNOTATIONS_RRNA_RM( ch_meta, PARSE_ANNOTATIONS_TABLE.out.gene_annotations_url, DGE_DESEQ2_RRNA_RM.out.dge_table, dge_annotations_script, "_rRNArm" )
        
        // MultiQC
        ch_multiqc_config = params.multiqc_config ? Channel.fromPath( params.multiqc_config ) : Channel.fromPath("NO_FILE")
        RAW_READS_MULTIQC(samples_txt, raw_fastqc_zip, ch_multiqc_config, "raw_")
        TRIMMED_READS_MULTIQC(samples_txt, trimmed_fastqc_zip | concat( TRIMGALORE.out.reports ) | collect, ch_multiqc_config, "trimmed_")
        ALIGN_MULTIQC(samples_txt, bowtie2_alignment_logs, ch_multiqc_config, "align_")
        INFER_EXPERIMENT_MULTIQC(samples_txt, INFER_EXPERIMENT.out.log | map { it[1] } | collect, ch_multiqc_config, "infer_exp_")
        GENEBODY_COVERAGE_MULTIQC(samples_txt, GENEBODY_COVERAGE.out.log | map { it[1] } | collect, ch_multiqc_config, "geneBody_cov_")
        INNER_DISTANCE_MULTIQC(samples_txt, INNER_DISTANCE.out.log | map { it[1] } | collect, ch_multiqc_config, "inner_dist_")
        READ_DISTRIBUTION_MULTIQC(samples_txt, READ_DISTRIBUTION.out.log | map { it[1] } | collect, ch_multiqc_config, "read_dist_")
        COUNT_MULTIQC(samples_txt, FEATURECOUNTS.out.summary, ch_multiqc_config, "FeatureCounts_")
        
        all_multiqc_input = raw_fastqc_zip
            | concat( TRIMGALORE.out.reports )
            | concat( trimmed_fastqc_zip )
            | concat( bowtie2_alignment_logs )
            | concat( GENEBODY_COVERAGE.out.log | map { it[1] } | collect )
            | concat( INFER_EXPERIMENT.out.log | map { it[1] } | collect )
            | concat( INNER_DISTANCE.out.log | map { it[1] } | collect )
            | concat( READ_DISTRIBUTION.out.log | map { it[1] } | collect )
            | concat( FEATURECOUNTS.out.summary )
            | collect

        ALL_MULTIQC(samples_txt, all_multiqc_input, ch_multiqc_config, "all_")

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
        PARSE_QC_METRICS(
            ch_outdir,
            osd_accession,
            ch_meta,
            isa_archive.ifEmpty(file("ISA.zip")),  // Use a placeholder if isa_archive is empty
            all_multiqc_output,
            DGE_DESEQ2.out.norm_counts | map{ it -> it[1] }
        )

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
            TRIMGALORE.out.reports | collect
        )
        VV_BOWTIE2_ALIGNMENT(
            dp_tools_plugin,
            ch_outdir,
            ch_meta,
            runsheet_path,
            ALIGN_BOWTIE2.out.alignment_logs | collect,
            ALIGN_BOWTIE2.out.unmapped_reads | collect,
            SORT_AND_INDEX_BAM.out.sorted_bam | map{ it -> it[1] } | collect,
            SORT_AND_INDEX_BAM.out.sorted_bam | map{ it -> it[2] } | collect,
            ALIGN_MULTIQC.out.zipped_data,
            ALIGN_MULTIQC.out.html
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
            REMOVE_RRNA_FEATURECOUNTS.out.counts_rrnarm,
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
            ADD_GENE_ANNOTATIONS_RRNA_RM.out.annotated_dge_table
        )
        VV_CONCAT_FILTER( ch_outdir, VV_RAW_READS.out.log | mix( VV_TRIMMED_READS.out.log, // Concatenate and filter V&V logs
                                                    VV_BOWTIE2_ALIGNMENT.out.log,
                                                    VV_RSEQC.out.log,
                                                    VV_FEATURECOUNTS.out.log,
                                                    VV_DGE_DESEQ2.out.log,
                                                    ) | collect )

        // Software Version Capturing
        nf_version = '"NEXTFLOW":\n    nextflow: '.concat("${nextflow.version}\n")
        ch_nextflow_version = Channel.value(nf_version)
        ch_software_versions = Channel.empty()
        // Mix in versions from each process
        ch_software_versions = ch_software_versions
            | mix(GTF_TO_PRED.out.versions)
            | mix(PRED_TO_BED.out.versions)
            | mix(RAW_FASTQC.out.versions)
            | mix(TRIMGALORE.out.versions)
            | mix(ALIGN_BOWTIE2.out.versions)
            | mix(INFER_EXPERIMENT.out.versions)
            | mix(GENEBODY_COVERAGE.out.versions)
            | mix(INNER_DISTANCE.out.versions)
            | mix(READ_DISTRIBUTION.out.versions)
            | mix(FEATURECOUNTS.out.versions)
            | mix(RAW_READS_MULTIQC.out.versions)
            | mix(DGE_DESEQ2.out.versions)
            | mix(VV_RAW_READS.out.versions)
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

        GENERATE_PROTOCOL(ch_outdir,
            ch_meta,
            strandedness,
            SOFTWARE_VERSIONS.out.software_versions_yaml,
            reference_source,
            reference_version,
            genome_references_pre_ercc,
            runsheet_path
        )

    emit:
        SOFTWARE_VERSIONS.out.software_versions
}