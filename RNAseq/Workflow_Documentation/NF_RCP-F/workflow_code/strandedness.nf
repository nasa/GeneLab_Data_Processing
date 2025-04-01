// Workflow that determines the strandedness of reads compared to a reference genome bed file
include { SORT_INDEX_BAM } from './modules/rseqc.nf' addParams(PublishTo: "02-STAR_Alignment")
include { GENEBODY_COVERAGE } from './modules/rseqc.nf' addParams(PublishTo: "RSeQC_Analyses/02_geneBody_coverage")
include { INFER_EXPERIMENT } from './modules/rseqc.nf' addParams(PublishTo: "RSeQC_Analyses/03_infer_experiment")
include { INNER_DISTANCE } from './modules/rseqc.nf' addParams(PublishTo: "RSeQC_Analyses/04_inner_distance")
include { READ_DISTRIBUTION } from './modules/rseqc.nf' addParams(PublishTo: "RSeQC_Analyses/05_read_distribution")
include { ASSESS_STRANDEDNESS } from './modules/rseqc.nf'

include { MULTIQC as GENEBODY_COVERAGE_MULTIQC } from './modules/quality.nf' addParams(PublishTo: "RSeQC_Analyses/02_geneBody_coverage", MQCLabel:"geneBody_cov")
include { MULTIQC as INFER_EXPERIMENT_MULTIQC } from './modules/quality.nf' addParams(PublishTo: "RSeQC_Analyses/03_infer_experiment", MQCLabel:"infer_exp")
include { MULTIQC as INNER_DISTANCE_MULTIQC } from './modules/quality.nf' addParams(PublishTo: "RSeQC_Analyses/04_inner_distance", MQCLabel:"inner_dist")
include { MULTIQC as READ_DISTRIBUTION_MULTIQC } from './modules/quality.nf' addParams(PublishTo: "RSeQC_Analyses/05_read_distribution", MQCLabel:"read_dist")


workflow strandedness{
  take:
    bam_array // array: sample-wise tuples (meta, bam_file)
    genome_bed
    samples_ch
 
  main:
     
     bam_array | SORT_INDEX_BAM 
     SORT_INDEX_BAM.out.bam | combine( genome_bed ) |  set { ch_bam_bed }

     INFER_EXPERIMENT( ch_bam_bed )
     GENEBODY_COVERAGE( ch_bam_bed )
     INNER_DISTANCE( ch_bam_bed )
     READ_DISTRIBUTION( ch_bam_bed )
    
     // duplicated in each subworkflow, could use refactoring 
     ch_multiqc_config = params.multiqcConfig ? Channel.fromPath( params.multiqcConfig ) : Channel.fromPath("NO_FILE")
 
     INFER_EXPERIMENT_MULTIQC( samples_ch, INFER_EXPERIMENT.out.log | map { it[1] } | collect, ch_multiqc_config )
     GENEBODY_COVERAGE_MULTIQC( samples_ch, GENEBODY_COVERAGE.out.log | map { it[1] } | collect, ch_multiqc_config )
     INNER_DISTANCE_MULTIQC( samples_ch, INNER_DISTANCE.out.log | map { it[1] } | collect, ch_multiqc_config )
     READ_DISTRIBUTION_MULTIQC( samples_ch, READ_DISTRIBUTION.out.log | map { it[1] } | collect, ch_multiqc_config )
     
     ch_software_versions = Channel.empty()
     ch_software_versions | mix(INFER_EXPERIMENT.out.version,
                                GENEBODY_COVERAGE.out.version,
                                SORT_INDEX_BAM.out.version,
                                INNER_DISTANCE.out.version,
                                READ_DISTRIBUTION.out.version)
                          | set{ ch_software_versions }

     INFER_EXPERIMENT.out.log | map { it[1] }
                               | collect
                               | set { ch_infer_expt }
    
     ch_infer_expt | ASSESS_STRANDEDNESS
 
     ch_rseqc_mqc_reports = Channel.empty()
     ch_rseqc_mqc_reports | mix(INFER_EXPERIMENT_MULTIQC.out.zipped_report,
                                INNER_DISTANCE_MULTIQC.out.zipped_report,
                                READ_DISTRIBUTION_MULTIQC.out.zipped_report,
                                GENEBODY_COVERAGE_MULTIQC.out.zipped_report)
                          | collect
                          | set{ ch_rseqc_mqc_reports }
 
     ch_rseqc_logs = Channel.empty()
     ch_rseqc_logs | mix(INFER_EXPERIMENT.out.log_only,
                         GENEBODY_COVERAGE.out.all_output,
                         INNER_DISTANCE.out.all_output,
                         READ_DISTRIBUTION.out.log_only)
                   | collect
                   | set{ ch_rseqc_logs }
     

  emit:
     strandedness = ASSESS_STRANDEDNESS.out 
     infer_expt = ch_infer_expt
     versions = ch_software_versions 
     rseqc_logs = ch_rseqc_logs
     infer_expt_mqc = INFER_EXPERIMENT_MULTIQC.out.data
     mqc_reports = ch_rseqc_mqc_reports
     bam_bed = SORT_INDEX_BAM.out.bam_only_files
     genebody_coverage_multiqc = Channel.empty() | mix(GENEBODY_COVERAGE_MULTIQC.out.zipped_report, GENEBODY_COVERAGE_MULTIQC.out.unzipped_report) | collect
     infer_experiment_multiqc = Channel.empty() | mix(INFER_EXPERIMENT_MULTIQC.out.zipped_report, INFER_EXPERIMENT_MULTIQC.out.unzipped_report) | collect
     inner_distance_multiqc = Channel.empty() | mix(INNER_DISTANCE_MULTIQC.out.zipped_report, INNER_DISTANCE_MULTIQC.out.unzipped_report) | collect | ifEmpty({ file("NO_FILES.placeholder") }) // Ensures this channel is populated as a placeholder for single end studies
     read_distribution_multiqc = Channel.empty() | mix(READ_DISTRIBUTION_MULTIQC.out.zipped_report, READ_DISTRIBUTION_MULTIQC.out.unzipped_report) | collect
}
