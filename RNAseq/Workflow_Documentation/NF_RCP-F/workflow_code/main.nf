nextflow.enable.dsl=2
// color defs
c_back_bright_red = "\u001b[41;1m";
c_bright_green = "\u001b[32;1m";
c_blue = "\033[0;34m";
c_reset = "\033[0m";

include { FASTQC as RAW_FASTQC } from './modules/quality.nf'
include { FASTQC as TRIMMED_FASTQC } from './modules/quality.nf'
include { MULTIQC as RAW_MULTIQC } from './modules/quality.nf' addParams(MQCLabel:"raw")
include { MULTIQC as TRIMMED_MULTIQC } from './modules/quality.nf' addParams(MQCLabel:"trimmed")
include { MULTIQC as TRIM_MULTIQC } from './modules/quality.nf' addParams(MQCLabel:"trimming")
include { MULTIQC as ALIGN_MULTIQC } from './modules/quality.nf' addParams(MQCLabel:"align")
include { MULTIQC as COUNT_MULTIQC } from './modules/quality.nf' addParams(MQCLabel:"RSEM_count")
include { MULTIQC as ALL_MULTIQC } from './modules/quality.nf' addParams(MQCLabel:"all")
include { TRIMGALORE } from './modules/quality.nf'
include { BUILD_STAR;
          ALIGN_STAR;
          BUILD_RSEM;
          COUNT_ALIGNED;
          SUBSAMPLE_GENOME;
          CONCAT_ERCC;
          QUANTIFY_STAR_GENES;
          QUANTIFY_RSEM_GENES } from './modules/genome.nf'
include { DGE_BY_DESEQ2 } from './modules/DGE_BY_DESEQ2'
include { VV_RAW_READS;
          VV_TRIMMED_READS;
          VV_STAR_ALIGNMENTS;
          VV_RSEQC;
          VV_RSEM_COUNTS;
          VV_DESEQ2_ANALYSIS;
          VV_CONCAT_FILTER } from './modules/vv.nf'
include { GET_MAX_READ_LENGTH } from './modules/fastqc.nf'
include { UPDATE_ISA_TABLES;
          GENERATE_MD5SUMS;
          SOFTWARE_VERSIONS } from './modules/genelab.nf'

/**************************************************
* HELP MENU  **************************************
**************************************************/
if (params.help) {
  println("┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅")
  println("┇ RNASeq Consensus Pipeline: $workflow.manifest.version  ┇")
  println("┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅")
  println("Usage example 1: Processing GLDS datasets using genome fasta and gtf from Ensembl")
  println("   > nextflow run ./main.nf --gldsAccession GLDS-194")
  println()
  println("Usage example 2: Processing GLDS datasets using local genome fasta and gtf")
  println("   Note: ensemblVersion and ref_source are used here to label subdirectories for derived reference files.")
  println("   > nextflow run ./main.nf --gldsAccession GLDS-194 --ensemblVersion 96 --ref_source <reference_label>  --ref_fasta </path/to/fasta> --ref_gtf </path/to/gtf>")
  println()
  println("Usage example 3: Processing Other datasets")
  println("   Note: This requires a user-created runsheet.")
  println("   > nextflow run ./main.nf --runsheetPath </path/to/runsheet>")
  println()
  println("arguments:")
  println("  --help                show this help message and exit")
  println("  --gldsAccession GLDS-000")
  println("                        the GLDS accession id to process through the RNASeq Concensus Pipeline.")
  println("  --runsheetPath        Use a local runsheet instead one automatically generated from a GLDS ISA archive.")
  println("  --ensemblVersion n    Specifies the ensembl Version to use for the reference genome. The default version is ")
  println("  --skipVV              Skip automated V&V. Default: false")
  println("  --outputDir           Directory to save staged raw files and processed files. Default: <launch directory>")
  println("  --limitSamplesTo n    limit the number of samples staged to a number.")
  println("  --genomeSubsample n   subsamples genome fasta and gtf files to the supplied chromosome.")
  println("  --truncateTo n        limit number of reads downloaded and processed to *n* reads , for paired end limits number of reverse and forward read files to *n* reads each.")
  println("  --force_single_end    forces analysis to use single end processing.  For paired end datasets, this means only R1 is used.  For single end studies, this should have no effect.")
  println("  --stageLocal          download the raw reads files for the supplied GLDS accession id.  Set to false to retrieve metadata and generate a runsheet for GLDS datasets to disable raw read download and processing.  Default: true")
  println("  --ref_fasta           specifies a reference fasta from a local path. This an is an alternative approach from the automatic retrieval of reference files from ensembl")  
  println("  --ref_gtf             specifies a reference gtf from a local path. This an is an alternative approach from the automatic retrieval of reference files from ensembl")  
  println("  --referenceStorePath  specifies the directory where fetched reference files are downloaded to")  
  println("  --derivedStorePath    specifies the directory where derivative reference files are saved. Examples of such files in this pipeline included BED and PRED files generated from the reference gtf")  
  println("  --ref_source          a string to label subdirectories in 'StorePath' paths. Examples include 'ensembl' or 'ensembl_plants'.")  
  println("  -stub-run             runs the workflow forcing 'unstranded' RSEM settings and using dummy gene counts in the differential gene expression (DGE) analysis. Useful when combined with the --truncateTo parameter this often leads to low gene counts and errors in the DGE analysis")  
  exit 0
  }

println "PARAMS: $params"
println "\n"
println "Storing any newly fetched primary references files here: ${params.referenceStorePath}"
println "Storing any newly generated derived reference files here: ${params.derivedStorePath}"

/**************************************************
* CHECK REQUIRED PARAMS AND LOAD  *****************
**************************************************/
// Get all params sourced data into channels
// Set up channel containing glds accession number
if ( params.gldsAccession ) {ch_glds_accession = Channel.from( params.gldsAccession )} else { exit 1, "Missing Required Parameter: gldsAccession. Example for setting on CLI: --gldsAccession GLDS-194"}

// Check conditionally required parameter (if using direct fasta, an ensemblVersion must also be supplied)
if ( params.ref_fasta ) {
  if ( !params.ensemblVersion ) { exit 1, "Missing Required Parameter: ensemblVersion. Example for setting on CLI: --ensemblVersion 96" }
}

if ( !params.outputDir ) {  params.outputDir = "$workflow.launchDir" }

ch_multiqc_config = params.multiqcConfig ? Channel.fromPath( params.multiqcConfig ) : Channel.fromPath("NO_FILE")

/**************************************************
* DEBUG WARNING  **********************************
**************************************************/
if ( params.limitSamplesTo || params.truncateTo || params.force_single_end || params.genomeSubsample) {
  println("${c_back_bright_red}WARNING WARNING: DEBUG OPTIONS ENABLED!")
  params.limitSamplesTo ? println("Samples limited to ${params.limitSamplesTo}") : println("No Sample Limit Set")
  params.truncateTo ? println("Truncating reads to first ${params.truncateTo} records") : println("No Truncation By Record Limit Set")
  params.genomeSubsample ? println("Subsampling reference genome to chromosome '${params.genomeSubsample}'") : println("No subsampling of reference genome")
  params.force_single_end ? println("Forcing analysis to used only forward reads if paired end (i.e. as though single ended") : println("No forcing single end analysis")
  println("WARNING WARNING: DEBUG OPTIONS ENABLED!${c_reset}")
} else {
  params.limitSamplesTo ? println("Samples limited to ${params.limitSamplesTo}") : println("No Sample Limit Set")
  params.truncateTo ? println("Truncating reads to first ${params.truncateTo} records") : println("No Truncation By Record Limit Set")
  params.genomeSubsample ? println("Subsampling reference genome to chromosome '${params.genomeSubsample}'") : println("No subsampling of reference genome")
  params.force_single_end ? println("Forcing analysis to used only forward reads if paired end (i.e. as though single ended") : println("No forcing single end analysis")
}

/**************************************************
* WORKFLOW SPECIFIC PRINTOUTS  ********************
**************************************************/
if ( params.stageLocal && params.truncateTo ) {
  // download truncated raw reads
  println("${c_bright_green}Staging truncated raw reads for ${params.gldsAccession}${c_reset}")
} else if ( params.stageLocal && !params.truncateTo ) {
  // download full raw reads
  println("${c_bright_green}Staging raw reads for ${params.gldsAccession}${c_reset}")
} else {
  // maybe print some nice data from the samplesheet
  println("${c_bright_green}No Staging of raw reads.  Only getting Metadata for ${params.gldsAccession}${c_reset}")
}

include { staging as STAGING } from './stage_analysis.nf'
include { references as REFERENCES } from './references.nf'
include { strandedness as STRANDEDNESS } from './strandedness.nf'

workflow {
	main:
    STAGING( ch_glds_accession, params.stageLocal )
    // This process can use a single meta and a collection of read paths
    STAGING.out.raw_reads | first 
                          | map{it -> it[0]} 
                          | view { meta -> "${c_bright_green}Autodetected Processing Metadata:\n\t hasERCC: ${meta.has_ercc}\n\t pairedEND: ${meta.paired_end}\n\t organism: ${meta.organism_sci}${c_reset}"  }
                          | set { ch_meta }

    STAGING.out.raw_reads | map{ it -> it[1] } | collect | set { ch_all_raw_reads }
    STAGING.out.raw_reads | map { it[0].id }
                          | collectFile(name: "samples.txt", sort: true, newLine: true)
                          | set { ch_samples_txt }

    STAGING.out.raw_reads | RAW_FASTQC

    RAW_FASTQC.out.fastqc | map { it -> [ it[1], it[2] ] }
                          | flatten
                          | unique
                          | collect
                          | set { raw_mqc_ch }
    RAW_FASTQC.out.fastqc | map { it -> [ it[2] ] }
                          | flatten
                          | GET_MAX_READ_LENGTH

    GET_MAX_READ_LENGTH.out.length  | max { it.toInteger() }
                                    | set { max_read_length_ch }

    STAGING.out.raw_reads |  TRIMGALORE

    TRIMGALORE.out.reads | TRIMMED_FASTQC

    TRIMMED_FASTQC.out.fastqc | map { it -> [ it[1], it[2] ] } \
                              | flatten \
                              | unique \
                              | collect \
                              | set { trim_mqc_ch }

    REFERENCES( ch_meta | map { it.organism_sci }, ch_meta | map { it.has_ercc } )
    REFERENCES.out.genome_annotations | set { genome_annotations }      

    BUILD_STAR( 
      genome_annotations, 
      ch_meta, 
      max_read_length_ch,
      REFERENCES.out.reference_version_and_source
    )

    TRIMGALORE.out.reads | combine( BUILD_STAR.out.build ) | ALIGN_STAR


    STRANDEDNESS ( ALIGN_STAR.out.bam_by_coord, REFERENCES.out.genome_bed, ch_samples_txt ) 
    STRANDEDNESS.out.strandedness | map { it.text.split(":")[0] } | set { strandedness_ch }

    BUILD_RSEM( 
      genome_annotations, 
      ch_meta,
      REFERENCES.out.reference_version_and_source
      )

    ALIGN_STAR.out.bam_to_transcriptome | combine( BUILD_RSEM.out.build ) | set { aligned_ch }
    QUANTIFY_STAR_GENES( 
        ch_samples_txt, 
        ALIGN_STAR.out.read_per_gene | toSortedList,
        strandedness_ch
      )
      
    COUNT_ALIGNED( aligned_ch, strandedness_ch )

    ALIGN_STAR.out.alignment_logs       | collect 
                                        | set { align_mqc_ch }


    COUNT_ALIGNED.out.counts | map { it[1] } | collect | set { rsem_ch }

    QUANTIFY_RSEM_GENES( ch_samples_txt, rsem_ch )


    // Note: This is loaded as a zip file to ensure correct caching (directories don't seem to yield identical hases)
    ch_r_scripts = channel.fromPath( "${ projectDir }/bin/dge_annotation_R_scripts.zip" ) 

    DGE_BY_DESEQ2(STAGING.out.runsheet, 
                  COUNT_ALIGNED.out.gene_counts | toSortedList, 
                  ch_meta, 
                  REFERENCES.out.gene_annotations, 
                  ch_r_scripts
                  )

    // ALL MULTIQC
    RAW_MULTIQC( ch_samples_txt, raw_mqc_ch, ch_multiqc_config  )
    TRIMMED_MULTIQC( ch_samples_txt, trim_mqc_ch, ch_multiqc_config ) // refering to the trimmed reads
    TRIM_MULTIQC( ch_samples_txt, TRIMGALORE.out.reports | collect, ch_multiqc_config ) // refering to the trimming process
    ALIGN_MULTIQC( ch_samples_txt, align_mqc_ch, ch_multiqc_config )
    COUNT_MULTIQC( ch_samples_txt, rsem_ch, ch_multiqc_config )
    raw_mqc_ch | concat( trim_mqc_ch ) 
                | concat( ALIGN_STAR.out.alignment_logs ) 
                | concat( STRANDEDNESS.out.rseqc_logs )
                | concat( rsem_ch )
                | concat( TRIMGALORE.out.reports )
                | collect | set { all_mqc_ch }
    ALL_MULTIQC( ch_samples_txt, all_mqc_ch, ch_multiqc_config )

    VV_RAW_READS( ch_meta,
                  STAGING.out.runsheet,
                  ch_all_raw_reads,
                  RAW_FASTQC.out.fastqc | map { it -> [ it[1], it[2] ] } | flatten | collect,
                  RAW_MULTIQC.out.zipped_report,
                  RAW_MULTIQC.out.unzipped_report,
                  "${ projectDir }/bin/dp_tools__NF_RCP" // dp_tools plugin
                )
    VV_TRIMMED_READS( ch_meta,
                      STAGING.out.runsheet,
                      TRIMGALORE.out.reads | map { it -> it[1] } | flatten | collect,
                      TRIMMED_FASTQC.out.fastqc | map { it -> [ it[1], it[2] ] } | flatten | collect,
                      TRIMMED_MULTIQC.out.zipped_report,
                      TRIMMED_MULTIQC.out.unzipped_report,
                      TRIMGALORE.out.reports | collect,
                      TRIM_MULTIQC.out.zipped_report,
                      TRIM_MULTIQC.out.unzipped_report,
                      "${ projectDir }/bin/dp_tools__NF_RCP" // dp_tools plugin
                    )
    VV_STAR_ALIGNMENTS( STAGING.out.runsheet,
                        ALIGN_STAR.out.publishables | collect,
                        QUANTIFY_STAR_GENES.out.publishables | collect,
                        ALIGN_MULTIQC.out.zipped_report,
                        ALIGN_MULTIQC.out.unzipped_report,
                        STRANDEDNESS.out.bam_bed | collect,
                        "${ projectDir }/bin/dp_tools__NF_RCP" // dp_tools plugin
                      )
    VV_RSEQC( ch_meta,
              STAGING.out.runsheet,
              STRANDEDNESS.out.rseqc_logs,
              STRANDEDNESS.out.genebody_coverage_multiqc,
              STRANDEDNESS.out.infer_experiment_multiqc,
              STRANDEDNESS.out.inner_distance_multiqc,
              STRANDEDNESS.out.read_distribution_multiqc,
              "${ projectDir }/bin/dp_tools__NF_RCP" // dp_tools plugin
            )
    VV_RSEM_COUNTS( STAGING.out.runsheet,
                    COUNT_ALIGNED.out.only_counts | collect,
                    QUANTIFY_RSEM_GENES.out.publishables,
                    COUNT_MULTIQC.out.zipped_report,
                    COUNT_MULTIQC.out.unzipped_report,
                    "${ projectDir }/bin/dp_tools__NF_RCP" // dp_tools plugin
                  )
    VV_DESEQ2_ANALYSIS( ch_meta,
                        STAGING.out.runsheet,
                        QUANTIFY_RSEM_GENES.out.publishables,
                        COUNT_MULTIQC.out.zipped_report,
                        COUNT_MULTIQC.out.unzipped_report,
                        DGE_BY_DESEQ2.out.norm_counts,
                        DGE_BY_DESEQ2.out.dge,
                        DGE_BY_DESEQ2.out.norm_counts_ercc | ifEmpty( { file("NO_FILES.placeholder") }),
                        DGE_BY_DESEQ2.out.dge_ercc | ifEmpty( { file("NO_FILES.placeholder") }),
                        "${ projectDir }/bin/dp_tools__NF_RCP" // dp_tools plugin
                  )

    // Software Version Capturing
    nf_version = "Nextflow Version:".concat("${nextflow.version}\n<><><>\n")
    ch_nextflow_version = Channel.value(nf_version)
    ch_software_versions = Channel.empty()
    RAW_FASTQC.out.version | mix(ch_software_versions) | set{ch_software_versions}
    RAW_MULTIQC.out.version | mix(ch_software_versions) | set{ch_software_versions}
    TRIMGALORE.out.version | mix(ch_software_versions) | set{ch_software_versions}
    TRIMMED_FASTQC.out.version | mix(ch_software_versions) | set{ch_software_versions}
    TRIMMED_MULTIQC.out.version | mix(ch_software_versions) | set{ch_software_versions}
    ALIGN_STAR.out.version | mix(ch_software_versions) | set{ch_software_versions}
    COUNT_ALIGNED.out.version | mix(ch_software_versions) | set{ch_software_versions}
    DGE_BY_DESEQ2.out.version | mix(ch_software_versions) | set{ch_software_versions}
    STRANDEDNESS.out.versions | mix(ch_software_versions) | set{ch_software_versions}
    ch_software_versions | map { it.text + "\n<><><>\n"}
                          | unique
                          | mix(ch_nextflow_version)
                          | collectFile(name: "software_versions.txt", newLine: true, cache: false)
        | set{ch_final_software_versions}

    // VV processes
      
    VV_CONCAT_FILTER( VV_RAW_READS.out.log | mix( VV_TRIMMED_READS.out.log,
                                                  VV_STAR_ALIGNMENTS.out.log,
                                                  VV_RSEQC.out.log,
                                                  VV_RSEM_COUNTS.out.log,
                                                  VV_DESEQ2_ANALYSIS.out.log,
                                                  ) | collect )

    // Generate final versions output
    SOFTWARE_VERSIONS(ch_final_software_versions)

    emit:
      meta = ch_meta 
      runsheet = STAGING.out.runsheet 
      done = SOFTWARE_VERSIONS.out 
}

workflow.onComplete {
    println "${c_bright_green}Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    if ( workflow.success ) {
      println "Raw and Processed data location: ${ params.outputDir }/${ params.gldsAccession }"
      println "V&V logs location: ${ params.outputDir }/${ params.gldsAccession }/VV_Logs"
      println "Pipeline tracing/visualization files location:  ${ params.outputDir }/${ params.gldsAccession }/Resource_Usage${c_reset}"
    }
}

workflow STAGING_ONLY {
  main:
    STAGING( ch_glds_accession, false )
}

workflow POST_PROCESSING {
  main:
    ch_processed_directory = Channel.fromPath("${ params.outputDir }/${ params.gldsAccession }", checkIfExists: true)
    ch_runsheet = Channel.fromPath("${ params.outputDir }/${ params.gldsAccession }/Metadata/*_runsheet.csv", checkIfExists: true)
    GENERATE_MD5SUMS(ch_processed_directory, ch_runsheet, "${ projectDir }/bin/dp_tools__NF_RCP" )
    UPDATE_ISA_TABLES(ch_processed_directory, ch_runsheet, "${ projectDir }/bin/dp_tools__NF_RCP" )
}