nextflow.enable.dsl=2
// color defs
c_back_bright_red = "\u001b[41;1m";
c_bright_green = "\u001b[32;1m";
c_blue = "\033[0;34m";
c_reset = "\033[0m";

include { PARSE_ANNOTATION_TABLE } from './modules/PARSE_ANNOTATION_TABLE.nf'
include { VV_AGILE1CH } from './modules/VV_AGILE1CH.nf'
include { AGILE1CH } from './modules/AGILE1CH.nf'
include { RUNSHEET_FROM_GLDS } from './modules/RUNSHEET_FROM_GLDS.nf'

/**************************************************
* HELP MENU  **************************************
**************************************************/
if (params.help) {
  println("┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅")
  println("┇ Microarray Agilent 1 Channel Pipeline: $workflow.manifest.version  ┇")
  println("┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅")
  println("Usage example 1: Processing GLDS datasets using genome fasta and gtf from Ensembl")
  println("   > nextflow run ./main.nf --gldsAccession GLDS-194")
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
  println("  --skipVV              Skip automated V&V. Default: false")
  println("  --outputDir           Directory to save staged raw files and processed files. Default: <launch directory>")
  // println("  -stub-run             runs the workflow forcing 'unstranded' RSEM settings and using dummy gene counts in the differential gene expression (DGE) analysis. Useful when combined with the --truncateTo parameter this often leads to low gene counts and errors in the DGE analysis")  
  exit 0
  }

println "PARAMS: $params"
println "\n"

/**************************************************
* CHECK REQUIRED PARAMS AND LOAD  *****************
**************************************************/
// Get all params sourced data into channels
// Set up channel containing glds accession number
if ( !params.outputDir ) {  params.outputDir = "$workflow.launchDir" }

/**************************************************
* WORKFLOW SPECIFIC PRINTOUTS  ********************
**************************************************/

workflow {
	main:
    if ( !params.runsheetPath ) {
        RUNSHEET_FROM_GLDS( 
          params.osdAccession,
          params.gldsAccession,
        ) 
        RUNSHEET_FROM_GLDS.out.runsheet | set{ ch_runsheet }
    } else {
        ch_runsheet = channel.fromPath( params.runsheetPath )
    }


    ch_runsheet | splitCsv(header: true) | first | view |  set{ ch_meta }

    PARSE_ANNOTATION_TABLE(
        params.annotation_file_path,
        ch_meta | map { it.organism }
    )

    AGILE1CH(
      channel.fromPath( "${ projectDir }/bin/Agile1CMP.qmd" ),
      ch_runsheet,
      params.biomart_attribute ? params.biomart_attribute : false, // supply biomart_attribute if parameter exists
      PARSE_ANNOTATION_TABLE.out.annotations_db_url,
      ch_meta | map { it.organism },
      params.limit_biomart_query
    )

    VV_AGILE1CH( 
      ch_runsheet, 
      AGILE1CH.out.de,
      params.skipVV
      )
    /*

    // Software Version Capturing
    nf_version = "Nextflow Version:".concat("${nextflow.version}\n<><><>\n")
    ch_nextflow_version = Channel.value(nf_version)
    ch_software_versions = Channel.empty()
    RAW_FASTQC.out.version | mix(ch_software_versions) | set{ch_software_versions}
    RAW_MULTIQC.out.version | mix(ch_software_versions) | set{ch_software_versions}

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
    */

    emit:
      meta = ch_meta 
      runsheet = ch_runsheet
      // done = SOFTWARE_VERSIONS.out 
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
    GENERATE_MD5SUMS(ch_processed_directory, ch_runsheet )
    UPDATE_ISA_TABLES(ch_processed_directory, ch_runsheet )
}