nextflow.enable.dsl=2
// color defs
c_back_bright_red = "\u001b[41;1m";
c_bright_green = "\u001b[32;1m";
c_blue = "\033[0;34m";
c_reset = "\033[0m";

include { PARSE_ANNOTATION_TABLE } from './modules/PARSE_ANNOTATION_TABLE.nf'
include { VV_AFFYMETRIX } from './modules/VV_AFFYMETRIX.nf'
include { PROCESS_AFFYMETRIX } from './modules/PROCESS_AFFYMETRIX.nf'
include { RUNSHEET_FROM_GLDS } from './modules/RUNSHEET_FROM_GLDS.nf'
include { RUNSHEET_FROM_ISA } from './modules/RUNSHEET_FROM_ISA.nf'
include { GENERATE_SOFTWARE_TABLE } from './modules/GENERATE_SOFTWARE_TABLE'
include { DUMP_META } from './modules/DUMP_META'

/**************************************************
* HELP MENU  **************************************
**************************************************/
if (params.help) {
  println("┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅")
  println("┇ Affymetrix Microarray Pipeline: $workflow.manifest.version  ┇")
  println("┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅")
  println("Usage example 1: Processing GLDS datasets using genome fasta and gtf from Ensembl")
  println("   > nextflow run ./main.nf --osdAccession OSD-266 --gldsAccession GLDS-266")
  println()
  println("Usage example 2: Processing Other datasets")
  println("   Note: This requires a user-created runsheet.")
  println("   > nextflow run ./main.nf --runsheetPath </path/to/runsheet>")
  println()
  println("arguments:")
  println("  --help                show this help message and exit")
  println("  --osdAccession OSD-000")
  println("                        the OSD accession id to process through the Affymetrix Microarray Pipeline.")
  println("  --gldsAccession GLDS-000")
  println("                        the GLDS accession id to process through the Affymetrix Microarray Pipeline.")
  println("  --runsheetPath        Use a local runsheet instead one automatically generated from a GLDS ISA archive.")
  println("  --skipVV              Skip automated V&V. Default: false")
  println("  --resultsDir           Directory to save staged raw files and processed files. Default: <launch directory>")
  exit 0
  }

println "PARAMS: $params"
println "\n"

/**************************************************
* CHECK REQUIRED PARAMS AND LOAD  *****************
**************************************************/
println("Resolved output directory: ${ params.resultsDir }")

/**************************************************
* WORKFLOW SPECIFIC PRINTOUTS  ********************
**************************************************/

workflow {
	main:
    if ( !params.runsheetPath && !params.isaArchivePath) {
        RUNSHEET_FROM_GLDS( 
          params.osdAccession,
          params.gldsAccession,
          "${ projectDir }/bin/dp_tools__affymetrix" // dp_tools plugin
        ) 
        RUNSHEET_FROM_GLDS.out.runsheet | set{ ch_runsheet }
    } else if ( !params.runsheetPath && params.isaArchivePath ) {
        RUNSHEET_FROM_ISA( 
          params.osdAccession,
          params.gldsAccession,
          params.isaArchivePath,
          "${ projectDir }/bin/dp_tools__affymetrix" // dp_tools plugin
        )
        RUNSHEET_FROM_ISA.out.runsheet | set{ ch_runsheet }
    } else if ( params.runsheetPath && !params.isaArchivePath ) {
        ch_runsheet = channel.fromPath( params.runsheetPath )
    } else if ( params.runsheetPath && params.isaArchivePath ) {
        System.err.println("Error: User supplied both runsheetPath and isaArchivePath.  Only one or neither is allowed to be supplied!") // Print error message to System.err
        System.exit(1) // Exit with error code 1
    }


    ch_runsheet | splitCsv(header: true) | first | view |  set{ ch_meta }

    PARSE_ANNOTATION_TABLE(
        params.annotation_file_path,
        ch_meta | map { it.organism }
    )

    PROCESS_AFFYMETRIX(
      channel.fromPath( "${ projectDir }/bin/Affymetrix.qmd" ),
      ch_runsheet,
      PARSE_ANNOTATION_TABLE.out.annotations_db_url,
      ch_meta | map { it.organism },
      params.limit_biomart_query
    )

    VV_AFFYMETRIX( 
      ch_runsheet, 
      PROCESS_AFFYMETRIX.out.de,
      params.skipVV,
      "${ projectDir }/bin/dp_tools__affymetrix" // dp_tools plugin
      )

    // Software Version Capturing
    nf_version = "- name: nextflow\n  ".concat(
"""
  version: ${nextflow.version}
  homepage: https://www.nextflow.io
  workflow task: N/A
""")
    ch_software_versions = Channel.value(nf_version)
    PROCESS_AFFYMETRIX.out.versions | map{ it -> it.text } | mix(ch_software_versions) | set{ch_software_versions}
    VV_AFFYMETRIX.out.versions | map{ it -> it.text } | mix(ch_software_versions) | set{ch_software_versions}

    GENERATE_SOFTWARE_TABLE(
      ch_software_versions | unique | collectFile(newLine: true, sort: true, cache: false),
      ch_runsheet | splitCsv(header: true) | first | map{ row -> row['Array Data File Name'] }
    )

    // export meta for post processing usage
    ch_meta | DUMP_META
}

workflow.onComplete {
    println "${c_bright_green}Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    if ( workflow.success ) {
      println "Raw and Processed data location: ${ params.resultsDir }"
      println "V&V logs location: ${ params.resultsDir }/VV_Logs"
      println "Pipeline tracing/visualization files location:  ${ params.resultsDir }/Resource_Usage${c_reset}"
    }
}
