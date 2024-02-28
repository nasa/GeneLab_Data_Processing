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
include { GENERATE_SOFTWARE_TABLE } from './modules/GENERATE_SOFTWARE_TABLE'

/**************************************************
* HELP MENU  **************************************
**************************************************/
if (params.help) {
  println("┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅")
  println("┇ Microarray Agilent 1 Channel Pipeline: $workflow.manifest.version  ┇")
  println("┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅")
  println("Usage example 1: Processing GLDS datasets")
  println("   > nextflow run ./main.nf --osdAccession OSD-548 --gldsAccession GLDS-548")
  println()
  println("Usage example 2: Processing Other datasets")
  println("   Note: This requires a user-created runsheet.")
  println("   > nextflow run ./main.nf --runsheetPath </path/to/runsheet>")
  println()
  println("arguments:")
  println("  --help                show this help message and exit")
  println("  --osdAccession OSD-000")
  println("                        the OSD accession id to process through the NF_MAAgilent1ch.")
  println("  --gldsAccession GLDS-000")
  println("                        the GLDS accession id to process through the NF_MAAgilent1ch.")
  println("  --runsheetPath        Use a local runsheet instead one automatically generated from a GLDS ISA archive.")
  println("  --skipVV              Skip automated V&V. Default: false")
  println("  --outputDir           Directory to save staged raw files and processed files. Default: <launch directory>")
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
          "${ projectDir }/bin/dp_tools__agilent_1_channel" // dp_tools plugin
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
      PARSE_ANNOTATION_TABLE.out.annotations_db_url,
      ch_meta | map { it.organism },
      params.limit_biomart_query
    )

    VV_AGILE1CH( 
      ch_runsheet, 
      AGILE1CH.out.de,
      params.skipVV,
      "${ projectDir }/bin/dp_tools__agilent_1_channel" // dp_tools plugin
      )

    // Software Version Capturing
    nf_version = "- name: nextflow\n  ".concat(
"""
  version: ${nextflow.version}
  homepage: https://www.nextflow.io
  workflow task: N/A
""")
    ch_software_versions = Channel.value(nf_version)
    AGILE1CH.out.versions | map{ it -> it.text } | mix(ch_software_versions) | set{ch_software_versions}
    VV_AGILE1CH.out.versions | map{ it -> it.text } | mix(ch_software_versions) | set{ch_software_versions}
    ch_software_versions | unique 
                         | collectFile(
                            newLine: true, 
                            sort: true,
                            cache: false
                            )
                         | GENERATE_SOFTWARE_TABLE

    emit:
      meta = ch_meta 
      runsheet = ch_runsheet
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
