nextflow.enable.dsl=2
// color defs
c_back_bright_red = "\u001b[41;1m";
c_bright_green = "\u001b[32;1m";
c_blue = "\033[0;34m";
c_reset = "\033[0m";

include { GENERATE_MD5SUMS } from './modules/GENERATE_MD5SUMS.nf'
include { UPDATE_ISA_TABLES } from './modules/UPDATE_ISA_TABLES.nf'

/**************************************************
* HELP MENU  **************************************
**************************************************/
if (params.help) {
  println("┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅")
  println("┇ Microarray Affymetrix Post Processing Pipeline: $workflow.manifest.version  ┇")
  println("┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅")
  println("Post processing workflow. Generates md5sum of output files and updates ISA archive tables. Help menu refinements to come")
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
    ch_processed_directory = Channel.fromPath("${ params.outputDir }/${ params.gldsAccession }", checkIfExists: true)
    ch_runsheet = Channel.fromPath("${ params.outputDir }/${ params.gldsAccession }/Metadata/*_runsheet.csv", checkIfExists: true)
    GENERATE_MD5SUMS(      
      ch_processed_directory, 
      ch_runsheet,       
      "${ projectDir }/bin/dp_tools__affymetrix_channel" // dp_tools plugin
    )
    UPDATE_ISA_TABLES(
      ch_processed_directory, 
      ch_runsheet,       
      "${ projectDir }/bin/dp_tools__affymetrix_channel" // dp_tools plugin
    )
}