nextflow.enable.dsl=2
// color defs
c_back_bright_red = "\u001b[41;1m";
c_bright_green = "\u001b[32;1m";
c_blue = "\033[0;34m";
c_reset = "\033[0m";

include { GENERATE_MD5SUMS } from './modules/GENERATE_MD5SUMS.nf'
include { UPDATE_ISA_TABLES } from './modules/UPDATE_ISA_TABLES.nf'
include { GENERATE_PROTOCOL } from './modules/POST_PROCESSING/GENERATE_PROTOCOL'

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
println("Resolved output directory: ${ params.resultsDir }")

/**************************************************
* WORKFLOW SPECIFIC PRINTOUTS  ********************
**************************************************/

workflow {
  main:
    ch_processed_directory = Channel.fromPath("${ params.resultsDir }", checkIfExists: true)
    ch_runsheet = Channel.fromPath("${ params.resultsDir }/Metadata/*_runsheet.csv", checkIfExists: true)
    ch_software_versions = Channel.fromPath("${params.resultsDir}/GeneLab/software_versions_GLmicroarray.md", checkIfExists: true)
    ch_processing_meta = Channel.fromPath("${params.resultsDir}/GeneLab/meta.sh", checkIfExists: true)
    GENERATE_MD5SUMS(      
      ch_processed_directory, 
      ch_runsheet,       
      "${ projectDir }/bin/${ params.skipDE ? 'dp_tools__affymetrix_skipDE' : 'dp_tools__affymetrix' }" // dp_tools plugin
    )
    UPDATE_ISA_TABLES(
      ch_processed_directory, 
      ch_runsheet,       
      "${ projectDir }/bin/${ params.skipDE ? 'dp_tools__affymetrix_skipDE' : 'dp_tools__affymetrix' }" // dp_tools plugin
    )
    GENERATE_PROTOCOL(
      ch_software_versions,
      ch_processing_meta,
      params.skipDE
    )
}