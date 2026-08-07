nextflow.enable.dsl=2

include { validateParameters } from 'plugin/nf-schema'
include { paramsSummaryLog } from 'plugin/nf-schema'

include { PURGE_PROCESSING_INFO } from './modules/purge_processing_info.nf'
include { GENERATE_MD5SUMS } from './modules/GENERATE_MD5SUMS.nf'
include { UPDATE_ISA_TABLES } from './modules/UPDATE_ISA_TABLES.nf'

/**************************************************
* WORKFLOW SPECIFIC PRINTOUTS  ********************
**************************************************/

workflow {

  
  main:

    // color defs
    c_back_bright_red = "\u001b[41;1m";
    c_reset = "\033[0m";


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

    // validate parameters and print parameter summary log (includes only parameters that are not set to default values)
    validateParameters(cast_cli_params: true)
    log.info paramsSummaryLog(workflow)

    /**************************************************
    * CHECK REQUIRED PARAMS AND LOAD  *****************
    **************************************************/
    
    // Resolve processed directory (mirrors ch_outdir handling in main.nf)
    def matches = params.accession
      ? file("GLDS-*", type: 'dir')
      : file("results", type: 'dir')

    def matches_list = matches instanceof List ? matches : [matches]
    if (matches_list.size() == 0) {
      error "No matching output directory found (looked for ${params.accession ? 'GLDS-*' : 'results'})"
    }
    if (matches_list.size() > 1) {
      error "Expected exactly one output directory, but found ${matches_list.size()}: ${matches_list}"
    }

    def processed_dir = matches_list[0]
    println "Resolved output directory: ${processed_dir}"
    if (!processed_dir?.exists()) {
      error "No matching output directory found (looked for ${params.accession ? 'GLDS-*' : 'results'})"
    }

    // Extract GLDS accession from the directory name, if present
    def dir_name = processed_dir.name
    def accession_matcher = (dir_name =~ /^(GLDS-\d+)$/)
    def glds_accession = accession_matcher.matches() ? accession_matcher.group(1) : null

    if (params.accession && !glds_accession) {
      error "Expected a GLDS-# directory name since params.accession was set, but got: ${dir_name}"
    }

    println "Resolved GLDS accession: ${glds_accession ?: 'none (results dir)'}"
  
    ch_processed_directory = channel.fromPath(processed_dir, checkIfExists: true)
    ch_runsheet = channel.fromPath("${ processed_dir }/Metadata/*_runsheet.csv", checkIfExists: true)
    ch_processing_info = channel.fromPath("$launchDir/processing_scripts/nextflow_processing_info_GLmicroarray.txt", checkIfExists: true)
    ch_glds_accession = channel.value(glds_accession ?: 'NA')


    PURGE_PROCESSING_INFO(
      processed_dir,
      ch_processing_info
    )
    
    GENERATE_MD5SUMS(      
      ch_processed_directory,
      PURGE_PROCESSING_INFO.out.purged_processing_info
    )

    def isa_file = file("${ processed_dir }/Metadata/*ISA*.zip")
    if ( isa_file ) {
      ch_isa = channel.fromPath("${ processed_dir }/Metadata/*ISA*.zip")
      UPDATE_ISA_TABLES(
        processed_dir,
        ch_runsheet,
        ch_isa,
        ch_glds_accession
      )
    } else {
      println "${ c_back_bright_red }WARNING: No ISA archive found in ${ processed_dir }/Metadata/ -- skipping UPDATE_ISA_TABLES${ c_reset }"
    }
    
}