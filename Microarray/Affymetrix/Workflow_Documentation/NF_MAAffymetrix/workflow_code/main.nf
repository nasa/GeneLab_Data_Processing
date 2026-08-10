nextflow.enable.dsl=2

include { paramsHelp } from 'plugin/nf-schema'
include { validateParameters } from 'plugin/nf-schema'
include { paramsSummaryLog } from 'plugin/nf-schema'

include { STAGE_ANALYSIS } from './subworkflows/stage_analysis.nf'
include { PARSE_ANNOTATION_TABLE } from './modules/parse_annotation_table.nf'
include { VV_AFFYMETRIX } from './modules/vv_affymetrix.nf'
include { PROCESS_AFFYMETRIX } from './modules/process_affymetrix.nf'
include { GENERATE_SOFTWARE_TABLE } from './modules/generate_software_table.nf'
include { DUMP_META } from './modules/dump_meta.nf'
include { GENERATE_PROTOCOL } from './modules/generate_protocol.nf'

ch_dp_tools_plugin = params.dp_tools_plugin ? 
  channel.value(file(params.dp_tools_plugin)) 
  : params.skipDE ? channel.value(file("$projectDir/bin/dp_tools__affymetrix_skipDE")) : channel.value(file("$projectDir/bin/dp_tools__affymetrix"))
ch_isa_archive_path = params.isaArchivePath ? file(params.isaArchivePath) : null
ch_runsheet = params.runsheetPath ? channel.fromPath(params.runsheetPath) : null

ch_outdir = params.outdir ? channel.fromPath(params.outdir, checkIfExists: true) : null

workflow {
	main:

    // color defs
    c_back_bright_red = "\u001b[41;1m";
    c_bright_green = "\u001b[32;1m";
    c_blue = "\033[0;34m";
    c_reset = "\033[0m";

    /**************************************************
    * HELP MENU  **************************************
    **************************************************/
    if (params.help) {
      before_text = """
************************************************
* Microarray Affymetrix Pipeline: ${workflow.manifest.version} *
************************************************

Usage example 1: Processing OSDR datasets
    > nextflow run ./main.nf --osdAccession OSD-266 --gldsAccession GLDS-266

Usage example 2: Processing Other datasets (requires a user-created runsheet)
    > nextflow run ./main.nf --runsheetPath </path/to/runsheet>


"""
      log.info paramsHelp(
        beforeText: before_text, 
        afterText: "For more information, please see the README.md file in the workflow code directory.",
        fullHelp: true)
      exit(0)
      }

      // validate parameters and print parameter summary log (includes only parameters that are not set to default values)
      validateParameters(cast_cli_params: true)
      log.info paramsSummaryLog(workflow)

      //  ---------------------  Sanity Checks ------------------------------------- //
      // Test input requirement
      if (!params.accession &&  !params.runsheetPath){
        error("""${c_back_bright_red}INPUT ERROR! 
                Please supply either an accession (OSD or Genelab number) or an input CSV file
                by passing either to the --accession or --runsheetPath parameter, respectively.
                ${c_reset}""")
      } 

      // Test ISA archive and accession
      if (params.isaArchivePath && !params.accession) {
          error """${c_back_bright_red}INPUT ERROR!
              --isaArchivePath requires --accession to resolve OSD/GLDS accessions
              for the ISA-to-runsheet conversion.${c_reset}"""
      }

      // Stage analysis setup (directory structure, inputs, and raw reads)
      STAGE_ANALYSIS(
          ch_outdir,
          ch_dp_tools_plugin,
          params.accession,
          ch_isa_archive_path,
          ch_runsheet,
          params.api_url
      )
      ch_outdir = STAGE_ANALYSIS.out.ch_outdir
      samples = STAGE_ANALYSIS.out.samples
      array_data_files = STAGE_ANALYSIS.out.array_data_files | map { meta, f -> f } | collect
      runsheet_path = STAGE_ANALYSIS.out.runsheet_path
      isa_archive = STAGE_ANALYSIS.out.isa_archive
      osd_accession = STAGE_ANALYSIS.out.osd_accession
      glds_accession = STAGE_ANALYSIS.out.glds_accession
      dp_tools_version = STAGE_ANALYSIS.out.dp_tools_version

      // Get dataset-wide metadata
      samples | first 
              | map { meta, reads -> meta }
              | set { ch_meta }

      ch_meta | map { meta -> meta.organism_sci }
              | set { organism_sci }


    PARSE_ANNOTATION_TABLE(params.annotation_file_path, organism_sci)

    PROCESS_AFFYMETRIX(
      ch_outdir,
      channel.fromPath( "${ projectDir }/bin/Affymetrix.qmd" ),
      runsheet_path,
      array_data_files,
      PARSE_ANNOTATION_TABLE.out.annotations_db_url,
      PARSE_ANNOTATION_TABLE.out.reference_version_and_source,
      channel.fromPath( params.referenceStorePath ),
      channel.fromPath( params.array_annot_path ),
      params.skipDE
    )

    VV_AFFYMETRIX( 
      ch_outdir,
      runsheet_path,
      PROCESS_AFFYMETRIX.out.de,
      params.skipVV,
      ch_dp_tools_plugin
      )

    // Software Version Capturing
    ch_software_versions = channel.empty()
    nf_version = """\
    - name: nextflow
      version: ${nextflow.version}
      homepage: https://www.nextflow.io
      workflow task: N/A
    """.stripIndent()
    ch_nextflow_version = channel.value(nf_version)
    ch_process_AFFYMETRIX = PROCESS_AFFYMETRIX.out.versions | map{ it -> it.text }
    ch_dp_tools_version = params.skipVV
      ? ( dp_tools_version ? dp_tools_version.map { it -> it.text} : channel.empty() )
      : ( VV_AFFYMETRIX.out.versions | map{ it -> it.text } )

    ch_software_versions = ch_software_versions
      | mix(ch_process_AFFYMETRIX)
      | mix(ch_dp_tools_version)
      | mix(ch_nextflow_version)

    GENERATE_SOFTWARE_TABLE(
      ch_outdir,
      ch_software_versions | unique | collectFile(newLine: true, sort: true, cache: false),
      runsheet_path | splitCsv(header: true, quote: '"') | first | map{ row -> row['Array Data File Name'] },
      params.skipDE
    )

    // export meta for post processing usage
    DUMP_META(ch_outdir, ch_meta)

    GENERATE_PROTOCOL(
      ch_outdir,
      ch_meta,
      ch_software_versions | unique | collectFile(newLine: true, sort: true, cache: false),
      PARSE_ANNOTATION_TABLE.out.reference_version_and_source,
      PARSE_ANNOTATION_TABLE.out.bioconductor_annotations,
      PARSE_ANNOTATION_TABLE.out.annotations_db_info_url,
      params.skipDE
    )

    ch_outdir.subscribe { outdir_value = it }
    workflow.onComplete = { 
      println "${c_bright_green}Pipeline completed at: $workflow.complete"
      println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
      if ( workflow.success ) {
        println "Raw and Processed data location: ${ outdir_value }"
        println "V&V logs location: ${ outdir_value }/VV_Logs"
        println "Pipeline tracing/visualization files location:  ${ outdir_value }/Resource_Usage${c_reset}"
      }
    }
}
