// main.nf
nextflow.enable.dsl=2

def colorCodes = [
    c_line: "┅" * 70,
    c_back_bright_red: "\u001b[41;1m",
    c_bright_green: "\u001b[32;1m",
    c_blue: "\033[0;34m",
    c_yellow: "\u001b[33;1m",
    c_reset: "\033[0m"
]

// Command for: 'nextflow run main.nf --version'
if (params.version) {
    println """${workflow.manifest.name}
Workflow Version: ${workflow.manifest.version}"""
    exit 0
}

// Print the pipeline version on start
println """
${colorCodes.c_bright_green}
${workflow.manifest.name}
Workflow Version: ${workflow.manifest.version}
${colorCodes.c_reset}
""".stripIndent()

// Debug warning
println("${colorCodes.c_yellow}")
if (params.limit_samples_to || params.truncate_to || params.force_single_end || params.genome_subsample) {
    println("WARNING: Debugging options enabled!")
    println("Sample limit: ${params.limit_samples_to ?: 'Not set'}")
    println("Read truncation: ${params.truncate_to ? "First ${params.truncate_to} records" : 'Not set'}")
    println("Reference genome subsampling: ${params.genome_subsample ? "Region '${params.genome_subsample}'" : 'Not set'}")
    println("Force single-end analysis: ${params.force_single_end ? 'Yes' : 'No'}")
} else {
    println("No debugging options enabled")
}
println("${colorCodes.c_reset}")

include { RNASEQ } from './workflows/rnaseq.nf'
include { RNASEQ_MICROBES } from './workflows/rnaseq_microbes.nf'

include { GENERATE_MD5SUMS } from './modules/generate_md5sums.nf'
include { UPDATE_ASSAY_TABLE } from './modules/update_assay_table.nf'

ch_dp_tools_plugin = params.dp_tools_plugin ? 
    Channel.value(file(params.dp_tools_plugin)) : 
    Channel.value(file(params.mode == 'microbes' ? 
        "$projectDir/bin/dp_tools__NF_RCP_Bowtie2" : 
        "$projectDir/bin/dp_tools__NF_RCP"))

ch_accession = params.accession ? Channel.value(params.accession) : null
ch_runsheet = params.runsheet_path ? Channel.fromPath(params.runsheet_path) : null
ch_isa_archive = params.isa_archive_path ? Channel.fromPath(params.isa_archive_path) : null
ch_reference_table = Channel.value(params.reference_table)
ch_api_url = Channel.value(params.api_url)

ch_truncate_to = Channel.value(params.truncate_to)
ch_genome_subsample = Channel.value(params.genome_subsample)
ch_force_single_end = Channel.value(params.force_single_end)

ch_reference_store_path = Channel.value(params.reference_store_path)
ch_derived_store_path = Channel.value(params.derived_store_path)

// set reference params
ch_reference_source = params.reference_source ? Channel.value(params.reference_source) : null
ch_reference_version = params.reference_version ? Channel.value(params.reference_version) : null
ch_reference_fasta = params.reference_fasta ? Channel.fromPath(params.reference_fasta) : null
ch_reference_gtf = params.reference_gtf ? Channel.fromPath(params.reference_gtf) : null

// Set outdir based on the presence of an accession input.. currently not implemented as nextflow.config's outdir is only used for nextflow run info
//  and ch_outdir is used for all processes' publishDir in order to set either ./results or ./{accession}
ch_outdir = params.outdir ? channel.fromPath(params.outdir, checkIfExists: true) : null

// Main workflows
workflow {
    if (params.mode == 'microbes') {
        RNASEQ_MICROBES(
            ch_outdir,
            ch_dp_tools_plugin,
            ch_reference_table,
            ch_accession,
            ch_isa_archive,
            ch_runsheet,
            ch_api_url,
            ch_force_single_end,
            ch_truncate_to,
            ch_reference_source,
            ch_reference_version,
            ch_reference_fasta,
            ch_reference_gtf,
            ch_reference_store_path,
            ch_derived_store_path
        )
    } else {
        RNASEQ(
            ch_outdir,
            ch_dp_tools_plugin,
            ch_reference_table,
            ch_accession,
            ch_isa_archive,
            ch_runsheet,
            ch_api_url,
            ch_force_single_end,
            ch_truncate_to,
            ch_reference_source,
            ch_reference_version,
            ch_reference_fasta,
            ch_reference_gtf,
            ch_reference_store_path,
            ch_derived_store_path
        )
    }
}

workflow POST_PROCESSING {
  main:
    ch_processed_directory = Channel.fromPath("${ params.outdir }/${ params.accession }", checkIfExists: true)
    ch_runsheet = Channel.fromPath("${ params.outdir }/${ params.accession }/Metadata/*_runsheet.csv", checkIfExists: true)
    GENERATE_MD5SUMS(ch_processed_directory)
    UPDATE_ASSAY_TABLE(ch_processed_directory)
}