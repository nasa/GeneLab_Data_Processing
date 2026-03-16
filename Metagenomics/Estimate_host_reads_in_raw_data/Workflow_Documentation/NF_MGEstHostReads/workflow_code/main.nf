// main.nf
nextflow.enable.dsl=2

include { paramsHelp         } from 'plugin/nf-schema'
include { validateParameters } from 'plugin/nf-schema'

include { FETCH_ISA          } from './modules/fetch_isa.nf'
include { ISA_TO_RUNSHEET    } from './modules/isa_to_runsheet.nf'
include { PARSE_RUNSHEET     } from './parse_runsheet.nf'
include { COPY_READS         } from './modules/copy_reads.nf'

include { PARSE_HOSTS_TABLE  } from './modules/parse_hosts_table.nf'

include { KRAKEN2_DB         } from './modules/kraken2_db.nf'
include { KRAKEN_2           } from './modules/kraken2.nf'
include { SUMMARY            } from './modules/summary.nf'
include { COMPILE_SUMMARY    } from './modules/summary.nf'

include { UPDATE_ASSAY_TABLE } from './modules/update_assay_table.nf'

include { SOFTWARE_VERSIONS  } from './modules/software_versions.nf'
include { GENERATE_PROTOCOL  } from './modules/generate_protocol.nf'


dp_tools_plugin_ch = params.dp_tools_plugin ?
    channel.value(file(params.dp_tools_plugin)) : 
    channel.value(file("$projectDir/bin/dp_tools__metagenomics_estHost"))

runsheet_ch = params.runsheet_path ? channel.fromPath(params.runsheet_path) : null
isa_archive_ch = params.isa_archive_path ? channel.fromPath(params.isa_archive_path) : null


workflow {
    main:
    
    // Check input parameters
    validateParameters()

    // Capture software versions
    software_versions_ch = channel.empty()

    isa_archive_out = channel.empty()
    runsheet_out = channel.empty()
    if (params.osd) {
        if (runsheet_ch == null) { // if runsheet_path is not provided, set it up from ISA input
            if (isa_archive_ch == null) { // if isa_archive_path is not provided, fetch the ISA
                FETCH_ISA()
                isa_archive_out = FETCH_ISA.out.isa_archive
            } else { // isa_archive_path is already a channel, use it directly
                isa_archive_out = isa_archive_ch
            }
            ISA_TO_RUNSHEET( isa_archive_out, dp_tools_plugin_ch )
            runsheet_out = ISA_TO_RUNSHEET.out.runsheet
            ISA_TO_RUNSHEET.out.version | mix(software_versions_ch) | set{software_versions_ch}
        } else { // runsheet_path is already a channel, use it directly, and use ISA.zip path if provided
            runsheet_out = runsheet_ch
            isa_archive_out = isa_archive_ch ?: channel.empty() // fallback to empty if null
        }

        // Get sample metadata from runsheet
        PARSE_RUNSHEET( runsheet_out )

        samples = PARSE_RUNSHEET.out.samples

        COPY_READS( samples )
        generated_reads_ch = COPY_READS.out.raw_reads
    }
    else {
        channel
            .fromPath(params.sample_id_list, checkIfExists: true)
            .splitText()
            .map { it.trim() }
            .filter { it } //Ignore blank lines
            .map { sample_id ->
                def meta = [id: sample_id, paired_end: !params.is_single]
                
                def reads = params.is_single ?
                    file("$params.reads_dir/${sample_id}$params.single_suffix") :
                    [file("$params.reads_dir/${sample_id}$params.R1_suffix"), file("$params.reads_dir/${sample_id}$params.R2_suffix")]

                return tuple(meta, reads)
            }
            .set {generated_reads_ch}
        }

    def host_id = params.host.replaceAll(' ', '_').toLowerCase()
    def host_db = file("${params.ref_dbs_dir}/kraken2-${host_id}-db")

    if (host_db.exists()) {
        // Database already exists in storeDir - assembly name and accession might be used by GENERATE_PROTOCOL (if mainfest.txt doesn't exist in db folder)
        reference_fasta     = channel.value("")
        reference_genome    = params.assembly_name ? channel.value(params.assembly_name) : channel.value("")
        reference_accession = params.assembly_acc  ? channel.value(params.assembly_acc)  : channel.value("")

    } else if ( params.db_url ) {
        // Option 1: pre-built DB - only assembly_name and assembly_acc needed
        reference_fasta       = channel.value("")
        reference_genome      = params.assembly_name  ? channel.value(params.assembly_name) : channel.value("unknown")
        reference_accession   = params.assembly_acc   ? channel.value(params.assembly_acc)  : channel.value("unknown")

    } else if ( params.ref_fasta ) {
        // Option 2: custom FASTA - assembly_name and assembly_acc should also be set
        reference_fasta       = channel.value(params.ref_fasta)
        reference_genome      = params.assembly_name  ? channel.value(params.assembly_name) : channel.value("unknown")
        reference_accession   = params.assembly_acc   ? channel.value(params.assembly_acc)  : channel.value("unknown")

    } else {
        // Option 3: parse hosts.csv
        PARSE_HOSTS_TABLE(params.hosts_table, host_id)
        reference_fasta       = PARSE_HOSTS_TABLE.out.reference_fasta_url
        reference_genome      = PARSE_HOSTS_TABLE.out.reference_genome
        reference_accession   = PARSE_HOSTS_TABLE.out.reference_accession
        // Note: for hosts like human where genome/accession are not in hosts.csv,
        // PARSE_HOSTS_TABLE will emit null - protocol will fall back to manifest.txt
    }

    KRAKEN2_DB(host_id, reference_fasta)
    database_ch = KRAKEN2_DB.out.build.first()

    KRAKEN_2(database_ch, generated_reads_ch)

    // Create sample IDs channel for generating summary depending on mode
    if (params.osd) {
        sample_ids_ch = PARSE_RUNSHEET.out.samples
                            .map { meta, reads -> meta.id }
                            .collectFile(name: "sample_ids.txt", newLine: true)
    } else {
        sample_ids_ch = channel.fromPath(params.sample_id_list)
    }
    
    // Generate summary and compile into one file
    SUMMARY(KRAKEN_2.out.output)
    COMPILE_SUMMARY(SUMMARY.out.sample_stats.collect(), sample_ids_ch, host_id)
    
    // Software Version Capturing - combining all captured software versions
    software_versions_ch = software_versions_ch
                            | mix(KRAKEN_2.out.version)
                            | mix(KRAKEN2_DB.out.version.ifEmpty(""))

    nf_version = "Nextflow Version ".concat("${nextflow.version}")
    nextflow_version_ch = channel.value(nf_version)

    //  Write software versions to file
    software_versions_ch | filter { it }
                         | map { it -> it.text.strip() }
                         | unique
                         | mix(nextflow_version_ch)
                         | collectFile({it -> it}, newLine: true, cache: false)
                         | SOFTWARE_VERSIONS
    
    protocol_out = GENERATE_PROTOCOL(
                        params.host, 
                        reference_genome, 
                        reference_accession,
                        SOFTWARE_VERSIONS.out,
                        database_ch)
    
    // Updated assay table (that includes an added "Parameter Value: Host Contamination" column) should be done only when running on datasets from OSDR
    assay_table_out = channel.empty()
    if (params.osd) {
        if (params.assay_table) {
            isa_input = channel.fromPath(params.assay_table)
        } else {
            isa_input = isa_archive_out
        }
        assay_table_out = UPDATE_ASSAY_TABLE(runsheet_out, isa_input, COMPILE_SUMMARY.out.summary_file)
    }
    

    publish:
    // Metadata
    isa_archive = isa_archive_out
    runsheet = runsheet_out

    // Kraken2 output
    kraken2_out = KRAKEN_2.out.output
    kraken2_report = KRAKEN_2.out.report
    summary_stats = COMPILE_SUMMARY.out.summary_file

    // Updated assay table
    assay_table = assay_table_out
    
    // Processing info
    protocol_out = protocol_out
    software_versions = SOFTWARE_VERSIONS.out
}

output {
    // Metadata
    isa_archive {
        path "Metadata"
    }

    runsheet {
        path "Metadata"
    }

    // Kraken2 output
    kraken2_out {
        path "results/kraken2-output"
    }

    kraken2_report {
        path "results/kraken2-output"
    }

    summary_stats {
        path "results"
    }

    assay_table {
        path "GeneLab"
    }

    // Processing info
    protocol_out {
        path "GeneLab"
    }

    software_versions {
        path "GeneLab"
    }
}