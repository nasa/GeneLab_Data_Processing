// main.nf
nextflow.enable.dsl=2

include { paramsHelp         } from 'plugin/nf-schema'
include { validateParameters } from 'plugin/nf-schema'

include { PARSE_HOSTS_TABLE  } from './modules/parse_hosts_table.nf'

include { KRAKEN2_DB         } from './modules/kraken2_db.nf'
include { KRAKEN_2           } from './modules/kraken2.nf'
include { SUMMARY            } from './modules/summary.nf'
include { COMPILE_SUMMARY    } from './modules/summary.nf'

include { SOFTWARE_VERSIONS  } from './modules/software_versions.nf'
include { GENERATE_PROTOCOL  } from './modules/generate_protocol.nf'

workflow {
    
    main:

    // check input parameters
    validateParameters()

    // Capture software versions
    software_versions_ch = channel.empty()

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

    channel
        .fromPath(params.sample_id_list)
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

    KRAKEN_2(database_ch, generated_reads_ch, params.out_suffix)
    
    // Generate summary and compile into one file
    SUMMARY(KRAKEN_2.out.output, KRAKEN_2.out.report)
    COMPILE_SUMMARY(SUMMARY.out.collect(), channel.fromPath(params.sample_id_list), host_id) 

    // Software Version Capturing - combining all captured software versions
    KRAKEN_2.out.version | mix(KRAKEN2_DB.out.version.ifEmpty(""))
                         | mix(software_versions_ch)
                         | set{software_versions_ch}

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
    
    publish:
    protocol_out = protocol_out
    software_versions = SOFTWARE_VERSIONS.out
    fastq_out = KRAKEN_2.out.host_removed
    kraken2_out = KRAKEN_2.out.output
    kraken2_report = KRAKEN_2.out.report
    summary_stats = COMPILE_SUMMARY.out.summary_file
    
}

output {
    protocol_out {
        path "GeneLab"
    }

    software_versions {
        path "GeneLab"
    }

    fastq_out {
        path "${params.reads_outdir}"
    }

    kraken2_out {
        path "results/kraken2-output"
    }

    kraken2_report {
        path "results/kraken2-output"
    }

    summary_stats {
        path "results"
    }
}