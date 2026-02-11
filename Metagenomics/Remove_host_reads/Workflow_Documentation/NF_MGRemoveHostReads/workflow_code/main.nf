// main.nf
nextflow.enable.dsl=2

// Terminal text color definitions
c_back_bright_red = "\u001b[41;1m";
c_reset           = "\033[0m";

include { KRAKEN2_DB } from './modules/kraken2_db.nf'
include { KRAKEN_2 } from './modules/kraken2.nf'
include { SUMMARY } from './modules/summary.nf'

include { SOFTWARE_VERSIONS } from './modules/utils.nf'
include { GENERATE_PROTOCOL } from './modules/generate_protocol.nf'

workflow {
    
    if (!params.ref_dbs_Dir){
       error("""${c_back_bright_red}INPUT ERROR! 
              Please supply the path to the directory storing kraken2 reference databases
              by passing --ref_dbs_Dir.
              ${c_reset}""")
    }

    // Capture software versions
    software_versions_ch = Channel.empty()

    // Get host info
    host_info = Channel
        .fromPath(params.hosts_table)
        .splitCsv(header:true)
        .filter { row -> row.name.toLowerCase() == params.host.toLowerCase() }  // match host
        .map { row ->
            def host_id   = row.name.replaceAll(' ', '_').toLowerCase() 
            tuple(row.name, host_id, row.species, row.refseq, row.genome, row.fasta) }

    host_info
    .ifEmpty { error("INPUT ERROR: Host '${params.host}' not found in hosts table '${params.hosts_table}'") }
    
    // Check if kraken2 database already exists or needs to be built
    def host_id = params.host.replaceAll(' ', '_').toLowerCase()
    def host_db = file("${params.ref_dbs_Dir}/kraken2-${host_id}-db")
    def db_exists = host_db.exists()

    if (db_exists) {
        database_ch = Channel.value(host_db)
    }
    else {
        build_ch = host_info.map { name, hostID, species, refseq, genome, fasta -> tuple(name, hostID, fasta) }
        KRAKEN2_DB(build_ch)
        database_ch = KRAKEN2_DB.out.first()
    }
    
    Channel
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

    KRAKEN_2(database_ch, generated_reads_ch)
    KRAKEN_2.out.version | mix(software_versions_ch) | set{software_versions_ch}
    
    // Generate summary and compile one file
    SUMMARY(KRAKEN_2.out.output, KRAKEN_2.out.report)
    SUMMARY.out
    .collect()
    .subscribe { summary_files ->
      def outfile = file("${params.outdir}/results/Host-read-removal-summary.tsv")
      def header = "Sample_ID\tTotal_fragments_before\tTotal_fragments_after\tPercent_host_reads_removed\n"
      outfile.text = header + summary_files.collect { it.text }.join()
      
      // summary.tmp cleanup
      summary_files.each { f ->
            def tmpFile = f.toFile()
            tmpFile.delete()
            } 
    }

    // Software Version Capturing - combining all captured software versions
    nf_version = "Nextflow Version ".concat("${nextflow.version}")
    nextflow_version_ch = Channel.value(nf_version)

    //  Write software versions to file
    software_versions_ch | map { it.text.strip() }
                         | unique
                         | mix(nextflow_version_ch)
                         | collectFile({it -> it}, newLine: true, cache: false)
                         | SOFTWARE_VERSIONS

    // Protocol always needs name, refseq ID, and genome build
    protocol_ch = host_info.map { name, hostID, species, refseq, genome, fasta -> tuple(name, refseq, genome) }
    
    GENERATE_PROTOCOL(protocol_ch, SOFTWARE_VERSIONS.out)

}