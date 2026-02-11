// main.nf
nextflow.enable.dsl=2

// Terminal text color definitions
c_back_bright_red = "\u001b[41;1m";
c_reset           = "\033[0m";

include { FETCH_ISA } from './modules/fetch_isa.nf'
include { ISA_TO_RUNSHEET } from './modules/isa_to_runsheet.nf'

include { PARSE_RUNSHEET } from './parse_runsheet.nf'
include { COPY_READS } from './modules/copy_reads.nf'

include { KRAKEN2_DB } from './modules/kraken2_db.nf'
include { KRAKEN_2 } from './modules/kraken2.nf'
include { SUMMARY } from './modules/summary.nf'

include { SOFTWARE_VERSIONS } from './modules/utils.nf'
include { GENERATE_PROTOCOL } from './modules/generate_protocol.nf'


dp_tools_plugin_ch = params.dp_tools_plugin ?
    Channel.value(file(params.dp_tools_plugin)) : 
    Channel.value(file("$projectDir/bin/dp_tools__metagenomics_estHost"))

runsheet_ch = params.runsheet_path ? Channel.fromPath(params.runsheet_path) : null
isa_archive_ch = params.isa_archive_path ? Channel.fromPath(params.isa_archive_path) : null


workflow {

    if (!params.osd && !params.glds &&  !params.sample_id_list){
       error("""${c_back_bright_red}INPUT ERROR! 
              Please supply either an accession (OSD and Genelab number) or an input CSV file
              by passing either to the --osd and --glds parameters or to --sample_id_list parameter, respectively.
              ${c_reset}""")
    }

    if (params.osd && !params.glds ||  !params.osd && params.glds){
       error("""${c_back_bright_red}INPUT ERROR! 
              Please supply both accession numbers (OSD and GLDS)
              by passing to both --osd and --glds parameters.
              ${c_reset}""")
    }
    
    if (!params.ref_dbs_Dir){
       error("""${c_back_bright_red}INPUT ERROR! 
              Please supply the path to the directory storing kraken2 reference databases
              by passing to the --ref_dbs_Dir parameter.
              ${c_reset}""")
    }

    // Capture software versions
    software_versions_ch = Channel.empty()

    if (params.osd) {
        //If no runsheet path, fetch ISA archive from OSDR (if needed) and convert to runsheet
        if (runsheet_ch == null) {
            if (isa_archive_ch == null) {
                FETCH_ISA()
                isa_archive_ch = FETCH_ISA.out.isa_archive
            }
            ISA_TO_RUNSHEET( isa_archive_ch, dp_tools_plugin_ch )
            runsheet_ch = ISA_TO_RUNSHEET.out.runsheet
            ISA_TO_RUNSHEET.out.version | mix(software_versions_ch) | set{software_versions_ch}
        }

        // Get sample metadata from runsheet
        PARSE_RUNSHEET( runsheet_ch )

        samples = PARSE_RUNSHEET.out.samples

        COPY_READS( samples )
        generated_reads_ch = COPY_READS.out.raw_reads
    }
    else {
        Channel
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
    }

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

    KRAKEN_2(database_ch, generated_reads_ch)
    KRAKEN_2.out.version | mix(software_versions_ch) | set{software_versions_ch}
    
    // Generate summary and compile one file
    SUMMARY(KRAKEN_2.out.output)
    SUMMARY.out
    .collect()
    .subscribe { summary_files ->
      def outfile = file("${params.outdir}/results/Host-read-count-summary.tsv")
      def header = "Sample_ID\tTotal_fragments\tTotal_host_fragments\tPercent_host\n"
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