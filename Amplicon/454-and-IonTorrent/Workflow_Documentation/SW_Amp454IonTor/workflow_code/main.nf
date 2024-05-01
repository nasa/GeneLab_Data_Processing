#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// color defs
c_back_bright_red = "\u001b[41;1m";
c_bright_green = "\u001b[32;1m";
c_blue = "\033[0;34m";
c_reset = "\033[0m";


// Read quality check and filtering
include { FASTQC as RAW_FASTQC ; MULTIQC as RAW_MULTIQC  } from './modules/quality_assessment.nf'
include { CUTADAPT; COMBINE_CUTADAPT_LOGS_AND_SUMMARIZE } from './modules/quality_assessment.nf'
include { BBDUK; COMBINE_BBDUK_LOGS_AND_SUMMARIZE } from './modules/quality_assessment.nf'
include { FASTQC as FILTERED_FASTQC ; MULTIQC as FILTERED_MULTIQC  } from './modules/quality_assessment.nf'
include { pick_otus } from './modules/vsearch.nf'
include { RUN_R} from './modules/assign_taxonomy.nf'
include { ZIP_BIOM } from './modules/zip_biom.nf'


workflow {

    Channel.fromPath(params.csv_file, checkIfExists: true)
           .splitCsv(header:true)
           .map{row -> tuple( "${row.sample_id}", [file("${row.read}")] )}
          .set{reads_ch} 

    // Read quality check and trimming
    raw_fastqc_files = RAW_FASTQC(reads_ch).flatten().collect()
    RAW_MULTIQC("raw", raw_fastqc_files)
    
    // Trim reads
    CUTADAPT(reads_ch)
    trim_counts = CUTADAPT.out.trim_counts.map{ sample_id, count -> file("${count}")}.collect()
    trim_logs = CUTADAPT.out.logs.map{ sample_id, log -> file("${log}")}.collect()
    COMBINE_CUTADAPT_LOGS_AND_SUMMARIZE(trim_counts, trim_logs)
    
    // Filter reads
    BBDUK(CUTADAPT.out.reads)
    filter_counts = BBDUK.out.filter_counts.map{ sample_id, count -> file("${count}")}.collect()
    filter_logs = BBDUK.out.logs.map{ sample_id, log -> file("${log}")}.collect()
    COMBINE_BBDUK_LOGS_AND_SUMMARIZE(filter_counts, filter_logs)

    filtered_fastqc_files = FILTERED_FASTQC(BBDUK.out.reads).flatten().collect()
    FILTERED_MULTIQC("filtered", filtered_fastqc_files)
  
    // Pick outs with vsearch
    pick_otus(BBDUK.out.reads)

    // Assign taxonomy
    RUN_R(pick_otus.out.otus, pick_otus.out.counts, 
          COMBINE_CUTADAPT_LOGS_AND_SUMMARIZE.out.counts,
           COMBINE_BBDUK_LOGS_AND_SUMMARIZE.out.counts)

    // Zip biom file
    ZIP_BIOM(RUN_R.out.biom)

}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Workflow completed without any error\n" : "Oops .. something went wrong" )
}
