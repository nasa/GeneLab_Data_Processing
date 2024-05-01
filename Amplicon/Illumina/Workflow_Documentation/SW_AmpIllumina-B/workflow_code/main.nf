#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// color defs
c_back_bright_red = "\u001b[41;1m";
c_bright_green = "\u001b[32;1m";
c_blue = "\033[0;34m";
c_reset = "\033[0m";



// create GLD runsheet

include { GET_RUNSHEET } from "./modules/create_runsheet.nf"

// Read quality check and filtering
include { FASTQC as RAW_FASTQC ; MULTIQC as RAW_MULTIQC  } from './modules/quality_assessment.nf'
include { CUTADAPT; COMBINE_CUTADAPT_LOGS_AND_SUMMARIZE } from './modules/quality_assessment.nf'
include { FASTQC as TRIMMED_FASTQC ; MULTIQC as TRIMMED_MULTIQC  } from './modules/quality_assessment.nf'

// Cluster ASvs
include { RUN_R_TRIM; RUN_R_NOTRIM } from './modules/run_dada.nf'
include { ZIP_BIOM } from './modules/zip_biom.nf'
include { R_VISUALIZATION } from './modules/visualization.nf'



// A function to delete white spaces from an input string and covert it to lower case
def deleteWS(string){

    return string.replaceAll(/\s+/, '').toLowerCase()

}



workflow {

   if(params.GLDS_accession){

       GET_RUNSHEET()
       GET_RUNSHEET.out.input_file
           .splitCsv(header:true)
           .set{file_ch}

       GET_RUNSHEET.out.params_file
                     .splitCsv(header:true)
                     .set{params_ch}

       target_region = params_ch.map{row -> "${row.target_region}"}.first()
       primers_ch =  params_ch.map{
                           row -> ["${row.f_primer}", "${row.r_primer}"] 
                           }.first()                 
       raw_read_suffix_ch =  params_ch.map{
                           row -> "${row.data_type}" == "PE" ? ["${row.raw_R1_suffix}", "${row.raw_R2_suffix}"] : ["${row.raw_R1_suffix}"] 
                           }.first() 


   }else{

        Channel.fromPath(params.csv_file, checkIfExists: true)
           .splitCsv(header:true)
           .set{file_ch}
   }

    file_ch.map{
                     row -> deleteWS(row.paired)  == 'true' ? tuple( "${row.sample_id}", [file("${row.forward}"), file("${row.reverse}")], deleteWS(row.paired)) : 
                                         tuple( "${row.sample_id}", [file("${row.forward}")], deleteWS(row.paired))
                }.set{reads_ch} 

     //reads_ch.view()
     //return
    // Generating a file with sample ids on a new line
    file_ch.map{row -> "${row.sample_id}"}
              .collectFile(name: "${baseDir}/unique-sample-IDs.txt", newLine: true)
              .set{sample_ids_ch}

    // Read quality check and trimming
    raw_fastqc_files = RAW_FASTQC(reads_ch).flatten().collect()
    RAW_MULTIQC("raw", params.multiqc_config,raw_fastqc_files)

    if(params.trim_primers){

        if(!params.GLDS_accession) primers_ch = Channel.value([params.F_primer, params.R_primer])
        CUTADAPT(reads_ch, primers_ch)
        logs = CUTADAPT.out.logs.map{ sample_id, log -> file("${log}")}.collect()
        counts = CUTADAPT.out.trim_counts.map{ sample_id, count -> file("${count}")}.collect()
        trimmed_reads = CUTADAPT.out.reads.map{ 
                                              sample_id, reads, isPaired -> reads instanceof List ? reads.each{file("${it}")}: file("${reads}")
                                              }.flatten().collect()

        COMBINE_CUTADAPT_LOGS_AND_SUMMARIZE(counts, logs)
        trimmed_fastqc_files = TRIMMED_FASTQC(CUTADAPT.out.reads).flatten().collect()
        TRIMMED_MULTIQC("filtered", params.multiqc_config, trimmed_fastqc_files)

        isPaired_ch = CUTADAPT.out.reads.map{ 
                                              sample_id, reads, isPaired -> isPaired
                                              }.first()

        samples_ch = sample_ids_ch.first()
                     .concat(isPaired_ch)
                     .collate(2)
        // Run dada2
        RUN_R_TRIM(samples_ch, trimmed_reads, COMBINE_CUTADAPT_LOGS_AND_SUMMARIZE.out.counts)

        dada_counts = RUN_R_TRIM.out.counts
        dada_taxonomy = RUN_R_TRIM.out.taxonomy
        dada_biom = RUN_R_TRIM.out.biom

    }else{

        raw_reads_ch = reads_ch.map{
                          sample_id, reads, isPaired -> reads instanceof List ? reads.each{file("${it}")}: file("${reads}")
                          }.flatten().collect()

        if(!params.GLDS_accession) {
            raw_read_suffix_ch =  reads_ch.map{
                                      sample_id, reads, isPaired -> isPaired  == 'true' ? [params.raw_R1_suffix, params.raw_R2_suffix] : [params.raw_R1_suffix]
                                         }
        }

        isPaired_ch = reads_ch.map{sample_id, reads, isPaired -> isPaired}.first()
        samples_ch = sample_ids_ch.first()
                     .concat(isPaired_ch)
                     .collate(2)
        RUN_R_NOTRIM(samples_ch, raw_reads_ch, raw_read_suffix_ch)

        dada_counts = RUN_R_NOTRIM.out.counts
        dada_taxonomy = RUN_R_NOTRIM.out.taxonomy
        dada_biom = RUN_R_NOTRIM.out.biom

    }


    // Zip biom file
    ZIP_BIOM(dada_biom)

    if(params.enable_visualizations){
        // Visualize
        runsheet = params.GLDS_accession ? GET_RUNSHEET.out.runsheet : params.runsheet
        R_VISUALIZATION(runsheet, sample_ids_ch, dada_counts, dada_taxonomy)

    }

}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Workflow completed without any error\n" : "Oops .. something went wrong" )
}
