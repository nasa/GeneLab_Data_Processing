nextflow.enable.dsl=2

/**************************************************
* HELP MENU  **************************************
**************************************************
if (params.help) {
  println("┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅")
  println("┇ RNASeq Consensus Pipeline: $workflow.manifest.version  ┇")
  println("┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅")
  println("Usage example 1: Processing GLDS datasets using genome fasta and gtf from Ensembl")
  println("   > nextflow run ./main.nf --gldsAccession GLDS-194 -resume -profile conda --paired true")
  println()
  println("Usage example 2: Processing GLDS datasets using local genome fasta and gtf")
  println("   Note: ensemblVersion and ref_source are used here to label subdirectories for derived reference files.")
  println("   > nextflow run ./main.nf --gldsAccession GLDS-194 --ensemblVersion 96 --ref_source <reference_label>  --ref_fasta </path/to/fasta> --ref_gtf </path/to/gtf>")
  println()
  println("Usage example 3: Processing Other datasets")
  println("   Note: This requires a user-created runsheet.")
  println("   > nextflow run ./main.nf --runsheetPath </path/to/runsheet>")
  println()
  println("arguments:")
  println("  --help                show this help message and exit")
  println("  --gldsAccession GLDS-000")
  println("                        the GLDS accession id to process through the RNASeq Concensus Pipeline.")
  println("  --runsheetPath        Use a local runsheet instead one automatically generated from a GLDS ISA archive.")
  println("  --ensemblVersion n    Specifies the ensembl Version to use for the reference genome. The default version is ")
  println("  --skipVV              Skip automated V&V. Default: false")
  println("  --paired              Are the input reads paired-end. Default: true. set to false if single-end")
  println("  --outputDir           Directory to save staged raw files and processed files. Default: <launch directory>")
  exit 0
  }

println "PARAMS: $params"
println "\n"
println "Storing any newly fetched primary references files here: ${params.referenceStorePath}"
println "Storing any newly generated derived reference files here: ${params.derivedStorePath}"

/**************************************************
* CHECK REQUIRED PARAMS AND LOAD  *****************
**************************************************
// Get all params sourced data into channels
// Set up channel containing glds accession number
if ( params.gldsAccession ){ 
    ch_glds_accession = Channel.from( params.gldsAccession )
    } else { 
      exit 1, "Missing Required Parameter: gldsAccession. Example for setting on CLI: --gldsAccession GLDS-194"
  }

// Check conditionally required parameter (if using direct fasta, an ensemblVersion must also be supplied)
if ( params.ref_fasta ) {
  if ( !params.ensemblVersion ) { exit 1, "Missing Required Parameter: ensemblVersion. Example for setting on CLI: --ensemblVersion 96" }
}

if ( !params.outputDir ) {  params.outputDir = "$workflow.launchDir" }

ch_multiqc_config = params.multiqcConfig ? Channel.fromPath( params.multiqcConfig ) : Channel.fromPath("NO_FILE")



*/

// color defs
c_back_bright_red = "\u001b[41;1m";
c_bright_green = "\u001b[32;1m";
c_blue = "\033[0;34m";
c_reset = "\033[0m";

// Processes to create the required database(s) if not provided
include { SETUP_CAT_DB; SETUP_KOFAMSCAN_DB; SETUP_GTDBTK_DB; 
          SETUP_CHOCOPHLAN; SETUP_UNIREF; SETUP_UTILITY_MAPPING;
          SETUP_METAPHLAN } from "./modules/database_creation.nf"
include { make_humann_db } from "./modules/database_creation.nf"

// Read quality check and filtering
include { quality_check as raw_qc; BBDUK } from "./modules/quality_assessment.nf"
include { quality_check as filtered_qc } from "./modules/quality_assessment.nf"

// Read-based workflow
include { read_based } from "./modules/read_based_processing.nf"

// Assembly-based workflow
include { assembly_based } from "./modules/assembly_based_processing.nf"


// Workflow to perform read-based analysis
workflow run_read_based_analysis {


    take:
        filtered_ch

    main:

     if(!params.database.chocophlan_dir ||!params.database.uniref_dir || 
        !params.database.metaphlan_db_dir || !params.database.utilities_dir) {
              
         make_humann_db()
         read_based(filtered_ch, 
                    make_humann_db.out.chocophlan_dir,
                    make_humann_db.out.uniref_dir,
                    make_humann_db.out.metaphlan_db_dir,
                    make_humann_db.out.utilities_dir)
       }else{

         read_based(filtered_ch, 
                    params.database.chocophlan_dir,
                    params.database.uniref_dir,
                    params.database.metaphlan_db_dir,
                    params.database.utilities_dir)
      }

}

// Workflow to perform assembly-based analysis
workflow run_assembly_based_analysis {

    take:
        file_ch
        filtered_ch


    main:
        kofam_db = params.database.ko_db_dir
         if(!params.database.ko_db_dir) {
             SETUP_KOFAMSCAN_DB()
             kofam_db = SETUP_KOFAMSCAN_DB.out.ko_db_dir
         }

         cat_db = params.database.cat_db
         if(!params.database.cat_db){

            SETUP_CAT_DB(params.database.CAT_DB_LINK)
            cat_db = SETUP_CAT_DB.out.cat_db
         }

         gtdbtk_db_dir = params.database.gtdbtk_db_dir
         if(!params.database.gtdbtk_db_dir){
              SETUP_GTDBTK_DB()
              gtdbtk_db_dir = SETUP_GTDBTK_DB.out.gtdbtk_db_dir
         }

        // Run assembly based workflow 
        assembly_based(file_ch, filtered_ch, kofam_db, 
                        cat_db, gtdbtk_db_dir, params.use_gtdbtk_scratch_location)

}



// A function to delete white spaces from an input string and covert it to lower case 
def deleteWS(string){

    return string.replaceAll(/\s+/, '').toLowerCase()

}

// Main workflow
workflow {
        
      // Parse file input
       if(params.GLDS_accession){

       GET_RUNSHEET()
       GET_RUNSHEET.out.input_file
           .splitCsv(header:true)
           .set{file_ch}

       GET_RUNSHEET.out.params_file
                     .splitCsv(header:true)
                     .set{params_ch} 


      }else{
 
       Channel.fromPath(params.csv_file, checkIfExists: true)
           .splitCsv(header:true)
           .set{file_ch}
      }


    file_ch.map{
                     row -> deleteWS(row.paired) == 'true'  ? tuple( "${row.sample_id}", [file("${row.forward}"), file("${row.reverse}")], deleteWS(row.paired)) : 
                                         tuple( "${row.sample_id}", [file("${row.forward}")], deleteWS(row.paired))
                }.set{reads_ch}
    //reads_ch.view()
    //return 

    // Qality check and trim the input reads
    raw_qc(Channel.of("raw"), params.multiqc_config,reads_ch)
    filtered_ch = BBDUK(reads_ch, params.adapters)
    filtered_qc(Channel.of("filtered"), params.multiqc_config, filtered_ch)

   // Run the analysis based on selection i.e, read-based, assembly-based or both
    // it will run both by default
    if(params.workflow == 'read-based'){
          run_read_based_analysis(filtered_ch)
    }else if(params.workflow == 'assembly-based') {
          run_assembly_based_analysis(file_ch,filtered_ch)
    }else{
          run_read_based_analysis(filtered_ch)
          run_assembly_based_analysis(file_ch, filtered_ch)
    }

}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Workflow completed without any error\n" : "Oops .. something went wrong" )
}
