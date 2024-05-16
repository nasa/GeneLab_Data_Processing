#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// color defs
c_back_bright_red = "\u001b[41;1m";
c_bright_green = "\u001b[32;1m";
c_blue = "\033[0;34m";
c_reset = "\033[0m";


/**************************************************
* HELP MENU  **************************************
**************************************************/
if (params.help) {
  println()
  println("Nextflow Metagenomics Illumina Consensus Pipeline: $workflow.manifest.version")
  println("USAGE:")
  println("Example 1: Submit and run jobs with slurm in singularity containers.")
  println("   > nextflow run main.nf -resume -profile slurm_sing --csv_file PE_file.csv")
  println()
  println("Example 2: : Submit and run jobs with slurm in conda environments.")
  println("   > nextflow run main.nf -resume -profile slurm_conda --csv_file SE_file.csv")
  println()
  println("Example 3: Run jobs locally in conda environments, supply a GLDS accession, and specify the path to an existing conda environment.")
  println("   > nextflow run main.nf -resume -profile conda --GLDS_accession OSD-456 --conda.qc <path/to/existing/conda/environment>")
  println()
  println("Required arguments:")
  println("""-profile [STRING] What profile should be used to run the workflow. Options are [singularity, docker, conda, slurm_sing, slurm_conda].
	         singularity, docker and conda will run the pipeline locally using singularity, docker, and conda, respectively.
             slurm_sing and slurm_conda will submit and run jobs using slurm in singularity containers and conda environments, respectively. """)			 
  println("--csv_file  [PATH] A 3-column (single-end) or 4-column (paired-end) input file (sample_id, forward, [reverse,] paired). Mandatory if a GLDS accession is not provided.")
  println("   Please see the files: SE_file.csv and PE_file.csv for single-end and paired-end examples, respectively.")
  println("   The sample_id column should contain unique sample ids.")
  println("   The forward and reverse columns should contain the absolute or relative path to the sample's forward and reverse reads.")
  println("   The paired column should be true for paired-end or anything else for single-end reads.")

  println("Optional arguments:")  
  println("  --help  Print this help message and exit")
  println("  --workflow [STRING] Which workflow should be run. Options are one of [read-based, assembly-based, both]. Default: both")
  println("  --publishDir_mode [STRING]  How should nextflow publish file outputs. Options can be found here https://www.nextflow.io/docs/latest/process.html#publishdir Default: link.")
  println("  --errorStrategy [STRING] How should nextflow handle errors. Options can be found here https://www.nextflow.io/docs/latest/process.html#errorstrategy. Default: ignore")
  println("  --swift_1S [BOOLEAN] Setting for trimming recommended when working with Swift 1S libraries.")
  println("    adds `swift=t` setting to bbduk quality trimming/filtering command. For info on this, see example, ")
  println("    https://swiftbiosci.com/wp-content/uploads/2019/03/16-0853-Tail-Trim-Final-442019.pdf.")
  println("    Set to true if data was generated with Swift 1S library prep. Default: false")
  println("  --adapters [PATH] Path to BBtools adapters for reads filtering. Default: config/bbtools_adapters.fa") 
  println("  --multiqc_config [PATH] Path to a custom multiqc config file. Default: config/multiqc.config")
  println("  --use_gtdbtk_scratch_location [BOOLEAN] Should a scratch location be used to store GTDBTK temp files? true or false.")
  println("    Scratch directory for gtdb-tk, if wanting to use disk space instead of RAM, can be memory intensive;")
  println("    see https://ecogenomics.github.io/GTDBTk/faq.html#gtdb-tk-reaches-the-memory-limit-pplacer-crashes")
  println("    leave empty if wanting to use memory, the default, put in quotes the path to a directory that")
  println("    already exists if wanting to use disk space. Default: false")

  println("MAG parameters: MAG filtering cutoffs based on checkm quality assessments (in percent); see https://github.com/Ecogenomics/CheckM/wiki/Reported-Statistics")
  println("	 --min_est_comp [INT] Minimum estimated completion. Default: 90") 
  println("	 --max_est_redund [INT] Minimum estimated redundancy. Default: 10") 
  println("	 --max_est_strain_het [INT] Minimum estimated strain heterogeneity. Default: 50")
  println("	 --reduced_tree [STRING] reduced_tree option for checkm, limits the RAM usage to 16GB; https://github.com/Ecogenomics/CheckM/wiki/Genome-Quality-Commands#tree.")
  println("    'True' for yes, anything else will be considered 'False' and the default full tree will be used. Default: 'True' ")
  println("	 --max_mem [INT] Maximum memory allowed passed to megahit assembler. Can be set either by proportion of available on system, e.g. 0.5")
  println("    or by absolute value in bytes, e.g. 100e9 would be 100 GB. Default: 100e9")
  
  println("	 --pileup_mem [STRING] pileup.sh paramater for calculating contig coverage and depth. Memory used by bbmap's pileup.sh (within the GET_COV_AND_DET process). ")
  println("	   passed as the -Xmx parameter, 20g means 20 gigs of RAM, 20m means 20 megabytes.")
  println("	   5g should be sufficient for most assemblies, but if that rule is failing, this may need to be increased.Default: '5g' ")
  println("	 --block_size [int] Block size variable for CAT/diamond, lower value means less RAM usage; see https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options#memory--performance-options. Default: 4")

  println("File Suffixes:")
  println("      --filtered_suffix [STRING]  Suffix to use for naming your quality filtered reads. Only applicable when input reads are single-end. Default: _filtered.fastq.gz")  
  println("      --filtered_R1_suffix [STRING]  Suffix to use for naming your quality filtered forward reads. Default: _R1_filtered.fastq.gz")
  println("      --filtered_R2_suffix [STRING]  Suffix to use for naming your quality filtered reverse reads. Default: _R2_filtered.fastq.gz")
  println("Output directories:")
  println("      --raw_reads_dir [PATH] Where should the fastqc report of the raw reads be stored. Default: Raw_Sequence_Data/")
  println("      --fastqc_out_dir [PATH] Where should multiqc outputs be stored. Default: FastQC_Outputs/")
  println("      --filtered_reads_dir [PATH] Where should your filtered reads be stored.  Default: Filtered_Sequence_Data/")
  println("      --assembly_based_dir [PATH] Where should the results of assembly-based analysis be stored. Default: Assembly-based_Processing/")
  println("      --assemblies_dir [PATH] Where should your assemblies be stored. Default: Assembly-based_Processing/assemblies/")
  println("      --genes_dir [PATH] Where should the predicted genes from your assemblies be stored. Default: Assembly-based_Processing/predicted-genes/")
  println("      --annotations_and_tax_dir [PATH] Contig taxonomy and annotation directory.  Default: Assembly-based_Processing/annotations-and-taxonomy/")
  println("      --mapping_dir [PATH] Read mapping to assembly directory.  Default: Assembly-based_Processing/read-mapping/")
  println("      --combined_output_dir [PATH] Assembly summuries and reports across samples directory.  Default: Assembly-based_Processing/combined-outputs/")
  println("      --bins_dir [PATH] Assembly bins directory.  Default: Assembly-based_Processing/bins/")
  println("      --MAGs_dir [PATH] Meta assembled genomes (MAGs) directory.  Default: Assembly-based_Processing/MAGs/")
  println("      --read_based_dir [PATH] Read-based analysis outputs directory.  Default: Read-based_Processing/")

  println("Genelab specific arguements:")
  println("      --GLDS_accession [STRING]  A Genelab accession number if the --csv_file parameter is not set. If this parameter is set, it will ignore the --csv_file parameter.")
  println("      --assay_suffix [STRING]  Genelabs assay suffix. Default: _GLmetagenomics.")
  println("      --additional_filename_prefix [STRING] additional prefix to add to output files that describe more than one sample (to make them unique compared to other datasets)")
  println("      include separator at end if adding one, e.g. Swift1S_ if wanted. Default: '' ")

  println("Paths to existing databases and database links.")
  println("CAT database directory strings:")
  println("The strings below will be added to the end of the --database.cat_db path arguement provided below")
  println("      --cat_taxonomy_dir [PATH] CAT taxonomy database directory. Default: 2021-01-07_taxonomy/")
  println("      --cat_db_sub_dir [PATH] CAT database sub directory. Default: 2021-01-07_CAT_database/")
  println("      --CAT_DB_LINK [URL] CAT database online download link. Default: https://tbb.bio.uu.nl/bastiaan/CAT_prepare/CAT_prepare_20210107.tar.gz.")
  println("      --database.cat_db [PATH] Path to CAT database. Example, /path/to/Reference_DBs/CAT_prepare_20210107/. Default: null.")
  println("      --database.ko_db_dir  [PATH] Path to kofam scan database. Example, /path/to/Reference_DBs/kofamscan_db/. Default: null.")
  println("      --database.metaphlan_db_dir [PATH] Path to metaphlan database. Example, /path/to/Reference_DBs/metaphlan4-db/. Default: null.")
  println("      --database.chocophlan_dir [PATH] Path to Humann's chocophlan nucleotide database. Example, /path/to/Reference_DBs/humann3-db/chocophlan/. Default: null.")
  println("      --database.uniref_dir [PATH] Path to Humann's Uniref protein database. Example, /path/to/Reference_DBs/humann3-db/uniref/. Default: null.")
  println("      --database.utilities_dir [PATH] Path to Humann's untilities database. Example, /path/to/Reference_DBs/humann3-db/utility_mapping/.  Default: null.")
  println("      --database.gtdbtk_db_dir  [PATH] Path to GTDBTK database. Example, /path/Reference_DBs/GTDB-tk-ref-db/. Default: null.")
 
  println("Paths to existing conda environments to use otherwise a new one will be created using the yaml file in envs/.")
  println("      --conda.qc [PATH] Path to a conda environment containing fastqc, multiqc, zip and python. Default: null.")
  println("      --conda.humann3 [PATH] Path to a conda environment with humann3 installed. Default: null.")
  println("      --conda.cat  [PATH] Path to a conda environment containing CAT (Contig annotation tool). Default: null.")
  println("      --conda.prodigal [PATH] Path to a conda environment with prodigal installed. Default: null.")
  println("      --conda.metabat [PATH] Path to a conda environment containing metabat. Default: null.")
  println("      --conda.gtdbtk [PATH] Path to a conda environment containing gtdbtk. Default: null.")
  println("      --conda.kegg_decoder [PATH] Path to a conda environment with kegg_decoder installed. Default: null.")
  println("      --conda.megahit  [PATH] Path to a conda environment containing megahit. Default: null.")
  println("      --conda.bit [PATH] Path to a conda environment with bit installed. Default: null.")
  println("      --conda.kofamscan [PATH] Path to a conda environment containing KOFAM SCAN . Default: null.")
  println("      --conda.mapping [PATH] Path to a conda environment with bowtie and samtools installed. Default: null.")
  println("      --conda.checkm [PATH] Path to a conda environment with checkm installed. Default: null.")
  print("Advanced users can edit the nextflow.config file for more control over default settings such container choice, number cpus, memory per task etc.")
  exit 0
  }

log.info """
         Nextflow Metagenomics Illumina Consensus Pipeline: $workflow.manifest.version
         
         You have set the following parameters:
         Profile: ${workflow.profile} 
         Input csv file : ${params.csv_file}
         GLDS Accession : ${params.GLDS_accession}
         Workflow : ${params.workflow}
         Nextflow Directory publishing mode: ${params.publishDir_mode}
         Swift 1S Libraries: ${params.swift_1S}
         Nextflow Error strategy: ${params.errorStrategy}
         BBDUK Adapters: ${params.adapters}
         Use GTDBTK Scratch Location: ${params.use_gtdbtk_scratch_location}
         MultiQC configuration file: ${params.multiqc_config}
         Megahit Maximum Memory: ${params.max_mem}
         Pile-up Memory: ${params.pileup_mem}
         CAT block size: ${params.block_size}

         File Suffixes:
         Filtered Reads Suffix (if single-end): ${params.filtered_suffix}
         Filtered Forward Reads Suffix: ${params.filtered_R1_suffix}
         Filtered Reverse Reads Suffix: ${params.filtered_R2_suffix}

         MAG Parameters:
         Minimum completion: ${params.min_est_comp}
         Maximum redundancy: ${params.max_est_redund}
         Maximum strain heterogeneity: ${params.max_est_strain_het}
         Use Reduced Tree: ${params.reduced_tree}
 
         Output Directories:
         Raw reads: ${params.raw_reads_dir}
         FastQC: ${params.fastqc_out_dir}
         Filtered Reads: ${params.filtered_reads_dir}
         Assembly-based Analysis: ${params.assembly_based_dir}
         Assemblies: ${params.assemblies_dir}
         Predicted Genes: ${params.genes_dir}
         Contigs Taxonomy and Annotation: ${params.annotations_and_tax_dir}
         Read mapping: ${params.mapping_dir}
         Assemblies Summary: ${params.combined_output_dir}
         Bins: ${params.bins_dir}
         Meta Assembled Genomes (MAGs): ${params.MAGs_dir}
         Read-based Analysis: ${params.read_based_dir}

         Genelab Assay Suffix: ${params.assay_suffix}
         Additional Filename Prefix: ${params.additional_filename_prefix}

         Conda Environments:
         qc: ${params.conda.qc}
         humann3: ${params.conda.humann3}
         CAT: ${params.conda.cat}
         prodigal: ${params.conda.prodigal}
         metabat: ${params.conda.metabat}
         gtdbtk: ${params.conda.gtdbtk}
         kegg decoder: ${params.conda.kegg_decoder}
         megahit: ${params.conda.megahit}
         bit: ${params.conda.bit}
         kofamscan: ${params.conda.kofamscan}
         mapping: ${params.conda.mapping}
         checkm: ${params.conda.checkm}         

         Databases:
         CAT Taxonomy: ${params.cat_taxonomy_dir}
         CAT DB sub directory: ${params.cat_db_sub_dir}
         CAT URL: ${params.database.CAT_DB_LINK}
         CAT DB: ${params.database.cat_db}
         KOFAM Scan: ${params.database.ko_db_dir}
         Metaphlan: ${params.database.metaphlan_db_dir}
         Chocophlan: ${params.database.chocophlan_dir}
         Uniref: ${params.database.uniref_dir}
         Utilities: ${params.database.utilities_dir}
         GTDBTK: ${params.database.gtdbtk_db_dir}
         """.stripIndent()

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

       chocophlanDirExists = params.database.chocophlan_dir != null
       unirefDirExists = params.database.uniref_dir != null
       metaphlanDirExists = params.database.metaphlan_db_dir != null
       utilitiesDirExists = params.database.utilities_dir != null

     // if any of the four databases  
     if(!chocophlanDirExists ||!unirefDirExists ||
       !metaphlanDirExists || !utilitiesDirExists) {
              
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
         if(params.database.ko_db_dir == null) {
             SETUP_KOFAMSCAN_DB()
             kofam_db = SETUP_KOFAMSCAN_DB.out.ko_db_dir
         }

         cat_db = params.database.cat_db
         if(params.database.cat_db == null){

            SETUP_CAT_DB(params.dataase.CAT_DB_LINK)
            cat_db = SETUP_CAT_DB.out.cat_db
         }

         gtdbtk_db_dir = params.database.gtdbtk_db_dir
         if(params.database.gtdbtk_db_dir == null){
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
