#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// color defs
c_back_bright_red = "\u001b[41;1m";
c_bright_green = "\u001b[32;1m";
c_blue = "\033[0;34m";
c_reset = "\033[0m";

params.help = false
params.debug = false
/**************************************************
* HELP MENU  **************************************
**************************************************/
if (params.help) {
  println()
  println("Nextflow Metagenomics Illumina Consensus Pipeline: $workflow.manifest.version")
  println("USAGE:")
  println("Example 1: Submit and run jobs with slurm in singularity containers.")
  println("   > nextflow run main.nf -resume -profile slurm,singularity --input_file PE_file.csv")
  println()
  println("Example 2: : Submit and run jobs with slurm in conda environments.")
  println("   > nextflow run main.nf -resume -profile slurm,conda --input_file SE_file.csv")
  println()
  println("Example 3: Run jobs locally in conda environments, supply a GLDS accession, and specify the path to an existing conda environment.")
  println("   > nextflow run main.nf -resume -profile conda --accession OSD-574 --conda.qc <path/to/existing/conda/environment>")
  println()
  println("Required arguments:")
  println("""-profile [STRING] Specifies the profile to be used to run the workflow. Options are [slurm, singularity, docker, and  conda].
	                    singularity, docker and conda will run the pipeline locally using singularity, docker, and conda, respectively.
                      To combine profiles, separate two or more profiles with comma. For example, to combine slurm and singularity profiles, pass 'slurm,singularity' as argument. """)			 
  println("--input_file  [PATH] A 3-column (single-end) or 4-column (paired-end) csv input file (sample_id, forward, [reverse,] paired). Required only if a GLDS accession is not provided. Default : null")
  println("   Please see the files: SE_file.csv and PE_file.csv for single-end and paired-end examples, respectively.")
  println("   The sample_id column should contain unique sample ids.")
  println("   The forward and reverse columns should contain the absolute or relative path to the sample's forward and reverse reads.")
  println("   The paired column should be true for paired-end or anything else for single-end reads.")
  println()
  println("Optional arguments:")  
  println("  --help  Print this help message and exit")
  println("  --workflow [STRING] Specifies that workflow to be run. Options are one of [read-based, assembly-based, both]. Default: both.")
  println("  --publishDir_mode [STRING]  Specifies how nextflow handles output file publishing. Options can be found here https://www.nextflow.io/docs/latest/process.html#publishdir Default: link.")
  println("  --errorStrategy [STRING] Specifies how nextflow handles errors. Options can be found here https://www.nextflow.io/docs/latest/process.html#errorstrategy. Default: ignore")
  println("  --swift_1S [BOOLEAN] Setting for trimming recommended when working with Swift 1S libraries.")
  println("    adds `swift=t` setting to bbduk quality trimming/filtering command. For info on this, see example, ")
  println("    https://swiftbiosci.com/wp-content/uploads/2019/03/16-0853-Tail-Trim-Final-442019.pdf.")
  println("    Set to true if data was generated with Swift 1S library prep. Default: false.")
  println("  --adapters [PATH] Path to BBtools adapters for reads filtering. Default: config/bbtools_adapters.fa.") 
  println("  --multiqc_config [PATH] Path to a custom multiqc config file. Default: config/multiqc.config.")
  println("  --use_gtdbtk_scratch_location [BOOLEAN] Should a scratch location be used to store GTDBTK temp files? true or false.")
  println("    Scratch directory for gtdb-tk, if wanting to use disk space instead of RAM, can be memory intensive;")
  println("    see https://ecogenomics.github.io/GTDBTk/faq.html#gtdb-tk-reaches-the-memory-limit-pplacer-crashes")
  println("    leave empty if wanting to use memory, the default, put in quotes the path to a directory that")
  println("    already exists if wanting to use disk space. Default: false.")
  println()
  println("MAG parameters: MAG filtering cutoffs based on checkm quality assessments (in percent); see https://github.com/Ecogenomics/CheckM/wiki/Reported-Statistics.")
  println("	 --min_est_comp [INT] Minimum estimated completion. Default: 90.") 
  println("	 --max_est_redund [INT] Minimum estimated redundancy. Default: 10.") 
  println("	 --max_est_strain_het [INT] Minimum estimated strain heterogeneity. Default: 50.")
  println("	 --reduced_tree [STRING] reduced_tree option for checkm, limits the RAM usage to 16GB; https://github.com/Ecogenomics/CheckM/wiki/Genome-Quality-Commands#tree.")
  println("    'True' for yes, anything else will be considered 'False' and the default full tree will be used. Default: 'True'. ")
  println("	 --max_mem [INT] Maximum memory allowed, passed to megahit assembler. Can be set either by proportion of available on system, e.g. 0.5")
  println("    or by absolute value in bytes, e.g. 100e9 would be 100 GB. Default: 100e9.")
  println()
  println("	 --pileup_mem [STRING] pileup.sh paramater for calculating contig coverage and depth. Memory used by bbmap's pileup.sh (within the GET_COV_AND_DET process). ")
  println("	   passed as the -Xmx parameter, 20g means 20 gigs of RAM, 20m means 20 megabytes.")
  println("	   5g should be sufficient for most assemblies, but if that rule is failing, this may need to be increased.Default: '5g' .")
  println("	 --block_size [int] Block size variable for CAT/diamond, lower value means less RAM usage; see https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options#memory--performance-options. Default: 4.")
  println()
  println("File Suffixes:")
  println("      --filtered_suffix [STRING]  Specifies the suffix for naming quality filtered reads. Only applicable when input reads are single-end. Default: _filtered.fastq.gz.")  
  println("      --filtered_R1_suffix [STRING]  Specifies the suffix for naming quality filtered forward reads. Default: _R1_filtered.fastq.gz.")
  println("      --filtered_R2_suffix [STRING]  Specifies the suffix for naming quality filtered reverse reads. Default: _R2_filtered.fastq.gz.")
  println()
  println("Output directories:")
  println("      --raw_reads_dir [PATH] Specifies where the fastqc report of the raw reads will be published. Default: ../Raw_Sequence_Data/.")
  println("      --fastqc_out_dir [PATH] Specifies where multiqc outputs will be published. Default: ../FastQC_Outputs/.")
  println("      --filtered_reads_dir [PATH] Specifies where filtered reads will be published.  Default: ../Filtered_Sequence_Data/.")
  println("      --assembly_based_dir [PATH] Specifies where the results of assembly-based analysis will be published. Default: ../Assembly-based_Processing/.")
  println("      --assemblies_dir [PATH] Specifies where assemblies will be published. Default: ../Assembly-based_Processing/assemblies/.")
  println("      --genes_dir [PATH] Specifies where predicted genes from the assemblies will be published. Default: ../Assembly-based_Processing/predicted-genes/.")
  println("      --annotations_and_tax_dir [PATH] Contig taxonomy and annotation directory.  Default: ../Assembly-based_Processing/annotations-and-taxonomy/.")
  println("      --mapping_dir [PATH] Read mapping to assembly directory.  Default: ../Assembly-based_Processing/read-mapping/.")
  println("      --combined_output_dir [PATH] Assembly summaries and reports across samples directory.  Default: ../Assembly-based_Processing/combined-outputs/.")
  println("      --bins_dir [PATH] Assembly bins directory.  Default: ../Assembly-based_Processing/bins/.")
  println("      --MAGs_dir [PATH] Meta assembled genomes (MAGs) directory.  Default: ../Assembly-based_Processing/MAGs/.")
  println("      --read_based_dir [PATH] Read-based analysis outputs directory.  Default: ../Read-based_Processing/.")
  println()
  println("Genelab specific arguements:")
  println("      --accession [STRING]  A Genelab accession number if the --input_file parameter is not set. If this parameter is set, it will ignore the --input_file parameter.")
  println("      --RawFilePattern [STRING]  If we do not want to download all files (which we often won't), we can specify a pattern here to subset the total files.")
  println("                                 For example, if we know we want to download just the fastq.gz files, we can say 'fastq.gz'. We can also provide multiple patterns")
  println("                                 as a comma-separated list. For example, If we want to download the fastq.gz files that also have 'NxtaFlex', 'metagenomics', and 'raw' in") 
  println("                                 their filenames, we can provide '-p fastq.gz,NxtaFlex,metagenomics,raw'. Default: null.")
  println("      --assay_suffix [STRING]  Genelab's assay suffix. Default: _GLmetagenomics.")
  println("      --additional_filename_prefix [STRING] additional prefix to add to output files that describe more than one sample (to make them unique compared to other datasets).")
  println("      include separator at end if adding one, e.g. Swift1S_ if wanted. Default: empty string .")
  println()
  println("Paths to existing databases and database links.")
  println("        --DB_ROOT [PATH]   FULL PATH to root directory where the databases will be downloaded if they don't exist.") 
  println("                  Relative paths such as '~/' and '../' will fail, please don't use them. Default: ../Reference_DBs/ ")
  println("CAT database directory strings:")
  println("    The strings below will be added to the end of the --database.cat_db path arguement provided below.")
  println("         --cat_taxonomy_dir [PATH] CAT taxonomy database directory. Default: 2021-01-07_taxonomy/.")
  println("         --cat_db_sub_dir [PATH] CAT database sub directory. Default: 2021-01-07_CAT_database/.")
  println("         --database.CAT_DB_LINK [URL] CAT database online download link. Default: https://tbb.bio.uu.nl/bastiaan/CAT_prepare/CAT_prepare_20210107.tar.gz.")
  println("CAT database ")
  println("         --database.cat_db [PATH] Path to CAT database. Example, /path/to/Reference_DBs/CAT_prepare_20210107/. Default: null.")
  println("Humann database:")
  println("      --database.metaphlan_db_dir [PATH] Path to metaphlan database. Example, /path/to/Reference_DBs/metaphlan4-db/. Default: null.")
  println("      --database.chocophlan_dir [PATH] Path to Humann's chocophlan nucleotide database. Example, /path/to/Reference_DBs/humann3-db/chocophlan/. Default: null.")
  println("      --database.uniref_dir [PATH] Path to Humann's Uniref protein database. Example, /path/to/Reference_DBs/humann3-db/uniref/. Default: null.")
  println("      --database.utilities_dir [PATH] Path to Humann's untilities database. Example, /path/to/Reference_DBs/humann3-db/utility_mapping/.  Default: null.")
  println("GTDBTK database:")
  println("      --database.GTDBTK_LINK [URL] GTDBTK database online download link. Default: https://data.gtdb.ecogenomic.org/releases/release220/220.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r220_data.tar.gz.")
  println("      --database.gtdbtk_db_dir  [PATH] Path to GTDBTK database. Example, /path/Reference_DBs/GTDB-tk-ref-db/. Default: null.")
  println("kofam scan database database:")
  println("      --database.ko_db_dir  [PATH] Path to kofam scan database. Example, /path/to/Reference_DBs/kofamscan_db/. Default: null.")
  println()
  println("Paths to existing conda environments to use, otherwise, new ones will be created using the yaml files in envs/.")
  println("      --conda.qc [PATH] Path to a conda environment containing fastqc, multiqc, zip and python. Default: null.")
  println("      --conda.humann3 [PATH] Path to a conda environment with humann3 installed. Default: null.")
  println("      --conda.cat  [PATH] Path to a conda environment containing CAT (Contig annotation tool). Default: null.")
  println("      --conda.prodigal [PATH] Path to a conda environment with prodigal installed. Default: null.")
  println("      --conda.metabat [PATH] Path to a conda environment containing metabat. Default: null.")
  println("      --conda.gtdbtk [PATH] Path to a conda environment containing gtdbtk. Default: null.")
  println("      --conda.kegg_decoder [PATH] Path to a conda environment with kegg_decoder installed. Default: null.")
  println("      --conda.megahit  [PATH] Path to a conda environment containing megahit. Default: null.")
  println("      --conda.bit [PATH] Path to a conda environment with bit installed. Default: null.")
  println("      --conda.kofamscan [PATH] Path to a conda environment containing KOFAM SCAN. Default: null.")
  println("      --conda.mapping [PATH] Path to a conda environment with bowtie and samtools installed. Default: null.")
  println("      --conda.checkm [PATH] Path to a conda environment with checkm installed. Default: null.")
  println()
  print("Advanced users can edit the nextflow.config file for more control over default settings such container choice, number of cpus, memory per task etc.")
  exit 0
  }

/************************************************
*********** Show pipeline parameters ************
*************************************************/

if (params.debug) {
log.info """
         Nextflow Metagenomics Illumina Consensus Pipeline: $workflow.manifest.version
         
         You have set the following parameters:
         Profile: ${workflow.profile} 
         Input csv file : ${params.input_file}
         GLDS or OSD Accession : ${params.accession}
         GLDS Raw File Pattern: ${params.RawFilePattern}         
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
         GTDBTK URL: ${params.database.GTDBTK_LINK}
         GTDBTK DB: ${params.database.gtdbtk_db_dir}
         """.stripIndent()
}

// Create GLDS runsheet
include { GET_RUNSHEET } from "./modules/create_runsheet.nf"

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

         software_versions_ch = Channel.empty()
         make_humann_db.out.versions | mix(software_versions_ch) | set{software_versions_ch}
         read_based.out.versions | mix(software_versions_ch) | set{software_versions_ch}

       }else{

         read_based(filtered_ch, 
                    params.database.chocophlan_dir,
                    params.database.uniref_dir,
                    params.database.metaphlan_db_dir,
                    params.database.utilities_dir)

         software_versions_ch =  read_based.out.versions
      }

    emit:
        versions =  software_versions_ch

}

// Workflow to perform assembly-based analysis
workflow run_assembly_based_analysis {

    take:
        file_ch
        filtered_ch


    main:
        software_versions_ch = Channel.empty()

        kofam_db = params.database.ko_db_dir
        cat_db = params.database.cat_db
        gtdbtk_db_dir = params.database.gtdbtk_db_dir

        // Run assembly based workflow 
        assembly_based(file_ch, filtered_ch, kofam_db, 
                        cat_db, gtdbtk_db_dir, params.use_gtdbtk_scratch_location)


        assembly_based.out.versions | mix(software_versions_ch) | set{software_versions_ch}


    emit:
        versions =  software_versions_ch

}

// A function to delete white spaces from an input string and covert it to lower case 
def deleteWS(string){

    return string.replaceAll(/\s+/, '').toLowerCase()

}

// Main workflow
workflow {

    // Sanity check : Test input requirement
    if (!params.accession &&  !params.input_file){
     
       error("""
              Please supply either an accession (OSD or Genelab number) or an input CSV file
              by passing either to the --accession or --input_file parameter, respectively.
              """)
    } 
        
     // Software Version Capturing - runsheet
     software_versions_ch = Channel.empty()
     // Parse file input
       if(params.accession){

       GET_RUNSHEET(params.accession)
       GET_RUNSHEET.out.input_file
           .splitCsv(header:true)
           .set{file_ch}

       GET_RUNSHEET.out.version | mix(software_versions_ch) | set{software_versions_ch}
      }else{
 
       Channel.fromPath(params.input_file, checkIfExists: true)
           .splitCsv(header:true)
           .set{file_ch}
      }


    file_ch.map{
                row -> deleteWS(row.paired) == 'true'  ? tuple( "${row.sample_id}", [file("${row.forward}", checkIfExists: true), file("${row.reverse}", checkIfExists: true)], deleteWS(row.paired)) : 
                                                         tuple( "${row.sample_id}", [file("${row.forward}", checkIfExists: true)], deleteWS(row.paired))
                }.set{reads_ch}



    // Quality check and trim the input reads
    raw_qc(Channel.of("raw"), params.multiqc_config,reads_ch)
    BBDUK(reads_ch, params.adapters)
    filtered_ch = BBDUK.out.reads
    filtered_qc(Channel.of("filtered"), params.multiqc_config, filtered_ch)

    // Quality check software capturing
    raw_qc.out.versions | mix(software_versions_ch) | set{software_versions_ch}
    BBDUK.out.version | mix(software_versions_ch) | set{software_versions_ch}
    filtered_qc.out.versions | mix(software_versions_ch) | set{software_versions_ch}

    // Run the analysis based on selection i.e, read-based, assembly-based or both
    // it will run both by default
    if(params.workflow == 'read-based'){

          run_read_based_analysis(filtered_ch)
          run_read_based_analysis.out.versions | mix(software_versions_ch) | set{software_versions_ch}
          
    }else if(params.workflow == 'assembly-based') {

          run_assembly_based_analysis(file_ch,filtered_ch)
          run_assembly_based_analysis.out.versions | mix(software_versions_ch) | set{software_versions_ch}

    }else{

          run_read_based_analysis(filtered_ch)
          run_assembly_based_analysis(file_ch, filtered_ch)

          run_read_based_analysis.out.versions | mix(software_versions_ch) | set{software_versions_ch}
          run_assembly_based_analysis.out.versions | mix(software_versions_ch) | set{software_versions_ch}
    }


     // Software Version Capturing - combining all captured sofware versions
     nf_version = "Nextflow Version ".concat("${nextflow.version}")
     nextflow_version_ch = Channel.value(nf_version)

     //  Write software versions to file
     software_versions_ch | map { it.text.strip() }
                          | unique
                          | mix(nextflow_version_ch)
                          | collectFile(name: "${params.metadata_dir}/software_versions.txt", newLine: true, cache: false)
                          | set{final_software_versions_ch}

}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Workflow completed without any error\n" : "Oops .. something went wrong" )
}
