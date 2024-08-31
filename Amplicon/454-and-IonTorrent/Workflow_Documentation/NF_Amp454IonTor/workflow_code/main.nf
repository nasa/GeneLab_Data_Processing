#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

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
  println("Nextflow Amp454IonTor Consensus Pipeline: $workflow.manifest.version")
  println("USAGE:")
  println("Example 1: Submit and run jobs with slurm in singularity containers with OSD accession as input.")
  println("   > nextflow run main.nf -resume -profile slurm,singularity --GLDS_accession OSD-72 --target_region 16S --min_bbduk_len 50 --min_bbduk_avg_quality 15")
  println()
  println("Example 2: : Submit and run jobs with slurm in conda environments.")
  println("   > nextflow run main.nf -resume -profile slurm,singularity --csv_file file.csv --target_region ITS --F_primer TCCGTAGGTGAACCTGCGG --R_primer GCTGCGTTCTTCATCGATGC --min_bbduk_len 50 --min_bbduk_avg_quality 15")
  println()
  println("Example 3: Run jobs locally in conda environments and specify the path to an existing conda environment")
  println("   > nextflow run main.nf -resume -profile conda --csv_file file.csv --target_region 16S --F_primer AGAGTTTGATCCTGGCTCAG --R_primer CTGCCTCCCGTAGGAGT --min_bbduk_len 50 --min_bbduk_avg_quality 15 --conda.qc <path/to/existing/conda/environment>")
  println()
  println()
  println("Required arguments:")
  println()  
  println("""-profile [STRING] Specifies the profile to be used to run the workflow. Options are [singularity, docker, conda, slurm].
	           singularity, docker and conda will run the pipeline locally using singularity, docker, and conda, respectively.
               You can combine profiles by separating them with comma. For example, to submit and run jobs using slurm in singularity containers pass 'slurm,singularity' as argument. """)			 
  println("--csv_file  [PATH]  Required only if --GLDS_accession is not set. A 2-column input file with these headers [sample_id, read] e.g. file.csv. The sample_id column should contain unique sample ids while the read column should contain the absolute or relative path to the sample's reads.")
  println("--target_region [STRING] Specifies the amplicon target region to be analyzed. options are one of [16S, ITS]. Default: 16S")
  println()
  println("Cutadapt (trimming) parameters:")
  println("	    --F_primer [STRING] Required only if --GLDS_accession is not set. Forward primer sequence e.g. AGAGTTTGATCCTGGCTCAG")
  println("	    --R_primer [STRING] Required only if --GLDS_accession is not set. Reverse primer sequence e.g. CTGCCTCCCGTAGGAGT")
  println()
  println("BBDUK (filtering) parameters:")
  println("	    --min_bbduk_len [INT] Specifies the minimum read length threshold for bbduk. Default: 50") 
  println("	    --min_bbduk_avg_quality [INT] Specifies the minimum average quality for bbduk. Default: 15")  
  println()
  println()	
  println("Optional arguments:")  
  println()
  println("  --help  Print this help message and exit")
  println("  --debug   Show a detailed log of the parameters set by the user when the workflow runs.")
  println("  --publishDir_mode [STRING]  Specifies how nextflow handles file outputs. Options can be found here https://www.nextflow.io/docs/latest/process.html#publishdir Default: link.")
  println("  --errorStrategy [STRING] Specifies how nextflow handles errors. Options can be found here https://www.nextflow.io/docs/latest/process.html#errorstrategy. Default: terminate")
  println()
  println("File Suffixes:")
  println("      --primer_trimmed_suffix [STRING] Suffix to use for naming your primer trimmed reads. Default: _trimmed.fastq.gz")
  println("      --filtered_suffix [STRING]  Suffix to use for naming your quality filtered reads. Default: _filtered.fastq.gz")
  println()  
  println("Output directories:")
  println("      --raw_reads_dir [PATH] Specifies where processed raw reads will be published. Default: ../Raw_Sequence_Data/")
  println("      --fastqc_out_dir [PATH] Specifies where fastqc and multiqc outputs will be published. Default: ../FastQC_Outputs/")
  println("      --trimmed_reads_dir [PATH] Specifies where cutadapt trimmed reads will be published. Default: ../Trimmed_Sequence_Data/")
  println("      --filtered_reads_dir [PATH] Specifies where BBDUK filtered reads will be published.  Default: ../Filtered_Sequence_Data/")
  println()
  println("Genelab specific arguements:")
  println("      --GLDS_accession [STRING]  A Genelab GLDS or OSD accession number if the --csv_file parameter is not set. If this parameter is set, it will ignore the --csv_file parameter.")
  println("      --RawFilePattern [STRING]  If we do not want to download all files (which we often won't), we can specify a pattern here to subset the total files.")
  println("                                 For example, if we know we want to download just the fastq.gz files, we can say 'fastq.gz'. We can also provide multiple patterns")
  println("                                 as a comma-separated list. For example, If we want to download the fastq.gz files that also have 'Amplicon', and 'raw' in") 
  println("                                 their filenames, we can provide '-p fastq.gz,Amplicon,raw'. Default: null.")
  println("      --assay_suffix [STRING]  Genelabs assay suffix. Default: _GLAmpSeq.")
  println("      --output_prefix [STRING] Unique name to tag on to output files. Default: ''")
  println()
  println("Paths to existing conda environments to use otherwise a new one will be created using the yaml file in envs/.")
  println("      --conda.qc [PATH] Path to a conda environment containing fastqc, multiqc, zip and python. Default: null.")
  println("      --conda.R [PATH] Path to a conda environment containing R along with the packages decipher and biomformat installed. Default: null.")
  println("      --conda.bbmap [PATH] Path to a conda environment containing bbmap. Default: null.")
  println("      --conda.cutadapt [PATH] Path to a conda environment containing cutadapt. Default: null.")
  println("      --conda.vsearch [PATH] Path to a conda environment containing vsearch and bit. Default: null.")
  println()
  print("Advanced users can edit the nextflow.config file for more control over default settings such as container choice, number of cpus, memory per task etc.")
  exit 0
  }


/************************************************
*********** Show pipeline parameters ************
*************************************************/

if(params.debug){

log.info """
         Nextflow Amp454IonTor Consensus Pipeline: $workflow.manifest.version
         You have set the following parameters:
         Profile: ${workflow.profile}
         GLDS_accession : ${params.GLDS_accession}
         Input csv file : ${params.csv_file}
         Amplicon target region : ${params.target_region}
         Foward Primer: ${params.F_primer}
         Reverse Primer: ${params.R_primer}
         Minimum read length: ${params.min_bbduk_len}bp
         Minimum read quality: ${params.min_bbduk_avg_quality}
         Directory publishing mode: ${params.publishDir_mode}
         Nextflow Error strategy: ${params.errorStrategy}

         File Suffixes:
         Primers Trimmed Reads Suffix: ${params.primer_trimmed_suffix}
         Filtered Reads Suffix: ${params.filtered_suffix} 

         Output directories:
         Raw reads: ${params.raw_reads_dir}
         FastQC: ${params.fastqc_out_dir}
         Trimmed Reads: ${params.trimmed_reads_dir}
         Filtered Reads: ${params.filtered_reads_dir}

         Genelab parameters:
         Genelab Assay Suffix: ${params.assay_suffix}
         Output prefix: ${params.output_prefix}

         Conda environments:
         qc: ${params.conda.qc}
         R: ${params.conda.R}
         bbmap: ${params.conda.bbmap}
         cutadapt: ${params.conda.cutadapt}
         vsearch: ${params.conda.vsearch}
         """.stripIndent()
}


// Create GLDS runsheet
include { GET_RUNSHEET } from "./modules/create_runsheet.nf"

// Read quality check and filtering
include { FASTQC as RAW_FASTQC ; MULTIQC as RAW_MULTIQC  } from './modules/quality_assessment.nf'
include { CUTADAPT; COMBINE_CUTADAPT_LOGS_AND_SUMMARIZE } from './modules/quality_assessment.nf'
include { BBDUK; COMBINE_BBDUK_LOGS_AND_SUMMARIZE } from './modules/quality_assessment.nf'
include { FASTQC as FILTERED_FASTQC ; MULTIQC as FILTERED_MULTIQC  } from './modules/quality_assessment.nf'
include { pick_otus } from './modules/vsearch.nf'
include { RUN_R} from './modules/assign_taxonomy.nf'
include { ZIP_BIOM } from './modules/zip_biom.nf'


workflow {

    // Capture software versions
    software_versions_ch = Channel.empty()

    if(params.GLDS_accession){

       GET_RUNSHEET(params.GLDS_accession, params.target_region)
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

       GET_RUNSHEET.out.version | mix(software_versions_ch) | set{software_versions_ch}

   }else{

       Channel.fromPath(params.csv_file, checkIfExists: true)
           .splitCsv(header:true)
           .set{file_ch} 

   }


    file_ch.map{row -> tuple( "${row.sample_id}", [file("${row.read}", checkIfExists: true)] )}
          .set{reads_ch}

    // Read quality check and trimming
    RAW_FASTQC(reads_ch)
    raw_fastqc_files = RAW_FASTQC.out.html.flatten().collect()
    RAW_MULTIQC("raw", raw_fastqc_files)
    
    // Trim reads
    if(!params.GLDS_accession) primers_ch = Channel.value([params.F_primer, params.R_primer])
    CUTADAPT(reads_ch, primers_ch)
    trim_counts = CUTADAPT.out.trim_counts.map{ sample_id, count -> file("${count}")}.collect()
    trim_logs = CUTADAPT.out.logs.map{ sample_id, log -> file("${log}")}.collect()
    COMBINE_CUTADAPT_LOGS_AND_SUMMARIZE(trim_counts, trim_logs)
    
    // Filter reads
    BBDUK(CUTADAPT.out.reads)
    filter_counts = BBDUK.out.filter_counts.map{ sample_id, count -> file("${count}")}.collect()
    filter_logs = BBDUK.out.logs.map{ sample_id, log -> file("${log}")}.collect()
    COMBINE_BBDUK_LOGS_AND_SUMMARIZE(filter_counts, filter_logs)

    FILTERED_FASTQC(BBDUK.out.reads)
    filtered_fastqc_files = FILTERED_FASTQC.out.html.flatten().collect()
    FILTERED_MULTIQC("filtered", filtered_fastqc_files)
  
    // Pick OTUs with vsearch
    pick_otus(BBDUK.out.reads)

    // Assign taxonomy
    RUN_R(pick_otus.out.otus, pick_otus.out.counts, 
          COMBINE_CUTADAPT_LOGS_AND_SUMMARIZE.out.counts,
           COMBINE_BBDUK_LOGS_AND_SUMMARIZE.out.counts)

    // Zip biom file
    ZIP_BIOM(RUN_R.out.biom)



    // Software Version Capturing - combining all captured sofware versions
    RAW_FASTQC.out.version | mix(software_versions_ch) | set{software_versions_ch}
    RAW_MULTIQC.out.version | mix(software_versions_ch) | set{software_versions_ch}
    CUTADAPT.out.version | mix(software_versions_ch) | set{software_versions_ch}
    BBDUK.out.version | mix(software_versions_ch) | set{software_versions_ch}
    FILTERED_FASTQC.out.version | mix(software_versions_ch) | set{software_versions_ch}
    FILTERED_MULTIQC.out.version | mix(software_versions_ch) | set{software_versions_ch}
    pick_otus.out.versions | mix(software_versions_ch) | set{software_versions_ch}
    RUN_R.out.version | mix(software_versions_ch) | set{software_versions_ch}


    nf_version = "Nextflow Version: ".concat("${nextflow.version}\n<><><>\n")
    nextflow_version_ch = Channel.value(nf_version)

     //  Write software versions to file
     software_versions_ch | map { it.text + "\n<><><>\n"}
                          | unique
                          | mix(nextflow_version_ch)
                          | collectFile(name: "${params.metadata_dir}/software_versions.txt", newLine: true, cache: false)
                          | set{final_software_versions_ch}
}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Workflow completed without any error\n" : "Oops .. something went wrong" )
}
