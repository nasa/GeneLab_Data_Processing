#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

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
  println("Nextflow Amp454IonTor Consensus Pipeline: $workflow.manifest.version")
  println("USAGE:")
  println("Example 1: Submit and run jobs with slurm in singularity containers.")
  println("   > nextflow run main.nf -resume -profile slurm_sing --csv_file file.csv --target_region 16S --F_primer AGAGTTTGATCCTGGCTCAG --R_primer CTGCCTCCCGTAGGAGT --min_bbduk_len 50 --min_bbduk_avg_quality 15")
  println()
  println("Example 2: : Submit and run jobs with slurm in conda environments.")
  println("   > nextflow run main.nf -resume -profile slurm_conda --csv_file file.csv --target_region 1TS --F_primer AGAGTTTGATCCTGGCTCAG --R_primer CTGCCTCCCGTAGGAGT --min_bbduk_len 50 --min_bbduk_avg_quality 15")
  println()
  println("Example 3: Run jobs locally in conda environments and specify the path to an existing conda environment")
  println("   > nextflow run main.nf -resume -profile conda --csv_file file.csv --target_region 16S --F_primer AGAGTTTGATCCTGGCTCAG --R_primer CTGCCTCCCGTAGGAGT --min_bbduk_len 50 --min_bbduk_avg_quality 15 --conda.qc <path/to/existing/conda/environment>")
  println()
  println("Required arguments:")
  println("""-profile [STRING] What profile should be used be use to run the workflow. Options are [singularity, docker, conda, slurm_sing, slurm_conda].
	         singularity, docker and conda will run the pipelne locally using singularity, docker, and conda, respectively.
             slurm_sing and slurm_conda will submit and run jobs using slurm in singularity containers and conda environments, respectively. """)			 
  println("--csv_file  [PATH]  A 2-column input file with these headers [sample_id, read] e.g. file.csv. The sample_id column should contain unique sample ids while the read column should contain the absolute or relative path to the sample's reads.")
  println("--target_region [STRING] What the amplicon target region to be aanalyzed. options are one of [16S, 18S, ITS]. Default: 16S")
  println("Cutadapt (trimming) parameters:")
  println("	    --F_primer [STRING] Forward primer sequence e.g. AGAGTTTGATCCTGGCTCAG")
  println("	    --R_primer [STRING] Reverse primer sequence e.g. CTGCCTCCCGTAGGAGT")
  println("BBDUK (filtering) parameters:")
  println("	    --min_bbduk_len [INT] Minimum read length threshold for bbduk. Default: 50") 
  println("	    --min_bbduk_avg_quality [INT] BBduk minimum average quality. Default: 15")  
	
  println("Optional arguments:")  
  println("  --help  Print this help message and exit")
  println("  --publishDir_mode [STRING]  How should nextflow handle file outputs. Options can be found here https://www.nextflow.io/docs/latest/process.html#publishdir Default: link.")
  println("File Suffixes:")
  println("      --primer_trimmed_suffix [STRING] Suffix to use for naming your primer trimmed reads. Default: _trimmed.fastq.gz")
  println("      --filtered_suffix [STRING]  Suffix to use for naming your quality filtered reads. Default: _filtered.fastq.gz")
  println("Output directories:")
  println("      --raw_reads_dir [PATH] Where should your processed raw reads be stored. Default: Raw_Sequence_Data/")
  println("      --fastqc_out_dir [PATH] Where should fastqc and multiqc outputs be stored. Default: FastQC_Outputs/")
  println("      --trimmed_reads_dir [PATH] Where should your cutadapt trimmed reads be stored. Default: Trimmed_Sequence_Data/")
  println("      --filtered_reads_dir [PATH] Where should your BBDUK filtered reads be stored.  Default: Filtered_Sequence_Data/")
  println("Genelab specific arguements:")
  println("      --assay_suffix [STRING]  Genelabs assay suffix. Default: GLAmpSeq.")
  println("      --output_prefix [STRING] Unique name to tag on to output files. Default: ''")
  println("Paths to existing conda environments to use otherwise a new one will be created using the yaml file in envs/.")
  println("      --conda.qc [PATH] Path to a conda environment containing fastqc, multiqc, zip and python. Default: null.")
  println("      --conda.R [PATH] Path to a conda environment containing R along with the packages decipher and biomformat installed. Default: null.")
  println("      --conda.bbmap [PATH] Path to a conda environment containing bbmap. Default: null.")
  println("      --conda.cutadapt [PATH] Path to a conda environment containing cutadapt. Default: null.")
  println("      --conda.vsearch [PATH] Path to a conda environment containing vsearch and bit. Default: null.")
  print("Advanced users can edit the nextflow.config file for more control over default settings such container choice, number cpus, memory per task etc.")
  exit 0
  }

log.info """
         Nextflow Amp454IonTor Consensus Pipeline: $workflow.manifest.version
         You have set the following parameters:
         Input csv file : ${params.csv_file}
         Amplicon target region : ${params.target_region}
         Foward Primer: ${params.F_primer}
         Reverse Primer: ${params.R_primer}
         Minimum read length: ${params.min_bbduk_len}bp
         Minimum read quality: ${params.min_bbduk_avg_quality}
         Directory publishing mode: ${publishDir_mode}
         Raw reads Directory: ${params.raw_reads_dir}
         FastQC Directory: ${params.fastqc_out_dir}
         Trimmed Reads Directory: ${params.trimmed_reads_dir}
         Filtered Reads Directory: ${params.}
         Genelab Assay Suffix: ${params.}
         Output prefix: ${params.}
         Conda environments:
         qc: ${params.conda.qc}
         R: ${params.conda.R}
         bbmap: ${params.conda.bbmap}
         cutadapt: ${params.conda.cutadapt}
         vsearch: ${params.conda.vsearch}
         """.stripIndent()

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
