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
  println("Nextflow AmpIllumina Consensus Pipeline: $workflow.manifest.version")
  println("USAGE:")
  println("Example 1: Submit and run jobs with slurm in singularity containers.")
  println("   > nextflow run main.nf -resume -profile slurm_sing --csv_file PE_file.csv --target_region 16S --F_primer AGAGTTTGATCCTGGCTCAG --R_primer CTGCCTCCCGTAGGAGT")
  println()
  println("Example 2: : Submit and run jobs with slurm in conda environments.")
  println("   > nextflow run main.nf -resume -profile slurm_conda --csv_file SE_file.csv --target_region 1TS --F_primer AGAGTTTGATCCTGGCTCAG --R_primer CTGCCTCCCGTAGGAGT")
  println()
  println("Example 3: Run jobs locally in conda environments, supplying a GLDS accession, and specifying the path to an existing conda environment")
  println("   > nextflow run main.nf -resume -profile conda --GLDS_accession GLDS-487 --target_region 16S --F_primer AGAGTTTGATCCTGGCTCAG --R_primer CTGCCTCCCGTAGGAGT --conda.qc <path/to/existing/conda/environment>")
  println()
  println("Required arguments:")
  println("""-profile [STRING] What profile should be used to run the workflow. Options are [singularity, docker, conda, slurm_sing, slurm_conda].
	         singularity, docker and conda will run the pipelne locally using singularity, docker, and conda, respectively.
             slurm_sing and slurm_conda will submit and run jobs using slurm in singularity containers and conda environments, respectively. """)			 
  println("--csv_file  [PATH] A 3-column (single-end) or 4-column (paired-end) input file (sample_id, forward, [reverse,] paired). Mandatory if a GLDS accession is not provided.")
  println(" Please see the files: SE_file.csv and PE_file.csv for single-end and paired-end examples, respectively.")
  println(" The sample_id column should contain unique sample ids.")
  println(" The forward and reverse columns should contain the absolute or relative path to the sample's forward and reverse reads.")
  println(" The paired column should be true for paired-end or anything else for single-end reads.")
  println("--target_region [STRING] What is the amplicon target region to be analyzed. Options are one of [16S, 18S, ITS]. Default: 16S.")
  println("--trim_primers [BOOLEAN] Should primers be trimmed? true or false. Default: true.") 
  println("PLEASE NOTE: This workflow assumes that all your raw reads end with the same suffix. If they don't please modify your filenames to have the same suffix as shown below.")
  println("--raw_R1_suffix [STRING] Raw forward reads suffix (region following the unique part of the sample names). e.g. _R1_raw.fastq.gz.") 
  println("--raw_R2_suffix [STRING] Raw reverse reads suffix (region following the unique part of the sample names). e.g. _R2_raw.fastq.gz.") 
  println()
  println("Cutadapt (trimming) parameters:")
  println("	    --F_primer [STRING] Forward primer sequence e.g. AGAGTTTGATCCTGGCTCAG. Default: emptry string.")
  println("	    --R_primer [STRING] Reverse primer sequence e.g. CTGCCTCCCGTAGGAGT. Default: emptry string.")
  println("	    --min_cutadapt_len [INTEGER] What should be the minimum read length after quality trimming with cutadapt. Default: 130.")
  println("	    --primers_linked [STRING] Are the primers linked?. https://cutadapt.readthedocs.io/en/stable/recipes.html#trimming-amplicon-primers-from-paired-end-reads. Default: TRUE. ")
  println("	    --discard_untrimmed [STRING] Should untrimmed reads be discarded? Any supplied string except TRUE will not discard them. Default: TRUE.")
  println()	
  println("Optional arguments:")  
  println("  --help  Print this help message and exit.")
  println("  --publishDir_mode [STRING]  How should nextflow publish file outputs. Options can be found here https://www.nextflow.io/docs/latest/process.html#publishdir. Default: link.")
  println("  --errorStrategy [STRING] How should nextflow handle errors. Options can be found here https://www.nextflow.io/docs/latest/process.html#errorstrategy. Default: terminate")
  println("  --enable_visualizations [BOOLEAN] Should ASV plots be made? true or false. if true supply a path to the ruhnsheet for plotting to the --runsheet option. Default: false.")
  println("  --runsheet [PATH] A 4-column file with these exact headers [Sample Name, read1_path, raw_R1_suffix, groups] for plotting. Only relevant if --enable_visualizations is true. Default: null.") 
  println("  --multiqc_config [PATH] Path to a custome multiqc config file. Default: config/multiqc.config.")
  println()
  println("Dada2 parameters passed to filterAndTrim() function:")
  println("	    --left_trunc [INTEGER] truncate the sequences to the left by this number of bases. Default: 0.") 
  println("	    --right_trunc [INTEGER] truncate the sequences to the right by this number of bases. Default: 0.") 
  println("	    --left_maxEE [INTEGER] Maximum allowed errors to the left. Default: 1.")
  println("	    --right_maxEE [INTEGER] Maximum allowed errors to the right. Default: 1.")
  println("	    --concatenate_reads_only [STRING] Concatenate only with dada2 instead of merging paired reads if TRUE.")
  println("      This is typically used with primers like 515-926, that captured 18S fragments that are typically too long to merge.")
  println("      Note that 16S and 18S should have been separated already prior to running this workflow. This should likely be left as FALSE for any option other than 18S above.") 	    
  println("	     Values are TRUE or FALSE Default: FALSE.")
  println()
  println("File Suffixes:")
  println("      --primer_trimmed_R1_suffix [STRING] Suffix to use for naming your primer trimmed forward reads. Default: _R1_trimmed.fastq.gz.")
  println("      --primer_trimmed_R2_suffix [STRING] Suffix to use for naming your primer trimmed reverse reads. Default: _R2_trimmed.fastq.gz.")  
  println("      --filtered_R1_suffix [STRING]  Suffix to use for naming your quality filtered forward reads. Default: _R1_filtered.fastq.gz.")
  println("      --filtered_R2_suffix [STRING]  Suffix to use for naming your quality filtered reverse reads. Default: _R2_filtered.fastq.gz.")
  println()
  println("Output directories:")
  println("      --raw_reads_dir [PATH] Where should the fastqc report of the raw reads be stored. Default: Raw_Sequence_Data/.")
  println("      --fastqc_out_dir [PATH] Where should multiqc outputs be stored. Default: workflow_output/FastQC_Outputs/.")
  println("      --trimmed_reads_dir [PATH] Where should your cutadapt trimmed reads be stored. Default: workflow_output/Trimmed_Sequence_Data/.")
  println("      --filtered_reads_dir [PATH] Where should your filtered reads be stored.  Default: workflow_output/Filtered_Sequence_Data/.")
  println("      --info_out_dir [PATH] Where should output metadata be stored. Default: workflow_output/Metadata/.")
  println("      --plots_dir [PATH] Where should your plots be stored if visualization is enabled. Default: workflow_output/Final_Outputs/Plots/.")
  println("      --final_outputs_dir [PATH] Where should most outputs and summary reports be stored.  Default: workflow_output/Final_Outputs/.")
  println()
  println("Genelab specific arguements:")
  println("      --GLDS_accession [STRING]  A Genelab accession number if the --csv_file parameter is not set. If this parameter is set, it will ignore the --csv_file parameter.")
  println("      --assay_suffix [STRING]  Genelabs assay suffix. Default: GLAmpSeq.")
  println("      --output_prefix [STRING] Unique name to tag onto output files. Default: empty string.")
  println()
  println("Paths to existing conda environments to use otherwise a new one will be created using the yaml file in envs/.")
  println("      --conda.qc [PATH] Path to a conda environment containing fastqc, multiqc, zip and python. Default: null.")
  println("      --conda.R [PATH] Path to a conda environment containing R along with the packages decipher and biomformat installed. Default: null.")
  println("      --conda.genelab  [PATH] Path to a conda environment containing genlab-utils. Default: null.")
  println("      --conda.cutadapt [PATH] Path to a conda environment containing cutadapt. Default: null.")
  println("      --conda.R_visualizations [PATH] Path to a conda environment containing R packages required for plotting. Default: null.")
  print("Advanced users can edit the nextflow.config file for more control over default settings such container choice, number of cpus, memory per task etc.")
  exit 0
  }

log.info """
         Nextflow AmpIllumina Consensus Pipeline: $workflow.manifest.version
         
         You have set the following parameters:
         Input csv file : ${params.csv_file}
         GLDS_accession : ${params.GLDS_accession}
         Amplicon target region : ${params.target_region}
         Nextflow Directory publishing mode: ${params.publishDir_mode}
         Trim Primers: ${params.trim_primers}
         Nextflow Error strategy: ${params.errorStrategy}
         Enable visualization: ${params.enable_visualizations}
         Runsheet For plotting: ${params.runsheet}
         MultiQC configuration file: ${params.multiqc_config}

         File Suffixes:
         Raw Forward Reads Suffix: ${params.raw_R1_suffix}
         Raw Reverse Reads Suffix: ${params.raw_R2_suffix}
         Trimmed Forward Reads Suffix: ${params.primer_trimmed_R1_suffix}
         Trimmed Reverse Reads Suffix: ${params.primer_trimmed_R2_suffix}
         Filtered Forward Reads Suffix: ${params.filtered_R1_suffix}
         Filtered Reverse Reads Suffix: ${params.filtered_R2_suffix}

         Cutadapt Parameters:
         Forward Primer: ${params.F_primer}
         Reverse Primer: ${params.R_primer}
         Minimum Trimmed Reads length: ${params.min_cutadapt_len}
         Primers Are linked: ${params.primers_linked}
         Discard Untrimmed Reads: ${params.discard_untrimmed}

 
         Dada2 Parameters:
         Truncate left: ${params.left_trunc}bp
         Truncate right: ${params.right_trunc}bp
         Max error left: ${params.left_maxEE}
         Max error right: ${params.right_maxEE}
         Concatenate Reads: ${params.concatenate_reads_only}
 
         Output Directories:
         Raw reads: ${params.raw_reads_dir}
         FastQC: ${params.fastqc_out_dir}
         Trimmed Reads: ${params.trimmed_reads_dir}
         Filtered Reads: ${params.filtered_reads_dir}
         Metadata: ${params.info_out_dir}
         Plots: ${params.plots_dir}
         Reports: ${params.final_outputs_dir}

         Genelab Assay Suffix: ${params.assay_suffix}
         Output Prefix: ${params.output_prefix}

         Conda Environments:
         qc: ${params.conda.qc}
         R: ${params.conda.R}
         genelab: ${params.conda.genelab}
         cutadapt: ${params.conda.cutadapt}
         R_visualizations: ${params.conda.R_visualizations}
         """.stripIndent()


// Create GLDS runsheet
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
