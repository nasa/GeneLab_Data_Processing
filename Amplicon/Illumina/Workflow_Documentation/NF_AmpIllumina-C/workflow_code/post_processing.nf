#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Terminal text color definitions
c_back_bright_red = "\u001b[41;1m"
c_bright_green    = "\u001b[32;1m"
c_blue            = "\033[0;34m"
c_reset           = "\033[0m"

params.help = false


/**************************************************
* HELP MENU  **************************************
**************************************************/
if(params.help){

  println()
  println("GeneLab Post Processing Pipeline: $workflow.manifest.version")
  println("USAGE:")
  println("Example: Submit and run jobs with slurm in singularity containers.")
  println("   > nextflow -C post_processing.config run post_processing.nf -resume -profile slurm,singularity")
  println()
  println("Required Parameters:")
  println("""-profile [STRING] Specifies the profile to be used to run the workflow. Options are [slurm, singularity, docker, and  conda].
	                 singularity, docker and conda will run the workflow locally using singularity, docker, and conda, respectively.
                      To combine profiles, separate two or more profiles with a comma. 
                      For example, to combine slurm and singularity profiles, pass 'slurm,singularity' as argument. """)	
  println("  --publishDir_mode [STRING]  Specifies how nextflow handles output file publishing. Options can be found here https://www.nextflow.io/docs/latest/process.html#publishdir Default: link.")
  println("  --GLDS_accession [STRING]  A Genelab GLDS accession number. Example GLDS-487. Default: null")
  println("  --OSD_accession [STRING]  A Genelab OSD accession number. Example OSD-487. Default: null")
  println("  --name [STRING] The analyst's full name. E.g. 'FirstName A. LastName'.  Default: FirstName A. LastName")
  println("  --email [STRING] The analyst's email address. E.g. 'mail@nasa.gov'.  Default: mail@nasa.gov")
  println("  --assay_suffix [STRING]  Genelab's assay suffix. Default: _GLAmpSeq.")
  println("  --output_prefix [STRING] Unique name to tag onto output files. Default: empty string.")
  println("  --V_V_guidelines_link [URL] Genelab metagenomics data validation and verification guidelines link. Default: https://genelab-tools.arc.nasa.gov/confluence/pages/viewpage.action?pageId=2428598.")
  println("  --target_files [STRING] A comma separated list of target files and/or directories to find in processing_info.zip. Default: main.nf,nextflow.config,unique-sample-IDs.txt,envs/,bin/,config/,modules/.")
  println()
  println("File Suffixes:")
  println("      --raw_suffix [STRING]  Suffix used for the raw reads during processing. Only applicable when input reads are single-end. Default: _raw.fastq.gz.")  
  println("      --raw_R1_suffix [STRING]  Suffix used for the raw forward reads during processing. Default: _R1_raw.fastq.gz.")
  println("      --raw_R2_suffix [STRING]  Suffix used for the raw reverse reads during processing. Default: _R2_raw.fastq.gz.")
  println("      --primer_trimmed_suffix [STRING]  Suffix used for primer trimmed reads during processing. Only applicable when input reads are single-end. Default: _trimmed.fastq.gz.")
  println("      --primer_trimmed_R1_suffix [STRING]  Suffix to use for primer trimmed reads during processing. Default: _R1_trimmed.fastq.gz.")
  println("      --primer_trimmed_R2_suffix [STRING]  Suffix to use for primer trimmed reverse reads during processing. Default: _R2_trimmed.fastq.gz.")
  println("      --filtered_suffix [STRING]  Suffix used for quality filtered reads during processing. Only applicable when input reads are single-end. Default: _filtered.fastq.gz.")
  println("      --filtered_R1_suffix [STRING]  Suffix to use for quality filtered forward reads during processing. Default: _R1_filtered.fastq.gz.")
  println("      --filtered_R2_suffix [STRING]  Suffix to use for quality filtered reverse reads during processing. Default: _R2_filtered.fastq.gz.")
  println()
  println("Extra parameters to scripts:")
  println("      --readme_extra [STRING] Extra parameters and arguments to GL-gen-processed-amplicon-data-readme command. Run 'GL-gen-processed-amplicon-readme --help' for extra parameters that can be set. Example '--raw-reads-dir  ../Raw_Sequence_Data/'. Default: empty string")
  println("      --validation_extra [STRING] Extra parameters and arguments to GL-validate-processed-amplicon-data command. Run 'GL-validate-processed-amplicon-data --help' for extra parameters that can be set. Example '--single-ended --R1-used-as-single-ended-data --skip_raw_multiqc'. Default: empty string ")
  println("      --file_association_extra [STRING] Extra parameters and arguments to GL-gen-amplicon-file-associations-table command. Run 'GL-gen-amplicon-file-associations-table --help' for extra parameters that can be set. Example '--single-ended --R1-used-as-single-ended-data'. Default: empty string ")
  println()
  println("Files:")
  println("    --main  [PATH] The main workflow script used for processing. Default: ./main.nf")
  println("    --nextflow_config  [PATH] The main workflow configuration file used for processing. Default: ./nextflow.config")
  println("    --samples  [PATH] A single column file with sample ids on each line generated after running the processing pipeline. Default: ./unique-sample-IDs.txt")
  println("You only need to supply one of --assay_table or --isa_zip. If you supply both it will use only the assay_table.")
  println("""  --assay_table  [PATH] GLDS assay table generated after running the processing pipeline with accession number as input.
                               Example, ../GeneLab/a_GLDS-487_amplicon-sequencing_16s_illumina-1.txt. Default: null""")
  println("""  --isa_zip  [PATH] GLDS ISA ZIP file generated after running the processing pipeline with accession number as input.
                               Example, ../GeneLab/OSD-487_metadata_GLDS-487-ISA.zip. Default: null""")
  println("    --runsheet  [PATH] A 3-column (single-end) or 4-column (paired-end) input file (sample_id, forward, [reverse,] paired) used to run the processing pipeline. This is the value set to the paremater --csv_file when run the processing pipeline with a csv file as input otherwise it is the GLfile.csv in the GeneLab directory if --GLDS_accession was used as input. Example '../GeneLab/GLfile.csv'.  Default: null")
  println("    --software_versions  [PATH] A file generated after running the processing pipeline listing the software versions used. Default: ../Metadata/software_versions.txt")
  println()
  println("Directories:")
  println("    --config_dir  [PATH] A directory containing configuration files used in the processing pipeline. Only relevent in Metagenomics and AmpIllumina workflows. Default: ./config/")
  println("    --bin  [PATH] A directory containing scripts used by nextflow. Default: ./bin/")
  println("    --envs  [PATH] A directory containing conda yaml files. Default: ./envs/")
  println("    --modules  [PATH] A directory containing nextflow module scripts. Default: ./modules/")
  println("    --Raw_Sequence_Data [PATH] A directory containing raw sequence and raw sequence outputs. Default: ../Raw_Sequence_Data/")
  println("    --FastQC_Outputs [PATH] A directory containing fastqc and multiqc zip reports. Default: ../workflow_outputs/FastQC_Outputs/")
  println("    --Trimmed_Sequence_Data  [PATH] A directory containing the outputs of read trimming after running the processing pipeline. Default: ../workflow_outputs/Trimmed_Sequence_Data/")
  println("    --Filtered_Sequence_Data  [PATH] A directory containing the outputs of read filtering after running the processing pipeline. Default: ../Filtered_Sequence_Data/")  
  println("    --Final_Outputs  [PATH] A directory containing the outputs of assembly based processing after running the processing pipeline. Default: ../workflow_outputs/Final_Outputs/")
  println("    --Output_dir  [PATH] Specifies the directory where outputs of this post-processing workflow will be published. Default: ../Post_Processing/")
  println()
  println("Optional arguments:")  
  println("    --help  Print this help message and exit")
  println("    --debug [BOOLEAN] Set to true if you'd like to see the values of your set parameters printed to the terminal. Default: false.")
  println()
  println("Paths to existing conda environments to use otherwise a new one will be created using the yaml file in envs/.")
  println("      --conda_dp_tools [PATH] Path to a conda environment containing dp_tools. Default: null.")
  exit 0
}


/************************************************
*********** Show pipeline parameters ************
*************************************************/
if(params.debug){

log.info """${c_blue}
         GeneLab Post Processing Pipeline: $workflow.manifest.version
         
         You have set the following parameters:
         Profile: ${workflow.profile} 
         Analyst's Name : ${params.name}
         Analyst's Email : ${params.email}
         GLDS Accession : ${params.GLDS_accession}
         OSD Accession : ${params.OSD_accession}
         Assay Suffix: ${params.assay_suffix}
         Output Prefix: ${params.output_prefix}
         V & V Link: ${params.V_V_guidelines_link}
         Target Files: ${params.target_files}
         README File: ${params.readme}
         Processing Info Zip File: ${params.processing_zip_file} 
         Nextflow Directory publishing mode: ${params.publishDir_mode}
         
         Suffixes:
         Raw Suffix: ${params.raw_suffix}
         Raw R1 suffix: ${params.raw_R1_suffix}
         Raw R2 suffix: ${params.raw_R2_suffix} 
         Primer Trimmed Suffix: ${params.primer_trimmed_suffix}
         Primer Trimmed R1 suffix: ${params.primer_trimmed_R1_suffix}
         Primer Trimmed R2 suffix: ${params.primer_trimmed_R2_suffix}
         Filtered Suffix: ${params.filtered_suffix}
         Filtered R1 suffix: ${params.filtered_R1_suffix}
         Filtered R2 suffix: ${params.filtered_R2_suffix}

         Extra scripts parameters:
         Readme Script Extra: ${params.readme_extra}
         Validation Script Extra : ${params.validation_extra}
         File Association Script Extra: ${params.file_association_extra}

         Files:
         Main Workflow Script: ${params.main}
         Nextflow Config File: ${params.nextflow_config}
         Samples: ${params.samples}
         Assay Table: ${params.assay_table}
         ISA Zip: ${params.isa_zip}
         Input Runsheet: ${params.runsheet}
         Software Versions: ${params.software_versions}

         Directories:
         Config: ${params.config_dir}
         Bin: ${params.bin}
         Conda Environments: ${params.envs}
         Modules: ${params.modules}
         Raw Reads Directory: ${params.Raw_Sequence_Data}
         Trimmed Reads Directory: ${params.Trimmed_Sequence_Data}
         Filtered Sequence Data: ${params.Filtered_Sequence_Data}
         FastQC Outputs: ${params.FastQC_Outputs}
         Final Outputs: ${params.Final_Outputs}
         Pipeline Outputs: ${params.Output_dir}
         ${c_reset}"""

}


include { CLEAN_FASTQC_PATHS; PACKAGE_PROCESSING_INFO; GENERATE_README; VALIDATE_PROCESSING;
           GENERATE_CURATION_TABLE; GENERATE_MD5SUMS; GENERATE_PROTOCOL} from './modules/genelab.nf'

workflow {

        // Make sure accessions numbers are set
        if(!params.GLDS_accession || !params.OSD_accession){
           error("""${c_back_bright_red}ACCESSION ERROR!. 
                    Please supply both --GLDS_accession and --OSD_accession.
                    They can be any string you choose but they must be set.
                 ${c_reset}""")
        }

       // ---------------------- Input channels -------------------------------- //
       // Input files
       sample_ids_file     =  Channel.fromPath(params.samples, checkIfExists: true)
       software_versions   =  Channel.fromPath(params.software_versions, checkIfExists: true)

       // Runsheet used to execute the processing workflow
       runsheet_ch = Channel.fromPath(params.runsheet, checkIfExists: true)
        
       // Directories
       Filtered_Sequence_Data  =  Channel.fromPath(params.Filtered_Sequence_Data,  type: 'dir', checkIfExists: true)
       Trimmed_Sequence_Data   =  Channel.fromPath(params.Trimmed_Sequence_Data,  type: 'dir', checkIfExists: true)
       Final_Outputs           =  Channel.fromPath(params.Final_Outputs,  type: 'dir', checkIfExists: true)

       // Input Value channels
       OSD_ch    =  Channel.of([params.name, params.email,
                                params.OSD_accession, params.protocol_id,
                                params.FastQC_Outputs, 
                                params.Filtered_Sequence_Data,
                                params.Trimmed_Sequence_Data,
                                params.Final_Outputs])

       GLDS_ch   =  Channel.of([params.GLDS_accession, params.V_V_guidelines_link, params.output_prefix,
                                params.target_files, params.assay_suffix,
                                params.raw_suffix, params.raw_R1_suffix, params.raw_R2_suffix,
                                params.primer_trimmed_suffix, params.primer_trimmed_R1_suffix,
                                params.primer_trimmed_R2_suffix, params.filtered_suffix, 
                                params.filtered_R1_suffix, params.filtered_R2_suffix])

       suffix_ch =  Channel.of([params.GLDS_accession, params.output_prefix, params.assay_suffix,
                                params.raw_suffix, params.raw_R1_suffix, params.raw_R2_suffix,
                                params.primer_trimmed_suffix, params.primer_trimmed_R1_suffix,
                                params.primer_trimmed_R2_suffix, params.filtered_suffix, 
                                params.filtered_R1_suffix, params.filtered_R2_suffix])    

        file_label_ch = Channel.of([params.processing_zip_file, params.readme])
        // Processed as paths but utilized as labels in the generate curation association table script
        dir_label_ch = Channel.of([params.Raw_Sequence_Data, params.Filtered_Sequence_Data,
                                   params.Trimmed_Sequence_Data, params.Final_Outputs])
                                   .collect()
                                  .map{ Raw_Sequence_Data, Filtered_Sequence_Data, 
                                        Trimmed_Sequence_Data, Final_Outputs -> 
                                        tuple(     file(Raw_Sequence_Data, checkIfExists: true),
                                                   file(Filtered_Sequence_Data, checkIfExists: true),
                                                   file(Trimmed_Sequence_Data, checkIfExists: true),
                                                   file(Final_Outputs, checkIfExists: true)
                                              )
                                      }  


        // If an assay table is provided use it as the input table otherwise use the provided ISA zip file
        input_table_ch = Channel.fromPath( params.assay_table ? params.assay_table : params.isa_zip,
                                          checkIfExists: true)


        // Files and directories to be package in processing_info.zip
        files_and_dirs_ch = Channel.of(params.config_dir, params.bin, params.modules, 
                                       params.envs, params.main, params.nextflow_config, 
                                       params.samples)
                                       .collect()
                                       .map{ config_dir, bin, modules, envs, main, config_file, samples -> 
                                            tuple( file(config_dir, checkIfExists: true),
                                                   file(bin, checkIfExists: true),
                                                   file(modules, checkIfExists: true),
                                                   file(envs, checkIfExists: true),
                                                   file(main, checkIfExists: true),
                                                   file(config_file, checkIfExists: true),
                                                   file(samples, checkIfExists: true)
                                                 ) }


        // ---------------------- Post-processing begins -------------------------------------- //
        PACKAGE_PROCESSING_INFO(files_and_dirs_ch)


        GENERATE_README(OSD_ch, PACKAGE_PROCESSING_INFO.out.zip)

        
        FastQC_Outputs_dir  =  Channel.fromPath(params.FastQC_Outputs, type: 'dir',  checkIfExists: true)
        CLEAN_FASTQC_PATHS(FastQC_Outputs_dir)

       // Automatic verification and validation
        VALIDATE_PROCESSING(GLDS_ch, sample_ids_file, 
                            GENERATE_README.out.readme,
                            PACKAGE_PROCESSING_INFO.out.zip, 
                            CLEAN_FASTQC_PATHS.out.clean_dir,
                            Filtered_Sequence_Data,
                            Trimmed_Sequence_Data,
                            Final_Outputs) 
        // Generate md5sums
        dirs_ch =  Channel.of(params.Trimmed_Sequence_Data,
                                   params.Filtered_Sequence_Data,
                                   params.Final_Outputs)
                                .concat(CLEAN_FASTQC_PATHS.out.clean_dir)
                                .collect()
                                .map{ trimmed_sequence, filtered_sequence, final_outputs, fastqc ->
                                      tuple( file(trimmed_sequence, checkIfExists: true),
                                             file(filtered_sequence, checkIfExists: true),
                                             file(final_outputs, checkIfExists: true),
                                             file(fastqc, checkIfExists: true)
                                           ) }
   

          GENERATE_MD5SUMS(PACKAGE_PROCESSING_INFO.out.zip,
                              GENERATE_README.out.readme, dirs_ch)

          // Generate file association table for curation
          GENERATE_CURATION_TABLE(suffix_ch, file_label_ch, 
                                  dir_label_ch, input_table_ch,
                                  CLEAN_FASTQC_PATHS.out.clean_dir)

          // Write methods
          GENERATE_PROTOCOL(software_versions, params.protocol_id)


}

workflow.onComplete {
    println("${c_bright_green}Pipeline completed at: $workflow.complete")
    println("""Execution status: ${ workflow.success ? 'OK' : "${c_back_bright_red}failed" }""")
    log.info ( workflow.success ? "\nDone! Workflow completed without any error\n" : "Oops .. something went wrong${c_reset}" )
    if ( workflow.success ) {
    println("Post-processing Outputs: ${params.Output_dir} ${c_reset}")
    }
}

