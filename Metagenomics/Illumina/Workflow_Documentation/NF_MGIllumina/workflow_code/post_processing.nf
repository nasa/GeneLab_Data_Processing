#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Terminal text color definitions
c_back_bright_red = "\u001b[41;1m";
c_bright_green    = "\u001b[32;1m";
c_blue            = "\033[0;34m";
c_reset           = "\033[0m";

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
  println("  --GLDS_accession [STRING]  A Genelab GLDS accession number. Example GLDS-574. Default: empty string")
  println("  --OSD_accession [STRING]  A Genelab OSD accession number. Example OSD-574. Default: empty string")
  println("  --name [STRING] The analyst's full name. E.g. 'FirstName A. LastName'.  Default: FirstName A. LastName")
  println("  --email [STRING] The analyst's email address. E.g. 'mail@nasa.gov'.  Default: mail@nasa.gov")
  println("  --logs_dir_basename [STRING]  Base directory name of directory containing per sample logs from processing - should always end with '/'. E.g. 'Logs/'.  Default: Logs/")
  println("  --runsheet_basename  [String] The runsheets base name. Example 'GLfile.csv' or 'PE_file.csv'. Default: null")
  println("  --assay_suffix [STRING]  Genelab's assay suffix. Default: _GLmetagenomics.")
  println("  --output_prefix [STRING] Unique name to tag onto output files. Default: empty string.")
  println("  --V_V_guidelines_link [URL] Genelab metagenomics data validation and verification guidelines link. Default: https://genelab-tools.arc.nasa.gov/confluence/pages/viewpage.action?pageId=8225175.")
  println("  --target_files [STRING] A comma separated list of target files and/or directories to find in processing_info.zip. Default: nextflow_processing_info_GLmetagenomics.txt,command.txt,unique-sample-IDs.txt,software_versions.txt,GLfile.csv,<--logs>.")
  println("File Suffixes:")
  println("      --raw_suffix [STRING]  Suffix used for the raw reads during processing. Only applicable when input reads are single-end. Default: _HRremoved_raw.fastq.gz.")  
  println("      --raw_R1_suffix [STRING]  Suffix used for the raw forward reads during processing. Default: _R1_HRremoved_raw.fastq.gz.")
  println("      --raw_R2_suffix [STRING]  Suffix used for the raw reverse reads during processing. Default: _R2_HRremoved_raw.fastq.gz.")

  println("      --filtered_suffix [STRING]  Suffix used for quality filtered reads during processing. Only applicable when input reads are single-end. Default: _filtered.fastq.gz.")
  println("      --filtered_R1_suffix [STRING]  Suffix to use for quality filtered forward reads during processing. Default: _R1_filtered.fastq.gz.")
  println("      --filtered_R2_suffix [STRING]  Suffix to use for quality filtered reverse reads during processing. Default: _R2_filtered.fastq.gz.")
  println()
  println("Extra parameters to scripts:")
  println("      --readme_extra [STRING] Extra parameters and arguments to GL-gen-processed-metagenomics-data-readme command. Run 'GL-gen-processed-metagenomics-readme --help' for extra parameters that can be set. Example '--raw-reads-dir  ../Raw_Sequence_Data/'. Default: empty string")
  println("      --validation_extra [STRING] Extra parameters and arguments to GL-validate-processed-metagenomics-data command. Run 'GL-validate-processed-metagenomics-data --help' for extra parameters that can be set. Example '--single-ended --R1-used-as-single-ended-data --skip_raw_multiqc'. Default: '--skip_raw_multiqc' ")
  println("      --file_association_extra [STRING] Extra parameters and arguments to GL-gen-metagenomics-file-associations-table command. Run 'GL-gen-metagenomics-file-associations-table --help' for extra parameters that can be set. Example '--single-ended --R1-used-as-single-ended-data'. Default: '--use-sample-names-from-assay-table' ")
  println()
  println("Files:")
  println("    --run_command  [PATH] File containing the nextflow run command used in processing. Default: ./processing_scripts/command.txt")
  println("    --processing_commands  [PATH] File containing all the process names and scripts used during processing. Default: ./processing_scripts/nextflow_processing_info_GLAmpliseq.txt")
  println("    --samples  [PATH] A single column file with sample ids on each line generated after running the processing pipeline. Default: ./unique-sample-IDs.txt")
  println("    --assay_table  [PATH] GLDS assay table generated after running the processing pipeline with accession number as input.")
  println("                   Example, ../Genelab/a_OSD-574_metagenomic-sequencing_whole-genome-shotgun-sequencing_illumina.txt. Default: empty string")
  println("    --isa_zip  [PATH] Genelab ISA zip files containing an assay atable for the OSD accession. This is only required if --files.assay_table is not set.")
  println("                   Example, ../Genelab/OSD-574_metadata_OSD-574-ISA.zip. Default: empty string")
  println("    --runsheet  [PATH] A 3-column (single-end) or 4-column (paired-end) input file (sample_id, forward, [reverse,] paired) used to run the processing pipeline. This is the value set to the paremater --input_file when run the processing pipeline with a csv file as input otherwise it is the GLfile.csv in the GeneLab directory if --GLDS_accession was used as input. Example '../GeneLab/GLfile.csv'.  Default: null")


 println("    --software_versions  [PATH] A file generated after running the processing pipeline listing the software versions used. Default: ../Metadata/software_versions.txt")
  println()
  println("Directories:")
  println("  --logs_dir [PATH]  Full or relative path to directory name of directory containing per sample logs from processing - should always end with '/'. E.g. 'Logs/'.  Default: Logs/")
  println("    --Raw_Sequence_Data [PATH] A directory containing raw sequence and raw sequence outputs. Default: ../Raw_Sequence_Data/")
  println("    --FastQC_Outputs [PATH] A directory containing fastqc and multiqc zip reports. Default: ../FastQC_Outputs/")
  println("    --Filtered_Sequence_Data  [PATH] A directory containing the outputs of read filtering after running the processing pipeline. Default: ../Filtered_Sequence_Data/")  
  println("    --Read_based_Processing  [PATH] A directory containing the outputs of read based processing after running the processing pipeline. Default: ../Read_based_Processing/")
  println("    --Assembly_based_Processing  [PATH] A directory containing the outputs of assembly based processing after running the processing pipeline. Default: ../Assembly_based_Processing/")
  println("    --Assemblies  [PATH] A directory containing sample contig assemblies after running the processing pipeline. Default: ../Assembly_based_Processing/assemblies/")
  println("    --Genes  [PATH] A directory containing sample predicted genes after running the processing pipeline. Default: ../Assembly_based_Processing/predicted-genes/")
  println("    --Annotations_And_Tax  [PATH] A directory containing sample gene and contig annotations after running the processing pipeline. Default: ../Assembly_based_Processing/annotations-and-taxonomy/")
  println("    --Mapping  [PATH] A directory containing sample read mapping (bam) files after running the processing pipeline. Default: ../Assembly_based_Processing/read-mapping/")
  println("    --Combined_Output  [PATH] A directory containing assembly summaries and reports across samples after running the processing pipeline. Default: ../Assembly_based_Processing/combined-outputs/")
  println("    --Bins  [PATH] A directory containing metagenome bins after running the processing pipeline. Default: ../Assembly_based_Processing/bins/")
  println("    --MAGS  [PATH] A directory containing metagenome assembled genomes (MAGS) after running the processing pipeline. Default: ../Assembly_based_Processing/MAGs/")
  println("    --Output_dir  [PATH] Specifies the directory where outputs of this post-processing workflow will be published. Default: ../Post_Processing/")
  println()
  println("Optional arguments:")  
  println("    --help [BOOLEAN] Print this help message and exit")
  println("    --debug [BOOLEAN] Set to true if you'd like to see the values of your set parameters printed to the terminal. Default: false.")
  println()
  println("Paths to existing conda environments to use otherwise a new one will be created using the yaml file in envs/.")
  println("    --conda_genelab [PATH] Path to a conda environment containing genelab-utils. Default: null.")
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
         Analyst's Name: ${params.name}
         Analyst's Email: ${params.email}
         GLDS Accession: ${params.GLDS_accession}
         OSD Accession: ${params.OSD_accession}
         Assay Suffix: ${params.assay_suffix}
         Output Prefix: ${params.output_prefix}
         V & V Link: ${params.V_V_guidelines_link}
         Target Files: ${params.target_files} 
         Nextflow Directory publishing mode: ${params.publishDir_mode}
        
         Base Names:
         Run sheet: ${params.runsheet_basename}
         Log Directory: ${params.logs_dir_basename}

         Suffixes:
         Raw Suffix: ${params.raw_suffix}
         Raw R1 suffix: ${params.raw_R1_suffix}
         Raw R2 suffix: ${params.raw_R2_suffix}          
         Filtered Suffix: ${params.filtered_suffix}
         Filtered R1 suffix: ${params.filtered_R1_suffix}
         Filtered R2 suffix: ${params.filtered_R2_suffix}

         Extra scripts parameters:
         Readme Script Extra: ${params.readme_extra}
         Validation Script Extra: ${params.validation_extra}
         File association Script Extra: ${params.file_association_extra}

         Files:
         Nextflow Command: ${params.run_command}
         Processing Commands: ${params.processing_commands}
         Samples: ${params.samples}
         Assay Table: ${params.assay_table}
         ISA Zip: ${params.isa_zip}
         Input Runsheet: ${params.runsheet}
         Software Versions: ${params.software_versions}

         Directories:
         Logs directory: ${params.logs_dir}
         Raw Reads Directory: ${params.Raw_Sequence_Data}
         Filtered Sequence Data: ${params.Filtered_Sequence_Data}
         FastQC Outputs: ${params.FastQC_Outputs}
         Read-based Processing: ${params.Read_based_Processing}
         Assemblies: ${params.Assemblies}
         Genes: ${params.Genes}
         Annotations And Taxonomy: ${params.Annotations_And_Tax}
         Mapping: ${params.Mapping}
         Combined Output: ${params.Combined_Output}
         Bins: ${params.Bins}
         MAGS: ${params.MAGS}
         Pipeline Outputs: ${params.Output_dir}
         """

}


include { CLEAN_FASTQC_PATHS; PACKAGE_PROCESSING_INFO; GENERATE_README; VALIDATE_PROCESSING;
           GENERATE_CURATION_TABLE; GENERATE_MD5SUMS; GENERATE_PROTOCOL} from './modules/genelab.nf'

workflow {

        // Make sure accession numbers are set
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

       // Directories
       Bins                =  Channel.fromPath(params.Bins,  type: 'dir', checkIfExists: true)
       MAGS                =  Channel.fromPath(params.MAGS,  type: 'dir', checkIfExists: true)

       // Input Value channels
       OSD_ch    =  Channel.of([params.name, params.email, params.output_prefix,
                                params.OSD_accession, params.protocol_id,
                                params.FastQC_Outputs, 
                                params.Filtered_Sequence_Data,
                                params.Read_Based_Processing,
                                params.Assembly_Based_Processing,
                                params.Assemblies, 
                                params.Genes,
                                params.Annotations_And_Tax,
                                params.Mapping,
                                params.Combined_Output])
      
       GLDS_ch   =  Channel.of([params.GLDS_accession, params.V_V_guidelines_link, params.output_prefix,
                                params.target_files, params.assay_suffix, params.logs_dir_basename,
                                params.raw_suffix, params.raw_R1_suffix, params.raw_R2_suffix,
                                params.filtered_suffix, params.filtered_R1_suffix, params.filtered_R2_suffix])

       suffix_ch =  Channel.of([params.GLDS_accession, params.output_prefix, params.assay_suffix,
                                params.raw_suffix, params.raw_R1_suffix, params.raw_R2_suffix,
                                params.filtered_suffix, params.filtered_R1_suffix, params.filtered_R2_suffix]) 

        file_label_ch = Channel.of([params.processing_zip_file, params.readme])

        // processed as paths but utilized as labels in the genberate curation association table script 
        dir_label_ch = Channel.of([params.Raw_Sequence_Data,
                                   params.Filtered_Sequence_Data,
                                   params.Read_Based_Processing,
                                   params.Assembly_Based_Processing,
                                   params.Annotations_And_Tax,
                                   params.Combined_Output])
                                  .collect()
                                  .map{ Raw_Sequence_Data, Filtered_Sequence_Data, Read_Based_Processing,
                                        Assembly_Based_Processing, Annotations_And_Tax, Combined_Output -> 
                                        tuple(    file(Raw_Sequence_Data, checkIfExists: true),
                                                   file(Filtered_Sequence_Data, checkIfExists: true),
                                                   file(Read_Based_Processing, checkIfExists: true),
                                                   file(Assembly_Based_Processing, checkIfExists: true),
                                                   file(Annotations_And_Tax, checkIfExists: true),
                                                   file(Combined_Output, checkIfExists: true)
                                              )
                                      }  

        // If the assay table is provided use it as the input table otherwise use the isa_zip
        assay_table_ch = Channel.fromPath( params.assay_table ?  params.assay_table : params.isa_zip,
                                          checkIfExists: true)

        // Runsheet used to execute the processing workflow
        runsheet_ch = Channel.fromPath(params.runsheet, checkIfExists: true)



       // Files and directories to be packaged in processing_info.zip
        files_and_dirs_ch = Channel.of(params.logs_dir, params.run_command, params.processing_commands,
                                       params.software_versions, params.runsheet, params.samples)
                                       .collect()
                                       .map{ logs, run_command, processing_commands, software_versions, runsheet, samples ->
                                            tuple( file(logs, checkIfExists: true),
                                                   file(run_command, checkIfExists: true),
                                                   file(processing_commands, checkIfExists: true),
                                                   file(software_versions, checkIfExists: true),
                                                   file(runsheet, checkIfExists: true),
                                                   file(samples, checkIfExists: true)
                                                 ) }

        // ---------------------- Post-processing begins ---------------------------------//
        PACKAGE_PROCESSING_INFO(files_and_dirs_ch)


        GENERATE_README(OSD_ch, PACKAGE_PROCESSING_INFO.out.zip, Bins, MAGS)

        
        FastQC_Outputs_dir  =  Channel.fromPath(params.FastQC_Outputs, type: 'dir',  checkIfExists: true)
        CLEAN_FASTQC_PATHS(FastQC_Outputs_dir)

        validation_dirs_ch =  Channel.of(params.Filtered_Sequence_Data,
                                         params.Read_Based_Processing,
                                         params.Assembly_Based_Processing, 
                                         params.Assemblies,
                                         params.Mapping,
                                         params.Genes, 
                                         params.Annotations_And_Tax,
                                         params.Bins, 
                                         params.MAGS,
                                         params.Combined_Output)
                            .concat(CLEAN_FASTQC_PATHS.out.clean_dir)
                            .collect()
                            .map{ filtered_sequence, read_based, assembly_based, assemblies, 
                                  mapping, genes, annotation, bins, mags, combined_output, fastqc ->
                                    tuple( file(filtered_sequence, checkIfExists: true),
                                           file(read_based, checkIfExists: true),
                                           file(assembly_based, checkIfExists: true),
                                           file(assemblies, checkIfExists: true),
                                           file(mapping, checkIfExists: true),
                                           file(genes, checkIfExists: true),
                                           file(annotation, checkIfExists: true),
                                           file(bins, checkIfExists: true),
                                           file(mags, checkIfExists: true),
                                           file(combined_output, checkIfExists: true),
                                           file(fastqc, checkIfExists: true)
                                          ) }

        // Automatic verification and validation
        VALIDATE_PROCESSING(GLDS_ch, validation_dirs_ch,
                            sample_ids_file, 
                            GENERATE_README.out.readme,
                            PACKAGE_PROCESSING_INFO.out.zip) 

        // Generate md5sums
        dirs_ch = Channel.of(params.Read_Based_Processing,
                             params.Filtered_Sequence_Data,
                             params.Assembly_Based_Processing)
                              .concat(CLEAN_FASTQC_PATHS.out.clean_dir)
                              .collect()
                              .map{ read_based, filtered_sequence, assembly_based, fastqc ->
                                    tuple( file(read_based, checkIfExists: true),
                                           file(filtered_sequence, checkIfExists: true),
                                           file(assembly_based, checkIfExists: true),
                                           file(fastqc, checkIfExists: true)
                                          ) }

         GENERATE_MD5SUMS(PACKAGE_PROCESSING_INFO.out.zip,
                            GENERATE_README.out.readme, dirs_ch)


          // Generate curation file association table
          curation_dirs_ch =   Channel.of(params.Assemblies,
                                          params.Genes,
                                          params.Mapping,
                                          params.Bins,
                                          params.MAGS)
                              .concat(CLEAN_FASTQC_PATHS.out.clean_dir)
                              .collect()
                              .map{ assemblies, genes, mapping, bins, mags, fastqc ->
                                    tuple( file(assemblies, checkIfExists: true),
                                           file(genes, checkIfExists: true),
                                           file(mapping, checkIfExists: true),
                                           file(bins, checkIfExists: true),
                                           file(mags, checkIfExists: true),
                                           file(fastqc, checkIfExists: true)
                                          ) } 

          GENERATE_CURATION_TABLE(suffix_ch, file_label_ch, 
                                  dir_label_ch, 
                                  curation_dirs_ch, 
                                  assay_table_ch, runsheet_ch)


         // Write methods
         GENERATE_PROTOCOL(software_versions, params.protocol_id)
}


workflow.onComplete {

    println("${c_bright_green}Pipeline completed at: $workflow.complete")
    println("""Execution status: ${ workflow.success ? 'OK' : "${c_back_bright_red}failed" }""")
    log.info ( workflow.success ? "\nDone! Workflow completed without any error\n" : "Oops .. something went wrong${c_reset}" )

    if ( workflow.success ) {

    println("Post-processing Outputs: ${params.Output_dir} ${c_reset}")
    println()

    }
}
