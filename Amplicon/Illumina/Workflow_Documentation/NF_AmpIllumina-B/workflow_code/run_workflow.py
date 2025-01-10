#!/usr/bin/env python
import argparse
import subprocess
import os
import sys



def create_config(target_region, raw_R1_suffix, raw_R2_suffix, trim_primers, input_file, min_cutadapt_len, primers_linked,
                  discard_untrimmed, F_primer, R_primer, left_trunc, right_trunc, left_maxEE, right_maxEE, concatenate_reads_only,
                  conda_genelab, conda_qc, conda_R, conda_cutadapt, conda_diversity, accession, assay_suffix, output_prefix, publishDir_mode,
                  primer_trimmed_R1_suffix, primer_trimmed_R2_suffix, filtered_R1_suffix, filtered_R2_suffix, output_dir, diff_abund_method,
                  group, samples_column, rarefaction_depth, errorStrategy, queueSize, default_cpus, default_memory, cutadapt_memory, R_cpus,
                  R_memory, diversity_cpus, diversity_memory, container_genelab, container_fastqc, container_multiqc, container_cutadapt, 
                  container_dada, container_ancom, container_diversity, singularity_cacheDir, remove_rare, prevalence_cutoff, library_cutoff):
    """ A function to create nextflow.config file by string interploation using the supplied arguements"""

    config =  \
    f"""
//***************************************** Global parameters *******************************************//
params {{
    // Mandatory parameters 
    target_region = "{target_region}"
    raw_R1_suffix = "{raw_R1_suffix}"
    raw_R2_suffix = "{raw_R2_suffix}"
    trim_primers  = "{trim_primers}" == "TRUE" ? true : false
    

    // -------- Required only if --accession is false ---------------//
    // A 4-column (single-end) or 5-column (paired-end) input csv file with the following headers ( sample_id, forward, [reverse,] paired, groups)
    input_file = "{input_file}" == "null" ? null : "{input_file}"


    // Cutadapt parameters
    min_cutadapt_len    = {min_cutadapt_len}
    primers_linked      = "{primers_linked}"
    discard_untrimmed   = "{discard_untrimmed}"
    F_primer            = "{F_primer}"  == "null" ? null : "{F_primer}"
    R_primer            = "{R_primer}"  == "null" ? null : "{R_primer}"

    // Dada2 parameters
    left_trunc     = {left_trunc}
    right_trunc    = {right_trunc}
    left_maxEE     = {left_maxEE}
    right_maxEE    = {right_maxEE}
    concatenate_reads_only = "{concatenate_reads_only}"


    // If using conda environments specify their locations so new ones won't be created
    conda{{
          // Specify the paths to existing conda environments (/path/to/envs/genelab-utils)
          // leave as is if you want to create a new conda environment
          genelab          = "{conda_genelab}" == "null" ? null :  "{conda_genelab}"     // /path/to/envs/genelab-utils
          qc               = "{conda_qc}" == "null" ? null : "{conda_qc}"   // /path/to/envs/qc
          R                = "{conda_R}" == "null" ? null : "{conda_R}"     // /path/to/envs/R
          cutadapt         = "{conda_cutadapt}" == "null" ? null :  "{conda_cutadapt}"  // /path/to/envs/cutadapt
          diversity        = "{conda_diversity}" == "null" ? null :  "{conda_diversity}"   // /path/to/envs/R_diversity
      }}


    // Mandatory parameters if using GLDS or OSD accession as input
    accession =  "{accession}" == "null" ? null :  "{accession}"

    assay_suffix    = "{assay_suffix}"
    output_prefix   = "{output_prefix}"
    publishDir_mode = "{publishDir_mode}"

    // Suffixes
    primer_trimmed_R1_suffix = "{primer_trimmed_R1_suffix}"
    primer_trimmed_R2_suffix = "{primer_trimmed_R2_suffix}"
    filtered_R1_suffix       = "{filtered_R1_suffix}"
    filtered_R2_suffix       = "{filtered_R2_suffix}"


    // Directories
    
    fastqc_out_dir      = "{output_dir}/workflow_output/FastQC_Outputs/"
    trimmed_reads_dir   = "{output_dir}/workflow_output/Trimmed_Sequence_Data/"
    filtered_reads_dir  = "{output_dir}/workflow_output/Filtered_Sequence_Data/"
    final_outputs_dir   = "{output_dir}/workflow_output/Final_Outputs/"

    raw_reads_dir       = "{output_dir}/Raw_Sequence_Data/" 
    metadata_dir        = "{output_dir}/Metadata/"
    genelab_dir         = "{output_dir}/GeneLab/"

    // Multiqc
    multiqc_config = "${{projectDir}}/config/multiqc.config"

    // -------- Differential abundance parameters ----- //
    diff_abund_method = "{diff_abund_method}"
    group             = "{group}"
    samples_column    = "{samples_column}"
    // Should rare features and samples be discarded. Values are true or false. If set to true then set the cutoffs below
    remove_rare       = "{remove_rare}" == "TRUE" ? true : false
    prevalence_cutoff = {prevalence_cutoff}  // a fraction between 0 and 1 that represents the prevalance in percentage of taxa to be retained
    library_cutoff    = {library_cutoff} // Samples with library sizes less than this number will be excluded in the analysis

    // Minimum desired sample rarefaction depth for diversity analysis
    rarefaction_depth  = {rarefaction_depth}
 

    errorStrategy  = "{errorStrategy}"
    debug          = false // set to true if you'd like to see the values of your set parameters
}}

// Setting the default container engine as singularity
params.containerEngine = "singularity"
// Conda shouldn't be used by default except when using conda-based profiles
params.use_conda = false


/*******************************************************************************************************
*************************************** Workflow Profiles **********************************************
********************************************************************************************************/

profiles {{

    slurm {{  
        process.executor      = 'slurm'
    }}

    conda {{   
        conda.enabled          = true
        params.use_conda       = true               
    }}

    singularity {{
        singularity.enabled    = true
        singularity.autoMounts = true
        singularity.cacheDir   = "{singularity_cacheDir}"
        params.containerEngine = "singularity"
    }}

    docker {{
        docker.enabled         = true
        docker.runOptions      = '-u $(id -u):$(id -g)'
        docker.userEmulation   = true
        params.containerEngine = "docker"
    }}

}}

// Maximum number of jobs to submit in parallel
executor.queueSize = {queueSize}


/******************************************************************************************************************
***************** Tune process specific resources (cpu, container, memory etc.) ***********************************
*******************************************************************************************************************/

process {{

    //******************* Default process settings ************************//
    errorStrategy = {{ params.errorStrategy ? params.errorStrategy : "ignore" }} 
    maxRetries = 2
    cpus = {default_cpus}
    memory = "{default_memory}"
    cache = 'lenient'
  //debug = true  // uncomment to see what is being emitted to the standard output

//************************* Accession runsheet and input file retrieval  **************************************//
    withName: GET_RUNSHEET {{
                  conda = {{params.conda.genelab ? params.conda.genelab : "envs/genelab.yaml"}}
                  container = "{container_genelab}"
                  publishDir = [path: params.genelab_dir, mode: params.publishDir_mode]
            }}

//********************************** Read quality control and assesment ********************************************//
    withLabel: fastqc {{
                  conda = {{params.conda.qc ? params.conda.qc : "envs/qc.yaml"}}
                  container = "{container_fastqc}"
            }}

    withName: RAW_FASTQC {{                  
                  publishDir = [path: params.raw_reads_dir, mode: params.publishDir_mode]
            }}

    withName: "RAW_MULTIQC|TRIMMED_MULTIQC" {{
                  conda = {{params.conda.qc ? params.conda.qc : "envs/qc.yaml"}}
                  container = "{container_multiqc}"
                  publishDir = [path: params.fastqc_out_dir, mode: params.publishDir_mode]
            }}

    withName: "CUTADAPT|COMBINE_CUTADAPT_LOGS_AND_SUMMARIZE" {{
                  conda = {{params.conda.cutadapt ?  params.conda.cutadapt : "envs/cutadapt.yaml"}}
                  container = "{container_cutadapt}"
                  memory = "{cutadapt_memory}"
                  publishDir = [path: params.trimmed_reads_dir, mode: params.publishDir_mode]
            }}
           
    withName: TRIMMED_FASTQC {{
                  publishDir = [path: params.filtered_reads_dir, mode: params.publishDir_mode ]
            }} 

//********************************** ASV table creation********************************************//
    withName: "RUN_R_TRIM|RUN_R_NOTRIM" {{
                  conda = {{params.conda.R ?  params.conda.R : "envs/R.yaml"}}
                  container = "{container_dada}"
                  memory = "{R_memory}"
                  cpus = {R_cpus} 
                  publishDir = [[path: params.filtered_reads_dir, pattern: "Filtered_Sequence_Data/*",
                                mode: params.publishDir_mode, saveAs: {{ fn -> fn.substring(fn.lastIndexOf('/')+1) }} ],
                                [path: params.final_outputs_dir , pattern: "final_outputs/*.{{tsv,biom,fasta}}",
                                mode: params.publishDir_mode, saveAs: {{ fn -> fn.substring(fn.lastIndexOf('/')+1)}} ]] 
          }}

    withName: ZIP_BIOM {{
                  conda = {{params.conda.qc ? params.conda.qc : "envs/qc.yaml"}}
                  container = "{container_multiqc}"
                  publishDir = [path: "${{params.final_outputs_dir}}${{params.output_prefix}}", mode: params.publishDir_mode]
            }}

//********************************** Diversity and differential abundance testing ********************************************//
    withLabel: visualization {{
                  conda = {{params.conda.diversity ? params.conda.diversity : "envs/diversity.yaml"}}
                  container = "{container_diversity}"
                  cpus = {diversity_cpus}
                  memory = "{diversity_memory}"
                  publishDir = [path: "${{params.final_outputs_dir}}${{params.output_prefix}}", mode: params.publishDir_mode]
           }}


    withName: ANCOMBC {{

             container = "{container_ancom}"

             }}

    withName: SOFTWARE_VERSIONS {{
                  publishDir = [path: params.metadata_dir, mode: params.publishDir_mode]
            }}

}}


/*****************************************************************************
********************** Workflow Resource Usage Capturing *********************
******************************************************************************/

// Adapted from : https://github.com/nf-core/rnaseq/blob/master/nextflow.config
def trace_timestamp = new java.util.Date().format( 'yyyy-MM-dd_HH-mm-ss')
timeline {{
    enabled = true
    file    = "../Resource_Usage/execution_timeline_${{trace_timestamp}}.html"
}}
report {{
    enabled = true
    file    = "../Resource_Usage/execution_report_${{trace_timestamp}}.html"
}}
trace {{
    enabled = true
    file    = "../Resource_Usage/execution_trace_${{trace_timestamp}}.txt"
}}



/******************************************************************************
**************************** Workflow Metadata ********************************
*******************************************************************************/

manifest {{
    author = 'Olabiyi Aderemi Obayomi, Mike D. Lee'
    homePage = 'https://github.com/nasa/GeneLab_Data_Processing/blob/master/Amplicon/'
    description = 'Amplicon Illumina workflow for pipeline document GL-DPPD-7104-B'
    mainScript = 'main.nf'
    defaultBranch = 'main'
    nextflowVersion = '>=24.04.4'
    version = '1.0.0'
}}

    """

    return config



def main():
    # Argument parser setup with short argument names and an automatic help option
    parser = argparse.ArgumentParser(
        description='Run workflow for GeneLab data processing.',
        add_help=True,
        usage="""
        Example 1: Use an OSD or Genelab acession as input
        %(prog)s --run --target-region 16S --accession GLDS-487 --profile slurm,singularity [options]
        
        Example 2: Use a csv file as input to the workflow
        %(prog)s --run --target-region 16S --input-file PE_file.csv --F-primer AGAGTTTGATCCTGGCTCAG --R-primer CTGCCTCCCGTAGGAGT --profile singularity [options]

        Example 3: Use a csv file as input to the workflow and supply extra arguments to nextflow run.
                  Here were want to mintor our jobs with nextflow tower.
        export TOWER_ACCESS_TOKEN=<ACCESS TOKEN>
        export TOWER_WORKSPACE_ID=<WORKSPACE ID>
        %(prog)s --run --target-region 16S --input-file PE_file.csv --F-primer AGAGTTTGATCCTGGCTCAG --R-primer CTGCCTCCCGTAGGAGT --profile slurm,conda --extra 'with-tower' [options]

        Example 4: Dry run: Just create an edited nextflow.config file but don't run the workflow
        %(prog)s --target-region 16S --accession GLDS-487 --profile slurm,singularity [options]
        """
    )
    
    parser.add_argument('-x', '--run',
                        action='store_true',
                        help="Set this flag if would like to run the workflow")

    parser.add_argument('-a', '--accession',
                        metavar='osd_number',
                        default="null",
                        help="""
                             A Genelab or OSD accession number if the --input-file parameter is not set.
                             If this parameter is set, it will ignore the --input-file parameter. Default: null
                             """,
                        type=str)
    
    parser.add_argument('-r', '--input-file',
                        metavar='/path/to/input_file.csv',
                        default="null",
                        help=""" 
                        A 4-column (single-end) or 5-column (paired-end) input file (sample_id, forward, [reverse,] paired, groups). Mandatory if a GLDS or OSD accession is not provided.
                        Please see the files: SE_file.csv and PE_file.csv for single-end and paired-end examples, respectively. The sample_id column should contain unique sample ids.
                        The forward and reverse columns should contain the absolute or relative path to the sample's forward and reverse reads.
                        The paired column should be true for paired-end or anything else for single-end reads.
                        The groups column should contain group levels / treatments to be compared during diversity and differential abundance testing analysis.
                        """,
                        type=str)

    parser.add_argument('-p', '--profile',
                        metavar='profile',
                        default="null",
                        help=""" 
                        What profile(s) should be used to run the workflow? Options are [singularity, docker, conda, slurm].
                        singularity, docker and conda will run the pipelne locally using singularity, docker, and conda, respectively.
                        To combine profiles, pass them together separated by comma. For example, to run jobs using slurm in singularity containers use 'slurm,singularity'.
                        .Default: null """,
                        type=str)
    
    parser.add_argument('-t', '--target-region',
                        default="16S",
                        choices=['16S', '18S', 'ITS'],
                        help="""
                        Specify the genomic target for the assay. Options: 16S, 18S, ITS. This is used to select the appropriate 
                        dataset from an OSD study when multiple options are available and also determines the database to use for taxonomy assignment. Default: 16S
                        """,
                        type=str)

    parser.add_argument('--raw-R1-suffix',
                        default='_R1_raw.fastq.gz',
                        help='Raw forward reads suffix (region following the unique part of the sample names). Default: _R1_raw.fastq.gz',
                        metavar='raw_R1_suffix',
                        type=str)

    parser.add_argument('--raw-R2-suffix',
                        default='_R2_raw.fastq.gz',
                        help='Raw reverse reads suffix (region following the unique part of the sample names). Default: _R2_raw.fastq.gz',
                        metavar='raw_R2_suffix',
                        type=str)

    parser.add_argument('--primer-trimmed-R1-suffix',
                        default='_R1_trimmed.fastq.gz',
                        help='Suffix to use for naming your primer trimmed forward reads. Default: _R1_trimmed.fastq.gz',
                        metavar='primer_trimmed_R1_suffix',
                        type=str)

    parser.add_argument('--primer-trimmed-R2-suffix',
                        default='_R2_trimmed.fastq.gz',
                        help='Suffix to use for naming your primer trimmed reverse reads. Default: _R2_trimmed.fastq.gz',
                        metavar='primer_trimmed_R2_suffix',
                        type=str)

    parser.add_argument('--filtered-R1-suffix',
                        default='_R1_filtered.fastq.gz',
                        help='Suffix to use for naming your quality filtered forward reads. Default: _R1_filtered.fastq.gz',
                        metavar='filtered_R1_suffix',
                        type=str)

    parser.add_argument('--filtered-R2-suffix',
                        default='_R2_filtered.fastq.gz',
                        help='Suffix to use for naming your quality filtered reverse reads. Default: _R2_filtered.fastq.gz',
                        metavar='filtered_R2_suffix',
                        type=str)

    parser.add_argument('--trim-primers',
                        choices=["TRUE", "FALSE"],
                        default="TRUE",
                        help='Specifies to trim primers (TRUE) or not (FALSE) using cutadapt. Default: "TRUE" ',
                        type=str)

    parser.add_argument('--F-primer',
                        default='null',
                        help='Forward primer sequence e.g. AGAGTTTGATCCTGGCTCAG. Default: null."',
                        metavar='FORWARD PRIMER',
                        type=str)

    parser.add_argument('--R-primer',
                        default='null',
                        help='Reverse primer sequence e.g. AGAGTTTGATCCTGGCTCAG. Default: emptry string."',
                        metavar='REVERSE PRIMER',
                        type=str)

    parser.add_argument('--group-column',
                        default='groups',
                        help='Column in input csv file with treatments to be compared. Default: groups',
                        metavar='group column',
                        type=str)

    parser.add_argument('--samples-column',
                        default='sample_id',
                        help="Column in input csv file with sample names belonging to each treatment group. Default: sample_id",
                        metavar='samples column',
                        type=str)
        
    parser.add_argument('--rarefaction-depth',
                        metavar='depth',
                        default=500, 
                        help='The Minimum desired sample rarefaction depth for diversity analysis. Default: 500',
                        type=int)
    
    parser.add_argument('--diff-abund-method',
                        choices=["all", "ancombc1", "ancombc2", "deseq2"],
                        default="all",
                        help="The method to use for differential abundance testing. Either ancombc1, ancombc2, or deseq2. Default: ancombc2" ,
                        type=str)

    parser.add_argument('-d', '--output-dir',
                        metavar='/path/to/output_dir/',
                        default='..', 
                        help='Specifies the output directory for the output files generated by the workflow. Default: ..',
                        type=str)

    parser.add_argument('-m', '--min-cutadapt-len',
                        metavar='length',
                        default=130,
                        help='Specifies the MINIMUM length of trimmed reads to pass to cutadapt. For paired-end data: if one read gets filtered, both reads are discarded. Default: 130',
                        type=int)

    parser.add_argument('--primers-linked',
                        choices=["TRUE", "FALSE"],
                        default="TRUE",
                        help='If set to TRUE, instructs cutadapt to treat the primers as linked. Default: FALSE',
                        type=str)


    parser.add_argument('--discard-untrimmed',
                        choices=['TRUE', 'FALSE'],
                        default='TRUE',
                        help='If set to TRUE, instructs cutadapt to remove reads if the primers were not found in the expected location; if FALSE, these reads are kept. Default: TRUE',
                        type=str)

    parser.add_argument('--left-trunc',
                        default=0,
                        help='Specifies the length of the forwards reads, bases beyond this length will be truncated and reads shorter than this length are discarded. Default: 0 (no truncation)',
                        metavar='length',
                        type=int)

    parser.add_argument('--right-trunc',
                        default=0,
                        help='Specifies the length of the reverse reads, bases beyond this length will be truncated and reads shorter than this length are discarded. Default: 0 (no truncation)',
                        metavar='length',
                        type=int)

    parser.add_argument('--left-maxEE',
                        default=1,
                        help='Specifies the maximum expected error (maxEE) allowed for each forward read, reads with higher than maxEE will be discarded. Default: 1',
                        metavar='left_max_error',
                        type=int)

    parser.add_argument('--right-maxEE',
                        default=1,
                        help='Specifies the maximum expected error (maxEE) allowed for each forward read, reads with higher than maxEE will be discarded. Default: 1',
                        metavar='right_max_error',
                        type=int)

    parser.add_argument('--concatenate-reads-only',
                        choices=['TRUE', 'FALSE'],
                        default='FALSE',
                        help='If set to TRUE, specifies to concatenate forward and reverse reads only with dada2 instead of merging paired reads. Default: FALSE',
                        type=str)

    parser.add_argument('--remove-rare',
                        choices=['TRUE', 'FALSE'],
                        default='FALSE',
                        help='If set to TRUE, rare features and samples with library size less than cutoff will be removed. Default: FALSE',
                        type=str)

    parser.add_argument('--prevalence-cutoff',
                        default=0,
                        help='If --remove-rare is TRUE, a numerical fraction between 0 and 1. Taxa with prevalences(the proportion of samples in which the taxon is present) less than this will be excluded from diversity and differential abundance analysis. Default is 0 , i.e. do not exclude any taxa. For example, to exclude taxa that are not present in at least 15% of the samples set it to 0.15.',
                        metavar='prevalence_cutoff',
                        type=float)

    parser.add_argument('--library-cutoff',
                        default=0,
                        help='If --remove-rare is TRUE, a numerical threshold for filtering samples based on library sizes. Samples with library sizes less than this number will be excluded in the analysis. Default is 0 i.e do not remove any sample. For example, if you want to discard samples with library sizes less than 100, then set to 100.',
                        metavar='library_cutoff',
                        type=int)

    parser.add_argument('--output-prefix',
                        default='',
                        help='Specifies the prefix to use on all output files to distinguish multiple primer sets, leave as an empty string if only one primer set is being processed. Default: empty string',
                        metavar='prefix',
                        type=str)

    parser.add_argument('--assay-suffix',
                        default='_GLAmpSeq',
                        help='Genelabs assay suffix. Default: GLAmpSeq',
                        metavar='suffix',
                        type=str)
    
    parser.add_argument('--publishDir-mode',
                        choices=['link', 'copy', 'sym'],
                        default='link',
                        help="How should nextflow publish file outputs. Options can be found here https://www.nextflow.io/docs/latest/process.html#publishdir. Default: link",
                        metavar='publishDir_mode',
                        type=str)

    parser.add_argument('--errorStrategy',
                        default='terminate',
                        choices=['terminate', 'ignore'],
                        help="How should nextflow handle errors. Options can be found here https://www.nextflow.io/docs/latest/process.html#errorstrategy. Default: terminate",
                        metavar='errorStrategy',
                        type=str)
    
    parser.add_argument('--queueSize',
                        default=20,
                        help="Maximum number of jobs to submit in parallel. Default: 20",
                        metavar='queueSize',
                        type=int)

    parser.add_argument('--default-cpus',
                        default=2,
                        help="Default number of cpus for each job. Default: 2",
                        metavar='default_cpus',
                        type=int)
    
    parser.add_argument('--default-memory',
                        default="5 GB",
                        help="Default amount of memory to use for each job. Default: '5 GB' ",
                        metavar='default_memory',
                        type=str)

    parser.add_argument('--R-cpus',
                        default=10,
                        help="Number of cpus to run R dada2 and DECIPHER. Default: 10",
                        metavar='R_cpus',
                        type=int)
    
    parser.add_argument('--R-memory',
                        default="20 GB",
                        help="Amount of memory R uses to run Dada2 and DECIPHER. Default: '20 GB' ",
                        metavar='R_memory',
                        type=str)

    parser.add_argument('--diversity-cpus',
                        default=5,
                        help="Number of cpus to run differential abundance and diversity analysis. Default: 5",
                        metavar='ancom_cpus',
                        type=int)
    
    parser.add_argument('--diversity-memory',
                        default="10 GB",
                        help="Amount of memory to use for each job. Default: '10 GB' ",
                        metavar='ancom_memory',
                        type=str)

    parser.add_argument('--cutadapt-memory',
                        default="10 GB",
                        help="Amount of memory used by cutadapt for read trimming. Default: '10 GB' ",
                        metavar='cutadapt_memory',
                        type=str)
        
    parser.add_argument('--conda-genelab',
                        metavar='/path/to/envs/genelab-utils',
                        default='null', 
                        help="Path to a conda environment containing genlab-utils. Default: null",
                        type=str)
    
    parser.add_argument('--conda-qc',
                        metavar='/path/to/envs/qc',
                        default='null', 
                        help="Path to a conda environment containing fastqc, multiqc, zip and python. Default: null",
                        type=str)

    parser.add_argument('--conda-R',
                        metavar='/path/to/envs/R',
                        default='null', 
                        help="Path to a conda environment containing R along with the packages decipher and biomformat installed. Default: null",
                        type=str)   

    parser.add_argument('--conda-cutadapt',
                        metavar='/path/to/envs/cutadapt',
                        default='null', 
                        help="Path to a conda environment containing cutadapt. Default: null",
                        type=str)

    parser.add_argument('--conda-diversity',
                        metavar='/path/to/envs/R_diversity',
                        default='null', 
                        help="Path to a conda environment containing R packages required for diversity and differential abundance testing. Default: null",
                        type=str)  

    parser.add_argument('--container-genelab',
                        default="olabiyi/genelab-utils:1.3.22",
                        help="Genelab utils container to be used to download raw sequences and retrieve sample metadata. Default: olabiyi/genelab-utils:1.3.22",
                        metavar='container_genelab',
                        type=str)
    
    parser.add_argument('--container-fastqc',
                        default="staphb/fastqc:0.12.1",
                        help="A docker container to be used to run fastqc. Default: staphb/fastqc:0.12.1",
                        metavar='container_fastqc',
                        type=str)
    
    parser.add_argument('--container-multiqc',
                        default="staphb/multiqc:1.19",
                        help="A docker container to be used to run multiqc. Default: staphb/multiqc:1.19",
                        metavar='container_multiqc',
                        type=str)

    parser.add_argument('--container-cutadapt',
                        default="zavolab/cutadapt:1.16",
                        help="A docker container to be used to run cutadapt. Default: zavolab/cutadapt:1.16",
                        metavar='container_cutadapt',
                        type=str)
    
    parser.add_argument('--container-dada',
                        default="olabiyi/r-dada-decipher-biomformat:1.0",
                        help="A docker container to be used to run dada2 and DECIPHER. Default: olabiyi/r-dada-decipher-biomformat:1.0",
                        metavar='container_dada',
                        type=str)

    parser.add_argument('--container-ancom',
                        default="quay.io/nasa_genelab/ancombc:2.6.0",
                        help="A docker container containing ancombc v2.6.0 . Default: quay.io/nasa_genelab/ancombc:2.6.0",
                        metavar='container_ancom',
                        type=str)

    parser.add_argument('--container-diversity',
                        default="quay.io/nasa_genelab/r-diversity:1.0",
                        help="A docker container to be used to run diversity analysis. Default: quay.io/nasa_genelab/r-diversity:1.0 ",
                        metavar='container_diversity',
                        type=str)
    
    parser.add_argument('--singularity-cacheDir',
                        default="singularity/",
                        help="A directory to store and retrieve singularity images. Default: singularity/",
                        metavar='singularity_directory',
                        type=str)

    parser.add_argument('--extra',
                        default="",
                        help="A comma separated string of extra arguement(s) to nextflow run e.g 'with-tower,name test'. \
                             Please run nextflow run -h for a full list of available options. Default: '' ",
                        metavar='extra_args',
                        type=str)
        
    # Check if no arguments were provided
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    try:
        args = parser.parse_args()
    except SystemExit:
        parser.print_help()
        sys.exit(1)

    # Test input requirement
    if (args.accession == "null" and args.input_file == "null"):
        print("""
              Please supply either an accession (OSD or Genelab number) or an input CSV file
              by passing either to the --accession or --input-file parameter, respectively.
              """)
        sys.exit(1)
    
    # Test input csv file
    if(args.input_file != "null"):
        # Test if input file exists
        if(not os.path.exists(args.input_file)):
           print(f"{args.input_file} does not exist. Please provide a correct path.")
           sys.exit(1)
        # Test primers
        if(args.F_primer == "null" or args.R_primer == "null"):
            print("""
                  When using a csv file as input (--input-file) to this workflow you must provide 
                  foward and reverse primer sequences. Please provide your forward 
                  and reverse primer sequences as arguements to the --F-primer 
                  and --R-primer parameters, respectively. 
                  """)
            sys.exit(1) 


    # Test profile 
    if(args.profile == "null"):
        print("Please provide a valid combination of profiles (conda, slurm, docker and singularity) to the --profile parameter.")
        sys.exit(1)

    # Create nextflow.config
    config_file = create_config(args.target_region, args.raw_R1_suffix, args.raw_R2_suffix, args.trim_primers,
                                args.input_file, args.min_cutadapt_len, args.primers_linked, args.discard_untrimmed, args.F_primer,
                                args.R_primer, args.left_trunc, args.right_trunc, args.left_maxEE, args.right_maxEE, 
                                args.concatenate_reads_only, args.conda_genelab, args.conda_qc, args.conda_R, args.conda_cutadapt,
                                args.conda_diversity, args.accession, args.assay_suffix, args.output_prefix, args.publishDir_mode,
                                args.primer_trimmed_R1_suffix, args.primer_trimmed_R2_suffix, args.filtered_R1_suffix, 
                                args.filtered_R2_suffix, args.output_dir, args.diff_abund_method, args.group_column, args.samples_column,
                                args.rarefaction_depth, args.errorStrategy, args.queueSize, args.default_cpus, args.default_memory,
                                args.cutadapt_memory, args.R_cpus, args.R_memory, args.diversity_cpus, args.diversity_memory,
                                args.container_genelab, args.container_fastqc, args.container_multiqc, args.container_cutadapt, 
                                args.container_dada, args.container_ancom, args.container_diversity, args.singularity_cacheDir,
                                args.remove_rare, args.prevalence_cutoff, args.library_cutoff)

    with open("nextflow.config", "w") as file:
        print(config_file, file=file)
        print("Nextflow workflow setup is complete.")
    
    # Run the nextflow workflow if --run is used
    if args.run:
        if(args.extra != ""):
            # Get extra arguement(s) to nextflow run
            extra = ""
            for opt in args.extra.split(","):
                extra += f"-{opt} "
            command = f"nextflow run main.nf -resume -profile {args.profile} {extra}"
        else:
            command = f"nextflow run main.nf -resume -profile {args.profile}"
        print(f"Running this nextflow command: {command}")
        subprocess.run(command, shell=True, check=True)

if __name__ == "__main__":
    main()
