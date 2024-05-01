#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/**************************************************************************************** 
*********************  Sequence quality assessment and control processes ****************
****************************************************************************************/

// a 2-column (single-end) or 3-column (paired-end) file
params.prefix = "raw" // "filetered"
params.csv_file = "file.csv" 
params.swift_1S = false
params.adapters = "${baseDir}/config/bbtools_dapters.fa"
params.multiqc_config = "config/multiqc.config"

process FASTQC {
  // FastQC performed on reads
  tag "Running fastqc on ${sample_id}"

  input:
    tuple val(sample_id), path(reads), val(isPaired)
  output:
    tuple path("*.html"), path("*.zip")
  script:
    """
    fastqc -o . \\
     -t ${task.cpus} \\
      ${reads}
    """
}

process MULTIQC {

  tag "Running multiqc on the ${prefix} files.."

  input:
    val(prefix)   
    path(multiqc_config)
    path(files)
  output:
    path("${params.additional_filename_prefix}${prefix}_multiqc${params.assay_suffix}_report.zip")
  script:
    """
      multiqc -q --filename ${params.additional_filename_prefix}${prefix}_multiqc \\
              --force --cl-config 'max_table_rows: 99999999' \\
              --interactive --config ${multiqc_config} \\
              --outdir ${params.additional_filename_prefix}${prefix}_multiqc_report  ${files} > /dev/null 2>&1

            
      # zipping and removing unzipped dir
      zip -q -r \\
           ${params.additional_filename_prefix}${prefix}_multiqc${params.assay_suffix}_report.zip \\
           ${params.additional_filename_prefix}${prefix}_multiqc_report

    """
  }


//  This process runs quality filtering/trimming on input fastq files.
process BBDUK {


    tag "Quality filtering ${sample_id}-s reads.."
    beforeScript "chmod +x ${baseDir}/bin/*"

    input:
        tuple val(sample_id), path(reads), val(isPaired)
        path(adapters)
    output:
        tuple val(sample_id), path("*${params.filtered_suffix}"), val(isPaired)
    script:
    def isSwift = params.swift_1S ? 't' : 'f'
    """
    if [ ${isPaired} == true ];then

        bbduk.sh in=${reads[0]} in2=${reads[1]} \\
                 out1=${sample_id}${params.filtered_R1_suffix} \\
                 out2=${sample_id}${params.filtered_R2_suffix} \\
                 ref=${adapters} \\
                 ktrim=l k=17 ftm=5 qtrim=rl \\
                 trimq=10 mlf=0.5 maxns=0 swift=${isSwift} 
    else
  
        bbduk.sh in=${reads[0]} out1=${sample_id}${params.filtered_suffix} \\
                  ref=${adapters} \\
				  ktrim=l k=17 ftm=5 qtrim=rl \\
                  trimq=10 mlf=0.5 maxns=0 swift=${isSwift}

    fi    
    """
}


workflow quality_check {

    take:
    prefix_ch
    multiqc_config
    reads_ch
    

    main:
    fastqc_ch = FASTQC(reads_ch).flatten().collect()
    MULTIQC(prefix_ch, multiqc_config, fastqc_ch)
}

workflow {

        Channel.fromPath(params.csv_file)
               .splitCsv()
               .map{ row -> row.paired ? tuple( "${row.sample_id}", [file("${row.forward}"), file("${row.reverse}")], row.paired) : 
                                         tuple( "${row.sample_id}", [file("${row.forward}")], row.paired)}
               .set{reads_ch}   

    res_ch = quality_check(Channel.of(params.prefix), params.multiqc_config, reads_ch)
    BBDUK(reads_ch)
}
