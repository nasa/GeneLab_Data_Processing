process VV_RAW_READS {
  // Publish VV log
  publishDir "${ publishdir }",
    pattern:  "VV_log.csv" ,
    mode: params.publish_dir_mode,
    saveAs: { "VV_Logs/VV_log_${ task.process.replace(":","-") }${ params.assay_suffix }.csv" }
  // Publish VVed data
  publishDir "${ publishdir }",
    pattern: '00-RawData/**',
    mode: params.publish_dir_mode

  label 'VV'

  input:
    path(dp_tools__NF_RCP)
    val(publishdir)
    val(meta)
    path(runsheet)              // Runsheet
    path("INPUT/raw_fastq/*")   // Raw reads
    path("INPUT/raw_fastqc/*")  // Raw FastQC reports 
    path("INPUT/raw_multiqc/*") // Unzipped Raw MultiQC report

  output:
    path("00-RawData/Fastq"),                                                        emit: VVed_raw_reads
    path("00-RawData/FastQC_Reports/*{_fastqc.html,_fastqc.zip}"),                   emit: VVed_raw_fastqc
    //path("00-RawData/FastQC_Reports/raw_multiqc${params.assay_suffix}_report"),      emit: VVed_raw_unzipped_multiqc_report
    path("00-RawData/FastQC_Reports/raw_multiqc${params.assay_suffix}_report.zip"), emit: VVed_raw_zipped_multiqc_report
    path("VV_log.csv"),                                                              optional: params.skipVV, emit: log
    path("vv.log"),                                                                optional: params.skipVV, emit: log_txt

  script:
    """
    mv INPUT/* . || true
    vv.py --assay-type rnaseq \
    --assay-suffix ${params.assay_suffix} \
    --runsheet-path ${runsheet} \
    --outdir ${publishdir} \
    --paired-end ${meta.paired_end} \
    --mode microbes \
    --raw-fastq raw_fastq/ \
    --raw-fastqc raw_fastqc/ \
    --raw-multiqc raw_multiqc/ \
    --run-components raw_reads
    """
}

process VV_TRIMMED_READS {
  // Log publishing
  publishDir "${ publishdir }",
    pattern:  "VV_log.csv" ,
    mode: params.publish_dir_mode,
    saveAs: { "VV_Logs/VV_log_${ task.process.replace(":","-") }${ params.assay_suffix }.csv" }
  // V&V'ed data publishing
  publishDir "${ publishdir }",
    pattern: '01-TG_Preproc/**',
    mode: params.publish_dir_mode

  label 'VV'

  input:
    path(dp_tools__NF_RCP)
    val(publishdir)
    val(meta)
    path(runsheet)              // Runsheet
    path("INPUT/trimmed_fastq/*")   // Trimmed reads
    path("INPUT/trimmed_fastqc/*")  // Trimmed FastQC reports
    path("INPUT/trimmed_multiqc/*") // Trimmed reads multiqc zipped report
    path("INPUT/trimming_reports/*") // Trimming reports
    path("INPUT/trimming_multiqc/*") // Trimming multiqc zipped report

  output:
    path("01-TG_Preproc/Fastq"),                                                      emit: VVed_trimmed_reads
    path("01-TG_Preproc/FastQC_Reports/*{_fastqc.html,_fastqc.zip}"),                 emit: VVed_trimmed_fastqc
    path("01-TG_Preproc/FastQC_Reports/trimmed_multiqc${params.assay_suffix}_report.zip"), emit: VVed_trimmed_zipped_multiqc_report
    path("01-TG_Preproc/Trimming_Reports"),                                           emit: VVed_trimming_reports_all
    path("VV_log.csv"),                                                               optional: params.skipVV, emit: log
    path("vv.log"),                                                                   optional: params.skipVV, emit: log_txt

  script:
    """
    mv INPUT/* . || true
    vv.py --assay-type rnaseq \\
    --assay-suffix ${params.assay_suffix} \\
    --runsheet-path ${runsheet} \\
    --outdir ${publishdir} \\
    --paired-end ${meta.paired_end} \\
    --mode microbes \\
    --trimmed-fastq trimmed_fastq/ \\
    --trimmed-fastqc trimmed_fastqc/ \\
    --trimmed-multiqc trimmed_multiqc/ \\
    --trimming-reports trimming_reports/ \\
    --trimming-multiqc trimming_multiqc/ \\
    --run-components trimmed_reads
    """
}

process VV_BOWTIE2_ALIGNMENT {
  // Log publishing
  publishDir "${ publishdir }",
    pattern:  "VV_log.csv" ,
    mode: params.publish_dir_mode,
    saveAs: { "VV_Logs/VV_log_${ task.process.replace(":","-") }${ params.assay_suffix }.csv" }
  // V&V'ed data publishing
  publishDir "${ publishdir }",
    pattern: '02-Bowtie2_Alignment/',
    mode: params.publish_dir_mode

  label 'VV'

  input:
    path(dp_tools__NF_RCP)
    val(publishdir)
    val(meta)
    path(runsheet)                       // Runsheet
    path("INPUT/bowtie2-alignment-log/*") // (log files *.bowtie2.log)
    path("INPUT/bowtie2-alignment-unmapped/*") // (unmapped reads *.Unmapped.fastq.gz)
    path("INPUT/bowtie2-alignment-multiqc/*") // (zipped multiqc report)
    path("INPUT/bowtie2-alignment-sorted/*") // (sorted BAMs *_sorted.bam)
    path("INPUT/alignment-sorted-index/*") // (sorted BAM index files *_sorted.bam.bai)
    

  output:
    path("02-Bowtie2_Alignment/")
    path("VV_log.csv"), optional: params.skipVV, emit: log

  script:
    """
    mv INPUT/* . || true
    vv.py --assay-type rnaseq \
    --assay-suffix ${params.assay_suffix} \
    --runsheet-path ${runsheet} \
    --outdir ${publishdir} \
    --paired-end ${meta.paired_end} \
    --mode microbes \
    --bowtie2-alignment-log bowtie2-alignment-log/ \
    --bowtie2-alignment-unmapped bowtie2-alignment-unmapped/ \
    --bowtie2-alignment-multiqc bowtie2-alignment-multiqc/ \
    --bowtie2-alignment-sorted bowtie2-alignment-sorted/ \
    --bowtie2-alignment-sorted-index alignment-sorted-index/ \
    --run-components bowtie2_alignment
    """
} 

process VV_RSEQC {
  // Log publishing
  publishDir "${ publishdir }",
    pattern: "VV_log.csv",
    mode: params.publish_dir_mode,
    saveAs: { "VV_Logs/VV_log_${ task.process.replace(":","-") }${ params.assay_suffix }.csv" }
  // V&V'ed data publishing
  publishDir "${ publishdir }",
    pattern: 'RSeQC_Analyses/**',
    mode: params.publish_dir_mode

  label 'VV'

  input:
      path(dp_tools__NF_RCP)
      val(publishdir)
      val(meta)
      path(runsheet)
      path("INPUT/rseqc-logs/*") 
      path("INPUT/genebody-coverage/*")
      path("INPUT/infer-experiment/*")
      path("INPUT/inner-distance/*")
      path("INPUT/read-distribution/*")

  output:
      path("RSeQC_Analyses/**"), emit: rseqc_outputs
      path("VV_log.csv"), emit: log
      path("vv.log"), optional: params.skipVV, emit: log_txt

  script:
  def inner_distance = meta.paired_end ? "--inner-distance inner-distance/" : ""
  """
  # If single-end data, remove inner-distance directory since it's not applicable
  if [ "${meta.paired_end}" == "false" ]; then
    rm -rf INPUT/inner-distance
  fi
  
  mv INPUT/* . || true
  
  vv.py --assay-type rnaseq \
  --assay-suffix ${params.assay_suffix} \
  --runsheet-path ${runsheet} \
  --outdir ${publishdir} \
  --paired-end ${meta.paired_end} \
  --mode microbes \
  --run-components rseqc \
  --genebody-coverage genebody-coverage/ \
  --infer-experiment infer-experiment/ \
  ${inner_distance} \
  --read-distribution read-distribution/
  """
}

process VV_STAR_ALIGNMENTS {
  // Log publishing
  publishDir "${ publishdir }",
    pattern:  "VV_log.csv" ,
    mode: params.publish_dir_mode,
    saveAs: { "VV_Logs/VV_log_${ task.process.replace(":","-") }${ params.assay_suffix }.csv" }
  // V&V'ed data publishing
  publishDir "${ publishdir }",
    pattern: '02-STAR_Alignment/',
    mode: params.publish_dir_mode

  label 'VV'

  input:
    val(publishdir)
    path("VV_INPUT/Metadata/*")
    path("VV_INPUT/02-STAR_Alignment/*") // direct STAR alignment output
    path("VV_INPUT/02-STAR_Alignment/*") // STAR alignment counts tables
    path("VV_INPUT/02-STAR_Alignment/*") // zipped multiqc report 
    path("VV_INPUT/02-STAR_Alignment/*") // unzipped multiqc report
    path("VV_INPUT/02-STAR_Alignment/*") // reindexed, sorted bam/bed files
    path(dp_tools__NF_RCP)

  output:
    path("02-STAR_Alignment")
    path("VV_log.csv"), optional: params.skipVV, emit: log

  script:
    """
    # move from VV_INPUT to task directory
    # This allows detection as output files for publishing
    mv VV_INPUT/* . || true
    sort_into_subdirectories_by_sample.py 02-STAR_Alignment 02-STAR_Alignment '_*'

    # Run V&V unless user requests to skip V&V
    if ${ !params.skipVV } ; then
      dpt validation run ${dp_tools__NF_RCP} . Metadata/*_runsheet.csv \\
                          --data-asset-key-sets  \\
                            'STAR alignments' \\
                          --run-components \\
                            'STAR Alignments,STAR Alignments By Sample' \\
                          --max-flag-code ${ params.max_flag_code } \\
                          --output VV_log.csv
    fi
    """

}



process VV_RSEM_COUNTS {
  // Log publishing
  publishDir "${ publishdir }",
    pattern:  "VV_log.csv" ,
    mode: params.publish_dir_mode,
    saveAs: { "VV_Logs/VV_log_${ task.process.replace(":","-") }${ params.assay_suffix }.csv" }
  // V&V'ed data publishing
  publishDir "${ publishdir }",
    pattern: '03-RSEM_Counts/',
    mode: params.publish_dir_mode

  label 'VV'

  input:
    val(publishdir)
    path("VV_INPUT/Metadata/*")
    path("VV_INPUT/03-RSEM_Counts/*") // RSEM sample wise output
    path("VV_INPUT/03-RSEM_Counts/*") // RSEM dataset output
    path("VV_INPUT/03-RSEM_Counts/*") // zipped multiqc report 
    path("VV_INPUT/03-RSEM_Counts/*") // unzipped multiqc report
    path(dp_tools__NF_RCP)
    

  output:
    path("03-RSEM_Counts")
    path("VV_log.csv"), optional: params.skipVV, emit: log
  
  script:
    """
    # move from VV_INPUT to task directory
    # This allows detection as output files for publishing
    mv VV_INPUT/* . || true

    # Run V&V unless user requests to skip V&V
    if ${ !params.skipVV } ; then
      dpt validation run ${dp_tools__NF_RCP} . Metadata/*_runsheet.csv \\
                          --data-asset-key-sets  \\
                            'RSEM counts' \\
                          --run-components \\
                            'RSEM Counts' \\
                          --max-flag-code ${ params.max_flag_code } \\
                          --output VV_log.csv
    fi
    """
}

process VV_DESEQ2_ANALYSIS {
  // Log publishing
  publishDir "${ publishdir }",
    pattern:  "VV_log.csv" ,
    mode: params.publish_dir_mode,
    saveAs: { "VV_Logs/VV_log_${ task.process.replace(":","-") }${ params.assay_suffix }.csv" }
  // V&V'ed data publishing
  publishDir "${ publishdir }",
    pattern: '{04-DESeq2_NormCounts,05-DESeq2_DGE}',
    mode: params.publish_dir_mode

  label 'VV'

  input:
    val(publishdir)
    val(meta)
    path("VV_INPUT/Metadata/*")
    path("VV_INPUT/03-RSEM_Counts/*") // RSEM dataset output
    path("VV_INPUT/03-RSEM_Counts/*") // zipped multiqc report 
    path("VV_INPUT/03-RSEM_Counts/*") // unzipped multiqc report
    path("VV_INPUT/04-DESeq2_NormCounts/*") // norm counts files
    path("VV_INPUT/05-DESeq2_DGE/*") // dge files
    path("VV_INPUT/04-DESeq2_NormCounts/*") // ERCC norm counts files
    path("VV_INPUT/05-DESeq2_DGE/ERCC_NormDGE/*") // ERCC dge files
    path(dp_tools__NF_RCP)

  output:
    path("04-DESeq2_NormCounts")
    path("05-DESeq2_DGE")
    path("VV_log.csv"), optional: params.skipVV, emit: log
  
  script:
    """
    # move from VV_INPUT to task directory
    # This allows detection as output files for publishing
    mv VV_INPUT/* . || true

    # Run V&V unless user requests to skip V&V
    if ${ !params.skipVV } ; then
      dpt validation run ${dp_tools__NF_RCP} . Metadata/*_runsheet.csv \\
                          --data-asset-key-sets  \\
                            'RSEM Output,DGE Output${ meta.has_ercc ? ",ERCC DGE Output" : ''}' \\
                          --run-components \\
                            'DGE Metadata${ meta.has_ercc ? ",DGE Metadata ERCC" : '' },DGE Output${ meta.has_ercc ? ",DGE Output ERCC" : '' }' \\
                          --max-flag-code ${ params.max_flag_code } \\
                          --output VV_log.csv
    fi

    # Remove all placeholder files and empty directories to prevent publishing
    find . -type f,l -name *.placeholder -delete
    find . -empty -type d -delete
    """
}

process VV_CONCAT_FILTER {
  publishDir "${ params.outputDir }/${ params.gldsAccession }/VV_Logs",
    mode: params.publish_dir_mode

  label 'VV'

  input:
    path("VV_in.csv")

  output:
    tuple path("VV_log_final_GLbulkRNAseq.csv"), path("VV_log_final_only_issues_GLbulkRNAseq.csv")

  script:
    """
    concat_logs.py
    filter_to_only_issues.py
    """
}