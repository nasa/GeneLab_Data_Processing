process VV_RAW_READS {
  // Arrange Raw reads, raw reads FastQC reports, Raw reads MultiQC report into the expected folder
  //  structure under 00-RawData
  // Run vv_raw_reads.py to generate VV_log.csv and publish it to the VV_Logs directory
  label 'VV'

  // Publish VVed data
  publishDir "${ publishdir }",
    pattern: '00-RawData/**',
    mode: params.publish_dir_mode
  // Publish VV log
  publishDir "${ publishdir }",
    pattern:  "VV_log.csv" ,
    mode: params.publish_dir_mode,
    saveAs: { "VV_Logs/VV_log_${ task.process.tokenize(':').last() }${ params.assay_suffix }.csv" }

  input:
    path(dp_tools__NF_RCP)
    val(publishdir)
    val(meta)
    path(runsheet)                             // Runsheet
    path("INPUT/00-RawData/Fastq/*")           // Raw reads
    path("INPUT/00-RawData/FastQC_Reports/*")  // Raw FastQC reports 
    path("INPUT/00-RawData/FastQC_Reports/*")  // Zipped Raw MultiQC report

  output:
    path("00-RawData/Fastq"),                                                        emit: VVed_raw_reads
    path("00-RawData/FastQC_Reports/*{_fastqc.html,_fastqc.zip}"),                   emit: VVed_raw_fastqc
    path("00-RawData/FastQC_Reports/raw_multiqc${params.assay_suffix}_report.zip"),  emit: VVed_raw_zipped_multiqc_report
    path("VV_log.csv"),                                                              optional: params.skipVV, emit: log

  script:
    """
    # Move inputs into task directory to be detected as outputs
    mv INPUT/* . || true

    vv_raw_reads.py --runsheet ${runsheet} --outdir .
    """
}

process VV_TRIMMED_READS {
  label 'VV'

  // Log publishing
  publishDir "${ publishdir }",
    pattern:  "VV_log.csv" ,
    mode: params.publish_dir_mode,
    saveAs: { "VV_Logs/VV_log_${ task.process.tokenize(':').last() }${ params.assay_suffix }.csv" }
  // V&V'ed data publishing
  publishDir "${ publishdir }",
    pattern: '01-TG_Preproc/**',
    mode: params.publish_dir_mode

  input:
    path(dp_tools__NF_RCP)
    val(publishdir)
    val(meta)
    path(runsheet)              // Runsheet
    path("INPUT/01-TG_Preproc/Fastq/*")   // Trimmed reads
    path("INPUT/01-TG_Preproc/FastQC_Reports/*")  // Trimmed FastQC reports
    path("INPUT/01-TG_Preproc/FastQC_Reports/*") // Trimmed reads multiqc zipped report
    path("INPUT/01-TG_Preproc/Trimming_Reports/*") // Trimming reports
    path("INPUT/01-TG_Preproc/Trimming_Reports/*") // Trimming multiqc zipped report

  output:
    path("01-TG_Preproc/Fastq"),                                                      emit: VVed_trimmed_reads
    path("01-TG_Preproc/FastQC_Reports/*{_fastqc.html,_fastqc.zip}"),                 emit: VVed_trimmed_fastqc
    path("01-TG_Preproc/FastQC_Reports/trimmed_multiqc${params.assay_suffix}_report.zip"), emit: VVed_trimmed_zipped_multiqc_report
    path("01-TG_Preproc/Trimming_Reports/*trimming_report.txt"), emit: VVed_trimming_reports
    path("01-TG_Preproc/Trimming_Reports/trimming_multiqc${params.assay_suffix}_report.zip"), emit: VVed_trimming_zipped_multiqc_report
    path("VV_log.csv"),                                                               optional: params.skipVV, emit: log

  script:
    """
    mv INPUT/* . || true
    
    vv_trimmed_reads.py --runsheet ${runsheet} --outdir .
    """
}

process VV_BOWTIE2_ALIGNMENT {
  // Log publishing
  publishDir "${ publishdir }",
    pattern:  "VV_log.csv" ,
    mode: params.publish_dir_mode,
    saveAs: { "VV_Logs/VV_log_${ task.process.tokenize(':').last() }${ params.assay_suffix }.csv" }
  // V&V'ed data publishing
  publishDir "${ publishdir }",
    pattern: '02-Bowtie2_Alignment/**',
    mode: params.publish_dir_mode

  label 'VV'

  input:
    path(dp_tools__NF_RCP)
    val(publishdir)
    val(meta)
    path(runsheet)                       // Runsheet
    path("INPUT/02-Bowtie2_Alignment/*") // (log files *.bowtie2.log)
    path("INPUT/02-Bowtie2_Alignment/*") // (unmapped reads *.Unmapped.fastq.gz)
    path("INPUT/02-Bowtie2_Alignment/*") // (sorted BAMs *_sorted.bam)
    path("INPUT/02-Bowtie2_Alignment/*") // (sorted BAM index files *_sorted.bam.bai)
    path("INPUT/02-Bowtie2_Alignment/*") // (zipped multiqc report)
    
  output:
    path("02-Bowtie2_Alignment/**")
    path("VV_log.csv"), optional: params.skipVV, emit: log

  script:
    """
    mv INPUT/* . || true

    sort_into_subdirectories.py --from 02-Bowtie2_Alignment --to 02-Bowtie2_Alignment --runsheet ${runsheet}

    vv_bowtie2_alignment.py --runsheet ${runsheet} --outdir .
    """
} 

process VV_RSEQC {
  // Log publishing
  publishDir "${ publishdir }",
    pattern: "VV_log.csv",
    mode: params.publish_dir_mode,
    saveAs: { "VV_Logs/VV_log_${ task.process.tokenize(':').last() }${ params.assay_suffix }.csv" }
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
      path("INPUT/RSeQC_Analyses/*")                      // direct logs
      path("INPUT/RSeQC_Analyses/02_geneBody_coverage/*") // genebody coverage multiqc
      path("INPUT/RSeQC_Analyses/03_infer_experiment/*")  // infer experiment multiqc
      path("INPUT/RSeQC_Analyses/04_inner_distance/*")    // inner distance multiqc
      path("INPUT/RSeQC_Analyses/05_read_distribution/*") // read distribution multiqc

  output:
      path("RSeQC_Analyses/**"), emit: rseqc_outputs
      path("VV_log.csv"), emit: log

  script:
    """
    mv INPUT/* . || true

    sort_into_subdirectories.py --from RSeQC_Analyses --to RSeQC_Analyses/02_geneBody_coverage --runsheet ${runsheet} --glob '.geneBodyCoverage.txt'
    sort_into_subdirectories.py --from RSeQC_Analyses --to RSeQC_Analyses/02_geneBody_coverage --runsheet ${runsheet} --glob '.geneBodyCoverage.curves.pdf'
    sort_into_subdirectories.py --from RSeQC_Analyses --to RSeQC_Analyses/02_geneBody_coverage --runsheet ${runsheet} --glob '.geneBodyCoverage.r'
    
    # These are not in sub directories
    mv RSeQC_Analyses/*.infer_expt.out RSeQC_Analyses/03_infer_experiment
    
    ${ meta.paired_end ? '' : '# Only for Paired end datasets: '} sort_into_subdirectories.py --from RSeQC_Analyses --to RSeQC_Analyses/04_inner_distance --runsheet ${runsheet} --glob '.inner_distance_freq.txt'
    ${ meta.paired_end ? '' : '# Only for Paired end datasets: '} sort_into_subdirectories.py --from RSeQC_Analyses --to RSeQC_Analyses/04_inner_distance --runsheet ${runsheet} --glob '.inner_distance_plot.pdf'
    ${ meta.paired_end ? '' : '# Only for Paired end datasets: '} sort_into_subdirectories.py --from RSeQC_Analyses --to RSeQC_Analyses/04_inner_distance --runsheet ${runsheet} --glob '.inner_distance_plot.r'
    ${ meta.paired_end ? '' : '# Only for Paired end datasets: '} sort_into_subdirectories.py --from RSeQC_Analyses --to RSeQC_Analyses/04_inner_distance --runsheet ${runsheet} --glob '.inner_distance.txt'
    
    # These are not in sub directories
    mv RSeQC_Analyses/*.read_dist.out RSeQC_Analyses/05_read_distribution

    # Run V&V unless user requests to skip V&V
    touch VV_log.csv

    # Remove all placeholder files and empty directories to prevent publishing
    #  (Removes Inner Distance output placeholder files for single-end runs)
    find RSeQC_Analyses -type f,l -name PLACEHOLDER -delete
    find RSeQC_Analyses -empty -type d -delete
    """
}

process VV_FEATURECOUNTS {
  // Log publishing
  publishDir "${ publishdir }",
    pattern: "VV_log.csv",
    mode: params.publish_dir_mode,
    saveAs: { "VV_Logs/VV_log_${ task.process.tokenize(':').last() }${ params.assay_suffix }.csv" }
  publishDir "${ publishdir }",
    pattern: '03-FeatureCounts/**',
    mode: params.publish_dir_mode

  label 'VV'

  input:
    path(dp_tools__NF_RCP)
    val(publishdir)
    val(meta)
    path(runsheet)
    path("INPUT/fc-counts/*") // featurecounts counts
    path("INPUT/fc-summary/*") // featurecounts summary
    path("INPUT/fc-counts-rrnarm/*") // featurecounts counts_rrnarm
    path("INPUT/fc-multiqc/*") // featurecounts multiqc zipped report

  output:
    path("03-FeatureCounts/**")
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
  --run-components featurecounts \
  --featurecounts-summary fc-summary/ \
  --featurecounts-counts fc-counts/ \
  --featurecounts-counts-rrnarm fc-counts-rrnarm/ \
  --featurecounts-multiqc fc-multiqc/
  """
}

process VV_DGE_MICROBES {
  // Log publishing
  publishDir "${ publishdir }",
    pattern:  "VV_log.csv" ,
    mode: params.publish_dir_mode,
    saveAs: { "VV_Logs/VV_log_${ task.process.tokenize(':').last() }${ params.assay_suffix }.csv" }
  // V&V'ed data publishing
  publishDir "${ publishdir }",
    pattern: '{04-DESeq2_NormCounts,05-DESeq2_DGE,04-DESeq2_NormCounts_rRNArm,05-DESeq2_DGE_rRNArm}',
    mode: params.publish_dir_mode

  label 'VV'

  input:
    path(dp_tools__NF_RCP)
    val(publishdir)
    val(meta)
    path(runsheet)
    path("INPUT/counts-norm-microbes/*") // unnormed counts, normed counts
    path("INPUT/counts-norm-microbes/*") // vst norm counts
    path("INPUT/dge-microbes/*") // sample table
    path("INPUT/dge-microbes/*") // contrasts
    path("INPUT/dge-microbes/*") // annotated dge table
    path("INPUT/counts-norm-microbes-rrnarm/*") // (rrna rm) unnormed counts, normed counts
    path("INPUT/counts-norm-microbes-rrnarm/*") // (rrna rm) vst norm counts
    path("INPUT/dge-microbes-rrnarm/*") // (rrna rm) sample table
    path("INPUT/dge-microbes-rrnarm/*") // (rrna rm) contrasts
    path("INPUT/dge-microbes-rrnarm/*") // (rrna rm) annotated dge table

  output:
    path("04-DESeq2_NormCounts")
    path("05-DESeq2_DGE")
    path("04-DESeq2_NormCounts_rRNArm")
    path("05-DESeq2_DGE_rRNArm")
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
  --run-components dge_microbes \
  --counts-norm-microbes counts-norm-microbes/ \
  --dge-microbes dge-microbes/ \
  --counts-norm-microbes-rrnarm counts-norm-microbes-rrnarm/ \
  --dge-microbes-rrnarm dge-microbes-rrnarm/ \
  ${meta.has_ercc ? '--has-ercc' : ''}
  """
}

process VV_STAR_ALIGNMENTS {
  // Log publishing
  publishDir "${ publishdir }",
    pattern:  "VV_log.csv" ,
    mode: params.publish_dir_mode,
    saveAs: { "VV_Logs/VV_log_${ task.process.tokenize(':').last() }${ params.assay_suffix }.csv" }
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
    sort_into_subdirectories.py --from 02-STAR_Alignment --to 02-STAR_Alignment --runsheet Metadata/*_runsheet.csv --glob '_*'

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
    saveAs: { "VV_Logs/VV_log_${ task.process.tokenize(':').last() }${ params.assay_suffix }.csv" }
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
    saveAs: { "VV_Logs/VV_log_${ task.process.tokenize(':').last() }${ params.assay_suffix }.csv" }
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