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
    path("INPUT/00-RawData/MultiQC_Reports/*")  // Zipped Raw MultiQC data directory
    path("INPUT/00-RawData/MultiQC_Reports/*")  // Raw MultiQC HTML report

  output:
    path("00-RawData/Fastq"),                                                        emit: VVed_raw_reads
    path("00-RawData/FastQC_Reports/*{_fastqc.html,_fastqc.zip}"),                   emit: VVed_raw_fastqc
    path("00-RawData/MultiQC_Reports/raw_multiqc${params.assay_suffix}_data.zip"),  emit: VVed_raw_zipped_multiqc_data
    path("00-RawData/MultiQC_Reports/raw_multiqc${params.assay_suffix}.html"),       emit: VVed_raw_multiqc_html
    path("VV_log.csv"),                                                              optional: params.skip_vv, emit: log
    path("versions.yml"), emit: versions

  script:
    """
    # Move inputs into task directory to be detected as outputs
    mv INPUT/* . || true

    # Run V&V unless user requests to skip V&V
    if ${ !params.skip_vv } ; then
      vv_raw_reads.py --runsheet ${runsheet} --outdir .
    fi

    echo '"${task.process}":' > versions.yml
    echo "    dp_tools: \$(pip show dp_tools | grep Version | sed 's/Version: //')" >> versions.yml
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
    path("INPUT/01-TG_Preproc/MultiQC_Reports/*") // Trimmed reads multiqc zipped data directory
    path("INPUT/01-TG_Preproc/MultiQC_Reports/*") // Trimmed reads multiqc HTML report
    path("INPUT/01-TG_Preproc/Trimming_Reports/*") // Trimming reports

  output:
    path("01-TG_Preproc/Fastq"),                                                              emit: VVed_trimmed_reads
    path("01-TG_Preproc/FastQC_Reports/*{_fastqc.html,_fastqc.zip}"),                         emit: VVed_trimmed_fastqc
    path("01-TG_Preproc/MultiQC_Reports/trimmed_multiqc${params.assay_suffix}_data.zip"),    emit: VVed_trimmed_zipped_multiqc_data
    path("01-TG_Preproc/MultiQC_Reports/trimmed_multiqc${params.assay_suffix}.html"),       emit: VVed_trimmed_multiqc_html
    path("01-TG_Preproc/Trimming_Reports/*trimming_report.txt"), emit: VVed_trimming_reports

    path("VV_log.csv"),                                                                       optional: params.skip_vv, emit: log

  script:
    """
    mv INPUT/* . || true

    # Run V&V unless user requests to skip V&V
    if ${ !params.skip_vv } ; then
      vv_trimmed_reads.py --runsheet ${runsheet} --outdir .
    fi
    """
}

process VV_BOWTIE2_ALIGNMENT {
  // Log publishing
  publishDir "${ publishdir }",
    pattern:  "VV_log.csv" ,
    mode: params.publish_dir_mode,
    saveAs: { "VV_Logs/VV_log_ALIGNMENT${ params.assay_suffix }.csv" }
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
    path("INPUT/02-Bowtie2_Alignment/MultiQC_Reports/*") // (zipped multiqc data directory)
    path("INPUT/02-Bowtie2_Alignment/MultiQC_Reports/*") // (multiqc HTML report)
    
  output:
    path("02-Bowtie2_Alignment/**")
    path("VV_log.csv"), optional: params.skip_vv, emit: log

  script:
    """
    mv INPUT/* . || true

    sort_into_subdirectories.py --from 02-Bowtie2_Alignment --to 02-Bowtie2_Alignment --runsheet ${runsheet}

    # Run V&V unless user requests to skip V&V
    if ${ !params.skip_vv } ; then
      vv_bowtie2_alignment.py --runsheet ${runsheet} --outdir .
    fi
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
      path("INPUT/RSeQC_Analyses/MultiQC_Reports/*") // genebody coverage multiqc data zip
      path("INPUT/RSeQC_Analyses/MultiQC_Reports/*") // genebody coverage multiqc HTML report
      path("INPUT/RSeQC_Analyses/MultiQC_Reports/*") // infer experiment multiqc data zip
      path("INPUT/RSeQC_Analyses/MultiQC_Reports/*") // infer experiment multiqc HTML report
      path("INPUT/RSeQC_Analyses/MultiQC_Reports/*") // inner distance multiqc data zip
      path("INPUT/RSeQC_Analyses/MultiQC_Reports/*") // inner distance multiqc HTML report
      path("INPUT/RSeQC_Analyses/MultiQC_Reports/*") // read distribution multiqc data zip
      path("INPUT/RSeQC_Analyses/MultiQC_Reports/*") // read distribution multiqc HTML report

  output:
      path("RSeQC_Analyses/**"), emit: rseqc_outputs
      path("VV_log.csv"), optional: params.skip_vv, emit: log

  script:
    """
    mv INPUT/* . || true

    # Create all necessary RSeQC subdirectories
    mkdir -p RSeQC_Analyses/02_geneBody_coverage
    mkdir -p RSeQC_Analyses/03_infer_experiment
    mkdir -p RSeQC_Analyses/05_read_distribution
    
    # Create inner distance directory only for paired-end data
    ${ meta.paired_end ? 'mkdir -p RSeQC_Analyses/04_inner_distance' : '# Skip inner distance for single-end data' }

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
    if ${ !params.skip_vv } ; then
      vv_rseqc.py --runsheet ${runsheet} --outdir .
    fi

    # Remove all placeholder files and empty directories to prevent publishing
    #  (Removes Inner Distance output placeholder files for single-end runs)
    find RSeQC_Analyses -type f,l -name PLACEHOLDER* -delete
    find RSeQC_Analyses -empty -type d -delete
    """
}

process VV_FEATURECOUNTS {
  // Log publishing
  publishDir "${ publishdir }",
    pattern: "VV_log.csv",
    mode: params.publish_dir_mode,
    saveAs: { "VV_Logs/VV_log_COUNTS${ params.assay_suffix }.csv" }
  publishDir "${ publishdir }",
    pattern: '03-FeatureCounts/**',
    mode: params.publish_dir_mode

  label 'VV'

  input:
    path(dp_tools__NF_RCP)
    val(publishdir)
    val(meta)
    path(runsheet)
    path("INPUT/03-FeatureCounts/*") // featurecounts counts
    path("INPUT/03-FeatureCounts/*") // featurecounts counts (rRNArm)
    path("INPUT/03-FeatureCounts/*") // featurecounts summary
    path("INPUT/03-FeatureCounts/*") // featurecounts num non zero genes
    path("INPUT/03-FeatureCounts/MultiQC_Reports/*") // featurecounts multiqc zipped data directory
    path("INPUT/03-FeatureCounts/MultiQC_Reports/*") // featurecounts multiqc HTML report

  output:
    path("03-FeatureCounts/**")
    path("VV_log.csv"), optional: params.skip_vv, emit: log

  script:
  """
  mv INPUT/* . || true

  # Run V&V unless user requests to skip V&V
  if ${ !params.skip_vv } ; then
    vv_featurecounts.py --runsheet ${runsheet} --outdir .
  fi
  """
}

process VV_DGE_DESEQ2 {
  // Log publishing
  publishDir "${ publishdir }",
    pattern:  "VV_log.csv" ,
    mode: params.publish_dir_mode,
    saveAs: { "VV_Logs/VV_log_DESEQ2_ANALYSIS${ params.assay_suffix }.csv" }
  // V&V'ed data publishing
  publishDir "${ publishdir }",
    pattern: '{04-*,05-*}',
    mode: params.publish_dir_mode

  label 'VV'

  input:
    path(dp_tools__NF_RCP)
    val(publishdir)
    val(meta)
    path(runsheet)
    path("INPUT/04-DESeq2_NormCounts/*") // unnormed counts, normed counts
    path("INPUT/04-DESeq2_NormCounts/*") // vst norm counts
    path("INPUT/05-DESeq2_DGE/*") // sample table
    path("INPUT/05-DESeq2_DGE/*") // contrasts
    path("INPUT/05-DESeq2_DGE/*") // annotated dge table
    path("INPUT/04-DESeq2_NormCounts/*") // (rrna rm) unnormed counts, normed counts
    path("INPUT/04-DESeq2_NormCounts/*") // (rrna rm) vst norm counts
    path("INPUT/05-DESeq2_DGE/*") // (rrna rm) annotated dge table

  output:
    path("04-*")
    path("05-*")
    path("VV_log.csv"), optional: params.skip_vv, emit: log

  script:
  def mode_params = params.mode == "microbes" ? "--mode microbes" : ""
  """
  mv INPUT/* . || true

  if ${ !params.skip_vv } ; then
    vv_dge_deseq2.py --runsheet ${runsheet} --outdir . ${mode_params}
  fi
  """
}

process VV_STAR_ALIGNMENT {
  // Log publishing
  publishDir "${ publishdir }",
    pattern:  "VV_log.csv" ,
    mode: params.publish_dir_mode,
    saveAs: { "VV_Logs/VV_log_ALIGNMENT${ params.assay_suffix }.csv" }
  // V&V'ed data publishing
  publishDir "${ publishdir }",
    pattern: '02-STAR_Alignment/**',
    mode: params.publish_dir_mode

  label 'VV'

  input:
    path(dp_tools__NF_RCP)
    val(publishdir)
    val(meta)
    path(runsheet)
    path("INPUT/02-STAR_Alignment/*") // direct STAR alignment output
    path("INPUT/02-STAR_Alignment/*") // STAR alignment counts tables
    path("INPUT/02-STAR_Alignment/*") // reindexed, sorted bam/bed files
    path("INPUT/02-STAR_Alignment/MultiQC_Reports/*") // multiqc HTML report
    path("INPUT/02-STAR_Alignment/MultiQC_Reports/*") // zipped multiqc data directory

  output:
    path("02-STAR_Alignment/**")
    path("VV_log.csv"), optional: params.skip_vv, emit: log
  script:
    """
    mv INPUT/* . || true

    sort_into_subdirectories.py --from 02-STAR_Alignment --to 02-STAR_Alignment --runsheet ${runsheet} --glob '_*'

    # Run V&V unless user requests to skip V&V
    if ${ !params.skip_vv } ; then
      vv_star_alignment.py --runsheet ${runsheet} --outdir .
    fi
    """
}

process VV_RSEM_COUNTS {
  // Log publishing
  publishDir "${ publishdir }",
    pattern:  "VV_log.csv" ,
    mode: params.publish_dir_mode,
    saveAs: { "VV_Logs/VV_log_COUNTS${ params.assay_suffix }.csv" }
  // V&V'ed data publishing
  publishDir "${ publishdir }",
    pattern: '03-RSEM_Counts/**',
    mode: params.publish_dir_mode

  label 'VV'

  input:
    path(dp_tools__NF_RCP)
    val(publishdir)
    val(meta)
    path(runsheet)
    path("INPUT/03-RSEM_Counts/*") // RSEM sample wise output
    path("INPUT/03-RSEM_Counts/*") // RSEM sample.genes.results (rRNArm)
    path("INPUT/03-RSEM_Counts/*") // RSEM dataset output
    path("INPUT/03-RSEM_Counts/MultiQC_Reports/*") // zipped multiqc data directory
    path("INPUT/03-RSEM_Counts/MultiQC_Reports/*") // multiqc HTML report
    
  output:
    path("03-RSEM_Counts/**")
    path("VV_log.csv"), optional: params.skip_vv, emit: log
  
  script:
    """
    # move from VV_INPUT to task directory
    # This allows detection as output files for publishing
    mv INPUT/* . || true

    # Run the simplified sorting script to organize files into sample subdirectories
    sort_into_subdirectories.py --from 03-RSEM_Counts --to 03-RSEM_Counts --runsheet ${runsheet} --glob '.*'
    
    # Explicitly sort the rRNArm files into sample directories
    sort_into_subdirectories.py --from 03-RSEM_Counts --to 03-RSEM_Counts --runsheet ${runsheet} --glob '*_rRNArm*'

    # Run V&V unless user requests to skip V&V
    if ${ !params.skip_vv } ; then
      vv_rsem_counts.py --runsheet ${runsheet} --outdir .
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
    path("VV_INPUT/03-RSEM_Counts/*") // zipped multiqc data directory
    path("VV_INPUT/03-RSEM_Counts/*") // multiqc HTML report
    path("VV_INPUT/04-DESeq2_NormCounts/*") // norm counts files
    path("VV_INPUT/05-DESeq2_DGE/*") // dge files
    path("VV_INPUT/04-DESeq2_NormCounts/*") // ERCC norm counts files
    path("VV_INPUT/05-DESeq2_DGE/ERCC_NormDGE/*") // ERCC dge files
    path(dp_tools__NF_RCP)

  output:
    path("04-DESeq2_NormCounts")
    path("05-DESeq2_DGE")
    path("VV_log.csv"), optional: params.skip_vv, emit: log
  
  script:
    """
    # move from VV_INPUT to task directory
    # This allows detection as output files for publishing
    mv VV_INPUT/* . || true

    # Run V&V unless user requests to skip V&V
    if ${ !params.skip_vv } ; then
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
  publishDir "${publishdir}/VV_Logs",
    mode: params.publish_dir_mode

  label 'VV'

  input:
    val(publishdir)
    path("VV_in.csv")

  output:
    tuple path("VV_log_final${params.assay_suffix}.csv"), path("VV_log_final_only_issues${params.assay_suffix}.csv")

  script:
    """
    concat_logs.py
    filter_to_only_issues.py --assay_suffix ${params.assay_suffix}
    """
}