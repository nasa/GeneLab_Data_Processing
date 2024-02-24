/*
 * Processes related to sequence quality assessment,
 *   quality control (e.g. trimming).
 */

process FASTQC {
  // FastQC performed on reads
  tag "Sample: ${ meta.id }"

  input:
    tuple val(meta), path(reads)

  output:
    tuple val(meta), path("${ meta.id }*.html"), path("${ meta.id }*.zip"), emit: fastqc
    path("versions.txt"), emit: version

  script:
    """
    fastqc -o . \
     -t $task.cpus \
      $reads

    fastqc -v > versions.txt
    """
}

process MULTIQC {
  label "fastLocal"

  input:
    path("samples.txt")
    path("mqc_in/*") // any number of multiqc compatible files
    path(multiqc_config)

  output:
    path("${ params.MQCLabel }_multiqc_GLbulkRNAseq_report/${ params.MQCLabel }_multiqc_GLbulkRNAseq.html"), emit: html
    path("${ params.MQCLabel }_multiqc_GLbulkRNAseq_report/${ params.MQCLabel }_multiqc_GLbulkRNAseq_data"), emit: data
    path("${ params.MQCLabel }_multiqc_GLbulkRNAseq_report.zip"), emit: zipped_report
    path("${ params.MQCLabel }_multiqc_GLbulkRNAseq_report"), emit: unzipped_report
    path("versions.txt"), emit: version

  script:
    config_arg =  multiqc_config.name != "NO_FILE" ? "--config ${ multiqc_config }" : ""
    """
    multiqc --sample-names samples.txt  \
            --interactive -o ${ params.MQCLabel }_multiqc_GLbulkRNAseq_report \
            -n ${ params.MQCLabel }_multiqc_GLbulkRNAseq mqc_in \
            ${ config_arg }

    zip -r '${ params.MQCLabel }_multiqc_GLbulkRNAseq_report.zip' '${ params.MQCLabel }_multiqc_GLbulkRNAseq_report'

    multiqc --version > versions.txt
    """
}

process TRIMGALORE {
  tag "Sample: ${ meta.id }"

  input:
    tuple val(meta), path(reads)

  output:
    tuple val(meta), path("${ meta.id }*trimmed.fastq.gz"), emit: reads
    path("${ meta.id }*.txt"), emit: reports
    path("versions.txt"), emit: version

  script:

    """
    trim_galore --gzip \
    --cores $task.cpus \
    --phred33 \
    ${ meta.paired_end ? '--paired' : '' } \
    $reads \
    --output_dir .

    # rename with _trimmed suffix
    ${ meta.paired_end ? \
      "mv ${ meta.id }_R1_raw_val_1.fq.gz ${ meta.id }_R1_trimmed.fastq.gz; \
      mv ${ meta.id }_R2_raw_val_2.fq.gz ${ meta.id }_R2_trimmed.fastq.gz" : \
      "mv ${ meta.id }_raw_trimmed.fq.gz ${ meta.id }_trimmed.fastq.gz"}

    trim_galore -v > versions.txt
    echo cutadapt version:\$(cutadapt --version) >> versions.txt
    """
}

