#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/**************************************************************************************** 
*********************  Sequence quality assessment and control processes ****************
****************************************************************************************/

// A 3-column (single-end) or 4-column (paired-end) file
//params.csv_file = "${projectDir}/PE_file.csv" 
//params.prefix = "raw"
//params.multiqc_config = "${projectDir}/config/multiqc.config"

// FastQC performed on reads
process FASTQC {

    tag "Running fastqc on ${sample_id}..."
    beforeScript "chmod +x ${projectDir}/bin/*"
    label "fastqc"

    input:
        tuple val(sample_id), path(reads), val(isPaired)
    output:
        tuple path("*.html"), path("*.zip"), emit: html
        path("versions.txt"), emit: version
    script:
        """
        fastqc -o . \\
        -t ${task.cpus} \\
        ${reads}

        fastqc --version > versions.txt
        """
}

process MULTIQC {

    tag "Running multiqc on the ${prefix} files..."
    beforeScript "chmod +x ${projectDir}/bin/*"
    
    input:
        val(prefix)   
        path(multiqc_config)
        path(files)
    output:
        path("${prefix}_multiqc_report"), emit: report_dir
        path("versions.txt"), emit: version
    script:
        """
        multiqc -q --filename ${prefix}_multiqc \\
                --force --cl-config 'max_table_rows: 99999999' \\
                --interactive --config ${multiqc_config} \\
                --outdir ${prefix}_multiqc_report  ${files} > /dev/null 2>&1

        multiqc --version > versions.txt
        """
  }


process ZIP_MULTIQC {

    tag "Zipping ${prefix} multiqc.."
    label "zip"
 
    input:
        val(prefix)
        path(multiqc_dir)

    output:
        path("${prefix}_multiqc${params.assay_suffix}_report.zip"), emit: report
        path("versions.txt"), emit: version

    script:
        """
        # Zipping
        zip -q -r ${prefix}_multiqc${params.assay_suffix}_report.zip ${multiqc_dir}
        zip -h | grep "Zip" | sed -E 's/(Zip.+\\)).+/\\1/' > versions.txt
        """
}
