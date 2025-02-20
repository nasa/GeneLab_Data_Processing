process FEATURECOUNTS {

  input:
    val(meta)
    tuple path(genomeFasta), path(genomeGtf)
    val(gtf_features)
    val(strandedness)
    path(bam_files)

  output:
    path("FeatureCounts${params.assay_suffix}.tsv"),   emit: counts
    path("FeatureCounts${params.assay_suffix}.tsv.summary"), emit: summary
    path("versions.yml"),                           emit: versions

  script:
    def pairedOption = meta.paired_end ? "-p --countReadPairs -d 10 -D 1000 -P -B" : ""
    def strandOption = (strandedness == "unstranded") ? 0 : (strandedness == "sense") ? 1 : 2
    """
    # Rename each BAM file to remove the '_sorted' substring
    for bam in *.bam; do
      new_bam=\$(echo \$bam | sed 's/_sorted//g')
      mv \$bam \$new_bam
    done

    # Build a list of the renamed BAM files
    bam_list=\$(ls *.bam | tr '\n' ' ')

    # Get all unique feature types in the reference GTF file
    GTF_FEATURES=\$(grep -v '^#' "${genomeGtf}" | cut -f3 | sort | uniq | tr '\n' ',' | sed 's/,\$//')

    # Run featureCounts using the renamed files
    featureCounts ${pairedOption} \\
      -T ${task.cpus} \\
      -G ${genomeFasta} \\
      -a ${genomeGtf} \\
      -t "\${GTF_FEATURES}" \\
      -s ${strandOption} \\
      -o "FeatureCounts${params.assay_suffix}.tsv" \\
      \$bam_list

    echo '"${task.process}":' > versions.yml
    echo "    subread: \$(featureCounts -v 2>&1 | awk '{print \$2}' | sed 's/^v//' | tr -d '\n')" >> versions.yml
    """
}
