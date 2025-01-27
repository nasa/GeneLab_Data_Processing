process FEATURECOUNTS {

  input:
    val(meta)
    tuple path(genomeFasta), path(genomeGtf)
    val(strandedness)
    path(bam_files)

  output:
    tuple path("FeatureCounts_GLbulkRNAseq.csv"), path("FeatureCounts_GLbulkRNAseq.csv.summary"), emit: publishables
    path("versions.yml"), emit: versions
  script:
    def pairedOption = meta.paired_end ? "-p" : ""
    def strandOption = (strandedness == "unstranded") ? 0 : (strandedness == "sense") ? 1 : 2
    def bamList = bam_files.join(' ')
    """
    featureCounts ${pairedOption} \
    -T ${ task.cpus } \
    -a ${ genomeGtf } \
    -s ${strandOption} \
    -t exon \
    -g gene_id \
    -o "FeatureCounts_GLbulkRNAseq.csv" \
    ${bamList}


    echo '"${task.process}":' > versions.yml
    echo "    featurecounts: \$(echo \$(featureCounts -v 2>&1) | sed 's/^.*featureCounts v//; s/ .*\$//')" >> versions.yml
    """

}