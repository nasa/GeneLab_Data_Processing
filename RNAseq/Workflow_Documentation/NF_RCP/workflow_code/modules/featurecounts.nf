process FEATURECOUNTS {

  input:
    val(meta)
    tuple path(genomeFasta), path(genomeGtf)
    val(strandedness)
    path(bam_files)

  output:
    path("FeatureCounts_GLbulkRNAseq.csv"),   emit: counts
    path("FeatureCounts_GLbulkRNAseq.csv.summary"), emit: summary
    path("versions.yml"),                           emit: versions
  script:
    def pairedOption = meta.paired_end ? "-p --countReadPairs -d 10 -D 1000" : ""
    def strandOption = (strandedness == "unstranded") ? 0 : (strandedness == "sense") ? 1 : 2
    def bamList = bam_files.join(' ')
    """
    featureCounts ${pairedOption} \
    -T ${ task.cpus } \
    -G ${ genomeFasta } \
    -a ${ genomeGtf } \
    -s ${strandOption} \
    -o "FeatureCounts_GLbulkRNAseq.csv" \
    ${bamList}


    echo '"${task.process}":' > versions.yml
    echo "    featurecounts: \$(echo \$(featureCounts -v 2>&1) | sed 's/^.*featureCounts v//; s/ .*\$//')" >> versions.yml
    """

}