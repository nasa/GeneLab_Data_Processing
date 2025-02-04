// http://qualimap.conesalab.org/doc_html/command_line.html
process QUALIMAP_BAM_QC {
    tag "Sample: ${ meta.id }"

    input:
    tuple val(meta), path(bam_file), path(bai_file)
    path(feature_file) // gtf, gff, or bed
    val(strandedness)

    output:
    path("${meta.id}_qualimap_bam_qc/"), emit: results
    path  "versions.yml",                                  emit: versions

    script:
    strandedness_opt_map = ["sense":"strand-specific-forward","antisense":"strand-specific-reverse","unstranded":"non-strand-specific"]
    def collect_pairs = meta.paired_end ? '--collect-overlap-pairs' : ''
    def memory = (task.memory.mega*0.8).intValue() + 'M'

    """
    mkdir -p tmp
    qualimap \
        --java-mem-size=${memory} \
        bamqc \
        -bam ${bam_file} \\
        --feature-file ${feature_file} \
        -p ${strandedness_opt_map.get(strandedness)} \
        $collect_pairs \
        -outdir ${meta.id}_qualimap_bam_qc \
        -nt ${task.cpus}

    echo '"${task.process}":' > versions.yml
    echo "    qualimap: \$(qualimap 2>&1 | sed 's/^.*QualiMap v.//; s/Built.*\$//')" >> versions.yml
    """
}

process QUALIMAP_RNASEQ_QC {
    tag "Sample: ${ meta.id }"

    input:
    tuple val(meta), path(bam_file), path(bai_file)
    path(gtf) 
    val(strandedness)

    output:
    path("${meta.id}_qualimap_rnaseq_qc/"), emit: results
    path  "versions.yml",                                        emit: versions

    script:
    strandedness_opt_map = ["sense":"strand-specific-forward","antisense":"strand-specific-reverse","unstranded":"non-strand-specific"]
    def paired_end = meta.paired_end ? '--pe' : ''
    def memory = (task.memory.mega*0.8).intValue() + 'M'

    """
    mkdir -p tmp
    qualimap \\
        --java-mem-size=${memory} \
        rnaseq \
        -bam ${bam_file} \
        --gtf ${gtf} \
        -p ${strandedness_opt_map.get(strandedness)} \
        $paired_end \
        -outdir ${meta.id}_qualimap_rnaseq_qc \
        -s 

    echo '"${task.process}":' > versions.yml
    echo "    qualimap: \$(qualimap 2>&1 | sed 's/^.*QualiMap v.//; s/Built.*\$//')" >> versions.yml
    """
}