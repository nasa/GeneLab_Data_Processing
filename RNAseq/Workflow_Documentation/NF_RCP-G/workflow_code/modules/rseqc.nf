process GENEBODY_COVERAGE {
  tag "Sample: ${ meta.id }"

  input:
    tuple val(meta), path(bam_file), path(bai_file) // bam file sorted by coordinate
    path(genome_bed)

  output:
    path("${ meta.id }.geneBodyCoverage.txt"), emit: log_only
    path("${ meta.id }.geneBodyCoverage.*"), emit: all_output
    tuple val(meta), path("${ meta.id }.geneBodyCoverage.txt"), emit: log
    path("versions.yml"), emit: versions

  script:
    def log_fname = "${ meta.id }.geneBodyCoverage.txt" 
    """    
    geneBody_coverage.py -r ${ genome_bed} -i ${ bam_file } -o ${ meta.id }

    # VERSIONS
    echo '"${task.process}":' > versions.yml
    echo "    geneBody_coverage.py: \$(geneBody_coverage.py --version | sed -e "s/geneBody_coverage.py //g")" >> versions.yml
    """
}

process INFER_EXPERIMENT {
  tag "Sample: ${meta.id}"

  input:
    tuple val(meta), path(bam_file), path(bai_file) // bam file sorted by coordinate
    path(genome_bed)

  output:
    path("${meta.id}_infer_expt.out"), emit: log_only
    tuple val(meta), path("${meta.id}_infer_expt.out"), emit: log
    path("versions.yml"), emit: versions

  script:
    def log_fname = "${meta.id}_infer_expt.out"
    """    
    infer_experiment.py -r ${genome_bed} -i ${bam_file} -s ${params.rseqc_sample_count} > ${log_fname}

    # VERSIONS
    echo '"${task.process}":' > versions.yml
    echo "    infer_experiment.py: \$(infer_experiment.py --version | sed -e 's/infer_experiment.py //g')" >> versions.yml
    """
}

process INNER_DISTANCE {
  tag "Sample: ${ meta.id }"

  input:
    tuple val(meta), path(bam_file), path(bai_file) // bam file sorted by coordinate
    path(genome_bed)
    val(max_read_length)

  output:
    path("${ meta.id }.inner_distance_freq.txt"), emit: log_only
    path("${ meta.id }.inner_distance*"), emit: all_output
    tuple val(meta), path("${ meta.id }.inner_distance_freq.txt"), emit: log
    path("versions.yml"), emit: versions

  when:
    meta.paired_end

  script:
    def log_fname = "${ meta.id }.inner_distance_freq.txt" 
    """    
    inner_distance.py -r ${ genome_bed } -i ${ bam_file } -k ${ params.rseqc_sample_count } -l -${ max_read_length_ch } -u 350 -o ${ meta.id } 

    # VERSIONS
    echo '"${task.process}":' > versions.yml
    echo "    inner_distance.py: \$(inner_distance.py --version | sed -e 's/inner_distance.py //g')" >> versions.yml
    """
}

process READ_DISTRIBUTION {
  tag "Sample: ${ meta.id }"

  input:
    tuple val(meta), path(bam_file), path(bai_file)  // bam file sorted by coordinate
    path(genome_bed)

  output:
    path("${ meta.id }_read_dist.out"), emit: log_only
    tuple val(meta), path("${ meta.id }_read_dist.out"), emit: log
    path("versions.yml"), emit: versions

  script:
    def log_fname = "${ meta.id }_read_dist.out"
    """    
    read_distribution.py -r ${ genome_bed } -i ${ bam_file } > ${ meta.id }_read_dist.out

    # VERSIONS
    echo '"${task.process}":' > versions.yml
    echo "    read_distribution.py: \$(read_distribution.py --version | sed -e 's/read_distribution.py //g')" >> versions.yml
    """
}