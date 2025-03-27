process ASSESS_STRANDEDNESS {

  input:
    path("infer_out/*") // a collection of infer_experiment stdout files

  output:
    path("result.txt")

  stub:
    """
    assess_strandedness.py infer_out
    echo "unstranded:0.48595" > result.txt # override original results, this is because heavy truncation and genome subsampling can result in an ambigious strand assignment, which normally is an issue, but should be ignore for stubruns 
    """

  script:
    """
    assess_strandedness.py infer_out
    """
}


process SORT_INDEX_BAM {
  tag "Sample:${ meta.id }"
  label 'big_mem'

  input:
    tuple val(meta), path(bam_file)

  output:
    tuple val(meta), path(sorted_bam_fname), path("${ sorted_bam_fname }.bai"), emit: bam
    tuple path(sorted_bam_fname), path("${ sorted_bam_fname }.bai"), emit: bam_only_files
    path("versions.txt"), emit: version

  script:
    sorted_bam_fname = bam_file.name.replaceAll('.out.bam','_sorted.out.bam')
    mem_MB_per_thread = (task.memory.toMega().intValue() * 0.8 / task.cpus).intValue()
    """    
    samtools sort -m ${ mem_MB_per_thread  }M \
                  --threads ${ task.cpus } \
                  -o ${ sorted_bam_fname } \
                  ${ bam_file }

    samtools index -@ ${ task.cpus  } ${ sorted_bam_fname }

    # VERSIONS
    echo "samtools version:\$(samtools --version-only)" > versions.txt
    """
}


process INFER_EXPERIMENT {
  tag "Sample:${ meta.id }"
  label 'big_mem'

  input:
    tuple val(meta), path(bam_file), path(bai_file), path(bed_file) // bam file sorted by coordinate

  output:
    path("${ meta.id }_infer_expt.out"), emit: log_only
    tuple val(meta), path("${ meta.id }_infer_expt.out"), emit: log
    path("versions.txt"), emit: version

  script:
    def log_fname = "${ meta.id }_infer_expt.out"
    """    
    infer_experiment.py -r ${ bed_file } -i ${ bam_file } -s ${ params.quality.rseqc_sample_count } > ${ log_fname }

    # VERSIONS
    echo "RSeQC infer_experiment version below:\n" > versions.txt 
    infer_experiment.py --version >> versions.txt
    """
}

process GENEBODY_COVERAGE {
  tag "Sample:${ meta.id }"
  label 'big_mem'

  input:
    tuple val(meta), path(bam_file), path(bai_file), path(bed_file) // bam file sorted by coordinate

  output:
    path("${ meta.id }.geneBodyCoverage.txt"), emit: log_only
    path("${ meta.id }.geneBodyCoverage.*"), emit: all_output
    tuple val(meta), path("${ meta.id }.geneBodyCoverage.txt"), emit: log
    path("versions.txt"), emit: version

  script:
    def log_fname = "${ meta.id }.geneBodyCoverage.txt" 
    """    
    geneBody_coverage.py -r ${ bed_file} -i ${ bam_file } -o ${ meta.id }

    # VERSIONS
    echo "RSeQC genebody_coverage version below:\n" > versions.txt 
    geneBody_coverage.py --version >> versions.txt
    """
}

process INNER_DISTANCE {
  tag "Sample:${ meta.id }"
  label 'big_mem'

  input:
    tuple val(meta), path(bam_file), path(bai_file), path(bed_file) // bam file sorted by coordinate

  output:
    path("${ meta.id }.inner_distance_freq.txt"), emit: log_only
    path("${ meta.id }.inner_distance*"), emit: all_output
    tuple val(meta), path("${ meta.id }.inner_distance_freq.txt"), emit: log
    path("versions.txt"), emit: version

  when:
    meta.paired_end

  script:
    def log_fname = "${ meta.id }.inner_distance_freq.txt" 
    """    
    inner_distance.py -r ${ bed_file } -i ${ bam_file } -k ${ params.quality.rseqc_sample_count } -l -150 -u 350 -o ${ meta.id } 

    # VERSIONS
    echo "RSeQC inner_distance version below:\n" > versions.txt 
    inner_distance.py --version >> versions.txt
    """
}

process READ_DISTRIBUTION {
  tag "Sample:${ meta.id }"
  label 'big_mem'

  input:
    tuple val(meta), path(bam_file), path(bai_file), path(bed_file) // bam file sorted by coordinate

  output:
    path("${ meta.id }_read_dist.out"), emit: log_only
    tuple val(meta), path("${ meta.id }_read_dist.out"), emit: log
    path("versions.txt"), emit: version

  script:
    def log_fname = "${ meta.id }_read_dist.out"
    """    
    read_distribution.py -r ${ bed_file } -i ${ bam_file } > ${ meta.id }_read_dist.out

    # VERSIONS
    echo "RSeQC read_distribution version below:\n" > versions.txt 
    read_distribution.py --version >> versions.txt
    """
}
