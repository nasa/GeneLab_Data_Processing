process GET_GTF_FEATURES {

    input:
        tuple path(genome_fasta), path(genome_gtf)

    output:
        path("gtf_features.txt"), emit: gtf_features

    script:
    """
    grep -v '^#' ${genome_gtf} | cut -f3 | sort | uniq | tr '\\n' ',' | sed 's/,\$//' > gtf_features.txt
    """
}
