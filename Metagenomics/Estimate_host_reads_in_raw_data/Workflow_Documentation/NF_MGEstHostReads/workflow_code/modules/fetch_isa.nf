process FETCH_ISA {

    tag "${params.osd}_${params.glds}"

    output:
        path "*.zip", emit: isa_archive
        
    script:
    """
    dpt-get-isa-archive --accession ${ params.osd }
    """
}
