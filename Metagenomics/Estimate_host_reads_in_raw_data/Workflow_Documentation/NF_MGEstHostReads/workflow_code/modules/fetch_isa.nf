process FETCH_ISA {

    tag "${params.osd}_${params.glds}"

    publishDir "${params.outdir}/Metadata", mode: params.publishDir_mode

    output:
        path "*.zip", emit: isa_archive
        
    script:
    """
    dpt-get-isa-archive --accession ${ params.osd }
    """
}
