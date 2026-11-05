process FETCH_ISA {
    tag "${osd_accession}_${glds_accession}"

    publishDir "${publishdir}/Metadata",
        mode: params.publish_dir_mode

    input:
    val(publishdir)
    val(osd_accession)
    val(glds_accession)
    
    output:
    path "*.zip", emit: isa_archive

    script:
    """
    fetch_isa.py --osd ${osd_accession} --outdir .
    """
}
