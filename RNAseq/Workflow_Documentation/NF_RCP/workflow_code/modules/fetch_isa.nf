process FETCH_ISA {

    tag "${osd_accession}"

    publishDir "${glds_accession}/Metadata",
        mode: params.publish_dir_mode

    input:
    val(osd_accession)
    val(glds_accession)
    output:
    path "*.zip", emit: isa_archive

    script:
    """
    python $projectDir/bin/fetch_isa.py --osd ${osd_accession} --outdir .
    """
}