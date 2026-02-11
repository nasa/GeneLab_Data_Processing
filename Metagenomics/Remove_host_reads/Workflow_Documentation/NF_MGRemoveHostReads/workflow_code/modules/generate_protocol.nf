process GENERATE_PROTOCOL {

    beforeScript "chmod +x ${projectDir}/bin/*"
    tag "Generating your analysis protocol..."
    publishDir "${params.outdir}/processing_info"

    input:
        tuple val(host), val(refSeq_ID), val(genome)
        path(software_versions)

    output:
        path("protocol.txt")
    
    script:
        """
        generate_protocol.sh ${software_versions} ${host} "${refSeq_ID}" ${genome} > protocol.txt
        """
}