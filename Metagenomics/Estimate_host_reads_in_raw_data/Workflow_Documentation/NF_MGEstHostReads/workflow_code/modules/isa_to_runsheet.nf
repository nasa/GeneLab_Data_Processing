process ISA_TO_RUNSHEET {
    tag "${params.osd}_${params.glds}"

    input: 
        path isa_archive
        path dp_tools_plugin

    output:
        path "*.csv", emit: runsheet
        path("versions.txt"), emit: version

    script:
    """
    dpt-isa-to-runsheet --accession ${params.osd} --isa-archive ${isa_archive} --plugin-dir ${dp_tools_plugin}
    echo "dp_tools \$(pip show dp_tools | grep Version | sed 's/Version: //')" >> versions.txt
    """
}
