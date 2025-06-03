process ISA_TO_RUNSHEET {
    tag "${osd_accession}_${glds_accession}"

    publishDir "${ch_outdir}/Metadata",
        mode: params.publish_dir_mode,
        pattern: "*.csv"

    publishDir "${ch_outdir}/Metadata",
        mode: params.publish_dir_mode,
        pattern: "isa_archive/*",
        saveAs: { filename ->
            if (filename.startsWith("isa_archive/")) return filename.replace("isa_archive/", "")
            else return filename
        }

    input: 
    val(ch_outdir)
    val(osd_accession)
    val(glds_accession)
    path(isa_archive)
    path(dp_tools_plugin)

    output:
    path("*.csv"), emit: runsheet
    path("isa_archive/${isa_archive}")
    //path("versions.yml"), emit: versions

    script:
    """
    dpt-isa-to-runsheet --accession ${osd_accession} --isa-archive ${isa_archive} --plugin-dir ${dp_tools_plugin}

    # Copy the ISA archive to the output directory
    mkdir -p isa_archive
    cp ${isa_archive} isa_archive/

    #echo '"${task.process}":' > versions.yml
    #echo "    dp_tools: \$(pip show dp_tools | grep Version | sed 's/Version: //')" >> versions.yml
    """
}