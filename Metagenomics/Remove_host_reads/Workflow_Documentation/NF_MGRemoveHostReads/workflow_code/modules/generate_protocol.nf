/*process GENERATE_PROTOCOL {

    beforeScript "chmod +x ${projectDir}/bin/*"
    tag "Generating analysis protocol text..."

    input:
        tuple val(host), val(genome), val(assembly_acc)
        path software_versions
        path protocol

    output:
        path("protocol.txt")
    
    script:
        if (protocol.exists())
        """
        cp ${protocol} protocol.txt
        """
        else
        """
        generate_protocol.sh ${software_versions} ${host} "${assembly_acc}" ${genome} > protocol.txt
        """
}*/

process GENERATE_PROTOCOL {
    tag "Generating analysis protocol text for ${host}"
    label 'python'

    input:
        val(host)
        val(genome)
        val(assembly_acc)
        path(software_versions)
        path(db_dir)

    output:
    path("protocol.txt"), emit: protocol

    script:
    """
    generate_protocol.py \\
        --host ${host} \\
        --db-dir ${db_dir} \\
        --assembly-name "${genome}" \\
        --assembly-acc "${assembly_acc}" \\
        --software-versions ${software_versions}
    """
}