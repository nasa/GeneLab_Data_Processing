process SOFTWARE_VERSIONS {
    input:
        path versions_file
    
    output:
        path "software_versions_GLbulkRNAseq.txt", emit: software_versions

    script:
    """
    software_versions.py ${versions_file} software_versions_GLbulkRNAseq.txt --assay rnaseq
    """
}