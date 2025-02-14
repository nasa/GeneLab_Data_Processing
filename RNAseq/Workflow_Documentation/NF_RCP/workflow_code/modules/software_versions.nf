process SOFTWARE_VERSIONS {
    input:
        path versions_file
    
    output:
        path "software_versions_GLbulkRNAseq.txt", emit: software_versions
        path "software_versions_GLbulkRNAseq.yaml", emit: software_versions_yaml

    script:
    """
    software_versions.py ${versions_file} software_versions_GLbulkRNAseq.txt --assay rnaseq
    """
}