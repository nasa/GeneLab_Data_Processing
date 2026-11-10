process UPDATE_ASSAY_TABLE {

    input:
        path(runsheet)
        path("input/*")
        path(summary_file)
    
    output:
        path("a_*.txt"), emit: updated_assay_table

    script:
    def isa_or_assay = params.assay_table ? "--assay_table input/*" : "--isa_zip input/*"
    """
    update_assay_table.py --runsheet ${runsheet} --summary_file ${summary_file} ${isa_or_assay}
    """
}