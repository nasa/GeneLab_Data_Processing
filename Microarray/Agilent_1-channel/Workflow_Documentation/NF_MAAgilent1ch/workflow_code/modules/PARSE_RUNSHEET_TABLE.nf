process PARSE_RUNSHEET_TABLE {
  // Extracts dataset metadata from runsheet

  input:
    val(runsheet_row)
  
  output:
    val(meta) // includes 'primary_keytype' 'organism_sci' 'array_design_ref'
  
  exec:
    // 'Short names' as per: http://www.pantherdb.org/panther/summaryStats.jsp
    def ORGANISMS = ["mus_musculus":"MOUSE",
                     "danio_rerio":"ZEBRAFISH",
                     "rattus_norvegicus":"RAT",
                     "homo_sapiens":"HUMAN",
                     "drosophila_melanogaster":"FLY",
                     "caenorhabditis_elegans":"WORM",
                     "brachypodium_distachyon":"BRADI",
                     "arabidopsis_thaliana":"ARABIDOPSIS"]

    def PRIMARY_KEYS = ["mus_musculus":"ENSEMBL",
                        "danio_rerio":"ENSEMBL",
                        "rattus_norvegicus":"ENSEMBL",
                        "homo_sapiens":"ENSEMBL",
                        "drosophila_melanogaster":"ENSEMBL",
                        "caenorhabditis_elegans":"ENSEMBL",
                        "brachypodium_distachyon":"ENSEMBL",
                        "arabidopsis_thaliana":"TAIR"]

    def meta = [:]
    meta.organism_sci               = runsheet_row.organism.replaceAll(" ","_").toLowerCase()
    meta.organism_non_sci           = ORGANISMS[meta.organism_sci]
    meta.primary_keytype            = PRIMARY_KEYS[meta.organism_sci]
    meta.biomart_attribute          = runsheet_row.array_design_ref
}
