workflow PARSE_ANNOTATIONS_TABLE {
  // Extracts data from GeneLab Reference Annotations Table or similarly formatted table
  // https://github.com/nasa/GeneLab_Data_Processing/blob/master/GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv
  take:
    annotations_csv_url_string
    organism_sci
  
  main:
    def colorCodes = [
        c_line: "┅" * 70,
        c_back_bright_red: "\u001b[41;1m",
        c_bright_green: "\u001b[32;1m",
        c_blue: "\033[0;34m",
        c_yellow: "\u001b[33;1m",
        c_reset: "\033[0m"
    ]
    
    def organisms = [:]
    println "${colorCodes.c_yellow}Fetching table from ${annotations_csv_url_string}${colorCodes.c_reset}"
    
    // Check if input is a URL or a local file path
    if (annotations_csv_url_string.startsWith('http://') || annotations_csv_url_string.startsWith('https://')) {
      // For URLs: use toURL() method
      annotations_csv_url_string.toURL().splitEachLine(",") {fields ->
            organisms[fields[1]] = fields
      }
    } else {
      // For local files: use File class
      new File(annotations_csv_url_string).splitEachLine(",") {fields ->
            organisms[fields[1]] = fields
      }
    }
    
    // extract required fields
    organism_key = organism_sci.capitalize().replace("_"," ")
    def fasta_url_val = organisms[organism_key][5]
    def gtf_url_val = organisms[organism_key][6]
    def gene_annotations_url_val = organisms[organism_key][10]
    def reference_version_val = organisms[organism_key][3]
    def reference_source_val = organisms[organism_key][4]

    println "${colorCodes.c_blue}Annotation table values parsed for '${organism_key}':${colorCodes.c_bright_green}"
    println "            Reference Fasta URL: ${fasta_url_val}"
    println "            Reference GTF URL: ${gtf_url_val}" 
    println "            Gene Annotations URL: ${gene_annotations_url_val}"
    println "            Reference Source: ${reference_source_val}${colorCodes.c_reset}"
    if (reference_source_val.toLowerCase().contains('ensembl')) {
        println "${colorCodes.c_bright_green}            Reference Version: ${reference_version_val}${colorCodes.c_reset}"
    }

  emit:
    reference_fasta_url = fasta_url_val
    reference_gtf_url = gtf_url_val
    gene_annotations_url = gene_annotations_url_val
    reference_source = reference_source_val
    reference_version = reference_version_val
}