def colorCodes = [
    c_line: "┅" * 70,
    c_back_bright_red: "\u001b[41;1m",
    c_bright_green: "\u001b[32;1m",
    c_blue: "\033[0;34m",
    c_yellow: "\u001b[33;1m",
    c_reset: "\033[0m"
]

process PARSE_ANNOTATIONS_TABLE {
  // Extracts data from this kind of table: 
  // https://github.com/nasa/GeneLab_Data_Processing/blob/master/GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110/GL-DPPD-7110_annotations.csv

  input:
    val(annotations_csv_url_string)
    val(organism_sci)
  
  output:
    val(fasta_url), emit: reference_fasta_url
    val(gtf_url), emit: reference_gtf_url
    val(gene_annotations_url), emit: gene_annotations_url
    val(reference_source), emit: reference_source
    val(reference_version), emit: reference_version
  
  exec:
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
    fasta_url = organisms[organism_key][5]
    gtf_url = organisms[organism_key][6]
    gene_annotations_url = organisms[organism_key][10]
    reference_version = organisms[organism_key][3]
    reference_source = organisms[organism_key][4]
    println "${colorCodes.c_blue}Annotation table values parsed for '${organism_key}':${colorCodes.c_bright_green}"
    println "            Reference Fasta URL: ${fasta_url}"
    println "            Reference GTF URL: ${gtf_url}" 
    println "            Gene Annotations URL: ${gene_annotations_url}"
    println "            Reference Source: ${reference_source}${colorCodes.c_reset}"
    if (reference_source.toLowerCase().contains('ensembl')) {
        println "${colorCodes.c_bright_green}            Reference Version: ${reference_version}${colorCodes.c_reset}"
    }
}