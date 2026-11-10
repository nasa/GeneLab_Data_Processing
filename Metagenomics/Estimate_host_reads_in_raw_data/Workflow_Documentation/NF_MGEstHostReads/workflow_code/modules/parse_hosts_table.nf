process PARSE_HOSTS_TABLE {

  input:
    val(hosts_table)
    val(organism_key)
  
  output:
    val(fasta_url), emit: reference_fasta_url
    val(genome), emit: reference_genome
    val(accession), emit: reference_accession
  
  exec:
    def colorCodes = [
        c_line: "┅" * 70,
        c_back_bright_red: "\u001b[41;1m",
        c_bright_green: "\u001b[32;1m",
        c_blue: "\033[0;34m",
        c_yellow: "\u001b[33;1m",
        c_reset: "\033[0m"
    ]

    def organisms = [:]
    println "${colorCodes.c_yellow}Fetching table from ${hosts_table}${colorCodes.c_reset}"
    
    new File(hosts_table).splitEachLine(",") {fields ->
        organisms[fields[0]] = fields 
    }
    
    // Check if the organism exists in the table
    if (organisms.containsKey(organism_key)) {
      fasta_url  = organisms[organism_key].size() > 4 ? organisms[organism_key][4] : null
      genome     = organisms[organism_key].size() > 3 ? organisms[organism_key][3] : null
      accession  = organisms[organism_key].size() > 2 ? organisms[organism_key][2] : null
      
      println "${colorCodes.c_blue}Hosts table values parsed for '${organism_key}':${colorCodes.c_bright_green}"
      println "--------------------------------------------------"
      println "- fasta_url: ${fasta_url}"
      println "- genome assembly: ${genome}"
      println "- assembly accession: ${accession}${colorCodes.c_reset}"
      println "--------------------------------------------------"
    }
    else {
      fasta_url = null
      genome = null
      accession = null
      
      println "${colorCodes.c_back_bright_red}WARNING: Organism '${organism_key}' not found in hosts table.${colorCodes.c_reset}"
      println "${colorCodes.c_yellow}Returning null values for all outputs.${colorCodes.c_reset}"
        
    }
}