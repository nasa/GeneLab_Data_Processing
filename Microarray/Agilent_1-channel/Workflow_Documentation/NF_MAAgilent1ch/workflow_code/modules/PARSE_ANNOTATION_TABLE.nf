process PARSE_ANNOTATION_TABLE {
  // Extracts data from this kind of table: 
  // https://github.com/nasa/GeneLab_Data_Processing/blob/master/GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110/GL-DPPD-7110_annotations.csv

  input:
    val(annotations_csv_url_string)
    val(organism_sci)
  
  output:
    val(annotations_db_url), emit: annotations_db_url
    tuple val(ensemblVersion), val(ensemblSource), emit: reference_version_and_source
  
  exec:
    def organisms = [:]
    println "Fetching table from ${annotations_csv_url_string}"
    
    // download data to memory
    annotations_csv_url_string.toURL().splitEachLine(",") {fields ->
          organisms[fields[1]] = fields
    }
    // extract required fields
    organism_key = organism_sci.capitalize().replace("_"," ")
    // fasta_url = organisms[organism_key][5]
    // gtf_url = organisms[organism_key][6]
    annotations_db_url = organisms[organism_key][9]
    ensemblVersion = organisms[organism_key][3]
    ensemblSource = organisms[organism_key][4]

    println "PARSE_ANNOTATION_TABLE:"
    println "Values parsed for '${organism_key}' using process:"
    println "--------------------------------------------------"
    // println "- fasta_url: ${fasta_url}"
    // println "- gtf_url: ${gtf_url}"
    println "- annotations_db_url: ${annotations_db_url}"
    println "- ensemblVersion: ${ensemblVersion}"
    println "- ensemblSource: ${ensemblSource}"
    println "--------------------------------------------------"
}
