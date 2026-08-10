process PARSE_ANNOTATION_TABLE {
  // Extracts data from this kind of table: 
  // https://github.com/nasa/GeneLab_Data_Processing/blob/master/GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv

  input:
    val(annotations_csv_url_string)
    val(organism_sci)
  
  output:
    val(annotations_db_url), emit: annotations_db_url
    tuple val(ensemblVersion), val(ensemblSource), emit: reference_version_and_source
    val(bioconductor_annotations), emit: bioconductor_annotations
    val(annotations_db_info_url), emit: annotations_db_info_url
  
  exec:
    println "Fetching table from ${annotations_csv_url_string}"
    
    // download data to memory and parse table
    def lines = new URL(annotations_csv_url_string).text.split("\n")
    def headers = lines[0].strip().split(",")
    def records = lines[1..-1].collect { line ->
      def values = line.strip().split(",")
      [headers, values].transpose().collectEntries()
    }

    def organism_key = organism_sci.capitalize().replace("_"," ")

    def organism_record = records.find { rec -> rec['species'] == organism_key }
    if (organism_record == null) {
      throw new Exception("Organism '${organism_key}' not found in annotation table at ${annotations_csv_url_string}")
    } else {
      annotations_db_url = organism_record['genelab_annots_link']
      ensemblVersion = organism_record['ensemblVersion']
      ensemblSource = organism_record['ref_source']
      bioconductor_annotations = organism_record['bioconductor_annotations']
      annotations_db_info_url = organism_record['genelab_annots_info_link']
      
      // Convert figshare ndownloader URL to API endpoint
      if (annotations_db_url != null && annotations_db_url.contains('figshare.com/ndownloader/files/')) {
        file_id = (annotations_db_url =~ /.*\/files\/([a-zA-Z0-9]+).*/)[0][1]
        annotations_db_url = "https://api.figshare.com/v2/file/download/${file_id}"
      }

      // Convert figshare ndownloader URL to API endpoint
      if (annotations_db_info_url != null && annotations_db_info_url.contains('figshare.com/ndownloader/files/')) {
        file_id = (annotations_db_info_url =~ /.*\/files\/([a-zA-Z0-9]+).*/)[0][1]
        annotations_db_info_url = "https://api.figshare.com/v2/file/download/${file_id}"
      }

      println "PARSE_ANNOTATION_TABLE:"
      println "Values parsed for '${organism_key}' using process:"
      println "--------------------------------------------------"
      println "- annotations_db_url: ${annotations_db_url}"
      println "- annotations_db_info_url: ${annotations_db_info_url}"
      println "- ensemblVersion: ${ensemblVersion}"
      println "- ensemblSource: ${ensemblSource}"
      println "--------------------------------------------------"
    }
}
