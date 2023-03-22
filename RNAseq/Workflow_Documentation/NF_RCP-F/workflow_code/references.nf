// This ensures DSL2 syntax and process imports
nextflow.enable.dsl=2

include { DOWNLOAD_ERCC } from './modules/download.nf'
include { CONCAT_ERCC;
          SUBSAMPLE_GENOME;
          TO_PRED;
          TO_BED } from './modules/genome.nf'


include { DOWNLOAD_GUNZIP_REFERENCES } from './modules/download.nf'

process PARSE_ANNOTATIONS_TABLE {
  // Extracts data from this kind of table: 
  // https://github.com/nasa/GeneLab_Data_Processing/blob/master/GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110/GL-DPPD-7110_annotations.csv

  input:
    val(annotations_csv_url_string)
    val(organism_sci)
  
  output:
    tuple val(organism_sci), val(fasta_url), val(gtf_url), emit: reference_genome_urls
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
    fasta_url = organisms[organism_key][5]
    gtf_url = organisms[organism_key][6]
    annotations_db_url = organisms[organism_key][9]
    ensemblVersion = organisms[organism_key][3]
    ensemblSource = organisms[organism_key][4]

    println "PARSE_ANNOTATIONS_TABLE:"
    println "Values parsed for '${organism_key}' using process:"
    println "--------------------------------------------------"
    println "- fasta_url: ${fasta_url}"
    println "- gtf_url: ${gtf_url}"
    println "- annotations_db_url: ${annotations_db_url}"
    println "- ensemblVersion: ${ensemblVersion}"
    println "- ensemblSource: ${ensemblSource}"
    println "--------------------------------------------------"
}

/**************************************************
* ACTUAL WORKFLOW  ********************************
**************************************************/
workflow references{
  take:
    organism_sci
    has_ercc
  main:
      // Must run in any approach to find the appropriate annotations database table
      PARSE_ANNOTATIONS_TABLE( params.reference_table, organism_sci)

      if (params.ref_fasta && params.ref_gtf) {
        genome_annotations_pre_subsample = Channel.fromPath([params.ref_fasta, params.ref_gtf], checkIfExists: true).toList()
        genome_annotations_pre_subsample | view
        Channel.value( [params.ensemblVersion, params.ref_source] ) | set { ch_ref_source_version }
      } else {
        // use assets table to find current fasta and gtf urls and associated metadata about those reference files
        
        DOWNLOAD_GUNZIP_REFERENCES( 
          PARSE_ANNOTATIONS_TABLE.out.reference_genome_urls,
          PARSE_ANNOTATIONS_TABLE.out.reference_version_and_source,
          )
        DOWNLOAD_GUNZIP_REFERENCES.out | set{ genome_annotations_pre_subsample }
        PARSE_ANNOTATIONS_TABLE.out.reference_version_and_source | set { ch_ref_source_version }
      }
      // use assets table to find current annotations file
      PARSE_ANNOTATIONS_TABLE.out.annotations_db_url | set{ ch_gene_annotations_url }

      // SUBSAMPLING STEP : USED FOR DEBUG/TEST RUNS
      if ( params.genomeSubsample ) {
        SUBSAMPLE_GENOME( genome_annotations_pre_subsample, organism_sci, ch_ref_source_version )
        SUBSAMPLE_GENOME.out.build | flatten | toList | set { genome_annotations_pre_ercc }
      } else {
        genome_annotations_pre_subsample | flatten | toList | set { genome_annotations_pre_ercc }
      }

      // ERCC STEP : ADD ERCC Fasta and GTF to genome files
      DOWNLOAD_ERCC(has_ercc).ifEmpty([file("ERCC92.fa"), file("ERCC92.gtf")]) | set { ch_maybe_ercc_refs }
      CONCAT_ERCC( genome_annotations_pre_ercc, ch_maybe_ercc_refs, organism_sci, has_ercc, PARSE_ANNOTATIONS_TABLE.out.reference_version_and_source )
      .ifEmpty { genome_annotations_pre_ercc.value }  | set { genome_annotations }


      TO_PRED( 
        genome_annotations | map { it[1] }, 
        organism_sci,
        PARSE_ANNOTATIONS_TABLE.out.reference_version_and_source
        )
      TO_BED( 
        TO_PRED.out, 
        organism_sci,
        PARSE_ANNOTATIONS_TABLE.out.reference_version_and_source
        )

  emit:
      genome_annotations = genome_annotations
      genome_bed = TO_BED.out
      gene_annotations = ch_gene_annotations_url
      reference_version_and_source = ch_ref_source_version
}
