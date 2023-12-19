/* Processes dealing with retrieving data from GeneLab
*/
// TODO Migrate these CLI args to updated dp tools package API
process RNASEQ_RUNSHEET_FROM_GLDS {
  // Downloads isazip and creates run sheets using GeneLab API
  tag "${ glds_accession }"
  publishDir "${ params.outputDir }/${ params.gldsAccession }/Metadata",
    mode: params.publish_dir_mode

  input:
    val(glds_accession)

  output:
    path("${ glds_accession }_bulkRNASeq_v*.csv"), emit: runsheet
    path("*.zip"), emit: isazip

  script:
    """
    dpt-get-isa-archive --accession ${ glds_accession }

    dpt-isa-to-runsheet --accession ${ glds_accession } \
      --config-type bulkRNASeq --isa-archive *.zip
    """
}


process STAGE_RAW_READS {
  // Stages the raw reads into appropriate publish directory
  tag "${ meta.id }"

  input:
    tuple val(meta), path("?.gz")

  output:
    tuple val(meta), path("${meta.id}*.gz"), emit: raw_reads

  script:
    if ( meta.paired_end ) {
      """
      cp -P 1.gz ${meta.id}_R1_raw.fastq.gz
      cp -P 2.gz ${meta.id}_R2_raw.fastq.gz
      """
    } else {
      """
      cp -P 1.gz  ${meta.id}_raw.fastq.gz
      """
    }
}


process GENERATE_METASHEET {
  // Generates a metadata table, not used in further processing
  tag "${ params.gldsAccession }"
  publishDir "${ params.outputDir }/${ params.gldsAccession }/Metadata",
    mode: params.publish_dir_mode

  input:
    path("isa.zip")
    path(runsheet)

  output:
    path("${ params.gldsAccession }_metadata_table.txt"), emit: metasheet

  script:
    """
    create_table_v2.py --accession ${ params.gldsAccession }  \
                       --isa-zip isa.zip \
                       --output-dir . \
                       --runsheet ${ runsheet }
    """
}

process GENERATE_MD5SUMS {
  // Generates tabular data indicating genelab standard publishing files, md5sum generation, and tool version table formatting
  tag "${ params.gldsAccession }"
  publishDir "${ params.outputDir }/${ params.gldsAccession }/GeneLab",
    mode: params.publish_dir_mode

  input:
    path(data_dir)
    path(runsheet)
    path("dp_tools__NF_RCP")

  output:
    path("*md5sum*")
    path("Missing_md5sum_files_GLbulkRNAseq.txt"), optional: true

  script:
    """
    generate_md5sum_files.py  --root-path ${ data_dir } \\
                              --runsheet-path ${ runsheet } \\
                              --plug-in-dir "dp_tools__NF_RCP"
    """
}

process UPDATE_ISA_TABLES {
  // Generates tabular data indicating genelab standard publishing files, md5sum generation, and tool version table formatting
  tag "${ params.gldsAccession }"
  publishDir "${ params.outputDir }/${ params.gldsAccession }/GeneLab",
    mode: params.publish_dir_mode

  input:
    path(data_dir)
    path(runsheet)

  output:
    path("updated_curation_tables") // directory containing extended ISA tables

  script:
    """
    update_curation_table.py  --root-path ${ data_dir } \\
                              --runsheet-path ${ runsheet } \\
                              --plug-in-dir "dp_tools__NF_RCP"
    """
}

process SOFTWARE_VERSIONS {
  // Generates tabular data indicating genelab standard publishing files, md5sum generation, and tool version table formatting
  tag "${ params.gldsAccession }"
  publishDir "${ params.outputDir }/${ params.gldsAccession }/GeneLab",
    mode: params.publish_dir_mode

  input:
    path("software_versions.txt")

  output:
    path("software_versions_GLbulkRNAseq.md")

  script:
    """
    format_software_versions.py software_versions.txt
    mv software_versions.md software_versions_GLbulkRNAseq.md
    """
}

// Adapted from Function: https://github.com/nf-core/rnaseq/blob/master/modules/local/process/samplesheet_check.nf
// Original Function Credit: Dr. Harshil Patel
// Function to get list of [ meta, [ fastq_1_path, fastq_2_path ] ]
def get_runsheet_paths(LinkedHashMap row) {
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
    meta.id                         = row["Sample Name"]
    meta.organism_sci               = row.organism.replaceAll(" ","_").toLowerCase()
    meta.organism_non_sci           = ORGANISMS[meta.organism_sci]
    meta.primary_keytype            = PRIMARY_KEYS[meta.organism_sci]
    meta.paired_end                 = row.paired_end.toBoolean()
    meta.has_ercc                   = row.has_ERCC.toBoolean()

    def array = []
    def raw_reads = []
    raw_reads.add(file(row.read1_path))
    if (meta.paired_end) {
        raw_reads.add(file(row.read2_path))
      }
    array = [meta, raw_reads]
    return array
}
