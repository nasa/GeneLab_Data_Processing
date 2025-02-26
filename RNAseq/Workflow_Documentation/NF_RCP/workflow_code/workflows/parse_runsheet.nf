def colorCodes = [
    c_line: "┅" * 70,
    c_back_bright_red: "\u001b[41;1m",
    c_bright_green: "\u001b[32;1m",
    c_blue: "\033[0;34m",
    c_yellow: "\u001b[33;1m",
    c_reset: "\033[0m"
]
// Adapted from Function: https://github.com/nf-core/rnaseq/blob/master/modules/local/process/samplesheet_check.nf
// Function to get list of [ meta, [ fastq_1_path, fastq_2_path ] ]
def get_runsheet_paths(LinkedHashMap row) {
    def meta = [:]
    def GENE_ID_TYPES = [
        // Mammals
        "homo_sapiens": "ENSEMBL",
        "mus_musculus": "ENSEMBL",
        "rattus_norvegicus": "ENSEMBL",

        // Other Vertebrates
        "danio_rerio": "ENSEMBL",
        "oryzias_latipes": "ENSEMBL",

        // Invertebrates
        "caenorhabditis_elegans": "ENSEMBL",
        "drosophila_melanogaster": "ENSEMBL",

        // Plants
        "arabidopsis_thaliana": "TAIR",
        "brachypodium_distachyon": "ENSEMBL",
        "oryza_sativa": "ENSEMBL",

        // Microbes
        "bacillus_subtilis": "ENSEMBL",
        "escherichia_coli": "ENSEMBL",
        "lactobacillus_acidophilus": "LOCUS",
        "mycobacterium_marinum": "LOCUS",
        "pseudomonas_aeruginosa": "LOCUS",
        "salmonella_enterica": "ENSEMBL",
        "saccharomyces_cerevisiae": "ENSEMBL",
        "serratia_liquefaciens": "LOCUS",
        "staphylococcus_aureus": "LOCUS",
        "streptococcus_mutans": "LOCUS",
        "vibrio_fischeri": "LOCUS"
    ]

    meta.id = row["Sample Name"]
    meta.organism_sci = row.organism.replaceAll(" ","_").toLowerCase()
    meta.gene_id_type = GENE_ID_TYPES.get(meta.organism_sci, "gene_id")
    meta.paired_end = row.paired_end.toBoolean()
    meta.has_ercc = row.has_ERCC.toBoolean()

    // Extract factors
    meta.factors = row.findAll { key, value -> 
        key.startsWith("Factor Value[") && key.endsWith("]")
    }.collectEntries { key, value ->
        [(key[13..-2]): value] // Remove "Factor Value[" and "]"
    }

    def array = []
    def raw_reads = []
    raw_reads.add(file(row.read1_path))
    if (meta.paired_end) {
        raw_reads.add(file(row.read2_path))
      }
    array = [meta, raw_reads]
    return array
}

def mutate_to_single_end(it) {
    def new_meta = it[0].clone()  // Create a copy of the meta map
    new_meta.paired_end = false   // Set paired_end to false
    return [new_meta, [it[1][0]]] // Return only first read
}

workflow PARSE_RUNSHEET {
    take:
        runsheet_path
    
    main:
        sample_limit = params.limit_samples_to ? params.limit_samples_to : -1 // -1 in take means no limit

        ch_samples = runsheet_path 
            | splitCsv(header: true)
            | map { row -> get_runsheet_paths(row) }
            | map{ it -> params.force_single_end ? mutate_to_single_end(it) : it }
            | take( sample_limit )

        // Remove the redundant paired-end handling since mutate_to_single_end now handles it
        ch_samples | set { ch_samples }

        // ch_samples | view

        // Validate consistency across samples
        ch_samples
            .map { meta, reads -> [meta.has_ercc, meta.paired_end, meta.organism_sci] }
            .unique()
            .count()
            .subscribe { count ->
                if (count > 1) {
                    log.error "${colorCodes.c_back_bright_red}ERROR: Inconsistent metadata across samples. Please check the runsheet.${colorCodes.c_reset}"
                    exit 1
                } else {
                    println "${colorCodes.c_bright_green}Metadata consistency check passed.${colorCodes.c_reset}"
                }
            }

        // Print autodetected processing metadata for the first sample
        ch_samples.take(1) | view { meta, reads -> 
            """${colorCodes.c_blue}Autodetected Processing Metadata:${colorCodes.c_bright_green}
            Has ERCC: ${meta.has_ercc}
            Paired End: ${meta.paired_end}
            Organism: ${meta.organism_sci}
            Gene ID Type: ${meta.gene_id_type}${colorCodes.c_reset}"""
        }
        // Check that all read file paths are unique
        ch_samples
            .flatMap { meta, reads -> reads }
            .collect()
            .map { all_reads ->
                def total_count = all_reads.size()
                def unique_count = all_reads.toSet().size()
                
                if (unique_count != total_count) {
                    throw new RuntimeException("${colorCodes.c_back_bright_red}ERROR: Duplicate read file paths detected. Please check the runsheet.${colorCodes.c_reset}")
                } else {
                    println "${colorCodes.c_bright_green}All ${unique_count} read file paths are unique.${colorCodes.c_reset}"
                }
            }

    emit:
        samples = ch_samples
}