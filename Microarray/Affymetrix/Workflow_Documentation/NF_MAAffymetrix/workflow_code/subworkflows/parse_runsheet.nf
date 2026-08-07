def colorCodes = [
    c_line: "┅" * 70,
    c_back_bright_red: "\u001b[41;1m",
    c_bright_green: "\u001b[32;1m",
    c_blue: "\033[0;34m",
    c_yellow: "\u001b[33;1m",
    c_reset: "\033[0m"
]
// Adapted from Function: https://github.com/nf-core/rnaseq/blob/master/modules/local/process/samplesheet_check.nf
// Function to get list of [ meta, data_file ]
def get_runsheet_paths(LinkedHashMap row) {
    def meta = [:]
    meta.id = row["Sample Name"]
    meta.biomart_id = row["biomart_attribute"]
    meta.organism_sci = row.organism.replaceAll(" ","_").toLowerCase()
    
    // Extract factors
    meta.factors = row.findAll { key, value -> 
        key.startsWith("Factor Value[") && key.endsWith("]")
    }.collectEntries { key, value ->
        [(key[13..-2]): value] // Remove "Factor Value[" and "]"
    }

    meta.is_gz = row['Array Data File Path'].endsWith('.gz')
    // Local staged filename: always the manufacturer's original name, minus .gz
    // (decompress during staging so the QMD never has to deal with compression)
    meta.file_name = meta.is_gz
        ? row['Array Data File Name'].replaceAll(/\.gz$/, '')
        : row['Array Data File Name']

    return [meta, file(row['Array Data File Path'])]
}

workflow PARSE_RUNSHEET {
    take:
        runsheet_path
    
    main:
        // Process samples from the runsheet
        ch_samples = runsheet_path
            | splitCsv(header: true)
            | map { row -> get_runsheet_paths(row) }

        ch_samples | set { ch_samples }

        // Validate consistency across samples
        ch_samples
            .map { meta, data_file -> [meta.biomart_id, meta.organism_sci] }
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
        ch_samples.take(1) | view { meta, data_file -> 
            """${colorCodes.c_bright_green}Autodetected Processing Metadata:
            Biomart Attribute: ${meta.biomart_id}
            Organism: ${meta.organism_sci}${colorCodes.c_reset}"""
        }

        // Check that all read files are unique
        ch_samples
            .map { meta, data_file -> meta.file_name }
            .collect()
            .map { all_names ->
                if (all_names.toSet().size() != all_names.size()) {
                    throw new RuntimeException("${colorCodes.c_back_bright_red}ERROR: Duplicate staged filenames detected — two samples would collide in the analysis work directory.${colorCodes.c_reset}")
                } else {
                    println "${colorCodes.c_bright_green}All ${all_names.size()} staged filenames are unique.${colorCodes.c_reset}"
        }
            }
        ch_samples
            .map { meta, data_file -> data_file }
            .collect()
            .map { all_data_files ->                
                if (all_data_files.toSet().size() != all_data_files.size()) {
                    throw new RuntimeException("${colorCodes.c_back_bright_red}ERROR: Duplicate assay data files detected. Please check the runsheet.${colorCodes.c_reset}")
                } else {
                    println "${colorCodes.c_bright_green}All ${all_data_files.size()} assay data files are unique.${colorCodes.c_reset}"
                }
            }

    emit:
        samples = ch_samples
        runsheet = runsheet_path
}
