// Nextflow's path staging fetches http(s)/ftp URLs (or copies local paths) automatically;
// this process just normalizes the result to the manufacturer's expected filename,
// decompressing if the source was gzipped.
process COPY_ARRAY_DATA_FILES {
    tag "${ meta.id }"

    input:
        tuple val(meta), path(array_data_file, stageAs: "staged_input")

    output:
        tuple val(meta), path("${meta.file_name}")

    script:
    if (meta.is_gz) {
        """
        gunzip -c staged_input > "${meta.file_name}"
        """
    } else {
        """
        cp -P staged_input "${meta.file_name}"
        """
    }
}
