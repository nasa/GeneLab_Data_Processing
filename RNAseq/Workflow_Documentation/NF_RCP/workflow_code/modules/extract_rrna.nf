process EXTRACT_RRNA {
    /**
     * Extracts unique rRNA ENSEMBL gene IDs from a GTF file.
     */

    input:
        val organism_sci
        path gtf_file

    output:
        path "${organism_sci}_rrna_ids.txt", emit: rrna_ids

    script:
        """
        grep "rRNA" ${gtf_file} \
            | grep -o 'gene_id "[^"]*"' \
            | sed 's/gene_id "\\(.*\\)"/\\1/' \
            | sort -u > ${organism_sci}_rrna_ids.txt
        """
}