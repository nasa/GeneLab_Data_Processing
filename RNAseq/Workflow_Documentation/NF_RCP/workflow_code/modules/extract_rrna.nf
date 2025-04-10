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
        # Use grep with OR patterns to find all rRNA genes
        grep -E 'gene_biotype "rRNA"|gene_type "rRNA"|gbkey "rRNA"' ${gtf_file} \
            | grep -o 'gene_id "[^"]*"' \
            | sed 's/gene_id "\\(.*\\)"/\\1/' \
            | sort -u > ${organism_sci}_rrna_ids.txt || touch ${organism_sci}_rrna_ids.txt
        
        # Report count
        echo "Found \$(wc -l < ${organism_sci}_rrna_ids.txt) rRNA genes for ${organism_sci}" >&2
        """
}