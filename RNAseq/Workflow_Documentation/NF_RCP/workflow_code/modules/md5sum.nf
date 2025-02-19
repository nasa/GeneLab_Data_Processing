process MD5SUM {
    /**
     * Run recursive md5sum on files/directories and output md5sum files
     *
     * Inputs:
     *   path to file(s) or directory
     *   label for output file name
     *
     * Outputs:
     *   {MD5SumLabel}_md5sums.tsv - Raw md5sums
     */

    input:
        path files
        val(md5sum_label)

    output:
        path("${ md5sum_label }_md5sum_GLbulkRNAseq.tsv"), emit: md5sums

    script:
        """
        # Generate raw md5sums
        if [ -d "${ files }" ]; then
            find "${ files }" -type f -exec md5sum {} \\; > ${ md5sum_label }_md5sum_GLbulkRNAseq.tsv
        else
            md5sum ${ files } > ${ md5sum_label }_md5sum_GLbulkRNAseq.tsv
        fi

        """
}