process SUMMARY {

    tag "${kraken_output.simpleName.replaceFirst(/-kraken2-output$/, '')}"
    publishDir "${params.outdir}/results/kraken2-output"

    input:
    path kraken_output
    path kraken_report
    
    output:
    path "*-removal-info.tmp"

    script:
    """
    meta_id=\$(basename $kraken_output | sed 's/-kraken2-output.txt//')

    total_fragments=\$(wc -l $kraken_output | sed 's/^ *//' | cut -f 1 -d " ")
    fragments_retained=\$(grep -w -m 1 "unclassified" $kraken_report | cut -f 2)
    perc_removed=\$(printf "%.2f\\n" \$(echo "scale=4; 100 - \$fragments_retained / \$total_fragments * 100" | bc -l))

    echo -e "\$meta_id\\t\$total_fragments\\t\$fragments_retained\\t\$perc_removed\\n" > \$meta_id-removal-info.tmp
    """
}