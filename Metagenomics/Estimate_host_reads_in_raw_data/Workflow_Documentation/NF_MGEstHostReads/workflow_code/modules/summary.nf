process SUMMARY {

    tag "${kraken_output.simpleName.replaceFirst(/-kraken2-output$/, '')}"
    publishDir "${params.outdir}/results/kraken2-output"

    input:
    path kraken_output
    
    output:
    path "*-removal-info.tmp"

    script:
    """
    meta_id=\$(basename $kraken_output | sed 's/-kraken2-output.txt//')

    total_fragments=\$(wc -l $kraken_output | sed 's/^ *//' | cut -f 1 -d " ")
    fragments_classified=\$(grep -w -c "^C" $kraken_output)
    perc_host=\$(printf "%.2f\\n" \$(echo "scale=4; \$fragments_classified / \$total_fragments * 100" | bc -l))

    echo -e "\$meta_id\\t\$total_fragments\\t\$fragments_classified\\t\$perc_host\\n" > \$meta_id-removal-info.tmp
    """
}