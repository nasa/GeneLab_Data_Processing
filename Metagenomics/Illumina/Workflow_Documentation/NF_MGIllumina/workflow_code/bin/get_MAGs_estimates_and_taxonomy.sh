#!/usr/bin/env bash
MAGs_dir=$1
MAG_assembly_summaries=$2
MAGs_checkm_out=$3
gtdbtk_out=$4

# Making sure none of the intermediate files exist already
rm -rf checkm-estimates.tmp \
            gtdb-taxonomies.tmp \
            checkm-estimates-with-headers.tmp \
            gtdb-taxonomies-with-headers.tmp \
            MAGs-overview.tmp \
            MAGs-overview-header.tmp \
            MAGs-overview-sorted.tmp
for MAG in $(cut -f 1 ${MAG_assembly_summaries} | tail -n +2); do

	grep -w -m 1 "^${MAG}" ${MAGs_checkm_out} | \
        cut -f 12,13,14 >> checkm-estimates.tmp

        grep -w "^${MAG}" ${gtdbtk_out}/gtdbtk.*.summary.tsv | \
        cut -f 2 | sed 's/^.__//' | \
        sed 's/;.__/\t/g' | \
	awk 'BEGIN{ OFS=FS="\t" } { for (i=1; i<=NF; i++) if ( $i ~ /^ *$/ ) $i = "NA" }; 1' \
                                >> gtdb-taxonomies.tmp
done

