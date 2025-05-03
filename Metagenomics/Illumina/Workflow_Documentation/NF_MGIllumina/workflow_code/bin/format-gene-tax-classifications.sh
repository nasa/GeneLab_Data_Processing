#!/usr/bin/env bash

 awk -F $'\t' ' BEGIN { OFS = FS } { if ( $3 == "lineage" ) { print $1,$3,$5,$6,$7,$8,$9,$10,$11 } \
     else if ( $2 == "ORF has no hit to database" || $2 ~ /^no taxid found/ ) { print $1,"NA","NA","NA","NA","NA","NA","NA","NA" } \
     else { n=split($3,lineage,";"); print $1,lineage[n],$5,$6,$7,$8,$9,$10,$11 } } ' ${1} \
     | sed 's/no support/NA/g' | sed 's/superkingdom/domain/' | sed 's/# ORF/gene_ID/' | sed 's/lineage/taxid/' > ${2}
