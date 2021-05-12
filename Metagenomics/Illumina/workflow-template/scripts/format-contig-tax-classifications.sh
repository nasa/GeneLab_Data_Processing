#!/usr/bin/env bash

awk -F $'\t' ' BEGIN { OFS = FS } { if ( $2 == "classification" ) { print $1,$4,$6,$7,$8,$9,$10,$11,$12 } \
    else if ( $2 == "no taxid assigned" ) { print $1,"NA","NA","NA","NA","NA","NA","NA","NA" } \
    else { n=split($4,lineage,";"); print $1,lineage[n],$6,$7,$8,$9,$10,$11,$12 } } ' ${1} \
    | sed 's/no support/NA/g' | sed 's/superkingdom/domain/' | sed 's/^# contig/contig_ID/' | sed 's/lineage/taxid/' > ${2}
