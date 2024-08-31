#!/usr/bin/env bash

# A script to the the input csv file rtequired by this pipeline when a GLDS accession is provided
# Rather than the required input csv file

ASSAY_TABLE=$1 # target_region_assay_table.txt

cat ${ASSAY_TABLE} | \
sed 's/"//g' | \
awk -v PWD=$PWD 'BEGIN{FS="\t"; OFS=","; print "sample_id,read,f_primer,r_primer,target_region"} \
                  NR==1{ for(i=1; i<=NF; i++) header[$i]=i} \
		  NR>1{ split ($header["Parameter Value[Primer Info]"], primers, ","); \
		  printf "%s,%s/Raw_Sequence_Data/%s,%s,%s,%s\n", $header["Sample Name"], PWD, $header["Raw Data File"], primers[2], primers[4],$header["Parameter Value[Target Molecule]"] }' | \
		  sed -E "s/5'-([A-Z]+)-3'/\1/g" | \
		  sed -E 's/\s+//g'
