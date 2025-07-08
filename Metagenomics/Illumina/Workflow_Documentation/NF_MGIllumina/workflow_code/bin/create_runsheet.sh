#!/usr/bin/env bash

# A script to the the input csv file rtequired by this pipeline when a GLDS accession is provided
# Rather than the required input csv file

ASSAY_TABLE=$1 # a_GLDS-466_metagenome-sequencing_whole-genome-shotgun-sequencing_illumina.txt
awk -v PWD=$PWD -F "\t" '
                 BEGIN{print "sample_id,forward,reverse,paired"} \
		 NR==1{for (i=1; i<=NF; i++) {ix[$i] = i}} \
		 NR>1{gsub(" ", "", $ix["Raw Data File"]); \
	         split($ix["Raw Data File"], reads_path, ","); \
		 gsub("PAIRED","true",$ix["Parameter Value[Library Layout]"]); \
		 gsub("SINGLE","false",$ix["Parameter Value[Library Layout]"]); \
		 printf "%s,%s/Raw_Sequence_Data/%s,%s/Raw_Sequence_Data/%s,%s\n",$ix["Sample Name"], PWD, reads_path[1], PWD,reads_path[2],$ix["Parameter Value[Library Layout]"]}
		' ${ASSAY_TABLE}
