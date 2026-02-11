#!/usr/bin/env bash

KRAKEN2=$(grep -i 'kraken2' $1 | sort -u |awk '{print $2}' |sed -E 's/v//')

HOST=$2 # ex: mouse
REFSEQ_ID=$3
GENOME=$4

PROTOCOL="FASTQ files were screened using Kraken2 v"$KRAKEN2" against a "$HOST" genome reference database that was constructed from NCBI's RefSeq ("$REFSEQ_ID") "$GENOME", with --no-masking, and defaults of kmer-length 35 and minimizer-length of 31, and is fully compatible with Kraken2 v"$KRAKEN2". Reads classified as host were removed after being quantified and expressed as a percentage of total reads."

echo ${PROTOCOL}