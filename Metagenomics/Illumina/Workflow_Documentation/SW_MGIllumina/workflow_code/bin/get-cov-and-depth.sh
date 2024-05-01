#!/usr/bin/env bash
SAMPLE_ID=$1
ASSEMBLY=$2
NT=$3
BAM=$4
PILEUP_MEM=$5


# Only running if the assembly produced anything
if [ -s ${ASSEMBLY} ]; then

    # Only running on genes also if genes were identified
    if [ -s ${NT} ]; then

        pileup.sh -Xmx${PILEUP_MEM} -in ${BAM} \
                    fastaorf=${NT} outorf=${SAMPLE_ID}-gene-cov.tmp \
                    out=${SAMPLE_ID}-contig-cov-and-det.tmp 

        # Filtering coverages based on detection
            # Genes
        grep -v "#" ${SAMPLE_ID}-gene-cov-and-det.tmp | \
        awk -F $'\t' ' BEGIN {OFS=FS} { if ( $10 <= 0.5 ) $4 = 0 } { print \$1,\$4 } ' \
                > ${SAMPLE_ID}-gene-cov.tmp

        cat <( printf "gene_ID\tcoverage\n" ) ${SAMPLE_ID}-gene-cov.tmp \
             > ${SAMPLE_ID}-gene-coverages.tsv

            # Contigs
        grep -v "#" ${SAMPLE_ID}-contig-cov-and-det.tmp | \
        awk -F $'\t' ' BEGIN {OFS=FS} { if ( $5 <= 50 ) $2 = 0 } { print $1,$2 } ' \
            > ${SAMPLE_ID}-contig-cov.tmp

        cat <( printf "contig_ID\tcoverage\n" ) ${SAMPLE_ID}-contig-cov.tmp \
            > ${SAMPLE_ID}-contig-coverages.tsv

        # Removing intermediate files
        rm  ${SAMPLE_ID}-gene-cov-and-det.tmp ${SAMPLE_ID}-contig-cov-and-det.tmp \
            ${SAMPLE_ID}-gene-cov.tmp ${SAMPLE_ID}-contig-cov.tmp

    else

        pileup.sh -in ${BAM} out=${SAMPLE_ID}-contig-cov-and-det.tmp

        # Filtering coverages based on detection
            # Contigs
        grep -v "#" ${SAMPLE_ID}-contig-cov-and-det.tmp | \
        awk -F $'\t' ' BEGIN {OFS=FS} { if ( $5 <= 50 ) $2 = 0 } { print $1,$2 } ' \
            > ${SAMPLE_ID}-contig-cov.tmp
        cat <( printf "contig_ID\tcoverage\n" ) ${SAMPLE_ID}-contig-cov.tmp \
            > ${SAMPLE_ID}-contig-coverages.tsv

        # Writing out empty genes coverage file
        printf "gene_ID\tcoverage\n" > ${SAMPLE_ID}-gene-coverages.tsv
        printf "\n\nGene-level coverage info not recovered because the assembly didn't have any genes identified.\n"

        # Removing intermediate files
        rm ${SAMPLE_ID}-contig-cov-and-det.tmp ${SAMPLE_ID}-contig-cov.tmp

    fi

else

    printf "gene_ID\tcoverage\n" > ${SAMPLE_ID}-gene-coverages.tsv
    printf "contig_ID\tcoverage\n" > ${SAMPLE_ID}-contig-coverages.tsv
    printf "Coverage info not recovered because the assembly didn't produce anything.\n"

fi