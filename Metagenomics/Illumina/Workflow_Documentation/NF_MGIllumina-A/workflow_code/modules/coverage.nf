#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
//params.pileup_mem = "5g"

/*
    This process pulls out coverage and detection information for each sample, gene-level and contig-level,
    and filters the gene-level coverage information based on requiring at least 50% detection.
*/

process GET_COV_AND_DET {

    tag "Calculating gene and contig coverage for ${sample_id}..."

    input:
        tuple val(sample_id), path(bam), path(assembly), path(aa), path(nt)
    output:
        // Gene_covs and contig_covs
        tuple val(sample_id), path("${sample_id}-gene-coverages.tsv"), path("${sample_id}-contig-coverages.tsv"), emit: coverages
        path("versions.txt"), emit: version 
    script:
        """
        # Only running if the assembly produced anything
        if [ -s ${assembly} ]; then

            # Only running on genes also if genes were identified
            if [ -s ${nt} ]; then

                pileup.sh -Xmx${params.pileup_mem} -in ${bam} \\
                          fastaorf=${nt} outorf=${sample_id}-gene-cov-and-det.tmp \\
                          out=${sample_id}-contig-cov-and-det.tmp 

                # Filtering coverages based on detection
                    # Genes
                grep -v "#" ${sample_id}-gene-cov-and-det.tmp | \\
                awk -F \$'\\t' ' BEGIN {OFS=FS} { if ( \$10 <= 0.5 ) \$4 = 0 } { print \$1,\$4 } ' \\
                      > ${sample_id}-gene-cov.tmp

                cat <( printf "gene_ID\\tcoverage\\n" ) ${sample_id}-gene-cov.tmp > ${sample_id}-gene-coverages.tsv

                    # Contigs
                grep -v "#" ${sample_id}-contig-cov-and-det.tmp | \\
                awk -F \$'\\t' ' BEGIN {OFS=FS} { if ( \$5 <= 50 ) \$2 = 0 } { print \$1,\$2 } ' \\
                    > ${sample_id}-contig-cov.tmp

                cat <( printf "contig_ID\\tcoverage\\n" ) ${sample_id}-contig-cov.tmp > ${sample_id}-contig-coverages.tsv

                # Removing intermediate files
                rm  ${sample_id}-gene-cov-and-det.tmp ${sample_id}-contig-cov-and-det.tmp \\
                   ${sample_id}-gene-cov.tmp ${sample_id}-contig-cov.tmp

            else

                pileup.sh -in ${bam} out=${sample_id}-contig-cov-and-det.tmp

                # Filtering coverages based on detection
                    # Contigs
                grep -v "#" ${sample_id}-contig-cov-and-det.tmp | \\
                awk -F \$'\\t' ' BEGIN {OFS=FS} { if ( \$5 <= 50 ) \$2 = 0 } { print \$1,\$2 } ' \\
                    > ${sample_id}-contig-cov.tmp
                cat <( printf "contig_ID\\tcoverage\\n" ) ${sample_id}-contig-cov.tmp > ${sample_id}-contig-coverages.tsv

                # Writing out empty genes coverage file
                printf "gene_ID\\tcoverage\\n" > ${sample_id}-gene-coverages.tsv
                printf "\\n\\nGene-level coverage info not recovered because the assembly didn't have any genes identified.\\n"

                # Removing intermediate files
                rm ${sample_id}-contig-cov-and-det.tmp ${sample_id}-contig-cov.tmp

            fi

        else

            printf "gene_ID\\tcoverage\\n" > ${sample_id}-gene-coverages.tsv
            printf "contig_ID\\tcoverage\\n" > ${sample_id}-contig-coverages.tsv
            printf "Coverage info not recovered because the assembly didn't produce anything.\\n"

        fi
        VERSION=`bbversion.sh`
        echo "bbtools \${VERSION}" > versions.txt
        """
}
