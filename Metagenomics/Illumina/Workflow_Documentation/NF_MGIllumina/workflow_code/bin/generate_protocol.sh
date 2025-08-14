#!/usr/bin/env bash

# Generate protocol according to a pipeline document

# USAGE:
# generate_protocol.sh <software_versions> <protocol_id>
# EXAMPLE
# generate_protocol.sh ../Metadata/software_versions.txt GL-DPPD-7107-A

FASTQC=`grep -i 'fastqc' $1 | awk '{print $2}' |sed -E 's/v//'`
MULTIQC=`grep -i 'multiqc' $1 | awk '{print $3}'`
BBMAP=`grep -i 'bbtools' $1 | awk '{print $2}'`
HUMANN=`grep -i 'humann' $1 | awk '{print $2}'|sed -E 's/v//'`
MEGAHIT=`grep -i 'megahit' $1 | awk '{print $2}'|sed -E 's/v//'`
PRODIGAL=`grep -i 'prodigal' $1 | awk '{print $2}'|sed -E 's/[vV:]//g'`
CAT=`grep 'CAT' $1 | awk '{print $2}'|sed -E 's/v//'`
KOFAMSCAN=`grep 'exec_annotation' $1 | awk '{print $2}'`
BOWTIE2=`grep -i 'bowtie' $1 | awk '{print $3}'`
SAMTOOLS=`grep -i 'samtools' $1 | awk '{print $2}'`
METABAT2=`grep -i 'metabat' $1 | awk '{print $2}'`
BIT=`grep -i 'bioinformatics tools' $1 | awk '{print $3}' | sed 's/v//' | sed -E 's/.+([0-9]+.[0-9]+.[0-9]+).+/\1/'`
CHECKM=`grep -i 'checkm' $1 | awk '{print $2}' |sed -E 's/v//'`
GTDBTK=`grep -i '^GTDB' $1 | awk '{print $2}' |sed -E 's/v//' | head -n2` # If 2 versions are used, choose the second

PROTOCOL_ID=$2

PROTOCOL="Data were processed as described in ${PROTOCOL_ID} (https://github.com/nasa/GeneLab_Data_Processing/blob/master/Metagenomics/Illumina/Pipeline_GL-DPPD-7107_Versions/${PROTOCOL_ID}.md), using workflow NF_MGIllumina v1.0.0 (https://github.com/nasa/GeneLab_Data_Processing/tree/NF_MGIllumina_1.0.0/Metagenomics/Illumina/Workflow_Documentation/NF_MGIllumina). \
	  In breif, quality assessment of reads was performed with FastQC v${FASTQC} and reports were summarized with MultiQC v${MULTIQC}. \
          Quality trimming and filtering were performed with bbmap v${BBMAP}. Read-based processing was performed with humann3 v${HUMANN}. \
	  Individual samples were assembled with megahit v${MEGAHIT}. Genes were called with prodigal v${PRODIGAL}. \
	  Taxonomic classification of genes and contigs was performed with CAT v${CAT}. Functional annotation was done with KOFamScan v${KOFAMSCAN}. \
	  Reads were mapped to assemblies with bowtie2 v${BOWTIE2} and coverage information was extracted for reads and contigs with samtools v${SAMTOOLS} and bbmap v${BBMAP}. \
	  Binning of contigs was performed with metabat2 v${METABAT2}. Bins were summarized with bit v${BIT} and estimates of quality were generated with checkm v${CHECKM}. \
	  High-quality bins (> 90% est. completeness and < 10% est. redundancy) were taxonomically classified with gtdb-tk v${GTDBTK}."

echo ${PROTOCOL}
