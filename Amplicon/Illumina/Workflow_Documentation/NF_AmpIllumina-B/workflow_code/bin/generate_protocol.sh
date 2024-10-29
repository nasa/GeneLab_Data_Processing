#!/usr/bin/env bash

FASTQC=`grep -i 'fastqc' $1 | awk '{print $2}' |sed -E 's/v//'`
MULTIQC=`grep -i 'multiqc' $1 | awk '{print $3}' |sed -E 's/v//'`
DADA=`grep -i 'dada2' $1 | awk '{print $2}' |sed -E 's/v//'`
DECIPHER=`grep -i 'decipher' $1 | awk '{print $2}' |sed -E 's/v//'`
CUTADAPT=`grep -i 'cutadapt' $1 | awk '{print $2}' |sed -E 's/v//'`


PROTOCOL_ID=$2 # GL-DPPD-7104-A

PROTOCOL="Data were processed as described in ${PROTOCOL_ID} (https://github.com/nasa/GeneLab_Data_Processing/blob/master/Amplicon/Illumina/Pipeline_GL-DPPD-7104_Versions/${PROTOCOL_ID}.md), \
          using workflow NF_AmpIllumina v1.0.0 (https://github.com/nasa/GeneLab_Data_Processing/tree/NF_AmpIllumina_1.0.0/Amplicon/Illumina/Workflow_Documentation/NF_AmpIllumina). \
		  Quality assessment of reads was performed with FastQC v${FASTQC} and reports were summarized with MultiQC v${MULTIQC}. Primers were removed from raw reads using cutadapt v${CUTADAPT}. \
		  DADA2 v${DADA} was utilized for quality-trimming, filtering and inference of amplicon sequence variants (ASVs), and taxonomy was assigned with DECIPHER v${DECIPHER} against the SILVA r138 reference database."

echo ${PROTOCOL}