# Runsheet Specification

## Description

* The Runsheet is a csv file that contains the metadata required for processing bulk RNA sequence datasets through GeneLab's RNAseq consensus processing pipeline (RCP).


## Examples

1. [Runsheet for GLDS-48](single_end_runsheet/GLDS-48_bulkRNASeq_v1_runsheet.csv) (single end dataset)
2. [Runsheet for GLDS-194](paired_end_runsheet/GLDS-194_bulkRNASeq_v1_runsheet.csv) (paired end dataset)
3. [Runsheet for GLDS-605](paired_end_runsheet/GLDS-605_bulkRNASeq_v2_runsheet.csv) (paired end dataset with technical replicates)



## Required columns

| Column Name | Type | Description | Example |
|:------------|:-----|:------------|:--------|
| Sample Name | string | Sample Name, added as a prefix to sample-specific processed data output files. Should not include spaces or weird characters. | Mmus_BAL-TAL_LRTN_BSL_Rep1_B7 |
| has_ERCC | bool | Set to True if ERCC spike-ins are included in the samples. This ensures ERCC normalized DGE is performed in addition to standard DGE. | True |
| paired_end | bool | Set to True if the samples were sequenced as paired-end. If set to False, samples are assumed to be single-end. | False |
| organism | string | Species name used to map to the appropriate gene annotations file. Supported species can be found in the `species` column of the [GL-DPPD-7110-A_annotations.csv](../../../../../GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv) file. | Mus musculus |
| read1_path | string (url or local path) | Location of the raw reads file. For paired-end data, this specifies the forward reads fastq.gz file. | /my/data/sample_1.fastq.gz |
| read2_path | string (url or local path) | Location of the raw reads file. For paired-end data, this specifies the reverse reads fastq.gz file. For single-end data, this column should be omitted. | /my/data/sample_2.fastq.gz |
| Factor Value[<name, e.g. Spaceflight>] | string | A set of one or more columns specifying the experimental group the sample belongs to. In the simplest form, a column named 'Factor Value[group]' is sufficient. | Space Flight |
| Original Sample Name | string | Used to map the sample name that will be used for processing to the original sample name. This is often identical except in cases where the original name includes spaces or weird characters. | Mmus_BAL-TAL_LRTN_BSL_Rep1_B7 |

## Optional columns

| Column Name | Type | Description | Example |
|:------------|:-----|:------------|:--------|
| Source Name | string | Identifier linking samples. Used for handling technical replicates during differential gene expression analysis. Multiple samples with the same Source Name may be collapsed during analysis depending on the Has Tech Reps setting. | RR3_BSL_B7 |
| Has Tech Reps | bool | Indicates whether this sample is a technical replicate that should be collapsed with other samples sharing the same Source Name. Set to True for technical replicates that should be collapsed, False for distinct samples that should remain separate even if they share a Source Name. | False |
