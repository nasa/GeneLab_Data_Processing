# CSV Sample File Specification

## Description

* The sample file is a comma-separated file that contains the metadata required for processing 
metagenomics sequence datasets through GeneLab's GeneLab Illumina metagenomics sequencing data 
processing pipeline (MGIllumina).


## Examples

1. Samplefile for an example [paired-end dataset](paired_end_dataset/PE_file.csv)
2. Samplefile for an example [single-end dataset](single_end_dataset/SE_file.csv)


## Required columns

| Column Name | Type | Description | Example |
|:------------|:-----|:------------|:--------|
| sample_id | string | Unique Sample Name, added as a prefix to sample-specific processed data output files. Should not include spaces or weird characters. | RR23_FCS_FLT_F1 |
| forward | string (local path) | Location of the raw reads file. For paired-end data, this specifies the forward reads fastq.gz file. | /my/data/sample1_R1_HRremoved_raw.fastq.gz |
| reverse | string (local path) | Location of the raw reads file. For paired-end data, this specifies the reverse reads fastq.gz file. For single-end data, this column should be omitted. | /my/data/sample1_R2_HRremoved_raw.fastq.gz |
| paired | bool | Set to True if the samples were sequenced as paired-end. If set to False, samples are assumed to be single-end. | False |