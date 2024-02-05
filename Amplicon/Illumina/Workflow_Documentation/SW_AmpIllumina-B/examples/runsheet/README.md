# Runsheet Specification

## Description

* The Runsheet is a csv file that contains the metadata required for processing amplicon sequencing datasets through GeneLab's GeneLab Illumina amplicon sequencing data processing pipeline (AmpIllumina).


## Examples

1. [Runsheet for GLDS-487](paired_end_runsheet/GLDS-487_amplicon_v0_runsheet.csv) (paired end dataset)
2. [Sample runsheet for single end sequencing](single_end_runsheet/single_end_example.csv) (template for single end dataset)


## Required columns
> *Note: Avoid using sample and group names that begin with a number or weird character.*

| Column Name                      | Type                     | Description                                                                                     | Example                                                                                          |
|:---------------------------------|:-------------------------|:------------------------------------------------------------------------------------------------|:-------------------------------------------------------------------------------------------------|
| Sample Name                      | string                   | A unique identifier for each sample. Should not include spaces or weird characters.                                                           | Mmus_C57-6T_FCS_BSL_LAR_Rep1_B1                                                                                      |
| Parameter Value[Library Selection] | string                   | Library type, either '16S' or 'ITS'.                                                            | 16S                                                                                              |
| paired_end                       | bool                     | Set to True if the samples were sequenced as paired-end. If set to False, samples are assumed to be single-end.                  | True                                                                                             |
| F_Primer                         | string                   | Sequence of the forward primer. For single-end data, only this primer is required.              | GTGCCAGCMGCCGCGGTAA                                                                              |
| R_Primer                         | string                   | Sequence of the reverse primer. Required only for paired-end data.                              | GGACTACHVGGGTWTCTAAT                                                                             |
| read1_path                       | string (url or local path) | Location of the raw reads file. For paired-end data, this specifies the forward reads fastq.gz file. All local files should be in the same directory and should have unique identifiers to distinguish between different samples.                                                   | /my/data/sample1_R1_raw.fastq.gz                                |
| read2_path                       | string (url or local path) | Location of the raw reads file. For paired-end data, this specifies the reverse reads fastq.gz file. The names for forward and reverse read files should be identical except for their R1 and R2 suffixes, which distinguish between the two. For single-end data, this column should be omitted.                          | /my/data/sample1_R2_raw.fastq.gz                               |
| raw_R1_suffix                    | string                   | Suffix identifying the forward raw read files.                                                  | _R1_raw.fastq.gz                                                                                 |
| raw_R2_suffix                    | string                   | Suffix identifying the reverse raw read files. For single-end data, this column should be omitted.                   | _R2_raw.fastq.gz                                                                                 |
| Factor Value[<name, e.g. Spaceflight>]        | string                   | A set of one or more columns specifying the experimental group the sample belongs to. In the simplest form, a column named 'Factor Value[group]' is sufficient.                                         | Ground Control                                                                                     |
| Factor Value[<name, e.g. Dosage>]   | string                   | A set of one or more columns specifying the experimental group the sample belongs to. In the simplest form, a column named 'Factor Value[group]' is sufficient.                                      | Low                                                                             |
| groups                           | string                   | Define the experimental groups for each sample, using ' & ' to differentiate between factor levels. Use a consistent order of factors across all samples.     | Ground Control & Low                                                              |
