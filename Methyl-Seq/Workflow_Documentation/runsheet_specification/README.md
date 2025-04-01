# Runsheet Specification

## Description

* The Runsheet is a csv file that contains the metadata required for processing MethylSeq datasets through GeneLab's 
  Methyl-Seq consensus processing pipeline. The specified file paths may be either URLs pointing to file locations online
  (as in the single-end example) or absolute paths on the local file system (as in the paired-end example).

## Examples

1. [Runsheet for OSD-47](paired_end_runsheet/OSD-47_methylSeq_v2_runsheet.csv) (paired end DNA methylseq dataset)
2. [Runsheet for OSD-397](single_end_runsheet/OSD-397_methylSeq_v2_runsheet.csv) (single end DNA methylseq dataset)


## Required columns

| Column Name | Type | Description | Example |
|:------------|:-----|:------------|:--------|
| Sample Name | string | Sample Name, added as a prefix to sample-specific processed data output files. Should not include spaces or weird characters. | Mmus_C57-6T_LVR_BSL_Rep1_B1 |
| Study Assay Measurement Type | string | The assay measurement type. Used to determine correct processing mode. Should be one of "DNA Methylation Profiling" or "RNA Methylation Profiling" | DNA methylation profiling |
| Study Assay Technology Type | string | The assay technology type. Used to determine correct processing mode. {Reduced-Representation Bisulfite Sequencing, Whole Genome Bisulfite Sequencing} | Whole Genome Bisulfite Sequencing |
| organism | string | Species name used to map to the appropriate gene annotations file. Supported species can be found in the `species` column of the [GL-DPPD-7110-A_annotations.csv](../../../../GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv) file. | Mus musculus |
| paired_end | bool | Set to True if the samples were sequenced as paired-end. If set to False, samples are assumed to be single-end. | False |
| library_kit | string | Indicates which library kit was used to prepare the samples. Used to determine correct processing mode. Supported values can be found in the `Library Kit` column of the []() file. | EpiGnome |
| read1_path | string (url or local path) | Location of the raw reads file. For paired-end data, this specifies the forward reads fastq.gz file. | /my/data/sample_1.fastq.gz |
| read2_path | string (url or local path) | Location of the raw reads file. For paired-end data, this specifies the reverse reads fastq.gz file. For single-end data, this column should be omitted. | /my/data/sample_2.fastq.gz |
| Factor Value[<name, e.g. Space Flight>] | string | A set of one or more columns specifying the experimental group the sample belongs to. In the simplest form, a column named 'Factor Value[group]' is sufficient. | Space Flight |
| Original Sample Name | string | Used to map the sample name that will be used for processing to the original sample name. This is often identical except in cases where the original name includes spaces or weird characters. | Mmus_C57-6T_LVR_BSL_Rep1_B1 |
