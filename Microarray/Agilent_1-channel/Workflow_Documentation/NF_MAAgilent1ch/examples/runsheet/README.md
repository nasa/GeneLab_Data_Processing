# Runsheet Specification

## Description

* The Runsheet is a csv file that contains the metadata required for processing Agilent 1 Channel datasets through GeneLab's Agilent 1 Channel processing pipeline.


## Examples

1. [Runsheet for GLDS-367](GLDS-367_microarray_v0_runsheet.csv)
2. [Runsheet for GLDS-548](OSD-548_microarray_v0_runsheet.csv)


## Required columns

| Column Name | Type | Description | Example |
|:------------|:-----|:------------|:--------|
| Sample Name | string | Sample Name, added as a prefix to sample-specific processed data output files. Should not include spaces or weird characters. | Mmus_BAL-TAL_LRTN_BSL_Rep1_B7 |
| biomart_attribute | string | A bioMart attribute identifier denoting the microarray probe/probeset attribute used for annotation mapping. | agilent_020186 |
| organism | string | Species name used to map to the appropriate gene annotations file. | Mus musculus |
| Array Data File Path | string (url or local path) | Location of the raw data file for the sample. | /my/data/sample_1.txt.gz |
| Factor Value[<name, e.g. Spaceflight>] | string | A set of one or more columns specifying the experimental group the sample belongs to. In the simplest form, a column named 'Factor Value[group]' is sufficient. | Space Flight |
| Original Sample Name | string | Used to map the sample name that will be used for processing to the original sample name. This is often identical except in cases where the original name includes spaces or weird characters. | Mmus_BAL-TAL_LRTN_BSL_Rep1_B7 |
