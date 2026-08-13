# Runsheet Specification

## Description

* The Runsheet is a csv file that contains the metadata required for processing Affymetrix datasets through GeneLab's Affymetrix processing pipeline.


## Examples

1. [Runsheet for GLDS-3](OSD-3_microarray_v0_runsheet.csv)
2. [Runsheet for GLDS-213](OSD-213_microarray_v0_runsheet.csv)


## Required columns

| Column Name | Type | Description | Example |
|:------------|:-----|:------------|:--------|
| Sample Name | string | Sample Name, added as a prefix to sample-specific processed data output files. Should not include spaces or weird characters. | Dmel_Hml-GAL4-UAS-GFP_wo_FLT_3rd-Instar-Larva_Rep1 |
| biomart_attribute | string | A bioMart attribute identifier denoting the microarray probe/probeset attribute used for annotation mapping. | AFFY Drosophila 2 |
| organism | string | Species name used to map to the appropriate gene annotations file. | Drosophila melanogaster |
| Array Data File Path | string (url or local path) | Location of the raw data file for the sample. | /my/data/sample_1.CEL |
| Factor Value[<name, e.g. Spaceflight>] | string | A set of one or more columns specifying the experimental group the sample belongs to. In the simplest form, a column named 'Factor Value[group]' is sufficient. | Space Flight |
| Original Sample Name | string | Used to map the sample name that will be used for processing to the original sample name. This is often identical except in cases where the original name includes spaces or weird characters. | Dmel_Hml-GAL4-UAS-GFP_wo_FLT_3rd-Instar-Larva_Rep1 |
