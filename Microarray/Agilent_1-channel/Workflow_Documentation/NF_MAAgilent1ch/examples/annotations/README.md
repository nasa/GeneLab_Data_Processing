# Custom Annotations Specification

## Description

* If using custom gene annotations when processing Agilent 1-channel datasets through GeneLab's Agilent 1-channel processing pipeline, a csv config file must be provided as specified below.
* See [Agilent_array_annotations.csv](../Array_Annotations/Agilent_array_annotations.csv) for the latest config file used at GeneLab.


## Example

- [config.csv](config.csv)


## Required columns

| Column Name | Type | Description | Example |
|:------------|:-----|:------------|:--------|
| array_design | string | A bioMart attribute identifier denoting the microarray probe/probeset attribute used for annotation mapping. | AGILENT SurePrint G3 GE 8x60k v3 |
| annot_type | string | Used to determine how the custom annotations are parsed before merging to the data. Currently, only the below are supported: <ul><li>`agilent`: Annotations file is expected to be in the AA (All Annotations) format by [Agilent](https://earray.chem.agilent.com/earray/)</li><li>`custom`: Annotations file is merged as is, expected to have the following columns: `ProbesetID`, `ENTREZID`, `SYMBOL`, `GENENAME`, `ENSEMBL`, `REFSEQ`, `GOSLIM_IDS`, `STRING_id`, `count_gene_mappings`, `gene_mapping_source`</li></ul> | agilent |
| annot_filename | string | Name of the custom annotations file. | 072363_D_AA_20240521.txt |
