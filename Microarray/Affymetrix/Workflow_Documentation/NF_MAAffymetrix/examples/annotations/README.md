# Custom Annotations Specification

## Description

* If using custom gene annotations when processing Affymetrix datasets through GeneLab's Affymetrix processing pipeline, a csv file named `config.csv` must be provided as specified below.
* Both the `config.csv` and custom annotations files must be placed in the directory specified by `local_annotation_dir` in the pipeline.


## Example

- [config.csv](config.csv)


## Required columns

| Column Name | Type | Description | Example |
|:------------|:-----|:------------|:--------|
| array_design | string | A bioMart attribute identifier denoting the microarray probe/probeset attribute used for annotation mapping. | AFFY E coli Genome 2 0 |
| annot_type | string | Used to determine how the custom annotations are parsed before merging to the data. Currently, only the below are supported: <ul><li>`3prime-IVT`: Annotations file is expected to be in the format of the 3' IVT expression analysis arrays annotations by [Thermo Fisher](https://www.thermofisher.com/us/en/home/life-science/microarray-analysis/microarray-data-analysis/genechip-array-annotation-files.html)</li><li>`custom`: Annotations file is merged as is, expected to have the following columns: `ProbesetID`, `ENTREZID`, `SYMBOL`, `GENENAME`, `ENSEMBL`, `REFSEQ`, `GOSLIM_IDS`, `STRING_id`, `count_gene_mappings`, `gene_mapping_source`</li></ul> | 3prime-IVT |
| annot_filename | string | Name of the custom annotations file. | E_coli_2.na36.annot.csv |
