# GeneLab processing commands for GLDS-379
The GLDS version 3 of this dataset was processed with [GL-DPPD-7101-E](../../../Pipeline_GL-DPPD-7101_Versions/GL-DPPD-7101-E.md).

The raw and processed data are available from: [https://genelab-data.ndc.nasa.gov/genelab/accession/GLDS-379/](https://genelab-data.ndc.nasa.gov/genelab/accession/GLDS-379/)

These samples were processed in batch using the SLURM job scheduler. The following sub-directories contain the scripts used to generate each level of processed data.
  - [00-RawData](00-RawData): Scripts used to generate raw fastQC and multiQC reports
  - [01-TG_PreProc](01-TG_Preproc): Scripts used to generate trimmed reads, trimmed fastQC and multiQC reports
  - [02-STAR_Alignment](02-STAR_Alignment): Scripts used to generate STAR reference and STAR alignment files
  - [03-RSEM_Counts](03-RSEM_Counts): Scripts used to generate RSEM reference and RSEM count files
  - [04-05-DESeq2_NormCounts_DGE](04-05-DESeq2_NormCounts_DGE): Scripts used to generate raw and normlaized counts tables, differential gene expression (DGE), and DGE tables with annotations
  - [RSeQC_analyses](RSeQC_analyses): Scripts used to generate the RSeQC analyses: geneBody coverage, infer experiment, inner distance, and read distribution
  - [ERCC_Analysis](ERCC_Analysis): Jupyter Notebook used to evaluate the ERCC spike-in data

