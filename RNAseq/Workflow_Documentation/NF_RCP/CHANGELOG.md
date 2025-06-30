# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.1] - (https://github.com/nasa/GeneLab_Data_Processing/tree/NF_RCP_2.0.1/RNAseq/Workflow_Documentation/NF_RCP) - 2025-06-26

### Fixed

- Fixed fastqc metrics extraction in `parse_multiqc.py` script 
  - Added qc file validation output listing missing entries
  - Updated multiqc parsing for fastqc metrics

## [2.0.0](https://github.com/nasa/GeneLab_Data_Processing/tree/NF_RCP_2.0.0/RNAseq/Workflow_Documentation/NF_RCP) - 2025-04-10

### Added

- Prokaryotes pipeline support via `--microbes` parameter:  
  - Reads are aligned to a reference genome using Bowtie 2 rather than STAR, and gene counts are quantified using featureCounts instead of RSEM. Other steps remain unchanged.  
  - Added software versions:  
    - Bowtie 2 2.5.4  
    - subread 2.0.8  
- Read alignment now outputs unaligned reads as FASTQ files.  
- Added Variance-stabilizing transformation (VST) counts table.**  
- Incorporated rRNA removal into gene counts and differential gene expression (DGE) analysis.  
  - Separate results are generated for rRNA-removed DGE analysis, with new output directories:  
    - `04-DESeq2_NormCounts_rRNArm/`  
    - `05-DESeq2_DGE_rRNArm/`
- Added reference table support for Pseudomonas aeruginosa [#37](https://github.com/nasa/GeneLab_Data_Processing/issues/37)
- Added V&V check for adapter content removal using FastQC/MultiQC reports from trimmed reads [#42](https://github.com/nasa/GeneLab_Data_Processing/issues/42)
- Added generation of a CSV file summarizing parsed metrics from tool logs and MultiQC reports [#84](https://github.com/nasa/GeneLab_Data_Processing/issues/84)
- Added support for user-specified custom genomes that are not present in the provided GeneLab Annotation Reference table [#157](https://github.com/nasa/GeneLab_Data_Processing/issues/157)

### Changed

- Updated software versions:
  - FastQC 0.12.1
  - MultiQC 1.26
  - Cutadapt 4.2
  - TrimGalore! 0.6.10
  - STAR 2.7.11b
  - RSEM 1.3.3
  - Samtools 1.2.1
  - gtfToGenePred/genePredToBed 469
  - RSeQC tools 5.0.4
  - R 4.4.2
  - Bioconductor 3.20
  - BiocParallel 1.40.0
  - DESeq2 1.46.0
  - tximport 1.34.0
  - tidyverse 2.0.0
  - dplyr 1.1.4
  - knitr 1.49
  - stringr 1.5.1
  - yaml 2.3.10
  - dp_tools 1.3.8
  - pandas 2.2.3
  - seaborn 0.13.2
  - matplotlib 3.10.0
  - numpy 2.2.1
  - scipy 1.15.1
- Updated [Ensembl Reference Files](../../../GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv) now use: 
  - Animals: Ensembl release 112
  - Plants: Ensembl plants release 59
  - Bacteria: Ensembl bacteria release 59
- Added "_GLbulkRNAseq" suffix to output files
- RSeQC inner_distance minimum value now dynamically set based on read length
- DESeq2 analysis now handles technical replicates [#32](https://github.com/nasa/GeneLab_Data_Processing/issues/32)
- MultiQC reports replaced with separate data zip and html files
- Increased default memory allocation for the STAR alignment process to 40GB [#36](https://github.com/nasa/GeneLab_Data_Processing/issues/36)

### Fixed

- DGE validation script (`vv_dge_deseq2.py`) error with all-integer sample names [#112](https://github.com/nasa/GeneLab_Data_Processing/issues/112)
- The `--accession` parameter (formerly `--gldsAccession`) is now optional for runsheet-based workflows; if omitted, outputs default to the 'results' directory [#35](https://github.com/nasa/GeneLab_Data_Processing/issues/35)
- Metadata/ISA.zip is now optional when running post-processing workflow [#156](https://github.com/nasa/GeneLab_Data_Processing/issues/156)
- Add fallback for vst() function in case the representative subsampling implemented in vst() fails due to sparse matrix input as in testing runs. 


### Removed

- ERCC-normalized DGE analysis and associated output files
- GeneLab visualization output tables [#41](https://github.com/nasa/GeneLab_Data_Processing/issues/41)

## [1.0.4](https://github.com/nasa/GeneLab_Data_Processing/tree/NF_RCP-F_1.0.4/RNAseq/Workflow_Documentation/NF_RCP-F) - 2024-02-08

### Fixed

- Workflow usage files will all follow output directory set by workflow user
- ERCC Notebook:
  - Moved gene prefix definition to start of notebook
  - Added fallback for scenarios where every gene has zeros: use "poscounts" estimator to calculate a modified geometric mean
  - Reordered box-whisker plots from descending to ascending reference concentration order, ordered bar plots similarly
  
### Changed

- TrimGalore! will now use autodetect for adaptor type [#20](https://github.com/nasa/GeneLab_Data_Processing/issues/20)
- V&V migrated from dp_tools version 1.1.8 to 1.3.4 including:
  - Migration of V&V protocol code to this codebase instead of dp_tools
  - Fix for sample wise checks reusing same sample
- Added '_GLbulkRNAseq' to output file names

## [1.0.3](https://github.com/nasa/GeneLab_Data_Processing/tree/NF_RCP-F_1.0.3/RNAseq/Workflow_Documentation/NF_RCP-F) - 2023-01-25

### Added

- Test coverage using [nf-test](https://github.com/askimed/nf-test) approach

### Changed

- Updated software versions (via container update)
  - tximport == 1.27.1

### Fixed

- 'ERCC Non detection causes non-silent error' #65
- 'This function is not compatible with certain updated ISA archive metadata filenaming' #56
- 'Groups can become misassigned during group statistic calculation' #55
- 'sample to filename mapping fails when sample names are prefix substrings of other sample names' #60
- Fixed Singularity specific container issue related to DESeq2 steps

## [1.0.2](https://github.com/nasa/GeneLab_Data_Processing/tree/NF_RCP-F_1.0.2/RNAseq/Workflow_Documentation/NF_RCP-F) - 2022-11-30

### Added

- Manual tool version reporting functionality for [script](https://github.com/nasa/GeneLab_Data_Processing/tree/NF_RCP-F_1.0.2/RNAseq/Workflow_Documentation/NF_RCP-F/workflow_code/bin/format_software_versions.py) that consolidates tool versions for full workflow.
  - Currently includes manual version reporting for gtfToGenePred and genePredToBed

### Fixed

- Updated Cutadapt version in workflow from 3.4 to 3.7 in accordance with pipeline [specification](https://github.com/nasa/GeneLab_Data_Processing/tree/NF_RCP-F_1.0.2/RNAseq/Pipeline_GL-DPPD-7101_Versions/GL-DPPD-7101-F.md)

## [1.0.1](https://github.com/nasa/GeneLab_Data_Processing/tree/NF_RCP-F_1.0.1/RNAseq/Workflow_Documentation/NF_RCP-F) - 2022-11-17

### Changed

- Updated to dp_tools version 1.1.8 from 1.1.7: This addresses api changes from the release of the [OSDR](https://osdr.nasa.gov/bio/)

### Removed

- Docs: Recommendation to use Nextflow Version 21.10.6 removed as newer stable releases address original issue that had merited the recommendation

## [1.0.0](https://github.com/nasa/GeneLab_Data_Processing/tree/NF_RCP-F_1.0.0/RNAseq/Workflow_Documentation/NF_RCP-F) - 2022-11-04

### Added

- First internal production ready release of the RNASeq Consensus Pipeline Nextflow Workflow
