# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0](https://github.com/nasa/GeneLab_AmpliconSeq_Workflow/tree/NF_AmpIllumina_1.0.0) 

### Added

- ANCOMBC1 differential abundance analysis
- ANCOMBC2 differential abundance analysis
- Added workflow parameters for finer-grained control of diversity and differential abundance analyses
- Improved error handling and robustness in R scripts, plots
- Software version tracking
- Added software versions:
  - ANCOMBC 2.8.0
  - broom 1.0.7
  - DescTools 0.99.59
  - dp_tools 1.3.8
  - FSA 0.9.6
  - ggdendro 0.2.0
  - glue 1.8.0
  - hexbin 1.28.3
  - mia 1.14.0
  - taxize 0.10.0
  - vsn 3.74.0
- Added reference database support for PR2 v4.13
- Added persistent reference links to DECIPHER databases on Figshare

### Changed

- Workflow converted from Snakemake to Nextflow
- Updated software versions:
  - MultiQC 1.27.1
  - Cutadapt 5.0
  - R-base 4.4.2
  - DADA2 1.34.0
  - DECIPHER 3.20.0
  - biomformat 1.34.0
  - DESeq2 1.46.0
  - ggrepel 0.9.6
  - phyloseq 1.50.0
- Updated database versions:
  - SILVA SSU r138_2 2024
  - UNITE v2024 April2024

<br>

---

> ***Note:** All previous workflow changes associated with the previous version of the GeneLab Amplicon Pipeline
[GL-DPPD-7104-B](https://github.com/nasa/GeneLab_Data_Processing/blob/master/Amplicon/Illumina/Pipeline_GL-DPPD-7104_Versions/GL-DPPD-7104-B.md) and can be found in the
[change log of the Snakemake workflow (SW_AmpIllumina-B)](https://github.com/nasa/GeneLab_Data_Processing/blob/master/Amplicon/Illumina/Workflow_Documentation/SW_AmpIllumina-B/CHANGELOG.md).*
