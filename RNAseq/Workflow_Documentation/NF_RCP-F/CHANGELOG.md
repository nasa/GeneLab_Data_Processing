# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Github action support for CI testing (a83006ba91b1209e1857fefd96e9ff950ebb0cdc)

### Fixed
- Workflow usage files will all follow output directory set by workflow user (3e69f06432f62b7924d2e043ef4768c5d09bf614)
### Changed 
- TrimGalore! will now use autodetect for adaptor type (3b7e0bab4017e90481359c48f9cf7c8837ed54d2)
- V&V migrated from dp_tools version 1.1.8 to 1.3.2 including:
  - Migration of V&V protocol code to this codebase instead of dp_tools (b3684a4c1db5df06eab20916ef7e130c410c147c)
  - Fix for sample wise checks reusing same samples (dca4fdad7518ac9ead3ee2e4c5f57ac0fe25c715)

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