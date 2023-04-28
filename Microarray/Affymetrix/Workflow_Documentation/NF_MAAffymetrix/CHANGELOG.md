# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.1](https://github.com/asaravia-butler/GeneLab_Data_Processing/tree/NF_MAAffymetrix_1.0.1/Microarray/Affymetrix/Workflow_Documentation/NF_MAAffymetrix) - 2023-04-28

### Added

- Support for Arabidposis Thaliana datasets using the plants ensembl FTP server.
- Support for raw data FeatureSets (building on existing support for ExpressionSets)
- Better support for non-ascii characters in the runsheet, usually caused by such characters in the original ISA archive the runsheet is generated from.

### Fixed

- Typos related to shared code with Agilent 1 Channel platform.

### Changed

- Error message when encountering unique columns when reordering tables is now clearer about what unique columns were found.
- Post Processing Workflow: Assay Table Update now added '_array_' prefix to processed files instead of '_microarray_' prefix.

## [1.0.0](https://github.com/asaravia-butler/GeneLab_Data_Processing/tree/NF_MAAffymetrix_1.0.0/Microarray/Affymetrix/Workflow_Documentation/NF_MAAffymetrix) - 2023-04-24

### Added

- First internal production ready release of the Affymetrix Microarray Pipeline Nextflow Workflow