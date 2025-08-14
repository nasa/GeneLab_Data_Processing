# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.4](https://github.com/nasa/GeneLab_Data_Processing/tree/NF_MAAffymetrix_1.0.4/Microarray/Affymetrix/Workflow_Documentation/NF_MAAffymetrix) - 2024-05-17

### Changed

- Fix cache location issues that arose in `quarto render` when using Nextflow v.23.10.1 ([#82](https://github.com/nasa/GeneLab_Data_Processing/issues/82))
- Increase timeout that caused incomplete file download in `read.celfiles()` ([#86](https://github.com/nasa/GeneLab_Data_Processing/issues/86))

## [1.0.3](https://github.com/nasa/GeneLab_Data_Processing/tree/NF_MAAffymetrix_1.0.3/Microarray/Affymetrix/Workflow_Documentation/NF_MAAffymetrix) - 2024-02-26

### Added

- Retry wrapper for functions that utilize internet resources.  This is aimed to reduce failures due solely due to intermittent network issues. (ceb6d9a3)

### Fixed

- Missing Raw Data MA Plots when handling designs that loaded as `ExpressionFeatureSet` objects. (7af7192e)
  - Additionally, future unhandled raw data classes will raise an exception rather than fail to plot silently.

## [1.0.2](https://github.com/asaravia-butler/GeneLab_Data_Processing/tree/NF_MAAffymetrix_1.0.2/Microarray/Affymetrix/Workflow_Documentation/NF_MAAffymetrix) - 2023-05-24

### Added

- Workflow now produces a file called meta.sh (in the 'GeneLab' sub-directory) that contains information about the workflow run. This file is used by the post processing workflow to generate a protocol description. (5a8a255)
- POST_PROCESSING will now generate a protocol description using the contents of meta.sh and text templates. (801e2ad)
- Workflow can now be run using an ISA archive by supplying parameter: 'isaArchivePath' (as either a local path or public web uri) (8822069)

### Changed

- Update dp_tools from 1.3.2 to 1.3.4 (158ce5e)
  - This updates the POST_PROCESSING workflow assay table to join multiple files by ',' instead of ',<SPACE>' and enables max flag code setting.
- Slightly reduced stringency in V&V check for log2fc computation to account for rounding errors, specifically from 99.9% of rows within tolerance to 99.5%. (9fd2c11)
- Publish directory behavior reworked to use the OSD accession as part of the default name. Now uses `resultsDir` instead of `outputDir` as the parameter name when a user does control the published files directory. (97cba72)

### Fixed

- Halt level flags now properly trigger workflow halt. (0885175)
- Boxplots now show all y-axis labels when working with many samples. (7ec10d4s)
- Density plot legend cex (character expansion) now has a minimum of 0.35 (rather than raising an exception for very large numbers of samples) (9a54fdc)

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