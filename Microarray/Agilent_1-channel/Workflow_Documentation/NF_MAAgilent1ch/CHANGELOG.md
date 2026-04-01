# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.5](https://github.com/nasa/GeneLab_Data_Processing/tree/NF_MAAgilent1ch_1.0.5/Microarray/Agilent_1-channel/Workflow_Documentation/NF_MAAgilent1ch) - 2025-05-xx

### Added

- Support for custom annotations, see [specification](examples/annotations/README.md)
- Add option to skip differential expression analysis (`--skipDE`) ([#104](https://github.com/nasa/GeneLab_Data_Processing/issues/104))
- Add a retry wrapper for functions that utilize internet resources, syncing with [NF_MAAffymetrix_1.0.3](https://github.com/nasa/GeneLab_Data_Processing/tree/NF_MAAffymetrix_1.0.3/Microarray/Affymetrix/Workflow_Documentation/NF_MAAffymetrix)
- Workflow can now be run using an ISA archive by supplying parameter: 'isaArchivePath' (as either a local path or public web uri), syncing with [NF_MAAffymetrix_1.0.2](https://github.com/asaravia-butler/GeneLab_Data_Processing/tree/NF_MAAffymetrix_1.0.2/Microarray/Affymetrix/Workflow_Documentation/NF_MAAffymetrix)

### Changed

- Publish directory behavior reworked to use the OSD accession as part of the default name. Now uses `resultsDir` instead of `outputDir` as the parameter name when a user does control the published files directory. Syncing with [NF_MAAffymetrix_1.0.2](https://github.com/asaravia-butler/GeneLab_Data_Processing/tree/NF_MAAffymetrix_1.0.2/Microarray/Affymetrix/Workflow_Documentation/NF_MAAffymetrix)
- Small bug fixes in `Agile1CMP.qmd`
  - Update the custom `fetch_organism_specific_annotation_table()` function, used when loading organism-specific annotation metadata, to convert figshare ndownloader URLs to direct API endpoints, as ndownloader URLs require redirect handling that is not supported in all programmatic download contexts
  - Simplify group sample retrieval during differential expression group-wise statistics computation to use a more concise `filter/pull/sort` chain instead of `group_by/summarize/filter/pull`, addressing the deprecation warning in dplyr >= 1.1.0 where returning more than 1 row per `summarise()` group is deprecated

## [1.0.4](https://github.com/nasa/GeneLab_Data_Processing/tree/NF_MAAgilent1ch_1.0.4/Microarray/Agilent_1-channel/Workflow_Documentation/NF_MAAgilent1ch) - 2024-10-02

### Added

- Add automatic generation of processed data protocol ([#85](https://github.com/nasa/GeneLab_Data_Processing/issues/85))

### Changed

- Small bug fixes in `Agile1CMP.qmd`
  - Check if `getBM()` returned results before concatenating it to dataframe to avoid error in `bind_rows()` ([#96](https://github.com/nasa/GeneLab_Data_Processing/issues/96))
  - When renaming column names, specify which columns to rename to avoid unintentional renaming ([#97](https://github.com/nasa/GeneLab_Data_Processing/issues/97))
  - When renaming factor names, prevent cases where a factor is partially renamed because it contains a substring that is another factor ([#100](https://github.com/nasa/GeneLab_Data_Processing/issues/100))
- Update software table generation to exclude `R.utils` from table if data files are not compressed ([#99](https://github.com/nasa/GeneLab_Data_Processing/issues/99))

## [1.0.3](https://github.com/nasa/GeneLab_Data_Processing/tree/NF_MAAgilent1ch_1.0.3/Microarray/Agilent_1-channel/Workflow_Documentation/NF_MAAgilent1ch) - 2024-05-17

### Changed

- Fix cache location issues that arose in `quarto render` when using Nextflow v.23.10.1 ([#82](https://github.com/nasa/GeneLab_Data_Processing/issues/82))

## [1.0.2](https://github.com/nasa/GeneLab_Data_Processing/tree/NF_MAAgilent1ch_1.0.2/Microarray/Agilent_1-channel/Workflow_Documentation/NF_MAAgilent1ch) - 2023-04-28

### Added

- Support for Arabidposis Thaliana datasets using the plants ensembl FTP server.

### Changed

- When encountering error about column reordering, the expected order is saved for debugging purposes.
- Post Processing Workflow: Assay Table Update now added '_array_' prefix to processed files instead of '_microarray_' prefix.

## [1.0.1](https://github.com/asaravia-butler/GeneLab_Data_Processing/tree/NF_MAAgilent1ch_1.0.1/Microarray/Agilent_1-channel/Workflow_Documentation/NF_MAAgilent1ch) - 2023-03-31

### Removed

- Deprecated column renaming code (abcd380)

### Fixed

- Bumped dp_tools from 1.3.0 to 1.3.1 to address 'ISO-8859-1' encoded ISA archive files (example: OSD-271-v2) (d518f40)
- Added handling for raw data that lacks the ProbeUID column (example: OSD-271-v2) (efbc237)

### Changed

- Reordering error message is now more informative (007e36c)

## [1.0.0](https://github.com/asaravia-butler/GeneLab_Data_Processing/tree/NF_MAAgilent1ch_1.0.0/Microarray/Agilent_1-channel/Workflow_Documentation/NF_MAAgilent1ch) - 2023-03-22

### Added

- First internal production ready release of the Agilent 1 Channel Microarray Processing Workflow