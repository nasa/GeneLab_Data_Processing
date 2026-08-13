# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.5](https://github.com/nasa/GeneLab_Data_Processing/tree/NF_MAAffymetrix_1.0.5/Microarray/Affymetrix/Workflow_Documentation/NF_MAAffymetrix) - 2026-05-XX

### Added

- Support for custom annotations, see [specification](examples/annotations/README.md) ([#113](https://github.com/nasa/GeneLab_Data_Processing/issues/113))
- Add option to skip differential expression analysis (`--skipDE`) ([#104](https://github.com/nasa/GeneLab_Data_Processing/issues/104))
- Add nextflow schema support for parameter validation and help text generation
- Add conda support for easier local development and debugging

### Changed

- Replace `RUNSHEET_FROM_GLDS` and `RUNSHEET_FROM_ISA` processes and their associated workflow logic with a new staging analysis subworkflow supporting both accession-based and input-file-based execution modes
- Rework publish directory behavior as part of the staging analysis subworkflow where `outdir` is now the base directory for `GLDS-NNN/` output directory if `--accession` is provided, or the base directory for `results/` output directory if `--runsheet` is provided
- Bump gl-microarray image from version 1.0.0 to 1.1.0 to match R package updates in the [GL-DPPD-7114-A pipeline document](../../Pipeline_GL-DPPD-7114_Versions/GL-DPPD-7114-A.md)
- Rename `annotation_config_path` as `array_annot_path` and `config.csv` as `design_info.csv` throughout the workflow and documentation to better reflect the purpose of the file and its contents
- Convert `generate_protocol.sh` to a Python script for automated handling of reference/annotation parameters
- Add `create_date` to design_info.csv and parse it in the protocol
- Move protocol creation from post-processing to main nextflow script to make passing needed values easier and more robust
- Update software table generation to exclude `purrr` from table if custom annotations are not used
- Rename module files from UPPERCASE.nf to lowercase.nf following Nextflow community convention
- Flatten directory-based modules (PROCESS_NAME/ with scripts under resources/usr/bin/) to single lowercase process_name.nf files directly under modules/
- Move process scripts from modules/PROCESS_NAME/resources/usr/bin/ to the top-level bin/ directory
- Update processed data protocol to auto-populate workflow version from `nextflow.config` and add Caenorhabditis elegans, Saccharomyces cerevisiae, Escherichia coli, and Pseudomonas aeruginosa to supported organisms ([#98](https://github.com/nasa/GeneLab_Data_Processing/issues/98))
- Fixes in `Affymetrix.qmd`
  - Replace live biomaRt::getBM() queries in the QMD with direct downloads of Ensembl's FTP mart-dump tables; drop chunking/retry/Sys.sleep tied to those queries
  - When renaming column names, specify which columns to rename to avoid unintentional renaming ([#97](https://github.com/nasa/GeneLab_Data_Processing/issues/97))
  - When renaming factor names, prevent cases where a factor is partially renamed because it contains a substring that is another factor ([#100](https://github.com/nasa/GeneLab_Data_Processing/issues/100))
  - Update MA plot to support HTAFeatureSet ([#105](https://github.com/nasa/GeneLab_Data_Processing/issues/105))
  - Remove extra `.1` suffix in AFFY HTA 2 0 Probe IDs in the raw data to allow for merging to BioMart data ([#106](https://github.com/nasa/GeneLab_Data_Processing/issues/106))
  - Decrease legend size when sample names are long to prevent it from covering plot ([#107](https://github.com/nasa/GeneLab_Data_Processing/issues/107))
  - Simplify group sample retrieval during differential expression group-wise statistics computation to use a more concise `filter/pull/sort` chain instead of `group_by/summarize/filter/pull`, addressing the deprecation warning in dplyr >= 1.1.0 where returning more than 1 row per `summarise()` group is deprecated
- Changes to post-processing workflow
  - Resolve output directory `GLDS-NNN/` or `results/` to match main workflow behavior
  - Replace dp_tools dependency in assay table update and md5sum table generation with standalone scripts
  - Rename `UPDATE_ISA_TABLES` and `update_curation_table.py` to `UPDATE_ASSAY_TABLE` and `update_assay_table.py` to better reflect their purpose
  - Add new PURGE_PROCESSING_INFO Nextflow module to strip full paths in nextflow_processing_info_GLmicroarray.txt before publishing
  - Add parameter validation and summary log from nf-schema

### Removed

- Packages `R.utils`, and `biomaRt` are no longer used in the processing code, and have been removed from software table generation

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