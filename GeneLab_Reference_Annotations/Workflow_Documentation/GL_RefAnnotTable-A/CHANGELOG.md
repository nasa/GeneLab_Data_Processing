# Changelog  

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.0](https://github.com/nasa/GeneLab_Data_Processing/blob/DEV_GeneLab_Reference_Annotations_vGL-DPPD-7110-A/GeneLab_Reference_Annotations/Workflow_Documentation/GL_RefAnnotTable-A)  

### Added  

- Added AnnotationForge helper script to create local annotations databases if not available on Bioconductor

- Added support for BACSU, ECOLI, and ORYLA via install-annot-dbi.R

### Fixed  

- Fixed automated processing for ECOLI

### Changed  

- Updated Ensembl versions
    - Animals: Ensembl release 112
    - Plants: Ensembl plants release 59
    - Bacteria: Ensembl bacteria release 59
- Removed org.EcK12.eg.db and replaced with local annotations database creation since it is no longer on Bioconductor


## [1.0.0](https://github.com/nasa/GeneLab_Data_Processing/releases/tag/GL_RefAnnotTable_1.0.0)
