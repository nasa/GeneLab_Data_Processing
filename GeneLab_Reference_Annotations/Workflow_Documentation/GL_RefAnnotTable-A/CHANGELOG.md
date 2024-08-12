# Changelog  

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.0](https://github.com/nasa/GeneLab_Data_Processing/blob/DEV_GeneLab_Reference_Annotations_vGL-DPPD-7110-A/GeneLab_Reference_Annotations/Workflow_Documentation/GL_RefAnnotTable-A)  

### Added  

- Added support for:
    - Bacillus subtilis, subsp. subtilis 168  
    - Brachypodium distachyon  
    - Escherichia coli,str. K-12 substr. MG1655  
    - Oryzias latipes  
    - Lactobacillus acidophilus NCFM  
    - Mycobacterium marinum M  
    - Oryza sativa Japonica  
    - Pseudomonas aeruginosa UCBPP-PA14  
    - Salmonella enterica subsp. enterica serovar Typhimurium str. LT2  
    - Serratia liquefaciens ATCC 27592  
    - Staphylococcus aureus MRSA252  
    - Streptococcus mutans UA159  
    - Vibrio fischeri ES114  
- Added AnnotationForge helper script install-org-db.R to create organism-specific annotation packages (org.*.eg.db) in R if not available on Bioconductor. Used for:  
    - Bacillus subtilis, subsp. subtilis 168  
    - Brachypodium distachyon  
    - Escherichia coli,str. K-12 substr. MG1655  
    - Oryzias latipes  
    - Salmonella enterica subsp. enterica serovar Typhimurium str. LT2  
- Added NCBI as a source for FASTA and GTF files  

### Fixed  

- Fixed processing for ECOLI

### Changed  

- Updated Ensembl versions
    - Animals: Ensembl release 112
    - Plants: Ensembl plants release 59
    - Bacteria: Ensembl bacteria release 59
- Removed org.EcK12.eg.db and replaced it with a locally created annotations database, as it is no longer available on Bioconductor


## [1.0.0](https://github.com/nasa/GeneLab_Data_Processing/releases/tag/GL_RefAnnotTable_1.0.0)
