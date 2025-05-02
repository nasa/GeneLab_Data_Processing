# Workflow change log

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [1.0.0](https://github.com/nasa/GeneLab_Data_Processing/tree/NF_MGIllumina_1.0.0/Metagenomics/Illumina/Workflow_Documentation/NF_MGIllumina)

### Changed
- Update to the latest pipeline version [GL-DPPD-7101-A](../../Pipeline_GL-DPPD-7107_Versions/GL-DPPD-7107-A.md) 
of the GeneLab Metagenomics consensus processing pipeline.
- Pipeline implementation as a Nextflow workflow [NF_MGIllumina](./) rather than Snakemake as in 
previous workflow versions.
- Run checkm separately on each bin and combine results to improve performance

### Fixed
- Allow explicit specification of the humann3 database location ([#62](https://github.com/nasa/GeneLab_Data_Processing/issues/62))
- Package bin and MAGs fasta files into per sample zip archives ([#76](https://github.com/nasa/GeneLab_Data_Processing/issues/76))

<BR>

---

> ***Note:** All previous workflow changes were associated with the previous version of the GeneLab Metagenomics Pipeline
[GL-DPPD-7101](../../Pipeline_GL-DPPD-7107_Versions/GL-DPPD-7107.md) and can be found in the
[change log of the Snakemake workflow (SW_MGIllumina)](../SW_MGIllumina/CHANGELOG.md).*