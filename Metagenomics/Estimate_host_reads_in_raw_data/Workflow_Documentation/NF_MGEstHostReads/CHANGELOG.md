# Workflow change log

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [1.0.0](https://github.com/nasa/GeneLab_Data_Processing/tree/master/Metagenomics/Estimate_host_reads_in_raw_data/Workflow_Documentation/NF_MGEstHostReads)

### Changed
- Update to the latest pipeline version [GL-DPPD-7109-A](../../Pipeline_GL-DPPD-7107_Versions/GL-DPPD-7109-A.md) 
of the GeneLab Estimate-Host-Reads consensus processing pipeline.
- Pipeline implementation as a Nextflow workflow [NF_MGEstHostReads](./) rather than Snakemake as in 
previous workflow versions.

### Added
- Pull dataset from OSDR option using dp_tools
- Build kraken2 database from scratch using host organism's information pulled from [hosts.csv](workflow_code/assets/hosts.csv)
- Create protocol.txt as an output file describing workflow methods

### Removed
- kraken2-mouse-db/ no longer needed as part of the workflow files (can now be explicitly set or built from scratch in case it doesn't exist)

<BR>

---

> ***Note:** Change log of the Snakemake workflow (SW_MGEstHostReads) that is associated with the previous version of the GeneLab Estimate-Host-Reads Pipeline [GL-DPPD-7109](../../Pipeline_GL-DPPD-7107_Versions/GL-DPPD-7107.md) can be found [here](../SW_MGEstHostReads/CHANGELOG.md)*