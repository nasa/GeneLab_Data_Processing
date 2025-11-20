# Workflow change log

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [1.0.0](https://github.com/nasa/GeneLab_Data_Processing/tree/master/Metagenomics/Remove_host_reads/Workflow_Documentation/NF_MGRemoveHostReads)

### Changed
- Expand to support removal of host reads beyond human samples, forming the basis of the current MGRemoveHostReads workflow
- Update to the latest pipeline version [GL-DPPD-7105-B](../../Pipeline_GL-DPPD-7105_Versions/GL-DPPD-7105-B.md) 
of the GeneLab Remove-Host-Reads consensus processing pipeline.
- Pipeline implementation as a Nextflow workflow [NF_MGRemoveHostReads](./) rather than Snakemake as in 
previous workflow versions.

### Added
- Build kraken2 database from scratch using host organism's information pulled from [hosts.csv](workflow_code/assets/hosts.csv)
- Create protocol.txt as an output file describing workflow methods

### Removed
- kraken2-human-db/ no longer automatically downloaded to run with the workflow. It can now be explicitly set or built from scratch in case it doesn't exist.

<BR>

---

> ***Note:** Change log of the Snakemake workflow (SW_MGRemoveHumanReads) that is associated with the previous version of the GeneLab Remove-Host-Reads Pipeline [GL-DPPD-7105](../../Pipeline_GL-DPPD-7105_Versions/GL-DPPD-7105-A.md) can be found [here](../SW_MGRemoveHumanReads/CHANGELOG.md)*