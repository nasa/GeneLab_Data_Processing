# GeneLab Workflow Information for Removing Host Reads in Metagenomics Seq Data

> **GeneLab has wrapped each step of the removing human reads in metagenomics sequencing data pipeline (MGRemoveHumanReads), starting with pipeline version A, into a workflow. This workflow has since been expanded to support removal of host reads beyond human samples, forming the basis of the current MGRemoveHostReads workflow. The table below lists (and links to) the previous MGRemoveHumanReads version as well as the expanded MGRemoveHostReads versions and the corresponding workflow subdirectories, with the current implementation indicated. The workflow subdirectory contains information about the workflow along with instructions for installation and usage. Exact workflow run info and the MGRemoveHumanReads or MGRemoveHostReads version used to process specific datasets that have been released are provided with their processed data in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/).**  

## MGRemoveHumanReads Pipeline Version and Corresponding Workflow

|Pipeline Version|Current Workflow Version (for respective pipeline version)|Nextflow Version|
|:---------------|:---------------------------------------------------------|:---------------|
|*[GL-DPPD-7105-B.md](../Pipeline_GL-DPPD-7105_Versions/GL-DPPD-7105-B.md)|[NF_MGRemoveHostReads_1.0.0](NF_MGRemoveHostReads)|25.04.6|
|[GL-DPPD-7105-A.md](../Pipeline_GL-DPPD-7105_Versions/GL-DPPD-7105-A.md)|[SW_MGRemoveHumanReads_1.0.0](SW_MGRemoveHumanReads-A)|N/A (Snakemake v7.26.0)|

*Current GeneLab Pipeline/Workflow Implementation

> See the [workflow change log](NF_MGRemoveHostReads/CHANGELOG.md) to access previous workflow versions and view all changes associated with each version update.
