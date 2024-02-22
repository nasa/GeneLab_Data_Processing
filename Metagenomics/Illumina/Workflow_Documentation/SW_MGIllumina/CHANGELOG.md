# Workflow change log

## [2.0.4](https://github.com/nasa/GeneLab_Data_Processing/tree/SW_MGIllumina_2.0.4/Metagenomics/Illumina/Workflow_Documentation/SW_MGIllumina)
- assay-specific suffixes added for certain files as needed by OSDR system
- resource parameters added to specific rules as needed for when cluster limits are being enforced (only relevent when using a job scheduler)
- `scripts/slurm-status.py` added for when using slurm such that slurm jobs are tracked and reported through to snakemake output
  - this is set by adding `--cluster-status scripts/slurm-status.py` to the snakemake call
- fix for metaphlan auto-db setup by workflow
  - the metaphlan download command was pulling the wrong reference database, so the correct db for the utilized version had to be hard-coded (discussed here: https://forum.biobakery.org/t/metaphlan-v4-0-2-and-huma-3-6-metaphlan-taxonomic-profile-provided-was-not-generated-with-the-expected-database/4296/29)

## [2.0.3](https://github.com/nasa/GeneLab_Data_Processing/tree/SW_MGIllumina_2.0.3/Metagenomics/Illumina/Workflow_Documentation/SW_MGIllumina)
- added config file for multiqc to trim suffixes from sample names and not include paths in output report
- packaging multiqc html in with data dir in a zip to match what RNAseq workflow does
- restructured fastqc and multiqc rules using inheritance to reduce some redundancy
- added workflow version number to top of primary Snakefile

## [2.0.2](https://github.com/nasa/GeneLab_Data_Processing/tree/SW_MGIllumina_2.0.2/Metagenomics/Illumina/Workflow_Documentation/SW_MGIllumina)
- updating bit package from 1.8.47 to 1.8.53
    - some users were having trouble with conda installing an appropriate goatools (dependency of bit) on some systems, bumping to a later version of bit seemed to resolve the issue
        - goatools isn't used by this workflow anyway, and none of the code used here in bit changed


## [2.0.1](https://github.com/nasa/GeneLab_Data_Processing/tree/SW_MGIllumina_2.0.1/Metagenomics/Illumina/Workflow_Documentation/SW_MGIllumina)
- updated humann version from 3.5 to 3.6
    - in accordance with critical bug described here: 
[https://forum.biobakery.org/t/announcing-humann-3-6-critical-update/4155](https://forum.biobakery.org/t/announcing-humann-3-6-critical-update/4155)


## [2.0.0](https://github.com/nasa/GeneLab_Data_Processing/tree/SW_MGIllumina_2.0.0/Metagenomics/Illumina/Workflow_Documentation/SW_MGIllumina)
- humann version set to 3.5 with metaphlan 4.0.1 pinned to it
    - the humann conda installation was not pinned to a specific version of metaphlan, and new releases of metaphlan were incompatible
    - set as a major version update because the humann/metaphlan option for `--unknown_estimation` was changed to `--unclassified_estimation`, requring changes making it incompatible with the previous workflow


## [1.0.0](https://github.com/nasa/GeneLab_Data_Processing/tree/SW_MGIllumina_1.0.0/Metagenomics/Illumina/Workflow_Documentation/SW_MGIllumina)
- original workflow version
