# Workflow change log

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
