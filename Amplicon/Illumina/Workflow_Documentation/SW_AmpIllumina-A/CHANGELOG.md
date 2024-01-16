# Workflow change log

## [1.1.1](https://github.com/nasa/GeneLab_Data_Processing/tree/SW_AmpIllumina-A_1.1.1/Amplicon/Illumina/Workflow_Documentation/SW_AmpIllumina-A)
- assay-specific suffixes added for certain files as needed by OSDR system
- resource parameters added to specific rules as needed for when cluster limits are being enforced
- `scripts/slurm-status.py` added for when using slurm such that slurm jobs are tracked and reported through to snakemake output
  - this is set by adding `--cluster-status scripts/slurm-status.py` to the snakemake call
- small fix for where some ITS datasets would fail biome creation 

## [1.1.0](https://github.com/nasa/GeneLab_Data_Processing/tree/SW_AmpIllumina-A_1.1.0/Amplicon/Illumina/Workflow_Documentation/SW_AmpIllumina-A)
- 18S option added
  - will use the PR2 database for taxonomy provided by the DECIPHER developers here: http://www2.decipher.codes/Downloads.html
  - option added to be able to just concatenate reads instead of merge them (in dada2), which is useful in scenarios like when primers such as 515-926 are used which can capture 18S sequences that are not expected to overlap

## [1.0.1](https://github.com/nasa/GeneLab_Data_Processing/tree/SW_AmpIllumina-A_1.0.1/Amplicon/Illumina/Workflow_Documentation/SW_AmpIllumina-A)
- documentation links updated in config.yaml to match changes in host site
- added config file for multiqc to trim suffixes from sample names and not include paths in output report
- packaging multiqc html in with data dir in a zip to match what RNAseq workflow does
- restructured fastqc and multiqc rules using inheritance to reduce some redundancy
- added workflow version number to the top of the Snakefile
- fixed an issue with the taxonomy-and-counts.biom.zip file produced
  - while it could previously be unzipped at the command line, on some systems trying to unzip it through a browser was failing
  - it looks like this was due to there being relative paths involved during the initial zip
  - added the `-j` flag to drop them since there is no directory structure needed for this single file

## [1.0.0](https://github.com/nasa/GeneLab_Data_Processing/tree/SW_AmpIllumina-A_1.0.0/Amplicon/Illumina/Workflow_Documentation/SW_AmpIllumina-A)
- original workflow version
