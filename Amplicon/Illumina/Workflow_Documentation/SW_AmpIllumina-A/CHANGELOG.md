# Workflow change log

## [1.2.0](https://github.com/nasa/GeneLab_Data_Processing/tree/SW_AmpIllumina-A_1.2.0/Amplicon/Illumina/Workflow_Documentation/SW_AmpIllumina-A)
- Added runsheet dependency, runsheet definition in [config.yaml](workflow_code/config.yaml)
- Added [run_workflow.py](workflow_code/scripts/run_workflow.py)
  - Sets up runsheet for OSDR datasets. Uses a runsheet to set up [config.yaml](workflow_code/config.yaml) and [unique-sample-IDs.txt](workflow_code/unique-sample-IDs.txt). Runs the Snakemake workflow.
- Updated instructions in [README.md](README.md) to use [run_workflow.py](workflow_code/scripts/run_workflow.py)
- Added downstream analysis visualizations rule using [Illumina-R-Visualizations.R](workflow_code/scripts/Illumina-R-visualizations.R)
  - Volcano plots, dendrogram, PCoA, rarefaction, richness, taxonomy plots

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
