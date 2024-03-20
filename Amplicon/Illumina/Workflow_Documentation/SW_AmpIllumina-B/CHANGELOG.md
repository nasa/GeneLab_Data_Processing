# Workflow change log

## [1.2.2](https://github.com/nasa/GeneLab_Data_Processing/tree/SW_AmpIllumina-B_1.2.2/Amplicon/Illumina/Workflow_Documentation/SW_AmpIllumina-B)
- Visualizations are now optional with the default being off.
  - Enable with optional `run_workflow.py` argument `--visualizations TRUE` or setting `config.yaml` `enable_visualizations` to "TRUE"
- Moved visualization script and conda environment config file to `workflow_code/visualizations/`
- Added `workflow_code/visualizations/README.md` for instructions on running the visualization script manually
- Refactored Snakefile outputs

## [1.2.1](https://github.com/nasa/GeneLab_Data_Processing/tree/SW_AmpIllumina-B_1.2.1/Amplicon/Illumina/Workflow_Documentation/SW_AmpIllumina-B)
- Moved SW_AmpIllumina-A_1.2.1 to SW_AmpIllumina-B_1.2.1
- Workflow runs the [GL-DPPD-7104-B version](../../Pipeline_GL-DPPD-7104_Versions/GL-DPPD-7104-B.md) of the GeneLab standard pipeline, which includes data visualization outputs
- Removed wget SW_AmpIllumina-A_1.2.0.zip instructions and added GL-get-workflow instructions to SW_AmpIllumina usage instructions
- Removed dp-tools installation from SW_AmpIllumina usage instructions since it is now included in the genelab-utils installation
- Added a list of (edge-case) datasets that cannot be autoprocessed using run_workflow.py
- Set default for --anchor-primers to FALSE and set the default for --primers-linked to FALSE
- Visualization script changes: added samplewise taxonomy plots, renamed and moved plots related to alpha or beta diversity
- Fix where some ITS datasets would fail to create biom object
- assay-specific suffixes added for certain files as needed by OSDR system
- resource parameters added to specific rules as needed for when cluster limits are being enforced
- `scripts/slurm-status.py` added for when using slurm such that slurm jobs are tracked and reported through to snakemake output
  - this is set by adding `--cluster-status scripts/slurm-status.py` to the snakemake call 
 
## [1.2.0](https://github.com/nasa/GeneLab_Data_Processing/tree/SW_AmpIllumina-A_1.2.0/Amplicon/Illumina/Workflow_Documentation/SW_AmpIllumina-A)
- Added runsheet dependency, runsheet definition in [config.yaml](workflow_code/config.yaml)
- Added [run_workflow.py](workflow_code/scripts/run_workflow.py)
  - Sets up runsheet for OSDR datasets. Uses a runsheet to set up [config.yaml](workflow_code/config.yaml) and [unique-sample-IDs.txt](workflow_code/unique-sample-IDs.txt). Runs the Snakemake workflow.
- Updated instructions in [README.md](README.md) to use [run_workflow.py](workflow_code/scripts/run_workflow.py)
- Added downstream analysis visualizations rule using [Illumina-R-Visualizations.R](workflow_code/scripts/Illumina-R-visualizations.R)
  - Volcano plots, dendrogram, PCoA, rarefaction, richness, taxonomy plots

<br> 

---

<br> 

All previous workflow changes were associated with [version A of the GeneLab Amplicon Seq Illumina Pipeline](../../Pipeline_GL-DPPD-7104_Versions/GL-DPPD-7104-A.md), and can be found in the [change log of the SW_AmpIllumina-A workflow](../SW_AmpIllumina-A/CHANGELOG.md).