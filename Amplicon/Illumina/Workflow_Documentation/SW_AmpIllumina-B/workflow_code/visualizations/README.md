# SW_AmpIllumina-B Visualization Script Information and Usage Instructions<!-- omit in toc -->


## General  info <!-- omit in toc -->
The documentation for this script and its outputs can be found in [sections 6-10 of the GL-DPPD-7104-B.md processing protocol](/Amplicon/Illumina/Pipeline_GL-DPPD-7104_Versions/GL-DPPD-7104-B.md#6-amplicon-seq-data-analysis-set-up). This script is automatically executed as an optional step of the SW_AmpIllumina-B Snakemake workflow when the `run_workflow.py` argument `--visualizations TRUE` is used. Alternatively, the script can be executed manually as detailed below.

<br>

---

## Utilizing the script <!-- omit in toc -->


- [1. Setting up the execution environment](#1-run-the-workflow-using-run_workflowpy)  
- [2. Run the visualization script manually](#2-run-the-visualization-script-manually)  
- [3. Parameter definitions](#3-parameter-definitions)

<br>

___

### 1. Setting up the execution environment

To ensure that manual execution outputs that are consistent with those of the Snakemake workflow, the script should be executed from a Conda environment created using the [R_visualizations.yaml](/Amplicon/Illumina/Workflow_Documentation/SW_AmpIllumina-B/workflow_code/envs/R_visualizations.yaml/) environment file.

### 2. Run the visualization script manually  

To manually execute the script, the following variables: `runsheet_file`, `sample_info`, `counts`, `taxonomy`, `final_outputs_dir`, `assay_suffix`, and `output_prefix`, must be explicitly assigned with appropriate values. For manual execution, the lines within the Snakemake Configuration section should be commented out or removed. These lines are used for capturing variables during execution as part of the Snakemake workflow.

Additionally, the `RColorBrewer_Palette` variable can be modified.  This variable determines the color palette from the RColorBrewer package that is applied to the plots.

```R
####### Manual Configuration ########
# To run the script manually:
# 1. Uncomment and populate lines 23-30
# 2. Comment out the Snakemake configuration section (lines 35-42).

# runsheet_file <- "/path/to/SW_AmpIllumina-B/runsheet.csv"
# sample_info <- "/path/to/SW_AmpIllumina-B/unique-sample-IDs.txt"
# counts <- "/path/to/SW_AmpIllumina-B/workflow_output/Final_Outputs/counts_GLAmpSeq.tsv"
# taxonomy <- "/path/to/SW_AmpIllumina-B/workflow_output/Final_Outputs/taxonomy_GLAmpSeq.tsv"
# final_outputs_dir <- "/path/to/SW_AmpIllumina-B/workflow_output/Final_Outputs/" # Where visualization script outputs will be saved to

# assay_suffix <- "_GLAmpSeq" # Suffix appended to the end of output file names (Default: "_GLAmpSeq")
# output_prefix <- "" # Prefix appended to the start of output file names (Default: "")
#####################################

###### Snakemake Configuration ######
# Assigns arguments to variables based on the snakemake run command 
args <- commandArgs(trailingOnly = TRUE)
runsheet_file <- paste0(args[1])
sample_info <- paste0(args[2])
counts <- paste0(args[3])
taxonomy <- paste0(args[4])
assay_suffix <- paste(args[5])
final_outputs_dir <- paste0(args[6])
output_prefix <- paste0(args[7])
#####################################

RColorBrewer_Palette <- "Set1"

```
<br>

### 3. Parameter definitions 

**Parameter Definitions for Illumina-R-visualizations.R:**
* `runsheet_file` – Specifies the runsheet containing sample metadata required for processing (output from [GL-DPPD-7104-B step 6a](/Amplicon/Illumina/Pipeline_GL-DPPD-7104_Versions/GL-DPPD-7104-B.md#6a-create-sample-runsheet))
* `sample_info` – Specifies the text file containing the IDs of each sample used, required for running the SW_AmpIllumina workflow (output from [run_workflow.py](/Amplicon/Illumina/Workflow_Documentation/SW_AmpIllumina-B/README.md#5-additional-output-files))
* `counts` – Specifies the ASV counts table (output from [GL-DPPD-7104-B step 5g](/Amplicon/Illumina/Pipeline_GL-DPPD-7104_Versions/GL-DPPD-7104-B.md#5g-generating-and-writing-standard-outputs))
* `taxonomy` – Specifies the taxonomy table (output from [GL-DPPD-7104-B step 5g](/Amplicon/Illumina/Pipeline_GL-DPPD-7104_Versions/GL-DPPD-7104-B.md#5g-generating-and-writing-standard-outputs))
* `final_outputs_dir` – Specifies the path where output files will be saved.
* `assay_suffix` – Specifies a string that is prepended to the start of the output file names. Default: ""
* `output_prefix` – Specifies a string that is appended to the end of the output file names Default: "_GLAmpSeq"
* `RColorBrewer_Palette` – Specifies the RColorBrewer palette that will be used for coloring in the plots. Options include "Set1", "Accent", "Dark2", "Paired", "Pastel1", "Pastel2", "Set2", and "Set3". Default: "Set1"