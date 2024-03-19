# SW_AmpIllumina-B Visualization Script Information and Usage Instructions<!-- omit in toc -->


## General  info <!-- omit in toc -->
The documentation for this script and its outputs can be found in [sections 6-10 of the GL-DPPD-7104-B.md processing protocol](/Amplicon/Illumina/Pipeline_GL-DPPD-7104_Versions/GL-DPPD-7104-B.md#6-amplicon-seq-data-analysis-set-up). This script is automatically executed as an optional step of the SW_AmpIllumina-B Snakemake workflow when the `run_workflow.py` argument `--visualizations TRUE` is used. Alternatively, the script can be executed manually as detailed below.

<br>

---

## Utilizing the script <!-- omit in toc -->


- [1. Set up the execution environment](#1-run-the-workflow-using-run_workflowpy)  
- [2. Run the visualization script manually](#2-run-the-visualization-script-manually)  
- [3. Parameter definitions](#3-parameter-definitions)

<br>

___

### 1. Set up the execution environment

The script should be executed from a Conda environment created using the [R_visualizations.yaml](/Amplicon/Illumina/Workflow_Documentation/SW_AmpIllumina-B/workflow_code/visualizations/R_visualizations.yaml) environment file.

<br>

___

### 2. Run the visualization script manually  

To run the script, the variables `runsheet_file`, `sample_info`, `counts`, `taxonomy`, `assay_suffix`, `plots_dir`, and `output_prefix` must be specified. The script can be manually executed via the command line by providing  positional arguments.

Additionally, the `RColorBrewer_Palette` variable can be modified in the script.  This variable determines the color palette from the RColorBrewer package that is applied to the plots.

```R
# Store command line args as variables #
args <- commandArgs(trailingOnly = TRUE)
runsheet_file <- paste0(args[1])
sample_info <- paste0(args[2])
counts <- paste0(args[3])
taxonomy <- paste0(args[4])
assay_suffix <- paste(args[5])
plots_dir <- paste0(args[6])
output_prefix <- paste0(args[7])
########################################

RColorBrewer_Palette <- "Set1"

```

Example run command: 
```bash
Rscript /path/to/visualizations/Illumina-R-visualizations.R "{runsheet_file}" "{sample_info}" "{counts}" "{taxonomy}" "{assay_suffix}" "{plots_dir}" "{output_prefix}"
```

<br>

___

### 3. Parameter definitions 

**Parameter Definitions for Illumina-R-visualizations.R:**
* `runsheet_file` – specifies the runsheet containing sample metadata required for processing (output from [GL-DPPD-7104-B step 6a](/Amplicon/Illumina/Pipeline_GL-DPPD-7104_Versions/GL-DPPD-7104-B.md#6a-create-sample-runsheet))
* `sample_info` – specifies the text file containing the IDs of each sample used, required for running the SW_AmpIllumina workflow (output from [run_workflow.py](/Amplicon/Illumina/Workflow_Documentation/SW_AmpIllumina-B/README.md#5-additional-output-files))
* `counts` – specifies the ASV counts table (output from [GL-DPPD-7104-B step 5g](/Amplicon/Illumina/Pipeline_GL-DPPD-7104_Versions/GL-DPPD-7104-B.md#5g-generating-and-writing-standard-outputs))
* `taxonomy` – specifies the taxonomy table (output from [GL-DPPD-7104-B step 5g](/Amplicon/Illumina/Pipeline_GL-DPPD-7104_Versions/GL-DPPD-7104-B.md#5g-generating-and-writing-standard-outputs))
* `plots_dir` – specifies the path where output files will be saved
* `assay_suffix` – specifies a string that is prepended to the start of the output file names. Default: ""
* `output_prefix` – specifies a string that is appended to the end of the output file names. Default: "_GLAmpSeq"
* `RColorBrewer_Palette` – specifies the RColorBrewer palette that will be used for coloring in the plots. Options include "Set1", "Accent", "Dark2", "Paired", "Pastel1", "Pastel2", "Set2", and "Set3". Default: "Set1"