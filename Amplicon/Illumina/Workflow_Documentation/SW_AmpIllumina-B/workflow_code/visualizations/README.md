# SW_AmpIllumina-B Visualization Script Information and Usage Instructions<!-- omit in toc -->


## General info <!-- omit in toc -->
The documentation for this script and its outputs can be found in steps 6-10 of the [GL-DPPD-7104-B.md](https://github.com/nasa/GeneLab_Data_Processing/blob/master/Amplicon/Illumina/Pipeline_GL-DPPD-7104_Versions/GL-DPPD-7104-B.md) pipeline document. This script is automatically executed as an optional step of the [SW_AmpIllumina-B](../../) Snakemake workflow when the `run_workflow.py` argument, `--visualizations TRUE`, is used. Alternatively, the script can be executed independently as detailed below.

<br>

---

## Utilizing the script <!-- omit in toc -->


- [1. Set up the execution environment](#1-set-up-the-execution-environment)  
- [2. Run the visualization script manually](#2-run-the-visualization-script-manually)  

<br>

___

### 1. Set up the execution environment  

The script should be executed from a [conda](https://docs.conda.io/en/latest/) environment created using the [R_visualizations.yaml](R_visualizations.yaml) environment file.
> If you do not have conda installed, an introduction to conda with installation help and links to other resources can be found [here at Happy Belly Bioinformatics](https://astrobiomike.github.io/unix/conda-intro).

Download the [R_visualizations.yaml](R_visualizations.yaml) environment file and the [Illumina-R-visualizations.R](Illumina-R-visualizations.R) script by running the following commands:

```
curl -LO https://raw.githubusercontent.com/nasa/GeneLab_Data_Processing/dev2-amplicon-add-runsheet-visualizations/Amplicon/Illumina/Workflow_Documentation/SW_AmpIllumina-B/workflow_code/visualizations/R_visualizations.yaml

curl -LO https://raw.githubusercontent.com/nasa/GeneLab_Data_Processing/dev2-amplicon-add-runsheet-visualizations/Amplicon/Illumina/Workflow_Documentation/SW_AmpIllumina-B/workflow_code/visualizations/Illumina-R-visualizations.R
```

Next, create the AmpSeqVisualizations environment by running the following command:

```
conda env create -f R_visualizations.yaml -n AmpSeqVisualizations
```

Then activate the environment as follows:

```
conda activate AmpSeqVisualizations
```


<br>

___

### 2. Run the visualization script manually  

The [Illumina-R-visualizations.R](./Illumina-R-visualizations.R) script can be executed from the command line by providing `runsheet_file`, `sample_info`, `counts`, `taxonomy`, `plots_dir`, `output_prefix`, and `assay_suffix` as positional arguments, in their respecive order.

The example command below shows how to execute the script with the following parameters:
 * runsheet_file: /path/to/runsheet.csv  
 * sample_info: /path/to/unique-sample-IDs.txt
 * counts: /path/to/counts_GLAmpSeq.tsv
 * taxonomy: /path/to/taxonomy_GLAmpSeq.tsv
 * plots_dir: /path/to/Plots/
 * output_prefix: my_prefix_
 * assay_suffix: _GL_Ampseq

```bash
Rscript /path/to/visualizations/Illumina-R-visualizations.R /path/to/runsheet.csv /path/to/unique-sample-IDs.txt /path/to/counts_GLAmpSeq.tsv /path/to/taxonomy_GLAmpSeq.tsv /path/to/Plots/ "my_prefix_" "_GL_Ampseq"
```

Additionally, the `RColorBrewer_Palette` variable can be modified in the script. This variable determines the color palette from the RColorBrewer package that is applied to the plots.

**Parameter Definitions:**
* `runsheet_file` – specifies the table containing sample metadata required for processing 
* `sample_info` – specifies the text file containing the IDs of each sample used, required for running the SW_AmpIllumina workflow 
* `counts` – specifies the ASV counts table 
* `taxonomy` – specifies the taxonomy table 
* `plots_dir` – specifies the path where output files will be saved
* `output_prefix` – specifies a string that is prepended to the start of the output file names. Default: ""
* `assay_suffix` – specifies a string that is appended to the end of the output file names. Default: "_GLAmpSeq"
* `RColorBrewer_Palette` – specifies the RColorBrewer palette that will be used for coloring in the plots. Options include "Set1", "Accent", "Dark2", "Paired", "Pastel1", "Pastel2", "Set2", and "Set3". Default: "Set1"

**Input Data:**
* *runsheet.csv (output from [GL-DPPD-7104-B step 6a](https://github.com/nasa/GeneLab_Data_Processing/blob/master/Amplicon/Illumina/Pipeline_GL-DPPD-7104_Versions/GL-DPPD-7104-B.md#6a-create-sample-runsheet))
* unique-sample-IDs.txt (output from [run_workflow.py](/Amplicon/Illumina/Workflow_Documentation/SW_AmpIllumina-B/README.md#5-additional-output-files))
* counts_GLAmpSeq.tsv (output from [GL-DPPD-7104-B step 5g](https://github.com/nasa/GeneLab_Data_Processing/blob/master/Amplicon/Illumina/Pipeline_GL-DPPD-7104_Versions/GL-DPPD-7104-B.md#5g-generating-and-writing-standard-outputs))
* taxonomy_GLAmpSeq.tsv (output from [GL-DPPD-7104-B step 5g](https://github.com/nasa/GeneLab_Data_Processing/blob/master/Amplicon/Illumina/Pipeline_GL-DPPD-7104_Versions/GL-DPPD-7104-B.md#5g-generating-and-writing-standard-outputs))

**Output Data:**
* **{output_prefix}dendrogram_by_group{assay_suffix}.png** (dendrogram of euclidean distance - based hierarchical clustering of the samples, colored by experimental groups)
* **{output_prefix}rarefaction_curves{assay_suffix}.png** (Rarefaction curves plot for all samples)
* **{output_prefix}richness_and_diversity_estimates_by_sample{assay_suffix}.png** (Richness and diversity estimates plot for all samples)
* **{output_prefix}richness_and_diversity_estimates_by_group{assay_suffix}.png** (Richness and diversity estimates plot for all groups)
* **{output_prefix}relative_phyla{assay_suffix}.png** (taxonomic summaries plot based on phyla, for all samples)
* **{output_prefix}relative_classes{assay_suffix}.png** (taxonomic summaries plot based on class, for all samples)
* **{output_prefix}samplewise_phyla{assay_suffix}.png** (taxonomic summaries plot based on phyla, for all samples)
* **{output_prefix}samplewise_classes{assay_suffix}.png** (taxonomic summaries plot based on class, for all samples)
* **{output_prefix}PCoA_w_labels{assay_suffix}.png** (principle Coordinates Analysis plot of VST transformed ASV counts, with sample labels)
* **{output_prefix}PCoA_without_labels{assay_suffix}.png** (principle Coordinates Analysis plot of VST transformed ASV counts, without sample labels)
* **{output_prefix}normalized_counts{assay_suffix}.tsv** (size factor normalized ASV counts table)
* **{output_prefix}group1_vs_group2.csv** (differential abundance tables for all pairwise contrasts of groups)
* **{output_prefix}volcano_group1_vs_group2.png** (volcano plots for all pairwise contrasts of groups)
* {output_prefix}color_legend{assay_suffix}.png (color legend for all groups)