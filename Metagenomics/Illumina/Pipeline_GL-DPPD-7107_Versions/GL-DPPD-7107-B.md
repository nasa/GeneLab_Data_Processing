# Bioinformatics pipeline for Illumina metagenomics data  <!-- omit in toc -->

> **This document holds an overview and some example commands of how GeneLab processes Illumina metagenomics datasets. Exact processing commands for specific datasets that have been released are provided with their processed data in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/).**  

---

**Date:** April 3, 2026  
**Revision:** B  
**Document Number:** GL-DPPD-7107  

**Submitted by:**  
Olabiyi A. Obayomi (GeneLab Data Processing Team)  

**Approved by:**  
Jonathan Galazka (OSDR Project Manager)  
Danielle Lopez (OSDR Deputy Project Manager)  
Amanda Saravia-Butler (OSDR Subject Matter Expert)  
Barbara Novak (GeneLab Data Processing Lead)  


---

## Updates from previous version  <!-- omit in toc -->

Software Updates and Changes:

| Program      | Previous Version | New Version |
| :----------- | :--------------: | :---------: |
| MultiQC      |       1.19       |   1.27.1    |
| samtools     |       1.20       |   1.22.1    |
| Kaiju        |       N/A        |   1.10.1    |
| fastp        |       N/A        |    1.3.1    |
| Kaiju        |       N/A        |   1.10.1    |
| Kraken2      |       N/A        |    2.1.6    |
| KrakenTools  |       N/A        |     1.2     |
| Krona        |       N/A        |    2.8.1    |
| SPAdes       |       N/A        |    4.1.0    |
| R            |       N/A        |    4.5.3    |
| htmlwidgets  |       N/A        |    1.6.4    |
| pavian       |       N/A        |    1.2.0    |
| pheatmap     |       N/A        |   1.0.13    |
| phyloseq     |       N/A        |   1.54.0    |
| plotly       |       N/A        |   4.12.0    |
| dplyr        |       N/A        |    1.2.0    |
| ggplot2      |       N/A        |    4.0.2    |
| glue         |       N/A        |    1.8.0    |
| magrittr     |       N/A        |    2.0.5    |
| purrr        |       N/A        |    1.2.1    |
| readr        |       N/A        |    2.2.0    |
| scales       |       N/A        |    1.4.0    |
| stringr      |       N/A        |    1.6.0    |
| tibble       |       N/A        |    3.3.1    |
| tidyr        |       N/A        |    1.3.2    |
| htmlwidgets  |       N/A        |    1.6.4    |

- Synced this pipeline with the new low-biomass pipelines (update formatting and definitions)
- Added new processing steps for additional taxonomic profiling tools and downstream processed data outputs in R
  - Add additional read-based processing taxonomic profiling methods:
    - Kaiju taxonomic profiling ([Step 18](#18-taxonomic-profiling-using-kaiju))
    - Kraken2 taxonomic profiling ([Step 19](#19-taxonomic-profiling-using-kraken2))
  - summary plots for all taxonomic profiling and functional profiling for both read-based and assembly-based processing
    - barplots for read-based taxonomic profiling
    - heatmaps for read-based functional profiling
    - heatmaps for assembly-based taxonomy and functional profiling
  - Filtering (for rare taxa/features)
    - Read-based processing 
      - Kaiju and Kraken2 taxonomies (see [Step 18g](#18g-filter-kaiju-species-count-table) and [Step 19f](#19f-filter-kraken2-species-count-table))
      - HUMAnN/MetaPhlan taxonomies and functional profiling (see ([Step 20i](#20i-filter-metaphlan-species-count-table) and [Step 20k](#20k-filter-humann-output)))
    - Assembly-based processing (see [Step 16](#16-filtering-and-visualization-of-contig--and-gene-taxonomy-and-gene-function-outputs))
  - Added missing steps for generating assembly-based processing overview and failed assembly file
- replace bbduk with fastp for read quality filtering and adapter trimming

---

# Table of contents  <!-- omit in toc -->

- [Software used](#software-used)
- [General processing overview with example commands](#general-processing-overview-with-example-commands)
  - [Pre-processing](#pre-processing)
    - [1. Raw Data QC](#1-raw-data-qc)
      - [1a. Raw Data QC](#1a-raw-data-qc)
      - [1b. Compile Raw Data QC](#1b-compile-raw-data-qc)
    - [2. Quality filtering/trimming](#2-quality-filteringtrimming)
      - [2a. Filter Quality and Trim Adapters](#2a-filter-quality-and-trim-adapters)
      - [2b. Trim polyG](#2b-trim-polyg)
      - [2c. Filtered/Trimmed Data QC](#2c-filteredtrimmed-data-qc)
      - [2d. Compile Filtered/Trimmed Data QC](#2d-compile-filteredtrimmed-data-qc)
    - [3. R Environment Setup](#3-r-environment-setup)
      - [3a. Load libraries](#3a-load-libraries)
      - [3b. Define Custom Functions](#3b-define-custom-functions)
      - [3c. Set global variables](#3c-set-global-variables)
  - [Assembly-based Processing](#assembly-based-processing)
    - [4. Sample assembly](#4-sample-assembly)
    - [5. Rename Contigs and Summarize Assemblies](#5-rename-contigs-and-summarize-assemblies)
      - [5a. Rename Contig Headers](#5a-rename-contig-headers)
      - [5b. Summarize assemblies](#5b-summarize-assemblies)
    - [6. Gene prediction](#6-gene-prediction)
      - [6a. Generate Gene Predictions](#6a-generate-gene-predictions)
      - [6b. Remove Line Wraps In Gene Prediction Output](#6b-remove-line-wraps-in-gene-prediction-output)
    - [7. Functional annotation](#7-functional-annotation)
      - [7a. Download reference database of HMM models](#7a-download-reference-database-of-hmm-models)
      - [7b. Run KEGG annotation](#7b-run-kegg-annotation)
      - [7c. Filter KO Outputs](#7c-filter-ko-outputs)
    - [8. Taxonomic classification](#8-taxonomic-classification)
      - [8a. Pull and Unpack Pre-built Reference DB](#8a-pull-and-unpack-pre-built-reference-db)
      - [8b. Run Taxonomic Classification](#8b-run-taxonomic-classification)
      - [8c. Add taxonomy info from taxids to genes](#8c-add-taxonomy-info-from-taxids-to-genes)
      - [8d. Add Taxonomy Info From Taxids To Contigs](#8d-add-taxonomy-info-from-taxids-to-contigs)
      - [8e. Format Gene-level Output With awk and sed](#8e-format-gene-level-output-with-awk-and-sed)
      - [8f. Format Contig-level Output With awk and sed](#8f-format-contig-level-output-with-awk-and-sed)
    - [9. Read-Mapping](#9-read-mapping)
      - [9a. Build reference index](#9a-build-reference-index)
      - [9b. Align Reads to Sample Assembly](#9b-align-reads-to-sample-assembly)
      - [9c. Sort Assembly Alignments](#9c-sort-assembly-alignments)
    - [10. Get Coverage Information and Filter Based On Detection](#10-get-coverage-information-and-filter-based-on-detection)
      - [10a. Filter Coverage Levels Based On Detection](#10a-filter-coverage-levels-based-on-detection)
      - [10b. Filter Gene and Contig Coverage Based On Detection](#10b-filter-gene-and-contig-coverage-based-on-detection)
    - [11. Combine Gene-level Coverage, Taxonomy, and Functional Annotations For Each Sample](#11-combine-gene-level-coverage-taxonomy-and-functional-annotations-for-each-sample)
    - [12. Combine Contig-level Coverage and Taxonomy For Each Sample](#12-combine-contig-level-coverage-and-taxonomy-for-each-sample)
    - [13. Generating normalized, gene- and contig-level coverage summary tables of KO-annotations and taxonomy across samples](#13-generating-normalized-gene--and-contig-level-coverage-summary-tables-of-ko-annotations-and-taxonomy-across-samples)
      - [13a. Generate Gene-level Coverage Summary Tables](#13a-generate-gene-level-coverage-summary-tables)
      - [13b. Generate Contig-level Coverage Summary Tables](#13b-generate-contig-level-coverage-summary-tables)
    - [14. **M**etagenome-**A**ssembled **G**enome (MAG) recovery](#14-metagenome-assembled-genome-mag-recovery)
      - [14a. Bin contigs](#14a-bin-contigs)
      - [14b. Bin quality assessment](#14b-bin-quality-assessment)
      - [14c. Filter MAGs](#14c-filter-mags)
      - [14d. MAG taxonomic classification](#14d-mag-taxonomic-classification)
      - [14e. Generate Overview Table Of All MAGs](#14e-generate-overview-table-of-all-mags)
    - [15. Generate MAG-level functional summary overview](#15-generate-mag-level-functional-summary-overview)
      - [15a. Get KO annotations per MAG](#15a-get-ko-annotations-per-mag)
      - [15b. Summarize KO annotations with KEGG-Decoder](#15b-summarize-ko-annotations-with-kegg-decoder)
    - [16. Filtering and Visualization of Contig- and Gene-taxonomy and Gene-function Outputs](#16-filtering-and-visualization-of-contig--and-gene-taxonomy-and-gene-function-outputs)
      - [16a. Gene-level Taxonomy Heatmaps](#16a-gene-level-taxonomy-heatmaps)
      - [16b. Gene-level Taxonomy Feature Filtering](#16b-gene-level-taxonomy-feature-filtering)
      - [16c. Gene-level KO Functions Heatmaps](#16c-gene-level-ko-functions-heatmaps)
      - [16d. Gene-level KO Functions Feature Filtering](#16d-gene-level-ko-functions-feature-filtering)
      - [16e. Contig-level Heatmaps](#16e-contig-level-heatmaps)
      - [16f. Contig-level Feature Filtering](#16f-contig-level-feature-filtering)
    - [17. Generate Assembly-based Processing Overview](#17-generate-assembly-based-processing-overview)
  - [Read-based Processing](#read-based-processing)
    - [18. Taxonomic Profiling Using Kaiju](#18-taxonomic-profiling-using-kaiju)
      - [18a. Build Kaiju Database](#18a-build-kaiju-database)
      - [18b. Kaiju Taxonomic Classification](#18b-kaiju-taxonomic-classification)
      - [18c. Compile Kaiju Taxonomy Results](#18c-compile-kaiju-taxonomy-results)
      - [18d. Convert Kaiju Output To Krona Format](#18d-convert-kaiju-output-to-krona-format)
      - [18e. Compile Kaiju Krona Reports](#18e-compile-kaiju-krona-reports)
      - [18f. Create Kaiju Species Count Table](#18f-create-kaiju-species-count-table)
      - [18g. Filter Kaiju Species Count Table](#18g-filter-kaiju-species-count-table)
      - [18h. Kaiju Taxonomy Barplots](#18h-kaiju-taxonomy-barplots)
    - [19. Taxonomic Profiling Using Kraken2](#19-taxonomic-profiling-using-kraken2)
      - [19a. Download Kraken2 Database](#19a-download-kraken2-database)
      - [19b. Kraken2 Taxonomic Classification](#19b-kraken2-taxonomic-classification)
      - [19c. Compile Kraken2 Taxonomy Results](#19c-compile-kraken2-taxonomy-results)
        - [19ci. Create Merged Kraken2 Taxonomy Table](#19ci-create-merged-kraken2-taxonomy-table)
        - [19cii. Compile Kraken2 Taxonomy Reports](#19cii-compile-kraken2-taxonomy-reports)
      - [19d. Convert Kraken2 Output to Krona Format](#19d-convert-kraken2-output-to-krona-format)
      - [19e. Compile Kraken2 Krona Reports](#19e-compile-kraken2-krona-reports)
      - [19f. Filter Kraken2 Species Count Table](#19f-filter-kraken2-species-count-table)
      - [19g. Kraken2 Taxonomy Barplots](#19g-kraken2-taxonomy-barplots)
    - [20. Taxonomic Profiling Using HUMAnN/MetaPhlan](#20-taxonomic-profiling-using-humannmetaphlan)
      - [20a. Download and Install HUMAnN databases](#20a-download-and-install-humann-databases)
      - [20b. HUMAnN/MetaPhlAn Taxonomic Classification](#20b-humannmetaphlan-taxonomic-classification)
      - [20c. Merge Multiple Sample Functional Profiles](#20c-merge-multiple-sample-functional-profiles)
      - [20d. Split Results Tables](#20d-split-results-tables)
      - [20e. Normalize Gene Families and Pathway Abundances Tables](#20e-normalize-gene-families-and-pathway-abundances-tables)
      - [20f. Generate Normalized Gene-family Table Grouped by Kegg Orthologs (KOs)](#20f-generate-normalized-gene-family-table-grouped-by-kegg-orthologs-kos)
      - [20g. Combine MetaPhlan Taxonomy Tables](#20g-combine-metaphlan-taxonomy-tables)
      - [20h. Create MetaPhlan Species Count Table](#20h-create-metaphlan-species-count-table)
        - [20hi. Get Sample Read Counts](#20hi-get-sample-read-counts)
        - [20hii. Process MetaPhlan Taxonomy Table](#20hii-process-metaphlan-taxonomy-table)
      - [20i. Filter MetaPhlan Species Count Table](#20i-filter-metaphlan-species-count-table)
      - [20j. MetaPhlan Taxonomy Barplots](#20j-metaphlan-taxonomy-barplots)
      - [20k. Filter Humann Output](#20k-filter-humann-output)
      - [20l. Create Humann Function Heatmaps](#20l-create-humann-function-heatmaps)

---

# Software used

| Program      | Version | Relevant Links                                                                                                                                     |
| :----------- | :-----: | :------------------------------------------------------------------------------------------------------------------------------------------------- |
| BBTools      |  39.81  | [https://bbmap.org/](https://bbmap.org/)                                                                                                           |
| bit          | 1.13.15 | [https://github.com/AstrobioMike/bioinf_tools#bioinformatics-tools-bit](https://github.com/AstrobioMike/bioinf_tools#bioinformatics-tools-bit)     |
| bowtie2      |  2.5.5  | [https://bowtie-bio.sourceforge.net/bowtie2/index.shtml](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml)                                   |
| CAT          |   5.3   | [https://github.com/MGXlab/CAT_pack](https://github.com/MGXlab/CAT_pack)                                                                           |
| CheckM       |  1.2.5  | [https://github.com/Ecogenomics/CheckM](https://github.com/Ecogenomics/CheckM)                                                                     |
| fastp        |  1.3.1  | [https://github.com/OpenGene/fastp](https://github.com/OpenGene/fastp)                                                                             |
| FastQC       | 0.12.1  | [https://www.bioinformatics.babraham.ac.uk/projects/fastqc/](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)                           |
| GTDB-Tk      |  2.6.1  | [https://github.com/Ecogenomics/GTDBTk](https://github.com/Ecogenomics/GTDBTk)                                                                     |
| HUMAnN       |   3.9   | [https://github.com/biobakery/humann](https://github.com/biobakery/humann)                                                                         |
| Kaiju        | 1.10.1  | [https://bioinformatics-centre.github.io/kaiju/](https://bioinformatics-centre.github.io/kaiju/)                                                   |
| KEGG-Decoder |   1.3   | [https://github.com/bjtully/BioData/tree/master/KEGGDecoder#kegg-decoder](https://github.com/bjtully/BioData/tree/master/KEGGDecoder#kegg-decoder) |
| KOFamScan    |  1.3.0  | [https://github.com/takaram/kofam_scan#kofamscan](https://github.com/takaram/kofam_scan#kofamscan)                                                 |
| Kraken2      | 2.17.1  | [https://github.com/DerrickWood/kraken2](https://github.com/DerrickWood/kraken2)                                                                   |
| KrakenTools  |  1.2.1  | [https://ccb.jhu.edu/software/krakentools/](https://ccb.jhu.edu/software/krakentools/)                                                             |
| Krona        |  2.8.1  | [https://github.com/marbl/Krona/wiki](https://github.com/marbl/Krona/wiki)                                                                         |
| MEGAHIT      |  1.2.9  | [https://github.com/voutcn/megahit#megahit](https://github.com/voutcn/megahit#megahit)                                                             |
| MetaBAT      |  2.18   | [https://bitbucket.org/berkeleylab/metabat/src/master/](https://bitbucket.org/berkeleylab/metabat/src/master/)                                     |
| MetaPhlAn    |  4.1.0  | [https://github.com/biobakery/MetaPhlAn](https://github.com/biobakery/MetaPhlAn)                                                                   |
| MultiQC      | 1.27.1  | [https://multiqc.info/](https://multiqc.info/)                                                                                                     |
| Prodigal     |  2.6.3  | [https://github.com/hyattpd/Prodigal#prodigal](https://github.com/hyattpd/Prodigal#prodigal)                                                       |
| samtools     | 1.23.1  | [https://github.com/samtools/samtools#samtools](https://github.com/samtools/samtools#samtools)                                                     |
| SPAdes       |  4.2.0  | [https://github.com/ablab/spades](https://github.com/ablab/spades)                                                                                 |
| R            |  4.5.3  | [https://www.r-project.org](https://www.r-project.org)                                                                                             |
| dplyr        |  1.2.0  | [https://dplyr.tidyverse.org](https://dplyr.tidyverse.org)                                                                                         |
| ggplot2      |  4.0.2  | [https://ggplot2.tidyverse.org](https://ggplot2.tidyverse.org)                                                                                     |
| glue         |  1.8.0  | [https://glue.tidyverse.org](https://glue.tidyverse.org)                                                                                           |
| htmlwidgets  |  1.6.4  | [http://www.htmlwidgets.org](http://www.htmlwidgets.org)                                                                                           |
| magrittr     |  2.0.5  | [https://magrittr.tidyverse.org](https://magrittr.tidyverse.org)                                                                                   |
| pavian       | 1.2.0*  | [https://github.com/fbreitwieser/pavian](https://github.com/fbreitwieser/pavian)                                                                   |
| pheatmap     | 1.0.13  | [https://cran.r-project.org/package=pheatmap](https://cran.r-project.org/package=pheatmap)                                                         |
| phyloseq     | 1.54.0  | [https://www.bioconductor.org/packages/release/bioc/html/phyloseq.html](https://www.bioconductor.org/packages/release/bioc/html/phyloseq.html)     |
| plotly       | 4.12.0  | [https://plotly-r.com](https://plotly-r.com)                                                                                                       |
| purrr        |  1.2.1  | [https://purrr.tidyverse.org](https://purrr.tidyverse.org)                                                                                         |
| readr        |  2.2.0  | [https://readr.tidyverse.org](https://readr.tidyverse.org)                                                                                         |
| scales       |  1.4.0  | [https://scales.r-lib.org](https://scales.r-lib.org)                                                                                               |
| stringr      |  1.6.0  | [https://stringr.tidyverse.org](https://stringr.tidyverse.org)                                                                                     |
| tibble       |  3.3.1  | [https://tibble.tidyverse.org](https://tibble.tidyverse.org)                                                                                       |
| tidyr        |  1.3.2  | [https://tidyr.tidyverse.org](https://tidyr.tidyverse.org)                                                                                         |
> **Note:** pavian R package requires R version 4.0.5

---

# General processing overview with example commands

> Exact processing commands and output files listed in **bold** below are included with each Metagenomics Seq processed dataset in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/).  

## Pre-processing


### 1. Raw Data QC
> NOTE: It is NASA's policy that any human reads are to be removed from metagenomics datasets prior to being hosted in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/). As such this pipeline starts with fastq files that have had the human reads removed using the GeneLab Remove Human Reads pipeline ([GL-DPPD-7107-A](../../Remove_human_reads_from_raw_data/Pipeline_GL-DPPD-7105_Versions/GL-DPPD-7105-A.md))

#### 1a. Raw Data QC

```bash
fastqc -o HRrm_fastqc_output *HRrm_GLmetagenomics.fastq.gz
```

**Parameter Definitions:**

- `-o` – the output directory to store results
- `*HRrm_GLmetagenomics.fastq.gz` – the input reads are specified as a positional argument, and can be given all at once with wildcards like this, or as individual arguments with spaces in between them

**Input data:**

- *HRrm_GLmetagenomics.fastq.gz (raw reads, after human read removal)

**Output data:**

- *fastqc.html (FastQC output html summary)
- *fastqc.zip (FastQC output data)

#### 1b. Compile Raw Data QC

```bash
multiqc --zip-data-dir \
        --outdir raw_multiqc_report \
        --filename raw_multiqc_GLmetagenomics \
        --interactive 
        /path/to/raw_fastqc_output/
```

**Parameter Definitions:**

- `--zip-data-dir` - Compress the data directory.
- `--outdir` – Specifies the output directory to store results.
- `--filename` – Specifies the filename prefix of results.
- `--interactive` - Force multiqc to always create interactive javascript plots.
- `/path/to/raw_fastqc_output/` – The directory holding the output data from the FastQC run, provided as a positional argument.

**Input data:**

- /path/to/HRrm_fastqc_output/*fastqc.zip (FastQC output data, from [Step 1a](#1a-raw-data-qc))

**Output data:**

- **HRrm_multiqc_report/HRrm_multiqc_GLmetagenomics.html** (multiqc output html summary)
- **HRrm_multiqc_report/HRrm_multiqc_GLmetagenomics_data.zip** (directory containing multiqc output data)

<br>  

---

### 2. Quality filtering/trimming

#### 2a. Filter Quality and Trim Adapters

```bash
fastp --in1 sample_R1_HRrm_GLmetagenomics.fastq.gz --out1 temp_sample_R1_filtered.fastq.gz \
      --in2 sample_R2_HRrm_GLmetagenomics.fastq.gz --out2 temp_sample_R2_filtered.fastq.gz \
      --qualified_quality_phred  20 \
      --length_required 50 \
      --thread 2 \
      --detect_adapter_for_pe --disable_trim_poly_g \
      --json sample.fastp.json \
      --html sample.fastp.html 2> sample-fastp.log
```

**Parameter Definitions:**

- `--in1` - Specifies the forward input read file
- `--in2` - Specifies the reverse input read file
- `--in1` - Specifies the forward output read file
- `--in2` - Specifies the reverse output read file
- `--qualified_quality_phred` - the minimum quality value at which a base is qualified (default: 20)
- `--length_required` - the minimum read length. Shorter reads will be discarded (default: 50)
- `--thread` - number of worker threads (default: 2)
- `--detect_adapter_for_pe` - for paired end data, enable auto-detection of adapters
- `--disable_trim_poly_g` - explicitly disable automatic polyG trimming
- `--json` - Specifies the json format report file name
- `--html` - Specifies the html format report file name
- `2> sample-fastp.log` - Redirects the stderr output to a log file.

**Input Data:**

- *HRrm_GLmetagenomics.fastq.gz (raw sample reads with human reads removed)

**Output Data:**

- temp_*_filtered.fastq.gz (quality filtered and adapter trimmed reads)

#### 2b. Trim polyG

```bash
fastp --in1 temp_sample_R1_filtered.fastq.gz --out1 sample_R1_filtered_GLmetagenomics.fastq.gz \
      --in2 temp_sample_R2_filtered.fastq.gz --out2 sample_R2_filtered_GLmetagenomics.fastq.gz \
      --qualified_quality_phred  20 \
      --length_required 50 \
      --thread 2 \
      --detect_adapter_for_pe \
      --json sample.fastp.json \
      --html sample.fastp.html \
      --trim_poly_g 2> sample-fastp.log
```

**Parameter Definitions:**

- `--in1` - Specifies the forward input read file
- `--in2` - Specifies the reverse input read file
- `--in1` - Specifies the forward output read file
- `--in2` - Specifies the reverse output read file
- `--qualified_quality_phred` - the minimum quality value at which a base is qualified (default: 20)
- `--length_required` - the minimum read length. Shorter reads will be discarded (default: 50)
- `--thread` - number of worker threads (default: 2)
- `--detect_adapter_for_pe` - for paired end data, enable auto-detection of adapters
- `--json` - Specifies the json format report file name
- `--html` - Specifies the html format report file name
- `--trim_poly_g` - force polyG trimming
- `2> sample-fastp.log` - Redirects the stderr output to a log file.

**Input Data:**

- /path/to/filtered_data/temp_sample*.fastq.gz (round1 filtered/adapter trimmed reads, output from [Step 2a](#2a-filter-quality-and-trim-adapters)

**Output Data:**

- **\*filtered_GLmetagenomics.fastq.gz** (quality filtered and adapter trimmed, human removed reads)<br>

---

#### 2c. Filtered/Trimmed Data QC

```bash 
fastqc -o filtered_fastqc_output/ *filtered_GLmetagenomics.fastq.gz
```

**Parameter Definitions:**

-	`-o` – the output directory to store results  
-	`*filtered_GLmetagenomics.fastq.gz` – the input reads are specified as a positional argument, and can be given all at once with wildcards like this, or as individual arguments with spaces in between them  

**Input data:**

- *filtered_GLmetagenomics.fastq.gz (trimmed and filtered reads, from [Step 2b](#2b-trim-polyg))

**Output data:**

- *fastqc.html (FastQC output html summary)
- *fastqc.zip (FastQC output data)

#### 2d. Compile Filtered/Trimmed Data QC

```
multiqc --zip-data-dir \
        --outdir filtered_multiqc_report \
        --filename filtered_multiqc_GLmetagenomics \
        --interactive 
        /path/to/filtered_fastqc_output/
```

**Parameter Definitions:**

- `--zip-data-dir` - Compress the data directory.
- `--outdir` – Specifies the output directory to store results.
- `--filename` – Specifies the filename prefix of results.
- `--interactive` - Force multiqc to always create interactive javascript plots.
- `/path/to/filtered_fastqc_output/` – The directory holding the output data from the FastQC run, provided as a positional argument.

**Input Data:**

- /path/to/filtered_fastqc_output/*fastqc.zip (FastQC output data, from [Step 2c](#2c-filteredtrimmed-data-qc))

**Output Data:**

- **filtered_multiqc_report/filtered_multiqc_GLmetagenomics.html** (multiqc output html summary)
- **filtered_multiqc_report/filtered_multiqc_GLmetagenomics_data.zip** (zip archive containing multiqc output data)

<br>

### 3. R Environment Setup

> Taxonomy bar plots and heatmaps are performed in R.

#### 3a. Load libraries

```R
# load libraries
library(htmlwidgets)
library(pavian)
library(pheatmap)
library(phyloseq)

# load tidyverse libraries
library(dplyr)
library(ggplot2)
library(glue)
library(magrittr)
library(plotly)
library(purrr)
library(readr)
library(scales)
library(stringr)
library(tibble)
library(tidyr)
```

#### 3b. Define Custom Functions

#### get_last_assignment() <!-- omit in toc -->
<details>
  <summary>retrieves the last taxonomy assignment from a taxonomy string</summary>

  ```R
  get_last_assignment <- function(taxonomy_string, split_by = ';', remove_prefix = NULL) {

    # Spilt taxonomy string by the supplied delimiter 'split_by'
    # then convert the list of parts to a vector of parts
    split_names <- strsplit(x =  taxonomy_string , split = split_by) %>%
      unlist()
    # Get the last part of the split string
    level_name <- split_names[[length(split_names)]]
    
    if(level_name == "_"){
      return(taxonomy_string)
    }
    # remove an unwanted prefix if specified
    if(!is.null(remove_prefix)){
      level_name <- gsub(pattern = remove_prefix, replacement = "", x = level_name)
    }
    
    return(level_name)
  }
  ```

  **Function Parameter Definitions:**
  - `taxonomy_string` - a character string containing a list of taxonomy assignments separated by `split_by`
  - `split_by=` - a character string containing a regular expression used to split the `taxonomy_string`
  - `remove_prefix=` - a character string containing a regular expression to be matched and removed, default=`NULL`

  **Returns:** the last taxonomy assignment listed in the `taxonomy_string`

</details>

#### mutate_taxonomy() <!-- omit in toc -->
<details>
  <summary>mutate taxonomy column to contain the lowest taxonomy assignment</summary>

  ```R
  mutate_taxonomy <- function(df, taxonomy_column="taxonomy") {

    # make sure that the taxonomy column is always named taxonomy
    col_index <- which(colnames(df) == taxonomy_column)
    colnames(df)[col_index] <- "taxonomy"
    df <- df %>% dplyr::mutate(across(where(is.numeric), function(x) tidyr::replace_na(x, 0))) %>%
      dplyr::mutate(taxonomy=map_chr(taxonomy, .f = function(taxon_name = .x) {
        last_assignment <- get_last_assignment(taxon_name) 
        last_assignment  <- gsub(pattern = "\\[|\\]|'", replacement = "", x = last_assignment)
        trimws(last_assignment, which = "both")
      })) %>% 
      as.data.frame(check.names = FALSE, StringAsFactor = FALSE)
    # Ensure the taxonomy names are unique by aggregating duplicates
    df <- aggregate(.~taxonomy, data = df, FUN = sum)
    return(df)
  }
  ```

  **Custom Functions Used:**
  - [get_last_assignment()](#get_last_assignment)

  **Function Parameter Definitions:**
  - `df` - a dataframe containing the taxonomy assignments
  - `taxonomy_column=` - name of the column in `df` containing the taxonomy assignments, default="taxonomy"

  **Returns:** dataframe, `df`, with unique last taxonomy names stored in a column named "taxonomy"

</details>

#### process_kaiju_table() <!-- omit in toc -->
<details>
  <summary>reformat kaiju output table</summary>

  ```R
  process_kaiju_table <- function(file_path, taxon_col = "taxon_name") {
  
    # read input table
    kaiju_table <-  read_delim(file = file_path,
                               delim = "\t",
                               col_names = TRUE)

    # Create  a sample colname if the file column wasn't pre-edited
    if(colnames(kaiju_table)[1] ==  "file" ){
      kaiju_table <-  kaiju_table %>% rename(sample=file)
    }

    # filter out all kaiju database entries
    kaiju_table <- kaiju_table %>% 
      filter(!str_detect(sample, "dmp")) %>%
      mutate(sample=str_replace_all(sample, ".+/(.+)_kaiju.out", "\\1"))
 
    # keep only sample, reads, and taxonomy column (as defined by taxon_col argument) 
    # convert long dataframe to wide dataframe
    # mutate the taxonomy column such that it contains only lowest taxonomy assignment
    abs_abun_df <- kaiju_table %>%
      select(sample, reads, taxonomy=!!sym(taxon_col)) %>%
      pivot_wider(names_from = "sample", values_from = "reads", names_sort = TRUE) %>%
      mutate_taxonomy 
  
    # Set the taxon names as row names, drop the taxonomy column and convert to a matrix
    rownames(abs_abun_df) <- abs_abun_df[,"taxonomy"]
    abs_abun_df <- abs_abun_df[,-(which(colnames(abs_abun_df) == "taxonomy"))]
    abs_abun_matrix <- as.matrix(abs_abun_df)
    
    return(abs_abun_matrix)
  }
  ```

  **Custom Functions Used:**
  - [mutate_taxonomy()](#mutate_taxonomy)

  **Function Parameter Definitions:**
  - `file_path` - file path to the tab-delimited kaiju output table file
  - `taxon_col=`- name of the taxon column in the input data file, default="taxon_name"

  **Returns:** dataframe, `abs_abun_matrix`, with reformated kaiju output

</details>

#### merge_kraken_reports() <!-- omit in toc -->
<details>
  <summary>merge and process multiple kraken outputs to one species table</summary>

  ```R
  merge_kraken_reports <- function(reports_dir) {

    reports <- read_reports(reports_dir)

    # Retrieve sample names from file names
    samples <- names(reports) %>% str_split("-") %>% map_chr(function(x) pluck(x, 1))
    merged_reports  <- merge_reports2(reports, col_names = samples)
    taxonReads <- merged_reports$taxonReads
    cladeReads <- merged_reports$cladeReads
    tax_data <- merged_reports[["tax_data"]]

    species_table <- tax_data %>%
      bind_cols(cladeReads) %>%
      filter(taxRank %in% c("U", "S")) %>% # select unclassified and species rows 
      select(-contains("tax")) %>%
      zero_if_na() %>%
      filter(name != 0) %>% # drop unknown taxonomies
      group_by(name) %>%
      summarise(across(everything(), sum)) %>%
      ungroup() %>%
      as.data.frame %>%
      rename(species = name)

    # Set rownames as species name, drop species column
    # and convert table from dataframe to matrix
    species_names <- species_table[, "species"]
    rownames(species_table) <- species_names
    
    return(species_table)
  }
  ```

  **Function Parameter Definitions:**
  - `reports_dir` - path to a directory containing kraken2 reports 

  **Returns:** a kraken species count matrix, `species_table`, with samples and species as columns and rows, respectively.

</details>

#### get_abundant_features() <!-- omit in toc -->
<details>
  <summary>Find abundant features based on the sum of feature values</summary>
  
  ```R
  get_abundant_features <- function(mat, cpm_threshold = 1000){
  
    # Filtered out unassigned functions
    unassigned <- "UNMAPPED|UNGROUPED|UNINTEGRATED|Not annotated"
    mat <- mat %>%
      as.data.frame %>%
      rownames_to_column("Feature") %>%
      filter(str_detect(Feature, unassigned, negate = TRUE))
    rownames(mat) <- mat$Feature
    mat <- mat[, -1]

    features <- rowSums(mat, na.rm = TRUE) %>% sort()
    
    abund_features <- features[features > cpm_threshold] %>% names
    
    abund_features.m <- mat[abund_features, ]
    
    return(abund_features.m)
  }
  ```

  **Function Parameter Definitions:**
  - `mat` - a feature count matrix with features as rows and samples as columns
  - `cpm_threshold = 1000` - threshold to identify abundant features

  **Returns:** a matrix, `abund_features.m`, holding the features that pass the requested threshold
  
</details>

#### count_to_rel_abundance() <!-- omit in toc -->
<details>
  <summary>Convert species count matrix to relative abundance matrix</summary>

  ```R
  count_to_rel_abundance <- function(species_table) {

    # calculate species relative abundance per sample and
    # drop columns where none of the reads were classified or were non-microbial (NA)
    abund_table <- species_table %>%
      as.data.frame %>%
      mutate(across(everything(), function(x) (x/sum(x, na.rm = TRUE))*100)) %>%
        select(
          where( ~all(!is.na(.)))
        ) %>%
      rownames_to_column("Species")

    # Set rownames as species name and drop species column  
    rownames(abund_table) <- abund_table$Species
    abund_table <- abund_table[, -match(x = "Species", colnames(abund_table))] %>% t

    return(abund_table)
  }

  ```

  **Function Parameter Definitions:**
  - `species_table` - a species count matrix with samples and species as columns and rows, respectively.

  **Returns:** a species relative abundance matrix, `abund_table`, with samples and species as rows and columns, respectively.

</details>

#### filter_rare() <!-- omit in toc -->
<details>
  <summary>filter out rare and non_microbial taxonomy assignments based on relative abundance</summary>

  ```R
  filter_rare <- function(species_table, non_microbial, threshold=1) {
    
    # Drop species listed in 'non_microbial' regex
    clean_tab_count  <-  species_table %>% 
                         as.data.frame %>% 
                         rownames_to_column("Species") %>% 
                         filter(str_detect(Species, non_microbial, negate = TRUE))
    # Calculate species relative abundance
    clean_tab <- clean_tab_count %>%
      mutate(across(where(is.numeric), function(x) (x/sum(x, na.rm = TRUE))*100))
    # Set rownames as species name and drop species column
    rownames(clean_tab) <- clean_tab$Species
    clean_tab  <- clean_tab[, -1]
    
    # Get species with relative abundance less than `threshold` in all samples
    rare_species <- map(clean_tab, .f = function(col) rownames(clean_tab)[col < threshold])
    rare <- Reduce(intersect, rare_species)
    
    # Set rownames as species name and drop species column  
    rownames(clean_tab_count) <- clean_tab_count$Species
    clean_tab_count  <- clean_tab_count[,-1] 
    # Drop rare species
    abund_table <- clean_tab_count[!(rownames(clean_tab_count) %in% rare), ]
    
    return(abund_table)
  }
  ```

  **Function Parameter Definitions:**
  - `species_table` - the species matrix to filter with species and samples as rows and columns, respectively.
  - `non_microbial` - a regular expression denoting the names used to identify a species as non-microbial or unwanted
  - `threshold=` - abundance threshold used to determine if the relative abundance is rare, value denotes a percentage between 0 and 100.

  **Returns:** dataframe, `abund_table`, with rare and non_microbial/unwanted species removed

</details>

#### group_low_abund_taxa() <!-- omit in toc -->
<details>
  <summary>Group rare taxa or return a table with only rare taxa</summary>

  ```R
  group_low_abund_taxa <- function(abund_table, threshold = 0.05,
                                   rare_taxa = FALSE) {
    # If set to TRUE then a table with only the rare taxa will be returned 
    # initialize an empty vector that will contain the indices for the
    # low abundance columns/ taxa to group
    taxa_to_group <- c()
    # initialize the index variable of species with low abundance (taxa/columns)
    index <- 1
    
    #loop over every column or taxa check to see if the max abundance is less than the set threshold
    #if true save the index in the taxa_to_group vector variable
    for (column in ncol(abund_table)) {
      if(max(abund_table[,column], na.rm = TRUE) < threshold) {
        #print(column)
        taxa_to_group[index] <- column
        index = index + 1
      }
    }
    
    if(is.null(taxa_to_group)) {
      message(glue("Rare taxa were not grouped. please provide a higher 
                        threshold than {threshold} for grouping rare taxa, 
                        only numbers are allowed."))
      return(abund_table)
    }
    
    if(rare_taxa) {
      abund_table <- abund_table[,taxa_to_group,drop=FALSE]
    } else {
      #remove the low abundant taxa or columns
      abundant_taxa <-abund_table[,-(taxa_to_group), drop=FALSE]
      #get the rare taxa
      # rare_taxa <-abund_table[,taxa_to_group]
      rare_taxa <- subset(x = abund_table, select = taxa_to_group)
      #get the proportion of each sample that makes up the rare taxa
      rare <- rowSums(rare_taxa)
      #bind the abundant taxa to the rae taxa
      abund_table <- cbind(abundant_taxa,rare)
      #rename the columns i.e the taxa
      colnames(abund_table) <- c(colnames(abundant_taxa),"Rare")
    }
    
    return(abund_table)
  }
  ```

  **Function Parameter Definitions:**
  - `abund_table` - a relative abundance matrix with taxa as columns and  samples as rows
  - `rare_taxa` - a boolean specifying if only rare taxa should be returned
  - `threshold` - a max abundance threshold for defining taxa as rare

  **Returns:** a relative abundance matrix, `abund_table`, with rare taxa grouped or with non-rare taxa filtered out

</details>

#### make_plot() <!-- omit in toc -->
<details>
  <summary>Create stacked bar plots of relative abundance from input dataframes</summary>

  ```R
  # Make bar plot
  make_plot <- function(abund_table, metadata, custom_palette, publication_format,
                        samples_column="sample_id", prefix_to_remove="barcode"){
  
    abund_table_wide <- abund_table %>%
        as.data.frame() %>%
        rownames_to_column(samples_column) %>%
        inner_join(metadata) %>%
        select(!!!colnames(metadata), everything()) %>%
        mutate(!!samples_column := !!sym(samples_column) %>% str_remove(prefix_to_remove))
        
      
    abund_table_long <- abund_table_wide %>%
        pivot_longer(-colnames(metadata),
                     names_to = "Species",
                     values_to = "relative_abundance")
      
    p <- ggplot(abund_table_long, mapping = aes(x = !!sym(samples_column),
                                                y = relative_abundance, fill = Species)) +
         geom_col() +
         scale_fill_manual(values = custom_palette) +
         labs(x = NULL, y = "Relative Abundance (%)") +
         publication_format
    
    return(p)
  }
  ```

  **Function Parameter Definitions:**
  - `abund_table` - a relative abundance dataframe with rows summing to 100%
  - `metadata` - a metadata dataframe with samples as row and columns describing each sample
  - `custom_palette` - a vector of strings specifying a custom color palette for coloring plots
  - `publication_format` - a ggplot::theme object specifying a custom theme for plotting
  - `samples_column` - a character column specifying the column in `metadata` holding sample names, default is "Sample_ID"
  - `prefix_to_remove` - a string specifying a prefix or any character set to remove from sample names, default is "barcode"

  **Returns:** a relative abundance stacked bar plot, `p`

</details>

#### make_barplot()  <!-- omit in toc -->
<details>
  <summary>Parse Metadata and Feature table files in order to create stacked barplots of relative abundance.</summary>
  
  ```R
  make_barplot <- function(metadata_table_file, feature_table_file, 
                           feature_column = "species", samples_column = "sample_id", group_column = "group", 
                           output_prefix, assay_suffix = "_GLmetagenomics",
                           publication_format, custom_palette) {
    facet_by <- reformulate(group_column)
    # Prepare feature table
    feature_table <- read_delim(feature_table_file)
    rownames(feature_table) <- feature_table[[1]]
    feature_table <- feature_table[, -1]

    number_of_species <- nrow(feature_table)

    if (number_of_species > length(custom_palette)) {
      N <- number_of_species / length(custom_palette)
      custom_palette <- rep(custom_palette, times = N * 2)
    }

    # Prepare metadata
    metadata <- read_delim(metadata_table_file, delim = ",") %>% as.data.frame
    row.names(metadata) <- metadata[, samples_column]

    # compute abundances from counts
    abund_table <- count_to_rel_abundance(feature_table)

    metadata <- metadata %>%
                mutate(!!sym(group_column) := str_wrap(!!sym(group_column) %>%
                         str_replace_all("_", " "), width = 10)
                )
    
    # create plot
    p <- make_plot(abund_table, metadata, custom_palette, publication_format, samples_column) +
         facet_wrap(facet_by, nrow = 1, scales = "free_x", labeller = label_wrap_gen(width = 10)) +
         theme(axis.text.x = element_text(angle = 90))

    static_plot <- p
    number_of_species <- p$data$Species %>% unique() %>% length()
    # Don't save legend if the number of species to plot is greater than 30
    if(number_of_species > 30) {
      static_plot <- static_plot + theme(legend.position = "none")
    }
    
    width <- 2 * nrow(metadata) # 3.6 * number_of_samples
    if(width < 14) { width = 14 } # set minimum width to 14 inches
    if(width > 50) { width = 50 } # Cap plot with at 50 inches
    # Save Static
    ggsave(filename = glue("{output_prefix}_barplot{assay_suffix}.png"), 
           plot = static_plot,
           device = 'png', width = width,
           height = 10, units = 'in', dpi = 300 , limitsize = FALSE)

    # Save interactive
    htmlwidgets::saveWidget(ggplotly(p), glue("{output_prefix}_barplot{assay_suffix}.html"), selfcontained = TRUE)
  }
  ```

  **Custom Functions Used:**
  - [make_plot()](#make_plot)
  - [count_to_rel_abundance()](#count_to_rel_abundance)

  **Function Parameter Definitions:**
  - `metadata_table_file` - path to a file with samples as rows and columns describing each sample
  - `feature_table_file` - path to a tab separated samples feature table i.e. species/functions 
                           table with species/functions as the first column and samples as other columns.
  - `feature_column` - a character string containing the feature column name in the feature table ['Species', 'species', 'KO_ID'], default: "species".
  - `samples_column` - a character string specifying the column in `metadata` holding sample names, default: "sample_id"
  - `group_column` - a character string specifying the column in `metadata` used to facet/group plots, default: "group"
  - `output_prefix` - a character string specifying the unique name to add to the output file names 
                      used to denote the data type/source, for example "unfiltered-kaiju_species"
  - `assay_suffix` - a character string specifying the GeneLab assay suffix (default: "_GLmetagenomics")
  - `publication_format` - a ggplot::theme object specifying a custom theme for plotting, from [Step 3c](#3c-set-global-variables)
  - `custom_palette` - a vector of strings specifying a custom color palette for coloring plots, from [Step 3c](#3c-set-global-variables)

  **Output Data:** 2 barplot files, `{output_prefix}_barplot{assay_suffix}.png` and `{output_prefix}_barplot{assay_suffix}.html`, containing relative abundance stacked bar plot as output from [make_plot](#make_plot)
  
</details>

#### make_heatmap() <!-- omit in toc -->
<details>
  <summary>Creates heatmaps from a feature table file</summary>
  
  ```R
  make_heatmap <- function(metadata_table_file, feature_table_file, 
                           samples_column = "sample_id", group_column = "group", 
                           output_prefix, assay_suffix = "_GLmetagenomics",
                           custom_palette) {
    # Prepare feature table
    feature_table <- read_delim(feature_table_file) %>%  as.data.frame()
    rownames(feature_table) <- feature_table[[1]]
    feature_table <- feature_table[,-1] %>% as.matrix()
    colnames(feature_table) <-  colnames(feature_table) %>% str_remove_all("barcode")

    # Prepare metadata
    metadata <- read_delim(metadata_table_file) %>% as.data.frame()
    row.names(metadata) <- metadata[,samples_column] %>% str_remove_all("barcode")

    # GFet common samples and re-arrange feature table and metadata
    common_samples <- intersect(colnames(feature_table), rownames(metadata))
    feature_table <- feature_table[, common_samples]
    metadata <- metadata[common_samples,]
    metadata <- metadata %>% arrange(!!sym(group_column))

    # Create column annotation
    col_annotation <- as.data.frame(metadata)[, group_column, drop = FALSE]

    # Calculate output plot width and height
    number_of_samples <- ncol(feature_table)
    width <- 1 * number_of_samples
    if (width < 10) { width <- 10} # Set the minimum width to 10 inches
    if (width > 100) { width <- 100} # Set the maximum width to 100 inches
    number_of_features <- nrow(feature_table)
    height <- 0.2 * number_of_features
    if (height < 10) { height <- 10 } # Set the minimum height to 10 inches
    if (height > 100) { height <- 100 } # Set the maximum height to 100 inches (highest that won't generate an error)

    # Set colors by group
    groups <- metadata[[group_column]] %>%  unique()
    number_of_groups <-  length(groups)
    my_colors <- custom_palette[1:number_of_groups]
    names(my_colors) <- groups
    annotation_colors  <- list(my_colors)
    names(annotation_colors) <- group_column

    # create heatmap
    png(filename = glue("{output_prefix}_heatmap{assay_suffix}.png"), width = width,
        height = height, units = "in", res = 300)
    pheatmap(mat = feature_table[, rownames(col_annotation)],
             cluster_cols = FALSE,
             cluster_rows = FALSE,
             col = colorRampPalette(c('white','red'))(255), 
             angle_col = 0,
             display_numbers = TRUE,
             fontsize = 12,
             annotation_col = col_annotation,
             annotation_colors = annotation_colors,
             number_format = "%.0f")
    dev.off()

    sorted_features <- rowSums(feature_table) %>% sort(decreasing = TRUE)

    # Plot only top 50 features as it is often difficult to visualize all features at once
    if(length(sorted_features >= 50)) { 
      top50 <- sorted_features[1:50]

      png(filename = glue("{output_prefix}_top_50_heatmap{assay_suffix}.png"), width = width,
          height = 12, units = "in", res=300)
      pheatmap(mat = feature_table[names(top50), rownames(col_annotation)],
               cluster_cols = FALSE, 
               cluster_rows = FALSE,
               col = colorRampPalette(c('white','red'))(255), 
               angle_col = 90, 
               display_numbers = TRUE, 
               fontsize = 12, 
               annotation_col = col_annotation,
               annotation_colors = annotation_colors,
               number_format = "%.0f")
      dev.off()
    }
  }
  ```

  **Function Parameter Definitions:**
  - `metadata_table_file` - path to a file with samples as rows and columns describing each sample
  - `feature_table_file` - path to a tab separated samples feature table i.e. species/functions 
                           table with species/functions as the first column and samples as other columns.
  - `samples_column` - a character string specifying the column in `metadata` holding sample names, default: "sample_id"
  - `group_column` - a character string specifying the column in `metadata` used to facet/group plots, default: "group"
  - `output_prefix` - a character string specifying the unique name to add to the output file names 
                      used to denote the data type/source, for example "unfiltered-kaiju_species"
  - `assay_suffix` - a character string specifying the GeneLab assay suffix (default: "_GLmetagenomics")
  - `custom_palette` - a vector of strings specifying a custom color palette for coloring plots, from [Step 3c](#3c-set-global-variables)

  **Output Data:** 2 heatmap png files, `{output_prefix}_heatmap{assay_suffix}.png` and `{output_prefix}_top_50_heatmap{assay_suffix}.png`, of species/functions across samples from the input feature table
  
</details>

#### process_taxonomy() <!-- omit in toc -->
<details>
  <summary>process a taxonomy assignment table</summary>

  ```R
  process_taxonomy <- function(taxonomy, prefix='\\w__') { 
    
    taxonomy <- apply(X = taxonomy, MARGIN = 2, FUN = as.character) 

    # replace NAs and empty cells with "Other" and delete the `prefix` from taxonomy names
    for (rank in colnames(taxonomy)) {
      # Delete the taxonomy prefix
      taxonomy[,rank] <- gsub(pattern = prefix, x = taxonomy[, rank],
                              replacement = '')
      indices <- which(is.na(taxonomy[,rank]))
      taxonomy[indices, rank] <- rep(x = "Other", times=length(indices)) 
      # Replace empty cells with "Other"
      indices <- which(taxonomy[,rank] == "")
      taxonomy[indices,rank] <- rep(x = "Other", times=length(indices))
    }
    # Replace underscore with space
    taxonomy <- apply(X = taxonomy,MARGIN = 2,
                      FUN =  gsub,pattern = "_",replacement = " ") %>% 
      as.data.frame(stringAsfactor=FALSE)
    return(taxonomy)
  }
  ```
  **Function Parameter Definitions:**
  - `taxonomy` - is a taxonomy assignment dataframe with ranks [Phylum, Class .. Species] as columns and taxonomy assignments as rows
  - `prefix`  - is a regular expression specifying a character sequence to remove
                from taxon names

  **Returns:** dataframe, `taxonomy`, containing reformated taxonomy names
</details>

#### fix_names() <!-- omit in toc -->
<details>
  <summary>clean taxonomy names</summary>

  ```R
  fix_names<- function(taxonomy,stringToReplace="Other",suffix=";_"){
    
    for(index in seq_along(stringToReplace)){

      for (taxa_index in seq_along(taxonomy)) {    
        # Get the row indices of the current taxonomy columns
        # with rows matching the sting in `stringToReplace`
        indices <- grep(x = taxonomy[,taxa_index], pattern = stringToReplace)
        # Replace the value in that row with the value in the adjacent cell concatenated with `suffix`
        taxonomy[indices,taxa_index] <-
          paste0(taxonomy[indices,taxa_index-1],
                rep(x = suffix, times=length(indices)))
      }

    }
    return(taxonomy)
  }
  ```

  **Function Parameter Definitions:**
  - `taxonomy` -  taxonomy dataframe with taxonomy ranks as column names
  - `stringToReplace` - a regex string specifying what to replace
  - `suffix` - string specifying the replacement value

  **Returns:** dataframe, `taxonomy`, containing reformated/cleaned taxonomy names

</details>

#### read_taxonomy_table() <!-- omit in toc -->
<details>
  <summary>Read Assembly-based coverage annotation table</summary>

  ```R
  read_taxonomy_table <- function(df, sample_names){
  
    # Subset taxonomy portion (domain:species) of input table
    # and replace empty/Na domain assignments with "Unclassified"
    taxonomy_table <- df %>%
      select(domain:species) %>%
      mutate(domain=replace_na(domain, "Unclassified"))
    
    # Subset count table
    sample_names <- get_samples(df, sample_names)
    counts_table <- df %>% select(!!sample_names)

    # Mutate taxonomy names
    taxonomy_table  <- process_taxonomy(taxonomy_table)
    taxonomy_table <- fix_names(taxonomy_table, "Other", ";_")

    # Column bind taxonomy dataframe with species count dataframe
    df <- bind_cols(taxonomy_table, counts_table)
    
    return(df)
  }
  ```

  **Custom Functions Used:**
  [process_taxonomy](#process_taxonomy)
  [fix_names()](#fix_names)

  **Function Parameter Definitions:**
  - `df` - dataframe containing assembly-based coverage
  - `sample_names` - a character vector of sample names to keep in the final dataframe

  **Returns:** dataframe, `df`, containing cleaned taxonomy names and sample species count

</details>

#### get_samples() <!-- omit in toc -->
<details>
  <summary>retrieve sample names for which assemblies were generated</summary>

  ```R
  get_samples <- function(assembly_table_df, sample_names, end_col='species') {
    # Get common samples 
    cols <- colnames(df)
    index <- grep(end_col, cols)
    start <- grep(end_col, cols) + 1
    end <- (length(cols) - index)
    df_samples <- cols[start:end]
    sample_names <- intersect(df_samples, sample_names)

    return(sample_names)
  }
  ```

  **Function Parameter Definitions:**
  - `assembly_table_df` - dataframe containing assembly-based coverage
  - `sample_names` - a character vector of samples names to keep in the final dataframe
  - `end_col` - string containing the name of the last column

  **Returns:** a character vector, `sample_names`, of sample names that appear in both the assembly dataframe and the sample_names list

</details>

#### 3c. Set global variables

```R
# Define custom theme for plotting
publication_format <- theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.ticks.length=unit(-0.15, "cm"),
        axis.text.x=element_text(margin=ggplot2::margin(t=0.5,r=0,b=0,l=0,unit ="cm")),
        axis.text.y=element_text(margin=ggplot2::margin(t=0,r=0.5,b=0,l=0,unit ="cm")), 
        axis.title = element_text(size = 18,face ='bold.italic', color = 'black'), 
        axis.text = element_text(size = 16,face ='bold', color = 'black'),
        legend.position = 'right', legend.title = element_text(size = 15,face ='bold', color = 'black'),
        legend.text = element_text(size = 14,face ='bold', color = 'black'),
        strip.text =  element_text(size = 14,face ='bold', color = 'black'))

# Define custom palette for plotting
custom_palette <- c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F", "#FF7F00",
                    "#CAB2D6","#6A3D9A","#FF00FFFF","#B15928","#000000","#FFC0CBFF","#8B864EFF","#F0027F",
                    "#666666","#1B9E77", "#E6AB02","#A6761D","#FFFF00FF","#FFFF99","#00FFFFFF",
                    "#B2182B","#FDDBC7","#D1E5F0","#CC0033","#FF00CC","#330033",
                    "#999933","#FF9933","#FFFAFAFF",colors()) 
# Drop white colors
custom_palette <- custom_palette[-c(21:23,
                                    grep(pattern = "white|snow|azure|gray|#FFFAFAFF|aliceblue",
                                         x = custom_palette, 
                                         ignore.case = TRUE)
                                   )
                                ]                      
```

**Input Data:** 

*No input data required*

**Output Data:**

- `publication_format` (a ggplot::theme object specifying a custom theme for plotting)
- `custom_palette` (a vector of strings specifying a custom color palette for coloring plots)

<br>  

---

## Assembly-based Processing


### 4. Sample assembly
```
megahit -1 sample_R1_filtered_GLmetagenomics.fastq.gz -2 sample_R2_filtered_GLmetagenomics.fastq.gz \
        -o sample-assembly -t NumberOfThreads --min-contig-length 500 > sample-assembly.log 2>&1
```

**Parameter Definitions:**  

-	`-1 and -2` – specifies the input forward and reverse reads (if single-end data, then neither `-1` nor `-2` are used, instead single-end reads are passed to `-r`)
-	`-o` – specifies output directory
-	`-t` – specifies the number of threads to use
-	`--min-contig-length` – specifies the minimum contig length to write out
-	`> sample-assembly.log 2>&1` – sends stdout/stderr to log file

**Input data:**

- *_R[12]_filtered_GLmetagenomics.fastq.gz (filtered/trimmed reads from [Step 2b](#2b-trim-polyg) above)

**Output data:**

- sample-assembly/final.contigs.fa (assembly file)
- sample-assembly.log (log file)

<br>

---

### 5. Rename Contigs and Summarize Assemblies

#### 5a. Rename Contig Headers

```bash
bit-rename-fasta-headers -i sample/final.contigs.fa \
                         -w c_sample \
                         -o sample-assembly.fasta
```

**Parameter Definitions:**  

- `-i` – Specifies the input fasta file.
- `-w` – Specifies the wanted header prefix (a number will be appended for each contig), starts with a "c" to ensure they won't start with a number which can be problematic.
- `-o` – Specifies the output fasta file.


**Input data:**

- sample/final.contigs.fa (assembly file from [Step 4](#4-sample-assembly))

**Output files:**

- **sample-assembly_GLmetagenomics.fasta** (contig-renamed assembly file)

#### 5b. Summarize assemblies

```bash
bit-summarize-assembly -o assembly-summaries_GLmetagenomics.tsv \
                       *assembly_GLmetagenomics.fasta

# test assembly fasta files for absence of contigs
for assembly_file in *-assembly_GLmetagenomics.fasta; do 
  sample_id=${assembly_file%-assembly_GLmetagenomics.fasta} 
  if [ ! -s ${assembly_file} ]; then 
    printf "${sample_id}\tNo contigs assembled\n" >> Failed-assemblies_GLmetagenomics.tsv
  fi
done
```

**Parameter Definitions:**  

-	`-o` – Specifies the output summary table.
- `*-assembly_GLmetagenomics.fasta`	– Specifies the input assemblies to summarize, provided as positional arguments

**Input data:**

- *-assembly_GLmetagenomics.fasta (contig-renamed assembly files from [Step 5a](#5a-rename-contig-headers))

**Output files:**

- **assembly-summaries_GLmetagenomics.tsv** (table of assembly summary statistics)
- **Failed-assemblies_GLmetagenomics.tsv** (list of samples with no assembled contigs. Only present if no contigs were generated for at least one sample.)

<br>

---

### 6. Gene prediction

#### 6a. Generate Gene Predictions

```bash
prodigal -a sample-genes.faa \
         -d sample-genes.fasta \
         -f gff \
         -p meta \
         -c \
         -q \
         -o sample-genes_GLmetagenomics.gff \
         -i sample-assembly_GLmetagenomics.fasta
```
**Parameter Definitions:**

- `-a` – Specifies the output amino acid sequences file.
- `-d` – Specifies the output nucleotide sequences file.
- `-f` – Specifies the gene-calls output format, gff = GFF format.
- `-p` – Specifies which mode to run the gene-caller in. 
- `-c` – No incomplete genes reported. 
- `-q` – Run in quiet mode (don’t output process on each contig). 
- `-o` – Specifies the name of the output gene-calls file. 
- `-i` – Specifies the input assembly file.

**Input data:**

- sample-assembly.fasta (contig-renamed assembly file from [Step 5a](#5a-rename-contig-headers))

**Output data:**

- sample-genes.faa (gene-calls amino-acid fasta file)
- sample-genes.fasta (gene-calls nucleotide fasta file)
- **sample-genes.gff** (gene-calls in general feature format)

<br>

#### 6b. Remove Line Wraps In Gene Prediction Output

```bash
bit-remove-wraps sample-genes.faa > sample-genes.faa.tmp 2> /dev/null
mv sample-genes.faa.tmp sample-genes_GLmetagenomics.faa

bit-remove-wraps sample-genes.fasta > sample-genes.fasta.tmp 2> /dev/null
mv sample-genes.fasta.tmp sample-genes_GLmetagenomics.fasta
```

**Input Data:**

- sample-genes.faa (gene-calls amino-acid fasta file, output from [Step 6a](#6a-generate-gene-predictions))
- sample-genes.fasta (gene-calls nucleotide fasta file, output from [Step 6a](#6a-generate-gene-predictions))

**Output Data:**

- **sample-genes_GLmetagenomics.faa** (gene-calls amino-acid fasta file with line wraps removed)
- **sample-genes_GLmetagenomics.fasta** (gene-calls nucleotide fasta file with line wraps removed)

---

### 7. Functional annotation

> **Note:**  
> The annotation process overwrites the same temporary directory by default. When running multiple processes at a time, it is necessary to specify a specific temporary directory with the `--tmp-dir` argument as shown below.

#### 7a. Download reference database of HMM models

> **Note:** This step only needs to be done once.

```bash
curl -LO ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz
curl -LO ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz
tar -xzvf profiles.tar.gz
gunzip ko_list.gz 
```

#### 7b. Run KEGG annotation

```bash
exec_annotation -p profiles/ \
                -k ko_list \
                --cpu NumberOfThreads \
                -f detail-tsv \
                -o sample-KO-tab.tmp \
                --tmp-dir sample-tmp-KO \
                --report-unannotated \
                sample-genes_GLmetagenomics.faa 
```

**Parameter Definitions:**
- `-p` – Specifies the directory holding the downloaded reference HMMs.
- `-k` – Specifies the downloaded reference KO  (Kegg Orthology) terms. 
- `--cpu` – Specifies the number of searches to run in parallel.
- `-f` – Specifies the output format.
- `-o` – Specifies the output file name.
- `--tmp-dir` – Specifies the temporary directory to write to (needed if running more than one process concurrently, see Note above).
- `--report-unannotated` – Specifies to generate an output for each entry, event when no KO is assigned.
- `sample-genes_GLmetagenomics.faa` – Specifies the input file, provided as a positional argument. 

**Input data:**

- sample-genes_GLmetagenomics.faa (amino-acid fasta file, from [Step 6b](#6b-remove-line-wraps-in-gene-prediction-output))
- profiles/ (reference directory holding the KO HMMs, downloaded in [Step 7a](#7a-download-reference-database-of-hmm-models))
- ko_list (reference list of KOs to scan for, downloaded in [Step 7a](#7a-download-reference-database-of-hmm-models))

**Output data:**

- sample-KO-tab.tmp (table of KO annotations assigned to gene IDs)

#### 7c. Filter KO Outputs
*Filter KO outputs to retain only those passing the KO-specific score and top hits.*

```bash
bit-filter-KOFamScan-results -i sample-KO-tab.tmp \
                             -o sample-annotations.tsv

# removing temporary files
rm -rf sample-tmp-KO/ sample-KO-annots.tmp
```

**Parameter Definitions:**  

- `-i` – Specifies the input table.
- `-o` – Specifies the output table.

**Input data:**

- sample-KO-tab.tmp (table of KO annotations assigned to gene IDs from [Step 7b](#7b-run-kegg-annotation))

**Output data:**

- sample-annotations.tsv (table of KO annotations assigned to gene IDs)

<br>

---

### 8. Taxonomic classification

#### 8a. Pull and Unpack Pre-built Reference DB

```bash
wget tbb.bio.uu.nl/bastiaan/CAT_prepare/CAT_prepare_20200618.tar.gz
tar -xvzf CAT_prepare_20200618.tar.gz
```

#### 8b. Run Taxonomic Classification

```bash
CAT contigs -c sample-assembly_GLmetagenomics.fasta \
            -d CAT_prepare_20200618/2020-06-18_database/ \
            -t CAT_prepare_20200618/2020-06-18_taxonomy/ \
            -p sample-genes_GLmetagenomics.faa \
            -o sample-tax-out.tmp \
            -n NumberOfThreads -r 3 \
            --top 4 \
            --I_know_what_Im_doing \
            --no-stars
```

**Parameter Definitions:**  

- `-c` – Specifies the input assembly fasta file.
- `-d` – Specifies the CAT reference sequence database.
- `-t` – Specifies the CAT reference taxonomy database.
- `-p` – Specifies the input protein fasta file.
- `-o` – Specifies the output file prefix.
- `-n` – Specifies the number of CPU cores to use.
- `-r` – Specifies the number of top protein hits to consider in assigning taxonomy.
- `--top` – Specifies the number of protein alignments to store.
- `--I_know_what_Im_doing` – Allows us to alter the `--top` parameter.
- `--no-stars` - Suppress marking of suggestive taxonomic assignments.

**Input data:**

- CAT_prepare_20200618/2020-06-18_database/ (directory holding the CAT reference sequence database, output from [Step 8a](#8a-pull-and-unpack-pre-built-reference-db))
- CAT_prepare_20200618/2020-06-18_taxonomy/ (directory holding the CAT reference taxonomy database, output from [Step 8a](#8a-pull-and-unpack-pre-built-reference-db)
- sample-assembly_GLmetagenomics.fasta (assembly file from [Step 5a](#5a-rename-contig-headers))
- sample-genes_GLmetagenomics.faa (gene-calls amino-acid fasta file from [Step 6](#6b-remove-line-wraps-in-gene-prediction-output))

**Output data:**

- sample-tax-out.tmp.ORF2LCA.txt (gene-calls taxonomy file)
- sample-tax-out.tmp.contig2classification.txt (contig taxonomy file)

#### 8c. Add taxonomy info from taxids to genes

```bash
CAT add_names -i sample-tax-out.tmp.ORF2LCA.txt \
              -o sample-gene-tax-out.tmp \
              -t CAT_prepare_20200618/2020-06-18_taxonomy/ \
              --only_official \
              --exclude-scores
```

**Parameter Definitions:**  

- `-i` – Specifies the input taxonomy file.
- `-o` – Specifies the output file name.
- `-t` – Specifies the CAT reference taxonomy database.
- `--only_official` – Specifies to add only standard taxonomic ranks.
- `--exclude-scores` - Specifies to exclude bit-score support scores in the lineage.

**Input data:**

- sample-tax-out.tmp.ORF2LCA.txt (gene-calls taxonomy file from [Step 8b](#8b-run-taxonomic-classification))
- CAT_prepare_20200618/2020-06-18_taxonomy/ (directory holding the CAT reference taxonomy database, output from [Step 8a](#8a-pull-and-unpack-pre-built-reference-db)

**Output data:**

- sample-gene-tax-out.tmp (gene-calls taxonomy file with lineage info added)

#### 8d. Add Taxonomy Info From Taxids To Contigs

```bash
CAT add_names -i sample-tax-out.tmp.contig2classification.txt \
              -o sample-contig-tax-out.tmp \
              -t CAT_prepare_20200618/2020-06-18_taxonomy/ \
              --only_official \
              --exclude-scores
```

**Parameter Definitions:**  

- `-i` – Specifies the input taxonomy file.
- `-o` – Specifies the output file name.
- `-t` – Specifies the CAT reference taxonomy database.
- `--only_official` – Specifies to add only standard taxonomic ranks.
- `--exclude-scores` - Specifies to exclude bit-score support scores in the lineage.

**Input data:**

- sample-tax-out.tmp.contig2classification.txt (contig taxonomy file from [Step 8b](#8b-run-taxonomic-classification))
- CAT_prepare_20200618/2020-06-18_taxonomy/ (directory holding the CAT reference taxonomy database, output from [Step 8a](#8a-pull-and-unpack-pre-built-reference-db)

**Output data:**

- sample-contig-tax-out.tmp (contig taxonomy file with lineage info added)

#### 8e. Format Gene-level Output With awk and sed

```bash
awk -F $'\t' ' BEGIN { OFS=FS } { if ( $3 == "lineage" ) { print $1,$3,$5,$6,$7,$8,$9,$10,$11 } \
    else if ( $2 == "ORF has no hit to database" || $2 ~ /^no taxid found/ ) \
    { print $1,"NA","NA","NA","NA","NA","NA","NA","NA" } else { n=split($3,lineage,";"); \
    print $1,lineage[n],$5,$6,$7,$8,$9,$10,$11 } } ' sample-gene-tax-out.tmp | \
    sed 's/no support/NA/g' | sed 's/superkingdom/domain/' | sed 's/# ORF/gene_ID/' | \
    sed 's/lineage/taxid/'  > sample-gene-tax-out.tsv
```

**Input Data:**

- sample-gene-tax-out.tmp (gene-calls taxonomy file with lineage info added from [Step 8c](#8c-add-taxonomy-info-from-taxids-to-genes))

**Output Data:**

- sample-gene-tax.tsv (reformatted gene-calls taxonomy file with lineage info)

#### 8f. Format Contig-level Output With awk and sed

```bash
awk -F $'\t' ' BEGIN { OFS=FS } { if ( $2 == "classification" ) { print $1,$4,$6,$7,$8,$9,$10,$11,$12 } \
    else if ( $2 == "no taxid assigned" ) { print $1,"NA","NA","NA","NA","NA","NA","NA","NA" } \
    else { n=split($4,lineage,";"); print $1,lineage[n],$6,$7,$8,$9,$10,$11,$12 } } ' sample-contig-tax-out.tmp | \
    sed 's/no support/NA/g' | sed 's/superkingdom/domain/' | sed 's/^# contig/contig_ID/' | \
    sed 's/lineage/taxid/' > sample-contig-tax-out.tsv

  # clearing intermediate files
rm sample*.tmp*
```

**Input data:**

- sample-contig-tax-out.tmp (contig taxonomy file with lineage info added from [Step 8d](#8d-add-taxonomy-info-from-taxids-to-contigs))

**Output data:**

- sample-contig-tax-out.tsv (reformatted contig taxonomy file with lineage info)

<br>

---

### 9. Read-Mapping

#### 9a. Build reference index

```bash
bowtie2-build sample-assembly_GLmetagenomics.fasta sample-assembly-bt-index
```

**Parameter Definitions:**  

-	`sample-assembly_GLmetagenomics.fasta` - first positional argument specifies the input assembly
-	`sample-assembly-bt-index` - second positional argument specifies the prefix of the output index files

**Input Data:**

- `sample-assembly_GLmetagenomics.fasta` (contig-renamed assembly file, output from [Step 5a](#5a-rename-contig-headers))

**Output Data:**

- `sample-assembly-bt-index*` - the bowtie2 index files

#### 9b. Align Reads to Sample Assembly

```bash
bowtie2 --mm --quiet --threads ${task.cpus} \
        -x sample-index \
        -1 sample_R1_filtered_GLmetagenomics.fastq.gz \
        -2 sample_R2_filtered_GLmetagenomics.fastq.gz \
        --no-unal > sample.sam  2> sample-mapping-info_GLmetagenomics.txt 
```

**Parameter Definitions:**  

- `--mm` - Use memory-mapped I/O to load the index.
- `--quiet` - Print only error messages.
- `--threads` - Number of parallel processing threads.
- `-x` - specifies the prefix of the reference index files to map to, generated by bowtie2-build
-	`-1` - specifies the forward reads to map
- `-2` – specifies the reverse reads to map
- `--no-unal` - Suppress SAM records for reads that did not align.
- `> sample.sam` - Redirects the output of the map reads command to a SAM file.
- `2> sample-mapping-info_GLmetagenomics.txt` – capture the printed summary results in a log file

**Input Data**

- sample-assembly-bt-index (bowtie2 index files, output from [Step 9a](#9a-build-reference-index))
- *_R[12]_filtered_GLmetagenomics.fastq.gz (filtered and trimmed sample reads, output from [Step 2b](#2b-trim-polyg))

**Output Data**

- sample.sam (reads aligned to sample assembly in SAM format)
- **sample-mapping-info_GLmetagenomics.txt** (read mapping information)

#### 9c. Sort Assembly Alignments

```bash
# Sort Sam, convert to bam and create index
samtools sort --threads NumberOfThreads \
              -o sample_GLmetagenomics.bam \
              sample.sam > sample_sort.log 2>&1
```

**Parameter Definitions:**

*samtools sort*
- `--threads` - Number of parallel processing threads to use.
- `-o` - Specifies the output file for the sorted aligned reads.
- `sample.sam` - Positional argument specifying the input SAM file.
- `> sample_sort.log 2>&1` - Redirects the standard output and standard error to a separate file.

**Input Data:**

- sample.sam (reads aligned to sample assembly, output from [Step 9b](#9b-align-reads-to-sample-assembly))

**Output Data:**

- **sample_GLmetagenomics.bam** (sorted mapping to sample assembly, in BAM format)

<br>

---

### 10. Get Coverage Information and Filter Based On Detection
> **Note:**  
> “Detection” is a metric of what proportion of a reference sequence recruited reads (see [here](https://merenlab.org/2017/05/08/anvio-views/#detection)). Filtering based on detection is one way of helping to mitigate non-specific read-recruitment.

#### 10a. Filter Coverage Levels Based On Detection

```bash
# pileup.sh comes from the BBTools package
pileup.sh -in sample_GLmetagenomics.bam \
          fastaorf=sample-genes_GLmetagenomics.fasta \
          outorf=sample-gene-cov-and-det.tmp \
          out=sample-contig-cov-and-det.tmp
```

**Parameter Definitions:**  

- `-in` – Specifies the input BAM file.
- `fastaorf=` – Specifies the input gene-calls nucleotide fasta file.
- `outorf=` – Specifies the output gene-coverage tsv file name.
- `out=` – Specifies the output contig-coverage tsv file name.

**Input Data:**

- sample_GLmetagenomics.bam (sorted mapping to sample assembly BAM file, output from [Step 9c](#9c-sort-assembly-alignments))
- sample-genes_GLmetagenomics.fasta (gene-calls nucleotide fasta file, output from [Step 6b](#6b-remove-line-wraps-in-gene-prediction-output))

**Output Data:**

- sample-gene-cov-and-det.tmp (gene-coverage tsv file)
- sample-contig-cov-and-det.tmp (contig-coverage tsv file)

#### 10b. Filter Gene and Contig Coverage Based On Detection

> *The following commands filter gene and contig coverage tsv files to only keep genes and contigs with at least 50% detection (as defined above) then parse the tables to retain only gene IDs and respective coverage.*

```bash
# Filtering gene coverage
grep -v "#" sample-gene-cov-and-det.tmp | \
awk -F $'\t' ' BEGIN { OFS=FS } { if ( $10 <= 0.5 ) $4 = 0 } \
     { print $1,$4 } ' > sample-gene-cov.tmp

cat <( printf "gene_ID\tcoverage\n" ) sample-gene-cov.tmp > sample-gene-coverages.tsv

#Filtering contig coverage based on requiring 50% detection and parsing down to just contig ID and coverage:
grep -v "#" sample-contig-cov-and-det.tmp | awk -F $'\t' ' BEGIN { OFS=FS } { if ( $5 <= 50 ) $2 = 0 } \
     { print $1,$2 } ' > sample-contig-cov.tmp

cat <( printf "contig_ID\tcoverage\n" ) sample-contig-cov.tmp > sample-contig-coverages.tsv

# removing intermediate files
rm sample-*.tmp
```

**Input data:**

- sample-gene-cov-and-det.tmp (temporary gene-coverage tsv file, output from [Step 10a](#10a-filter-coverage-levels-based-on-detection))
- sample-contig-cov-and-det.tmp (temporary contig-coverage tsv file, output from [Step 10a](#10a-filter-coverage-levels-based-on-detection))

**Output data:**

- sample-gene-coverages.tsv (table with gene-level coverages)
- sample-contig-coverages.tsv (table with contig-level coverages)

<br>

---

### 11. Combine Gene-level Coverage, Taxonomy, and Functional Annotations For Each Sample
> **Notes**  
> Just uses `paste`, `sed`, and `awk` standard Unix commands to combine gene-level coverage, taxonomy, and functional annotations into one table for each sample. 

```bash
paste <( tail -n +2 sample-gene-coverages.tsv | sort -V -k 1 ) \
      <( tail -n +2 sample-annotations.tsv | sort -V -k 1 | cut -f 2- ) \
      <( tail -n +2 sample-gene-tax-out.tsv | sort -V -k 1 | cut -f 2- ) \
      > sample-gene-tab.tmp

paste <( head -n 1 sample-gene-coverages.tsv ) \
      <( head -n 1 sample-annotations.tsv | cut -f 2- ) \
      <( head -n 1 sample-gene-tax-out.tsv | cut -f 2- ) \
      > sample-header.tmp

cat sample-header.tmp sample-gene-tab.tmp > sample-gene-coverage-annotation-and-tax_GLmetagenomics.tsv

  # removing intermediate files
rm sample*tmp sample-gene-coverages.tsv sample-annotations.tsv sample-gene-tax-out.tsv
```

**Input data:**

- sample-gene-coverages.tsv (table with gene-level coverages from [Step 10b](#10b-filter-gene-and-contig-coverage-based-on-detection))
- sample-annotations.tsv (table of KO annotations assigned to gene IDs from [Step 7c](#7c-filter-ko-outputs))
- sample-gene-tax-out.tsv (gene-level taxonomic classifications from [Step 8f](#8f-format-contig-level-output-with-awk-and-sed))

**Output data:**

- **sample-gene-coverage-annotation-and-tax_GLmetagenomics.tsv** (table with combined gene coverage, annotation, and taxonomy info)

<br>

---

### 12. Combine Contig-level Coverage and Taxonomy For Each Sample
> **Note:**  
> Just uses `paste`, `sed`, and `awk` standard Unix commands to combine contig-level coverage and taxonomy into one table for each sample.

```bash
paste <( tail -n +2 sample-contig-coverages.tsv | sort -V -k 1 ) \
      <( tail -n +2 sample-contig-tax-out.tsv | sort -V -k 1 | cut -f 2- ) \
      > sample-contig.tmp

paste <( head -n 1 sample-contig-coverages.tsv ) \
      <( head -n 1 sample-contig-tax-out.tsv | cut -f 2- ) \
      > sample-contig-header.tmp
      
cat sample-contig-header.tmp sample-contig.tmp > sample-contig-coverage-and-tax_GLmetagenomics.tsv

# removing intermediate files
rm sample*tmp sample-contig-coverages.tsv sample-contig-tax-out.tsv
```

**Input data:**

- sample-contig-coverages.tsv (table with contig-level coverages from [Step 10b](#10b-filter-gene-and-contig-coverage-based-on-detection))
- sample-contig-tax-out.tsv (contig-level taxonomic classifications from [Step 8f](#8f-format-contig-level-output-with-awk-and-sed))

**Output data:**

- **sample-contig-coverage-and-tax_GLmetagenomics.tsv** (table with combined contig coverage and taxonomy info)

<br>

---

### 13. Generating normalized, gene- and contig-level coverage summary tables of KO-annotations and taxonomy across samples

> **Notes**  
> * To combine across samples to generate these summary tables, we need the same "units". This is done for annotations based on the assigned KO terms, and all non-annotated functions are included together as "Not annotated". It is done for taxonomic classifications based on taxids (full lineages included in the table), and any not classified are included together as "Not classified". 
> * The values we are working with are coverage per gene (so they are number of bases recruited to the gene normalized by the length of the gene). These have been normalized by making the total coverage of a sample 1,000,000 and setting each individual gene-level coverage its proportion of that 1,000,000 total. So basically percent, but out of 1,000,000 instead of 100 to make the numbers more friendly. 

#### 13a. Generate Gene-level Coverage Summary Tables

```bash
bit-GL-combine-KO-and-tax-tables *-gene-coverage-annotation-and-tax_GLmetagenomics.tsv \
                                 -o Combined
# add assay specific suffix
mv "Combined-gene-level-KO-function-coverages-CPM.tsv Combined-gene-level-KO-function-coverages-CPM_GLmetagenomics.tsv"
mv "Combined-gene-level-KO-function-coverages-CPM.tsv Combined-gene-level-KO-function-coverages-CPM_GLmetagenomics.tsv"
mv "Combined-gene-level-KO-function-coverages.tsv Combined-gene-level-KO-function-coverages_GLmetagenomics.tsv"
mv "Combined-gene-level-taxonomy-coverages.tsv Combined-gene-level-taxonomy-coverages_GLmetagenomics.tsv"
```

**Parameter Definitions:**  

- `*-gene-coverage-annotation-and-tax_GLmetagenomics.tsv` - Positional arguments specifying the input tsv files, can be provided as a space-delimited list of files, or with wildcards like above.
- `-o` – Specifies the output file prefix.

**Input data:**

- *-gene-coverage-annotation-and-tax_GLmetagenomics.tsv (tables with combined gene coverage, annotation, and taxonomy info generated for individual samples from [Step 11](#11-combine-gene-level-coverage-taxonomy-and-functional-annotations-for-each-sample))

**Output data:**

- **Combined-gene-level-KO-function-coverages-CPM_GLmetagenomics.tsv** (table with all samples combined based on KO annotations; normalized to coverage per million genes covered)
- **Combined-gene-level-taxonomy-coverages-CPM_GLmetagenomics.tsv** (table with all samples combined based on gene-level taxonomic classifications; normalized to coverage per million genes covered)
- **Combined-gene-level-KO-function-coverages_GLmetagenomics.tsv** (table with all samples combined based on KO annotations)
- **Combined-gene-level-taxonomy-coverages_GLmetagenomics.tsv** (table with all samples combined based on gene-level taxonomic classifications)

#### 13b. Generate Contig-level Coverage Summary Tables

```bash
bit-GL-combine-contig-tax-tables *-contig-coverage-and-tax_GLmetagenomics.tsv -o Combined
```

**Parameter Definitions:**  

- `*-contig-coverage-and-tax_GLmetagenomics.tsv` - Positional arguments specifying the input tsv files, can be provided as a space-delimited list of files, or with wildcards like above.
- `-o` – Specifies the output file prefix.

**Input data:**

- *-contig-coverage-annotation-and-tax_GLmetagenomics.tsv (tables with combined contig coverage, annotation, and taxonomy info generated for individual samples from [Step 12](#12-combine-contig-level-coverage-and-taxonomy-for-each-sample))

**Output data:**

- **Combined-contig-level-taxonomy-coverages-CPM_GLmetagenomics.tsv** (table with all samples combined based on contig-level taxonomic classifications; normalized to coverage per million genes covered)
- **Combined-contig-level-taxonomy-coverages_GLmetagenomics.tsv** (table with all samples combined based on contig-level taxonomic classifications)

<br>

---

### 14. **M**etagenome-**A**ssembled **G**enome (MAG) recovery

#### 14a. Bin contigs

```bash
jgi_summarize_bam_contig_depths --outputDepth sample-metabat-assembly-depth.tsv \
                                --percentIdentity 97 \
                                --minContigLength 1000 \
                                --minContigDepth 1.0  \
                                --referenceFasta sample-assembly_GLmetagenomics.fasta \
                                sample.bam

metabat2  --inFile sample-assembly_GLmetagenomics.fasta \
          --outFile sample \
          --abdFile sample-metabat-assembly-depth_GLmetagenomics.tsv \
          -t NumberOfThreads

mkdir sample-bins
mv sample*bin*.fasta sample-bins
zip -r sample-bins_GLmetagenomics.zip sample-bins
```

**Parameter Definitions:**  

*jgi_summarize_bam_contig_depths*
-  `--outputDepth` – Specifies the output depth file name.
-  `--percentIdentity` – Minimum end-to-end percent identity of a mapped read to be included.
-  `--minContigLength` – Minimum contig length to include.
-  `--minContigDepth` – Minimum contig depth to include.
-  `--referenceFasta` – Specifies the input assembly fasta file.
-  `sample_GLmetagenomics.bam` – Input alignment BAM file, specified as a positional argument.

*metabat2*
-  `--inFile` - Specifies the input assembly fasta file.
-  `--outFile` - Specifies the prefix of the identified bins output files.
-  `--abdFile` - The depth file generated by the previous `jgi_summarize_bam_contig_depths` command.
-  `-t` - Number of parallel processing threads to use.

**Input data:**

- sample-assembly_GLmetagenomics.fasta (assembly fasta file created in [Step 5a](#5a-rename-contig-headers))
- sample.bam (bam file created in [Step 9b](#9c-sort-assembly-alignments))

**Output data:**

- **sample-metabat-assembly-depth_GLmetagenomics.tsv** (tab-delimited summary of coverages)
- sample-bins/sample-bin\*.fasta (fasta files of recovered bins)
- **sample-bins_GLmetagenomics.zip** (zip file containing fasta files of recovered bins)

#### 14b. Bin quality assessment
> Utilizes the default `checkm` database available [checkm_data_2015_01_16.tar.gz](https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz).

```bash
checkm lineage_wf -f bins-overview_GLmetagenomics.tsv \
                  --tab_table \
                  -x fa \
                  ./ \
                  checkm-output-dir
```

**Parameter Definitions:**  

-  `lineage_wf` – Specifies the workflow being utilized.
-  `-f` – Specifies the output summary file name.
-  `--tab_table` – Specifies the output summary file should be a tab-delimited table.
-  `-x` – Specifies the extension that is on the bin fasta files that are being assessed.
-  `./` – Specifies the directory holding the bins, provided as a positional argument.
-  `checkm-output-dir` – Specifies the primary checkm output directory, provided as a positional argument.

**Input data:**

- sample-bins/sample-bin\*.fasta (bin fasta files generated in [Step 14a](#14a-bin-contigs))

**Output data:**

- **bins-overview_GLmetagenomics.tsv** (tab-delimited file with quality estimates per bin)
- checkm-output-dir/ (directory holding detailed checkm outputs)

#### 14c. Filter MAGs

```bash
cat <( head -n 1 bins-overview_GLmetagenomics.tsv ) \
    <( awk -F $'\t' ' $12 >= 90 && $13 <= 10 && $14 == 0 ' bins-overview_GLmetagenomics.tsv | sed 's/bin./MAG-/' ) \
    > checkm-MAGs-overview.tsv
    
# copying bins into a MAGs directory in order to run tax classification
awk -F $'\t' ' $12 >= 90 && $13 <= 10 && $14 == 0 ' bins-overview_GLmetagenomics.tsv | cut -f 1 > MAG-bin-IDs.tmp

mkdir MAGs
for ID in MAG-bin-IDs.tmp
do
    MAG_ID=$(echo $ID | sed 's/bin./MAG-/')
    cp ${ID}.fasta MAGs/${MAG_ID}.fasta
done

for SAMPLE in $(cat MAG-bin-IDs.tmp | sed 's/-bin.*//' | sort -u);
do
  mkdir ${SAMPLE}-MAGs
  mv ${SAMPLE}-*MAG*.fasta ${SAMPLE}-MAGs
  zip -r ${SAMPLE}-MAGs_GLmetagenomics.zip ${SAMPLE}-MAGs
done
```

**Input data:**

- bins-overview_GLmetagenomics.tsv (tab-delimited file with quality estimates per bin from [Step 14b](#14b-bin-quality-assessment))

**Output data:**

- checkm-MAGs-overview.tsv (tab-delimited file with quality estimates per MAG)
- MAGs/\*.fasta (directory holding high-quality MAGs)
- **\*-MAGs_GLmetagenomics.zip** (zip files containing directories of high-quality MAGs)

#### 14d. MAG taxonomic classification
> Uses default `gtdbtk` database setup with program's `download.sh` command.

```bash
gtdbtk classify_wf --genome_dir MAGs/ \
                   -x fa \
                   --out_dir gtdbtk-output-dir  \
                   --skip_ani_screen
```

**Parameter Definitions:**  

-  `classify_wf` – Specifies the workflow being utilized.
-  `--genome_dir` – Specifies the directory holding the MAGs to classify.
-  `-x` – Specifies the extension that is on the MAG fasta files that are being taxonomically classified.
-  `--out_dir` – Specifies the output directory name.
-  `--skip_ani_screen`  - Specifies to skip ani_screening step to classify genomes using mash and skani.

**Input data:**

- MAGs/\*.fasta (directory holding high-quality MAGs from [Step 14c](#14c-filter-mags))

**Output data:**

- gtdbtk-output-dir/gtdbtk.\*.summary.tsv (files with assigned taxonomy and info)

#### 14e. Generate Overview Table Of All MAGs

```bash
# combine summaries
for MAG in $(cut -f 1 assembly-summaries_GLmetagenomics.tsv | tail -n +2); do

    grep -w -m 1 "^${MAG}" checkm-MAGs-overview.tsv | cut -f 12,13,14 \
        >> checkm-estimates.tmp

    grep -w "^${MAG}" gtdbtk-output-dir/gtdbtk.*.summary.tsv | \
    cut -f 2 | sed 's/^.__//' | \
    sed 's/;.__/\t/g' | \
    awk 'BEGIN{ OFS=FS="\t" } { for (i=1; i<=NF; i++) if ( $i ~ /^ *$/ ) $i = "NA" }; 1' \
        >> gtdb-taxonomies.tmp

done

# Add headers
cat <(printf "est. completeness\test. redundancy\test. strain heterogeneity\n") checkm-estimates.tmp \
    > checkm-estimates-with-headers.tmp

cat <(printf "domain\tphylum\tclass\torder\tfamily\\tgenus\tspecies\n") gtdb-taxonomies.tmp \
    > gtdb-taxonomies-with-headers.tmp

paste assembly-summaries_GLmetagenomics.tsv \
checkm-estimates-with-headers.tmp \
gtdb-taxonomies-with-headers.tmp \
    > MAGs-overview.tmp

# Ordering by taxonomy
head -n 1 MAGs-overview.tmp > MAGs-overview-header.tmp

tail -n +2 MAGs-overview.tmp | sort -t \$'\t' -k 14,20 > MAGs-overview-sorted.tmp

cat MAGs-overview-header.tmp MAGs-overview-sorted.tmp \
    > MAGs-overview_GLmetagenomics.tsv
```

**Input data:**

- assembly-summaries_GLmetagenomics.tsv (table of assembly summary statistics from [Step 5b](#5b-summarize-assemblies))
- MAGs/\*.fasta (directory holding high-quality MAGs from [Step 14c](#14c-filter-mags))
- checkm-MAGs-overview.tsv (tab-delimited file with quality estimates per MAG from [Step 14c](#14c-filter-mags))
- gtdbtk-output-dir/gtdbtk.\*.summary.tsv (directory of files with assigned taxonomy and info from [Step 14d](#14d-mag-taxonomic-classification))

**Output data:**

- **MAGs-overview_GLmetagenomics.tsv** (a tab-delimited overview of all recovered MAGs)

<br>

---

### 15. Generate MAG-level functional summary overview

#### 15a. Get KO annotations per MAG
> This utilizes the helper script [`parse-MAG-annots.py`](https://github.com/nasa/GeneLab_Metagenomics_Workflow/blob/DEV/bin/parse-MAG-annots.py)

```bash
for file in $( ls MAGs/*.fasta )
do

    MAG_ID=$( echo ${file} | cut -f 2 -d "/" | sed 's/.fasta//' )
    sample_ID=$( echo ${MAG_ID} | sed 's/-MAG-[0-9]*$//' )

    grep "^>" ${file} | tr -d ">" > ${MAG_ID}-contigs.tmp

    python parse-MAG-annots.py -i annotations-and-taxonomy/${sample_ID}-gene-coverage-annotation-and-tax_GLmetagenomics.tsv \
                               -w ${MAG_ID}-contigs.tmp -M ${MAG_ID} \
                               -o MAG-level-KO-annotations_GLmetagenomics.tsv

    rm ${MAG_ID}-contigs.tmp

done
```

**Parameter Definitions:**  

- `-i` – Specifies the input sample TSV file containing sample coverage, annotation, and taxonomy info.
- `-w` – Specifies the appropriate temporary file holding all the contigs in the current MAG.
- `-M` – Specifies the current MAG unique identifier.
- `-o` – Specifies the output file name.

**Input data:**

- \*-gene-coverage-annotation-and-tax_GLmetagenomics.tsv (tables with combined gene coverage, annotation, and taxonomy info generated for individual samples, output from [Step 11](#11-combine-gene-level-coverage-taxonomy-and-functional-annotations-for-each-sample))
- MAGs/\*.fasta (directory holding high-quality MAGs from [Step 14c](#14c-filter-mags))

**Output data:**

- **MAG-level-KO-annotations_GLmetagenomics.tsv** (tab-delimited table holding MAGs and their KO annotations)

#### 15b. Summarize KO annotations with KEGG-Decoder

```bash
KEGG-decoder -v interactive \
             -i MAG-level-KO-annotations_GLmetagenomics.tsv \
             -o MAG-KEGG-Decoder-out_GLmetagenomics.tsv
```

**Parameter Definitions:**  

- `-v interactive` – Specifies to create an interactive html output.
- `-i` – Specifies the input tab-delimited table holding MAGs and their KO annotations.
- `-o` – Specifies the output table.

**Input data:**

- MAG-level-KO-annotations_GLmetagenomics.tsv (tab-delimited table holding MAGs and their KO annotations, output from [Step 15a](#15a-get-ko-annotations-per-mag))

**Output data:**

- **MAG-KEGG-Decoder-out_GLmetagenomics.tsv** (tab-delimited table holding MAGs and their proportions of genes held known to be required for specific pathways/metabolisms)
- **MAG-KEGG-Decoder-out_GLmetagenomics.html** (interactive heatmap html file of the above output table)

<br>

### 16. Filtering and Visualization of Contig- and Gene-taxonomy and Gene-function Outputs

#### 16a. Gene-level Taxonomy Heatmaps

```R
assembly_table <- "Combined-gene-level-taxonomy-coverages-CPM_GLmetagenomics.tsv"
assembly_summary <- "assembly-summaries_GLmetagenomics.tsv"
metadata_table <- "/path/to/sample/metadata"

# Read in assembly summary table
overview_table <- read_delim(assembly_summary, comment="#") %>%
  select(
    where(~all(!is.na(.)))
  )

col_names <- names(overview_table) %>% str_remove_all("-assembly")
sample_order <- col_names[-1] %>% sort()

# deduplicate rows by summing together species values
df <- read_delim(assembly_table, comment = "#")
sample_order <- get_samples(df, sample_order)

table2write <- read_taxonomy_table(df, sample_order) %>%
               select(species, !!sample_order) %>%
               group_by(species) %>%
               summarise(across(everything(), sum)) %>%
               filter(species != "Unclassified;_;_;_;_;_;_") %>%
               as.data.frame()

# Write out gene taxonomy table
write_tsv(x = table2write, file = "Combined-gene-level-taxonomy_unfiltered_GLmetagenomics.tsv")

make_heatmap(metadata_table_file = metadata_table, 
             feature_table_file = "Combined-gene-level-taxonomy_unfiltered_GLmetagenomics.tsv", 
             samples_column="sample_id", group_column = "group", 
             output_prefix = "Combined-gene-level-taxonomy_unfiltered", 
             assay_suffix = "_GLmetagenomics", 
             custom_palette = custom_palette)
```

**Custom Functions Used:**
- [get_samples()](#get_samples)
- [read_taxonomy_table()](#read_taxonomy_table)
- [make_heatmap()](#make_heatmap)

**Input data:**
- assembly-summaries_GLmetagenomics.tsv (table of assembly summary statistics, output from [Step 5b](#5b-summarize-assemblies))
- Combined-gene-level-taxonomy-coverages-CPM_GLmetagenomics.tsv (table with all samples combined based on gene-level 
  taxonomic classifications, output from [Step 13a](#13a-generate-gene-level-coverage-summary-tables)) 

**Output data:**
- Combined-gene-level-taxonomy_unfiltered_GLmetagenomics.tsv (aggregated gene-level taxonomy table with samples in columns and species in rows)
- **Combined-gene-level-taxonomy_unfiltered_heatmap_GLmetagenomics.png** (heatmap of all gene-level taxonomy assignments, output from [make_heatmap()](#make_heatmap))
- **Combined-gene-level-taxonomy_unfiltered_top_50_heatmap_GLmetagenomics.png** (heatmap of the top 50 gene-level taxonomy assignments, output from [make_heatmap()](#make_heatmap))

#### 16b. Gene-level Taxonomy Feature Filtering

```R
feature_table_file <- "Combined-gene-level-taxonomy_unfiltered_GLmetagenomics.tsv"
metadata_table <- "/path/to/sample/metadata"
threshold <- 1000

# read in feature table
feature_table <- read_delim(feature_table_file) %>%
                 mutate(across(where(is.numeric), function(col) replace_na(col, 0))) %>%
                 as.data.frame()
feature_name <- colnames(feature_table)[1]
rownames(feature_table) <- feature_table[,1]
feature_table <- feature_table[, -1]

table2write <- get_abundant_features(feature_table, cpm_threshold=threshold) %>%
               as.data.frame() %>%
               rownames_to_column(feature_name)

write_tsv(x = table2write, file = "Combined-gene-level-taxonomy_filtered_GLmetagenomics.tsv")

make_heatmap(metadata_table_file = metadata_table, 
             feature_table_file = "Combined-gene-level-taxonomy_filtered_GLmetagenomics.tsv", 
             samples_column="sample_id", group_column = "group", 
             output_prefix = "Combined-gene-level-taxonomy_filtered", 
             assay_suffix = "_GLmetagenomics", 
             custom_palette = custom_palette)

```

**Custom Functions Used:**
- [get_abundant_features()](#get_abundant_features)
- [make_heatmap()](#make_heatmap)

**Parameter Definitions:**

- `feature_table_file` - path to a tab separated samples feature table containing gene-level coverage data 
                         species/functions as the first column and samples as other columns.
- `metadata_table` - path to a file with samples as rows and columns describing each sample
- `threshold` - threshold to identify abundant features, default: 1000

**Input Data:**

- `Combined-gene-level-taxonomy_unfiltered_GLmetagenomics.tsv`(aggregated gene taxonomy table with samples in columns and species in rows, from [Step 16a](#16a-gene-level-taxonomy-heatmaps))
- `/path/to/sample/metadata` (a file containing sample-wise metadata, mapping sample names to group metadata)

**Output Data:**

- **Combined-gene-level-taxonomy_filtered_GLmetagenomics.tsv** (filtered gene-level taxonomy, output from [get_abundant_features()](#get_abundant_features))
- **Combined-gene-level-taxonomy_filtered_heatmap_GLmetagenomics.png** (heatmap of all gene-level taxonomy assignments after filtering out non-abundant features, output from [make_heatmap()](#make_heatmap))
- **Combined-gene-level-taxonomy_filtered_top_50_heatmap_GLmetagenomics.png** (heatmap of the top 50 gene taxonomy assignments after filtering out non-abundant features, output from [make_heatmap()](#make_heatmap))

#### 16c. Gene-level KO Functions Heatmaps

```R
assembly_table <- "Combined-gene-level-KO-function-coverages-CPM_GLmetagenomics.tsv"
assembly_summary <- "assembly-summaries_GLmetagenomics.tsv"
metadata_table <- "/path/to/sample/metadata"

# Read in assembly summary table and remove columns where the values are NA
overview_table <- read_delim(assembly_summary, comment="#") %>%
  select(
    where(~all(!is.na(.)))
  )

col_names <- names(overview_table) %>% str_remove_all("-assembly")
sample_order <- col_names[-1] %>% sort()

# deduplicate rows by summing together species values
df <- read_delim(assembly_table, comment = "#")
sample_order <- get_samples(df, sample_order, end_col="KO_function")

table2write <- df %>%
               select(KO_ID, !!sample_order)

# Write out gene taxonomy table
write_tsv(x = table2write, file = "Combined-gene-level-KO_unfiltered_GLmetagenomics.tsv")

make_heatmap(metadata_table_file = metadata_table, 
             feature_table_file = "Combined-gene-level-KO_unfiltered_GLmetagenomics.tsv",
             samples_column="sample_id", group_column = "group", 
             output_prefix = "Combined-gene-level-KO-function_unfiltered", 
             assay_suffix = "_GLmetagenomics", 
             custom_palette = custom_palette)

```

**Custom Functions Used:**
- [get_samples()](#get_samples)
- [make_heatmap()](#make_heatmap)

**Parameter Definitions:**

- `metadata_table` - path to a file with samples as rows and columns describing each sample
- `assembly_table` - path to a tab-separated table containing gene-level KO function coverage data with
                         species/functions as the first column and samples as other columns.
- `assembly_summary` - path to a tab-separated file containing statistics on assemblies created for each sample

**Input data:**

- assembly-summaries_GLmetagenomics.tsv (table of assembly summary statistics, output from [Step 5b](#5b-summarize-assemblies))
- Combined-gene-level-KO-function-coverages-CPM_GLmetagenomics.tsv (table with all samples combined based on KO annotations; 
  normalized to coverage per million genes covered, output from [Step 13a](#13a-generate-gene-level-coverage-summary-tables))
- `/path/to/sample/metadata` (a file containing sample-wise metadata, mapping sample names to group metadata)

**Output data:**

- Combined-gene-level-KO-function_unfiltered_GLmetagenomics.tsv (aggregated and subsetted gene-level KO function table)
- **Combined-gene-level-KO-function_unfiltered_heatmap_GLmetagenomics.png** (heatmap of all gene-level KO function assignments, output from [make_heatmap()](#make_heatmap))
- **Combined-gene-level-KO-function_unfiltered_top_50_heatmap_GLmetagenomics.png** (heatmap of the top 50 gene-level KO function assignments, output from [make_heatmap()](#make_heatmap))

#### 16d. Gene-level KO Functions Feature Filtering

```R
feature_table_file <- "Combined-gene-level-KO-function_unfiltered_GLmetagenomics.tsv"
metadata_table <- "/path/to/sample/metadata"
threshold <- 1000

# read in feature table
feature_table <- read_delim(feature_table_file) %>%
                 mutate(across(where(is.numeric), function(col) replace_na(col, 0))) %>%
                 as.data.frame()
feature_name <- colnames(feature_table)[1]
rownames(feature_table) <- feature_table[,1]
feature_table <- feature_table[, -1]

table2write <- get_abundant_features(feature_table, cpm_threshold=threshold) %>%
               as.data.frame() %>%
               rownames_to_column(feature_name)

write_tsv(x = table2write, file = "Combined-gene-level-KO_filtered_GLmetagenomics.tsv")

make_heatmap(metadata_table_file = metadata_table, 
             feature_table_file = "Combined-gene-level-KO_filtered_GLmetagenomics.tsv", 
             samples_column="sample_id", group_column = "group", 
             output_prefix = "Combined-gene-level-KO_filtered", 
             assay_suffix = "_GLmetagenomics", 
             custom_palette = custom_palette)

```

**Custom Functions Used:**
- [get_abundant_features()](#get_abundant_features)
- [make_heatmap()](#make_heatmap)

**Parameter Definitions:**

- `feature_table_file` - path to a tab separated samples feature table containing gene-level coverage data 
                         species/functions as the first column and samples as other columns.
- `metadata_table` - path to a file with samples as rows and columns describing each sample
- `threshold` - threshold to identify abundant features, default: 1000

**Input Data:**

- `Combined-gene-level-KO-function_unfiltered_GLmetagenomics.tsv`(aggregated gene taxonomy table with samples in columns and species in rows, from [Step 16c](#16c-gene-level-ko-functions-heatmaps))
- `/path/to/sample/metadata` (a file containing sample-wise metadata, mapping sample names to group metadata)

**Output Data:**

- **Combined-gene-level-KO-function_filtered_GLmetagenomics.tsv** (filtered gene-level KO function table, output from [get_abundant_features()](#get_abundant_features))
- **Combined-gene-level-KO-function_filtered_heatmap_GLmetagenomics.png** (heatmap of all gene-level KO function assignments after filtering out non-abundant features, output from [make_heatmap()](#make_heatmap))
- **Combined-gene-level-KO-function_filtered_top_50_heatmap_GLmetagenomics.png** (heatmap of the top 50 gene-level KO function assignments after filtering out non-abundant features, output from [make_heatmap()](#make_heatmap))

#### 16e. Contig-level Heatmaps

```R
assembly_table <- "Combined-contig-level-taxonomy-coverages-CPM_GLmetagenomics.tsv"
assembly_summary <- "assembly-summaries_GLmetagenomics.tsv"
metadata_table <- "/path/to/sample/metadata"

# Read in assembly summary table
overview_table <- read_delim(assembly_summary, comment="#") %>%
  select(
    where(~all(!is.na(.)))
  )

col_names <- names(overview_table) %>% str_remove_all("-assembly")
sample_order <- col_names[-1] %>% sort()

# deduplicate rows by summing together species values
df <- read_delim(assembly_table, comment = "#")
sample_order <- get_samples(df, sample_order)

table2write <- read_taxonomy_table(df, sample_order) %>%
               select(species, !!sample_order) %>%
               group_by(species) %>%
               summarise(across(everything(), sum)) %>%
               filter(species != "Unclassified;_;_;_;_;_;_") %>%
               as.data.frame()

# Write out contig taxonomy table
write_tsv(x = table2write, file = "Combined-contig-level-taxonomy_unfiltered_GLmetagenomics.tsv")

make_heatmap(metadata_table_file = metadata_table, 
             feature_table_file = "Combined-contig-level-taxonomy_unfiltered_GLmetagenomics.tsv", 
             samples_column="sample_id", group_column = "group", 
             output_prefix = "Combined-contig-level-taxonomy", 
             assay_suffix = "_GLmetagenomics", 
             custom_palette = custom_palette)
```

**Custom Functions Used:**
- [get_samples()](#get_samples)
- [read_taxonomy_table()](#read_taxonomy_table)
- [make_heatmap()](#make_heatmap)

**Parameter Definitions:**

- `metadata_table` - path to a file with samples as rows and columns describing each sample
- `assembly_table` - path to a tab-separated table containing gene-level KO function coverage data with
                         species/functions as the first column and samples as other columns.
- `assembly_summary` - path to a tab-separated file containing statistics on assemblies created for each sample

**Input data:**

- assembly-summaries_GLmetagenomics.tsv (table of assembly summary statistics, output from [Step 5b](#5b-summarize-assemblies))
- Combined-contig-level-taxonomy-coverages-CPM_GLmetagenomics.tsv (table with all samples combined based on contig-level 
  taxonomic classifications, output from [Step 13b](#13b-generate-contig-level-coverage-summary-tables)) 

**Output data:**

- Combined-contig-level-taxonomy_unfiltered_GLmetagenomics.tsv (aggregated contig-level taxonomy table with samples in columns and species in rows)
- **Combined-contig-level-taxonomy_unfiltered_heatmap_GLmetagenomics.png** (heatmap of all contig-level taxonomy assignments, output from [make_heatmap()](#make_heatmap))
- **Combined-contig-level-taxonomy_unfiltered_top_50_heatmap_GLmetagenomics.png** (heatmap of the top 50 contig-level taxonomy assignments, output from [make_heatmap()](#make_heatmap))

#### 16f. Contig-level Feature Filtering

```R
feature_table_file <- "Combined-contig-level-taxonomy_GLmetagenomics.tsv"
metadata_table <- "/path/to/sample/metadata"
threshold <- 1000

# read in feature table
feature_table <- read_delim(feature_table_file) %>%
                 mutate(across(where(is.numeric), function(col) replace_na(col, 0))) %>%
                 as.data.frame()
feature_name <- colnames(feature_table)[1]
rownames(feature_table) <- feature_table[,1]
feature_table <- feature_table[, -1]

table2write <- get_abundant_features(feature_table, cpm_threshold=threshold) %>%
               as.data.frame() %>%
               rownames_to_column(feature_name)

write_tsv(x = table2write, file = "Combined-contig-level-taxonomy_filtered_GLmetagenomics.tsv")

make_heatmap(metadata_table_file = metadata_table, 
             feature_table_file = "Combined-contig-level-taxonomy_filtered_GLmetagenomics.tsv", 
             samples_column="sample_id", group_column = "group", 
             output_prefix = "Combined-contig-level-taxonomy_filtered", 
             assay_suffix = "_GLmetagenomics", 
             custom_palette = custom_palette)
```

**Custom Functions Used:**
- [get_abundant_features()](#get_abundant_features)
- [make_heatmap()](#make_heatmap)

**Parameter Definitions:**

- `feature_table_file` - path to a tab separated samples feature table containing gene-level coverage data 
                         species/functions as the first column and samples as other columns.
- `metadata_table` - path to a file with samples as rows and columns describing each sample
- `threshold` - threshold to identify abundant features, default: 1000

**Input Data:**

- `Combined-contig-level-taxonomy_unfiltered_GLmetagenomics.tsv`(aggregated gene taxonomy table with samples in columns and species in rows, from [Step 16c](#16c-gene-level-ko-functions-heatmaps))
- `/path/to/sample/metadata` (a file containing sample-wise metadata, mapping sample names to group metadata)

**Output Data:**

- **Combined-contig-level-taxonomy_filtered_GLmetagenomics.tsv** (filtered contig-level taxonomy, output from [get_abundant_features()](#get_abundant_features))
- **Combined-contig-level-taxonomy_filtered_heatmap_GLmetagenomics.png** (heatmap of all contig-level taxonomy assignments after filtering out non-abundant features, output from [make_heatmap()](#make_heatmap))
- **Combined-contig-level-taxonomy_filtered_top_50_heatmap_GLmetagenomics.png** (heatmap of the top 50 contig-level taxonomy assignments after filtering out non-abundant features, output from [make_heatmap()](#make_heatmap))

### 17. Generate Assembly-based Processing Overview
> This utilizes the helper script [`generate-assembly-based-overview-table.sh`](https://github.com/nasa/GeneLab_Metagenomics_Workflow/blob/DEV/bin/generate-assembly-based-overview-table.sh) 

```bash
bash generate-assembly-based-overview-table.sh sample_ids_file.txt \
  assemblies/ predicted-genes/ read-mapping/ bins/ MAGs/ \
  Assembly-based-processing-overview_GLmetagenomics.tsv
```

**Parameter Definitions:**

- `sample_ids_file.txt` - A file listing the sample names, one on each row, provided as a positional argument.
- `assemblies/` - The directory holding the contig-renamed assembly files generated in [Step 5a](#5a-rename-contig-headers), provided as a positional argument.
- `predicted-genes/` - The directory holding the gene-calls ammino-acid fasta files generated in [Step 6a](#6a-generate-gene-predictions) and [Step 6b](#6b-remove-line-wraps-in-gene-prediction-output), provided as a positional argument.
- `read-mapping/` - The directory holding the sorted mapping to the sample assembly in BAM format generated in [Step 9c](#9c-sort-assembly-alignments), provided as a positional argument.
- `bins/` - The directory holding the recovered bins fasta files generated in [Step 14a](#14a-bin-contigs), provided as a positional argument.
- `MAGs/` - The directory holding the high-quality MAGs fasta files generated in [Step 14c](#14c-filter-mags), provided as a positional argument.
- `Assembly-based-processing-overview_GLmetagenomics.tsv` - name of the output file, provided as a positional argument.

**Input Data:**

- assemblies/\*.fasta (contig-renamed assembly files from [Step 5a](#5a-rename-contig-headers))
- predicted-genes/\*.faa (gene-calls amino-acid fasta file with line wraps removed, output from [Step 6b](#6b-remove-line-wraps-in-gene-prediction-output))
- read-mapping/\*.bam (sorted mapping to sample assembly, in BAM format, output from [Step 9c](#9c-sort-assembly-alignments))
- bins/\*.fasta (fasta files of recovered bins, output from [Step 14a](#14a-bin-contigs))
- MAGs/\*.fasta (directory holding high-quality MAGs, output from [Step 14c](#14c-filter-mags))

**Output Data:**

- **Assembly-based-processing-overview_GLmetagenomics.tsv** (Tab delimited text file providing a summary of assembly-based processing results for each sample)

<br>

---

## Read-based Processing


### 18. Taxonomic Profiling Using Kaiju

#### 18a. Build Kaiju Database

```bash
# Make a directory that will hold the downloaded kaiju database
mkdir kaiju-db/

# Download kaiju's reference database
kaiju-makedb -s kaiju_db/nr_euk -t NumberOfThreads

# Clean up
rm nr_euk/kaiju_db_nr_euk.bwt nr_euk/kaiju_db_nr_euk.sa
```

**Parameter Definitions:**

- `-s nr_euk` - Specifies to download the subset of the NCBI BLAST nr (non-redundant) database containing all proteins belonging to Archaea, bacteria, and viruses, and additionally include proteins from fungi and microbial eukaryotes.
- `-t` - Number of parallel processing threads to use.

**Input Data:**

*No input data required*

**Output Data:**

- kaiju-db/nr_euk/kaiju_db_nr_euk.fmi (FM-index file containing the main Kaiju database index)
- kaiju-db/nr_euk/kaiju_db_nr_euk.faa (FASTA amino acid file containing the protein sequences used to build the .fmi index file)
- kaiju-db/nodes.dmp (taxonomy hierarchy file from the NCBI Taxonomy database defining the parent-child relationships in the taxonomic tree)
- kaiju-db/names.dmp (taxonomy names file from the NCBI Taxonomy database that maps taxonomic IDs to their scientific names)
- kaiju-db/merged.dmp (merged taxonomy IDs file from the NCBI Taxonomy database that maps deprecated taxonomic IDs to current ones)

#### 18b. Kaiju Taxonomic Classification

```bash
kaiju -f kaiju-db/nr_euk/kaiju_db_nr_euk.fmi \
      -t kaiju-db/nodes.dmp \
      -z NumberOfThreads \
      -E 1e-05 \
      -i /path/to/sample_R1_filtered_GLmetagenomics.fastq.gz \
      -j /path/to/sample_R2_filtered_GLmetagenomics.fastq.gz \
      -o sample_kaiju.out
```

**Parameter Definitions:**

- `-f` - Specifies the path to the kaiju database index file (.fmi).
- `-t` - Specifies the path to the kaiju taxonomy hierarchy file (nodes.dmp).
- `-z` - Number of parallel processing threads to use.
- `-E` - Specifies the minimum E-value to use for filter matches (an E-value of 1e-05 means that there's a 0.001% chance that the matches identified occurred randomly).
- `-i` - Specifies path to the forward read input file.
- `-i` - Specifies path to the reverse read input file.
- `-o` - Specifies the name of the output file.

**Input Data:**

- kaiju-db/nr_euk/kaiju_db_nr_euk.fmi (FM-index file containing the main Kaiju database index, output from [Step 18a](#18a-build-kaiju-database))
- kaiju-db/nodes.dmp (kaiju taxonomy hierarchy nodes file, output from [Step 18a](#18a-build-kaiju-database))
- *_R[12]_filtered_GLmetagenomics.fastq.gz (filtered/trimmed reads from [Step 2b](#2b-trim-polyg) above)

**Output Data:**

- sample_kaiju.out (kaiju output file)

#### 18c. Compile Kaiju Taxonomy Results

```bash
# Merge kaiju reports to one table at the species level 
kaiju2table -t nodes.dmp \
            -n names.dmp \
            -p \
            -r "species" \
            -o merged_kaiju_summary_${TAXON_LEVEL}.tsv \
            *_kaiju.out

# Convert file names to sample names
sed -i -E 's/.+\/(.+)_kaiju\.out/\1/g' merged_kaiju_table.tsv && \
sed -i -E 's/file/sample/' merged_kaiju_table.tsv
```

**Parameter Definitions:**

- `-t` - Specifies the path to the kaiju taxonomy hierarchy file (nodes.dmp).
- `-n` - Specifies the path to the kaiju taxonomy names file (names.dmp).
- `-p` - Print the full taxon path instead of only the taxon name.
- `-r` - Specifies taxonomic rank to print the taxon path to, must be one of: phylum, class, order, family, genus, species. (Default: species).
- `-o` - Specifies the name of the kaiju taxon summary output file.
- `*_kaiju.out` - Positional argument specifying the path to the kaiju output files for each sample. 

**Input Data:**

- kaiju-db/nodes.dmp (kaiju taxonomy hierarchy nodes file, output from [Step 18a](#18a-build-kaiju-database))
- kaiju-db/names.dmp (kaiju taxonomy names file, output from [Step 18a](#18a-build-kaiju-database))
- *kaiju.out (kaiju output files, output from [Step 18b](#18b-kaiju-taxonomic-classification))

**Output Data:**

- merged_kaiju_table.tsv (compiled kaiju summary table at the species level)

#### 18d. Convert Kaiju Output To Krona Format

```bash
kaiju2krona -u \
            -n kaiju-db/names.dmp \
            -t kaiju-db/nodes.dmp \
            -i sample_kaiju.out \
            -o sample.krona
```

**Parameter Definitions:**

- `-u` - Include count for unclassified reads in output.
- `-n` - Specifies the path to the kaiju taxonomy names file (names.dmp).
- `-t` - Specifies the path to the kaiju taxonomy hierarchy file (nodes.dmp).
- `-i` - Specifies the path to the kaiju output file.
- `-o` - Specifies the name of krona formatted kaiju output file.

**Input Data:**
- kaiju-db/names.dmp (kaiju taxonomy names file, output from [Step 18a](#18a-build-kaiju-database))
- kaiju-db/nodes.dmp (kaiju taxonomy hierarchy nodes file, output from [Step 18a](#18a-build-kaiju-database))
- sample_kaiju.out (kaiju output file, output from [Step 18b](#18b-kaiju-taxonomic-classification))

**Output Data:**

- sample.krona (krona formatted kaiju output)

#### 18e. Compile Kaiju Krona Reports

```bash
# Create a file containing a sorted list of all .krona files 
find . -type f -name "*.krona" | sort -uV > krona_files.txt

# Create a file containing a sorted list of all sample names
FILES=($(find . -type f -name "*.krona"))
basename -a -s '.krona' ${FILES[*]} | sort -uV  > sample_names.txt

# Create ktImportText input format files
KTEXT_FILES=($(paste -d',' "krona_files.txt" "sample_names.txt"))

# Create html containing krona plot  
ktImportText  -o kaiju-report.html ${KTEXT_FILES[*]}
```

**Parameter Definitions:**

*find*
- `-type f` -  Specifies that the type of file to find is a regular file.
- `-name "*.krona"` - Specifies to find files ending with the .krona suffix.  

*sort*
- `-u` - Specifies to perform a unique sort.
- `-V` - Specifies to perform a mixed type of sorting with names containing numbers within text.
- `> krona_files.txt` - Redirects the sorted list to a separate text file.

*basename*
- `-a` - Support multiple arguments and treat each as a file name.
- `-s '.krona'` - Remove trailing '.krona' suffix.

*paste*
- `-d','` - Paste both krona and sample files together line by line delimited by comma ','.

*ktImportText*
- `-o` - Specifies the compiled output html file name.
- `${KTEXT_FILES[*]}` - An array positional argument with the following content: 
                        sample_1.krona,sample_1 sample_2.krona,sample_2 ... sample_n.krona,sample_n.

**Input Data:**

- *.krona (all sample .krona formatted files, output from [Step 18d](#18d-convert-kaiju-output-to-krona-format)) 
             
**Output Data:**

- krona_files.txt (sorted list of all *.krona files)
- sample_names.txt (sorted list of all sample names)
- **kaiju-report_GLmetagenomics.html** (compiled krona html report containing all samples)

#### 18f. Create Kaiju Species Count Table

```R
feature_table <- process_kaiju_table(file_path="merged_kaiju_table_GLmetagenomics.tsv")
table2write <- feature_table  %>%
               as.data.frame() %>%
               rownames_to_column("Species")
write_tsv(x = table2write, file = "kaiju_species_table_GLmetagenomics.tsv")
```

**Custom Functions Used:**
- [process_kaiju_table()](#process_kaiju_table)

**Parameter Definitions:**

- `file_path` - path to compiled kaiju table at the species taxon level
- `x`  - feature table dataframe to write to file
- `file` - path to where to write kaiju count table per sample

**Input Data:**

- merged_kaiju_table_GLmetagenomics.tsv (compiled kaiju table at the species taxon level, from [Step 18c](#18c-compile-kaiju-taxonomy-results))

**Output Data:**

- **kaiju_species_table_GLmetagenomics.tsv** (kaiju species count table in tsv format)

#### 18g. Filter Kaiju Species Count Table

```R
feature_table_file <- "kaiju_species_table_GLmetagenomics.tsv"
output_file <- "kaiju_filtered_species_table_GLmetagenomics.tsv"
threshold <- 0.5

# string used to define non-microbial taxa
non_microbial <- "UNCLASSIFIED|Unclassified|unclassified|Homo sapien|cannot|uncultured|unidentified"

# read in feature table
feature_table <- read_delim(feature_table_file) %>%
                 mutate(across(where(is.numeric), function(col) replace_na(col, 0))) %>%
                 as.data.frame()
feature_name <- colnames(feature_table)[1]
rownames(feature_table) <- feature_table[,1]
feature_table <- feature_table[, -1]

# convert count table to a relative abundance matrix
abund_table <- feature_table %>% rownames_to_column(feature_name) %>%
  mutate(across(where(is.numeric), function(x) (x / sum(x, na.rm = TRUE)) * 100)) %>%
  as.data.frame

rownames(abund_table) <- abund_table[,1]
abund_table <- abund_table[,-1] %>% t 

table2write <- group_low_abund_taxa(abund_table, threshold = threshold)  %>%
  t %>% as.data.frame %>%
  rownames_to_column(feature_name)

write_tsv(x = table2write, file = output_file)
```

**Custom Functions Used:**
- [group_low_abund_taxa()](#group_low_abund_taxa)

**Parameter Definitions:**

- `threshold` - threshold for filtering out rare taxa, a percentage between 0 and 100.
- `output_file` - output filename
- `input_file` - input filename

**Input Data:**

- kaiju_species_table_GLmetagenomics.tsv (path to kaiju species table from [Step 18f](#18f-create-kaiju-species-count-table))

**Output Data:**

- **kaiju_filtered_species_table_GLmetagenomics.tsv** (a file containing the filtered species table)

---

#### 18h. Kaiju Taxonomy Barplots

```R
species_table_file <- "kaiju_species_table_GLmetagenomics.tsv"
filtered_species_table_file <- "kaiju_filtered_species_table_GLmetagenomics.tsv"
metadata_file <- "/path/to/sample/metadata"

make_barplot(metadata_file = metadata_file, feature_table_file = species_table_file, 
             feature_column = "Species", samples_column = "sample_id", group_column = "group",
             output_prefix = "kaiju_unfiltered_species", assay_suffix = "_GLmetagenomics",
             publication_format = publication_format, custom_palette = custom_palette)

# Save static unfiltered plot
make_barplot(metadata_file = metadata_file, feature_table_file = filtered_species_table_file, 
             feature_column = "Species", samples_column = "sample_id", group_column = "group",
             output_prefix = "kaiju_filtered_species", assay_suffix = "_GLmetagenomics",
             publication_format = publication_format, custom_palette = custom_palette)
```

**Custom Functions Used:**
- [make_barplot](#make_barplot)

**Parameter Definitions:**

- `species_table_file` - a file containing the species count table
- `filtered_species_table_file` - a file containing the filtered species count table
- `metadata_file` - a file containing group information for each sample in the species count files

**Input Data:**

- `kaiju_species_table_GLmetagenomics.tsv` (a file containing the species count table, output from [Step 18f](#18f-create-kaiju-species-count-table))
- `kaiju_filtered_species_table_GLmetagenomics.tsv` (a file containing the filtered species count table, output from [Step 18g](#18g-filter-kaiju-species-count-table))
- `/path/to/sample/metadata` (a file containing sample-wise metadata, mapping sample names to group metadata)

**Output Data:**

- kaiju_unfiltered_species_barplot_GLmetagenomics.png (taxonomy barplot without filtering)
- **kaiju_unfiltered_species_barplot_GLmetagenomics.html** (interactive taxonomy barplot without filtering)
- kaiju_filtered_species_barplot_GLmetagenomics.png (taxonomy barplot after filtering rare and non-microbial taxa)
- **kaiju_filtered_species_barplot_GLmetagenomics.html** (interactive taxonomy barplot after filtering rare and non-microbial taxa)

<br>

---

### 19. Taxonomic Profiling Using Kraken2

#### 19a. Download Kraken2 Database

```bash 
## Download all microbial (including eukaryotes) - https://benlangmead.github.io/aws-indexes/k2

# Downloading and building kraken2's pluspfp database which contains the standard database (Refseq archaea, bacteria, viral, plasmid, human1, UniVec_Core) + plants + protists + fungi

mkdir kraken2-db/ && cd kraken2-db/

# Inspect file
INSPECT_URL=https://genome-idx.s3.amazonaws.com/kraken/pluspfp_20250714/inspect.txt
wget ${INSPECT_URL}

# Library report
LIBRARY_REPORT_URL=https://genome-idx.s3.amazonaws.com/kraken/pluspfp_20250714/library_report.tsv
wget ${LIBRARY_REPORT_URL}

# Md5sums
MD5_URL=https://genome-idx.s3.amazonaws.com/kraken/pluspfp_20250714/pluspfp.md5 
wget ${MD5_URL}

# Download and unzip the main database files
DB_URL=https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_20250714.tar.gz 
wget -O k2_pluspfp.tar.gz --timeout=3600 --tries=0 --continue ${DB_URL} && \
tar -xvzf k2_pluspfp.tar.gz
```

**Parameter Definitions:**

*wget*
- `O` - Name of file to download the url content to.
- `--timeout=3600` - Specifies the network timeout in seconds.
- `--tries=0` - Retry download infinitely.
- `--continue` -  Continue getting a partially-downloaded file.
- `*_URL` - Position argument specifying the url to download a particular resource from.

*tar*
- `-xvzf` - unpack the specified *tar.gz archive in verbose mode

**Input Data:**

- `INSPECT_URL=` (url specifying the location of kraken2 inspect file)
- `LIBRARY_REPORT_URL=` (url specifying the location of kraken2 library report file)
- `MD5_URL=` (url specifying the location of the md5 file of the kraken database)
- `DB_URL=` (url specifying the location of the main kraken database archive in .tar.gz format)

**Output Data:**

- kraken2-db/  (a directory containing kraken2 database files)

#### 19b. Kraken2 Taxonomic Classification

```bash
kraken2 --db kraken2-db/ \
        --gzip-compressed \
        --threads NumberOfThreads \
        --use-names \
        --output sample-kraken2-output.txt \
        --report sample-kraken2-report.tsv \
        /path/to/sample_R1_filtered_GLmetagenomics.fastq.gz /path/to/sample_R2_filtered_GLmetagenomics.fastq.gz
```

**Parameter Definitions:**

- `--db` - Specifies the directory holding the kraken2 database files. 
- `--gzip-compressed` - Specifies the input files are gzip-compressed.
- `--threads` - Number of parallel processing threads to use.
- `--use-names` - Specifies to add taxa names in addition to taxids.
- `--output` - Specifies the name of the kraken2 read-based output file.
- `--report` - Specifies the name of the kraken2 report output file.
- `sample_R1_filtered_GLmetagenomics.fastq.gz` - Positional argument specifying the forward read input file.
- `sample_R2_filtered_GLmetagenomics.fastq.gz` - Positional argument specifying the reverse read input file.

**Input Data:**

- kraken2-db/ (a directory containing kraken2 database files, output from [Step 19a](#19a-download-kraken2-database))
- *_R[12]_filtered_GLmetagenomics.fastq.gz (filtered/trimmed reads from [Step 2b](#2b-trim-polyg) above)

**Output Data:**

- sample-kraken2-output.txt (kraken2 read-based output file (one line per read))
- sample-kraken2-report.tsv (kraken2 report output file (one line per taxa, with number of reads assigned to it))

#### 19c. Compile Kraken2 Taxonomy Results

##### 19ci. Create Merged Kraken2 Taxonomy Table

```R
species_table <- merge_kraken_reports(reports-dir='/path/to/kraken2/reports')
write_tsv(x = species_table, file = "kraken2_species_table_GLmetagenomics.tsv")
```

**Custom Functions Used:**
- [merge_kraken_reports()](#merge_kraken_reports)

**Parameter Definitions:**

- `reports-dir` - path to compiled kraken reports
- `x`  - feature table dataframe to write to file
- `file` - path to where to write kraken2 species table table

**Input Data:**

- \*-kraken2-report.tsv (kraken report from each sample to compile, outputs from [Step 19b](#19b-kraken2-taxonomic-classification))

**Output Data:**

- **kraken2_species_table_GLmetagenomics.tsv** (kraken species count table in tsv format)

##### 19cii. Compile Kraken2 Taxonomy Reports

```bash
multiqc --zip-data-dir \ 
        --outdir kraken2_multiqc_report \
        --filename kraken2_multiqc_GLmetagenomics \
        --interactive \
        /path/to/*kraken2-report.tsv
```

**Parameter Definitions:**

- `--zip-data-dir` - Compress the data directory.
- `--outdir` - Specifies the output directory to store results.
- `--filename` - Specifies the filename prefix of results.
- `--interactive` - Force multiqc to always create interactive javascript plots.
- `/path/to/*kraken2-report.tsv` - The kraken2 output report files, provided as a positional argument.

**Input Data:**

- \*-kraken2-report.tsv (kraken report from each sample to compile, outputs from [Step 19b](#19b-kraken2-taxonomic-classification))

**Output Data:**

- **kraken2_multiqc_GLmetagenomics.html** (multiqc output html summary)
- **kraken2_multiqc_GLmetagenomics_data.zip** (zip archive containing multiqc output data)

#### 19d. Convert Kraken2 Output to Krona Format

```bash
kreport2krona.py --report-file sample-kraken2-report.tsv  \
                 --output sample.krona
```

**Parameter Definitions:**

- `--report-file` - Specifies the name of the input kraken2 report file.
- `--output` - Specifies the name of the krona output file.

**Input Data:**

- sample-kraken2-report.tsv (kraken report, output from [Step 19b](#19b-kraken2-taxonomic-classification))

**Output Data:**

- sample.krona (krona formatted kraken2 output)

#### 19e. Compile Kraken2 Krona Reports

```bash
# Find, list and write all .krona files to file 
find . -type f -name "*.krona" | sort -uV > krona_files.txt

FILES=($(find . -type f -name "*.krona"))
basename --multiple --suffix='.krona' ${FILES[*]} | sort -uV  > sample_names.txt

# Create ktImportText input format files
KTEXT_FILES=($(paste -d',' "krona_files.txt" "sample_names.txt"))

# Create html   
ktImportText -o kraken2-report_GLmetagenomics.html ${KTEXT_FILES[*]}
```

**Parameter Definitions:**

*find*
  - `-type f` -  Specifies that the type of file to find is a regular file.
  - `-name "*.krona"` - Specifies to find files ending with the .krona suffix. 

*sort*
  - `-u` - Specifies to perform a unique sort.
  - `-V` - Specifies to perform a mixed type of sorting with names containing numbers within text.
  - `> {}.txt` - Redirects the sorted list to a separate text file.

*basename*
  - `--multiple` - Support multiple arguments and treat each as a file name.
  - `--suffix='.krona'` - Remove a trailing '.krona' suffix.

*paste*
  - `-d','` - Paste both krona and sample files together line by line delimited by comma ','.

*ktImportText*
  - `-o` - Specifies the compiled output html file name.
  - `${KTEXT_FILES[*]}` - An array positional argument with the following content: sample_1.krona,sample_1 sample_2.krona,sample_2 .. sample_n.krona,sample_n.

**Input Data:**

- *.krona (all sample .krona formatted files, output from [Step 19d](#19d-convert-kraken2-output-to-krona-format)) 

                      
**Output Data:**

- krona_files.txt (sorted list of all *.krona files)
- sample_names.txt (sorted list of all sample names)
- **kraken2-report_GLmetagenomics.html** (compiled krona html report containing all samples)

#### 19f. Filter Kraken2 Species Count Table

```R
feature_table_file <- "kraken2_species_table_GLmetagenomics.tsv"
output_file <- "kraken2_filtered_species_table_GLmetagenomics.tsv"
threshold <- 0.5

# string used to define non-microbial taxa
non_microbial <- "UNCLASSIFIED|Unclassified|unclassified|Homo sapien|cannot|uncultured|unidentified"

# read in feature table
feature_table <- read_delim(feature_table_file) %>%
                 across(where(is.numeric), function(col) replace_na(col, 0))) %>%
                 as.data.frame()
feature_name <- colnames(feature_table)[1]
rownames(feature_table) <- feature_table[,1]
feature_table <- feature_table[, -1]

# read-based count table
table2write <- filter_rare(feature_table, non_microbial, threshold = threshold) %>%
  as.data.frame %>%
  rownames_to_column(feature_name)

write_tsv(x = table2write, file = output_file)
```

**Custom Functions Used:**
- [group_low_abund_taxa()](#group_low_abund_taxa)

**Parameter Definitions:**

- `threshold` - threshold for filtering out rare taxa, a percentage between 0 and 100.
- `output_file` - output filename
- `input_file` - input filename

**Input Data:**

- kraken2_species_table_GLmetagenomics.tsv (path to kaiju species table from [Step 19ci](#19ci-create-merged-kraken2-taxonomy-table))

**Output Data:**

- **kraken2_filtered_species_table_GLmetagenomics.tsv** (a file containing the filtered species table)

---

#### 19g. Kraken2 Taxonomy Barplots

```R
species_table_file <- "kraken2_species_table_GLmetagenomics.tsv"
filtered_species_table_file <- "kraken2_filtered_species_table_GLmetagenomics.tsv"
metadata_file <- "/path/to/sample/metadata"

make_barplot(metadata_file = metadata_file, feature_table_file = species_table_file, 
             feature_column = "species", samples_column = "sample_id", group_column = "group",
             output_prefix = "kraken2_unfiltered_species", assay_suffix = "_GLmetagenomics",
             publication_format = publication_format, custom_palette = custom_palette)

# Save static unfiltered plot
make_barplot(metadata_file = metadata_file, feature_table_file = filtered_species_table_file, 
             feature_column = "Species", samples_column = "sample_id", group_column = "group",
             output_prefix = "kraken2_filtered_species", assay_suffix = "_GLmetagenomics",
             publication_format = publication_format, custom_palette = custom_palette)
```

**Custom Functions Used:**
- [make_barplot()](#make_barplot)

**Parameter Definitions:**

- `species_table_file` - a file containing the species count table
- `filtered_species_table_file` - a file containing the filtered species count table
- `metadata_file` - a file containing group information for each sample in the species count files

**Input Data:**

- `kraken2_species_table_GLmetagenomics.tsv` (path to kaiju species table from [Step 19ci](#19ci-create-merged-kraken2-taxonomy-table))
- `kraken2_filtered_species_table_GLmetagenomics.tsv` (a file containing the filtered species count table, output from [Step 19f](#19f-filter-kraken2-species-count-table))
- `/path/to/sample/metadata` (a file containing sample-wise metadata, mapping sample names to group metadata)

**Output Data:**

- kraken2_unfiltered_species_barplot_GLmetagenomics.png (taxonomy barplot without filtering)
- **kraken2_unfiltered_species_barplot_GLmetagenomics.html** (interactive taxonomy barplot without filtering)
- kraken2_filtered_species_barplot_GLmetagenomics.png (taxonomy barplot after filtering rare and non-microbial taxa)
- **kraken2_filtered_species_barplot_GLmetagenomics.html** (interactive taxonomy barplot after filtering rare and non-microbial taxa)

<br>  

---

### 20. Taxonomic Profiling Using HUMAnN/MetaPhlan

#### 20a. Download and Install HUMAnN databases

```bash
mkdir -p /path/to/humann3-db
humann_databases --download chocophlan full /path/to/humann3-db/
humann_databases --download uniref uniref90_ec_filtered_diamond /path/to/humann3-db/
humann_databases --download utility_mapping full /path/to/humann3-db/
metaphlan --install
```

**Parameter Definition:**

*humann3_databases*
- `--download` - Specifies the databases to download:
  - `chocophlan full` - the full ChocoPhlAn pangenome database, which includes Archaea, Bacteria, Eukaryotes, and Viruses
  - `uniref uniref90_ec_filtered_diamond` - Download the EC-filtered UniRef90 translated search database
  - `utility_mapping full` - additional gene family to functional category mapping database
-`/path/to/humann3-db` - Specifies the database install location

*metaphlan*
`--install` - install the MetaPhlan clade markers and database locally

**Input Data**

*No input data required*

**Output Data**

- /path/to/humann3-db (the installed MetaPhlan databases)

#### 20b. HUMAnN/MetaPhlAn Taxonomic Classification
```bash
  # forward and reverse reads need to be provided combined if paired-end (if not paired-end, single-end reads are provided to the --input argument next)
cat sample_R1_filtered_GLmetagenomics.fastq.gz sample_R2_filtered_GLmetagenomics.fastq.gz > sample-combined.fastq.gz

humann --input sample-combined.fastq.gz \
       --output sample-humann3-out-dir \
       --threads NumberOfThreads \
       --output-basename sample \
       --metaphlan-options "--bowtie2db /path/to/humann3-db/ --unclassified_estimation --add_viruses --sample_id sample" \
       --nucleotide-database /path/to/humann3-db/ \
       --protein-database /path/to/humann3-db/ \
       --bowtie-options "--sensitive --mm"

mv sample-humann3-out-dir/sample_humann_temp/sample_metaphlan_bugs_list.tsv \
   sample-humann3-out-dir/sample_metaphlan_bugs_list.tsv
```

**Parameter Definitions:**  

-	`--input` – specifies the input (combined forward and reverse reads)
-	`--output` – specifies output directory
-	`--threads` – specifies the number of threads to use
-	`--output-basename` – specifies prefix of the output files
-	`--metaphlan-options` – options to be passed to metaphlan
	- `--bowtie2db` – path to bowtie2 indexes (stored in HUMAnN database folder)
  - `unclassified_estimation` - scale the relative abundance profile according to the percentage of reads mapping to a clade.
	- `--add_viruses` – include viruses in the reference database
	- `--sample_id` – specifies the sample identifier we want in the table (rather than full filename)

**Input Data:**

- `/path/to/humann3-db/` (HUMAnN databases installed in [Step 20a](#20a-download-and-install-humann-databases))
- *_R[12]_filtered_GLmetagenomics.fastq.gz (filtered/trimmed reads from [Step 2b](#2b-trim-polyg) above)

**Output Data:**

- sample-humann3-out-dir/ (humann output directories containing *genefamilies.tsv, *pathabundance.tsv, and *pathcoverage.tsv files)

#### 20c. Merge Multiple Sample Functional Profiles

```bash
  # they need to be in their own directories
mkdir genefamily-results/ pathabundance-results/ pathcoverage-results/

  # copying results from previous running humann3 step (16a) to get them all together in their own directories (as is needed)
cp *-humann3-out-dir/*genefamilies.tsv genefamily-results/
cp *-humann3-out-dir/*abundance.tsv pathabundance-results/
cp *-humann3-out-dir/*coverage.tsv pathcoverage-results/

humann_join_tables -i genefamily-results/ -o gene-families.tsv
humann_join_tables -i pathabundance-results/ -o path-abundances.tsv
humann_join_tables -i pathcoverage-results/ -o path-coverages.tsv
```

**Parameter Definitions:**  

- `-i` - the directory holding the input tables
- `-o` - the name of the output table holding combined data

**Input Data:**

- `sample-humann3-out-dir` (HUMAnN output directory, from [Step 20b](#20b-humannmetaphlan-taxonomic-classification))

**Output Data:**

- gene-families.tsv (Combined gene family table in tab-separated format.)
- pathway-abundances.tsv (Combined path abundances table in tab-separated format.)
- pathway-coverages.tsv (Combined path coverages table in tab-separated format.)

#### 20d. Split Results Tables

The read-based functional annotation tables have taxonomic info and non-taxonomic info mixed together initially. `humann` comes with a helper script to split these. Here we are using that to generate both non-taxonomically grouped functional info files and taxonomically grouped ones.

```bash
humann_split_stratified_table -i gene-families.tsv -o ./
mv gene-families_stratified.tsv Gene-families-grouped-by-taxa_GLmetagenomics.tsv
mv gene-families_unstratified.tsv Gene-families_GLmetagenomics.tsv

humann_split_stratified_table -i path-abundances.tsv -o ./
mv path-abundances_stratified.tsv Path-abundances-grouped-by-taxa_GLmetagenomics.tsv
mv path-abundances_unstratified.tsv Path-abundances_GLmetagenomics.tsv

humann2_split_stratified_table -i path-coverages.tsv -o ./
mv path-coverages_stratified.tsv Path-coverages-grouped-by-taxa_GLmetagenomics.tsv
mv path-coverages_unstratified.tsv Path-coverages_GLmetagenomics.tsv
```

**Parameter Definitions:**  

-	`-i` – the input combined table
-	`-o` – output directory (here specifying current directory)

**Input Data:**

- gene-families.tsv (Combined gene family table from [Step 20c](#20c-merge-multiple-sample-functional-profiles))
- pathway-abundances.tsv (Combined path abundances table from [Step 20c](#20c-merge-multiple-sample-functional-profiles))
- pathway-coverages.tsv (Combined path coverages table from [Step 20c](#20c-merge-multiple-sample-functional-profiles))

**Output Data:**

- **Gene-families_GLmetagenomics.tsv** (gene-family abundances)
- **Gene-families-grouped-by-taxa_GLmetagenomics.tsv** (gene-family abundances grouped by taxa)
- **Pathway-abundances_GLmetagenomics.tsv**  (pathway abundances)
- **Pathway-abundances-grouped-by-taxa_GLmetagenomics.tsv** (pathway abundances grouped by tax)
- **Pathway-coverages_GLmetagenomics.tsv** (pathway coverages)
- **Pathway-coverages-grouped-by-taxa_GLmetagenomics.tsv** (pathway coverages grouped by taxa)

#### 20e. Normalize Gene Families and Pathway Abundances Tables
Generates some normalized tables of the read-based functional outputs from humann that are more readily suitable for across sample comparisons.

```bash
humann_renorm_table -i Gene-families_GLmetagenomics.tsv -o Gene-families-cpm_GLmetagenomics.tsv --update-snames
humann_renorm_table -i Path-abundances_GLmetagenomics.tsv -o Path-abundances-cpm_GLmetagenomics.tsv --update-snames
```

**Parameter Definitions:**  

-	`-i` – the input combined table
-	`-o` – name of the output normalized table
-	`--update-snames` – change suffix of column names in tables to "-CPM"

**Input Data:**

- Gene-families_GLmetagenomics.tsv (gene-family abundances, from [Step 20d](#20d-split-results-tables))
- Pathway-abundances_GLmetagenomics.tsv (pathway abundances, from [Step 20d](#20d-split-results-tables))

**Output Data:**
- **Gene-families-cpm_GLmetagenomics.tsv** (gene-family abundances normalized to copies-per-million)
- **Pathway-abundances-cpm_GLmetagenomics.tsv** (pathway abundances normalized to copies-per-million)

#### 20f. Generate Normalized Gene-family Table Grouped by Kegg Orthologs (KOs)

```bash
humann_regroup_table -i Gene-families_GLmetagenomics.tsv -g uniref90_ko | \
humann_rename_table -n kegg-orthology | \
humann_renorm_table -o Gene-families-KO-cpm_GLmetagenomics.tsv --update-snames

```

**Parameter Definitions:**  

*humann_regroup_table*
-	`-i` – the input table
-	`-g` – the map to use to group uniref IDs into Kegg Orthologs
-	`|` – sending that output into the next humann command to add human-readable Kegg Orthology names

*humann_rename_table*
-	`-n` – specifying we are converting Kegg orthology IDs into Kegg orthology human-readable names
-	`|` – sending that output into the next humann command to normalize to copies-per-million

*humann_renorm_table*
-	`-o` – specifying the final output file name
-  `--update-snames` – change suffix of column names in tables to "-CPM"

**Input Data:**

- Gene-families_GLmetagenomics.tsv (Non-taxonomically grouped gene families, from [Step 20d](#20d-split-results-tables))

**Output Data:**

- **Gene-families-KO-cpm_GLmetagenomics.tsv** (KO term abundances normalized to copies-per-million)

#### 20g. Combine MetaPhlan Taxonomy Tables

```bash
merge_metaphlan_tables.py *-humann3-out-dir/*_humann_temp/*_metaphlan_bugs_list.tsv > Metaphlan-taxonomy_GLmetagenomics.tsv
```

**Parameter Definitions:**  

*merge_metaphlan_tables.py*
- positional argument specifying input files and output filename

*sed*
- `-i` - Perform the search/replace in-place on the input file.

**Input data:**

-	\*-humann3-out-dir/\*_humann_temp/\*_metaphlan_bugs_list.tsv (MetaPhlan bugs_list produced during humann3 run in [step 20b](#20b-humannmetaphlan-taxonomic-classification))

**Output data:**

- **Metaphlan-taxonomy_GLmetagenomics.tsv** (MetaPhlan estimated taxonomic relative abundances)

#### 20h. Create MetaPhlan Species Count Table

##### 20hi. Get Sample Read Counts

```bash
unzip filtered_multiqc_GLmetagenomics_data.zip

grep _R1_filtered multiqc_fastqc.txt | awk 'BEGIN{FS="\t"; OFS="\t"}{print $1,int($5)}' > reads_per_sample.tsv
```

**Input Data:**

- filtered_multiqc_GLmetagenomics_data.zip or HostRm_multiqc_GLmetagenomics_data.zip (multiqc data from [Step 2d](#2d-compile-filteredtrimmed-data-qc)
  
**Output Data:**

- reads_per_sample.txt (a 2-column tab delimited file with the sample names and read counts as column 1 and 2, respectively)

##### 20hii. Process MetaPhlan Taxonomy Table

```R
input_file <- "metaphlan-taxonomy_GLmetagenomics.tsv"
read_count_file <- "reads_per_sample.tsv"
output_file <- "metaphlan_species_table_GLmetagenomics.tsv"
threshold <- 0.5

taxon_levels <- c("Kingdom", "Phylum", "Class", "Order",
                  "Family", "Genus", "Species")

# read in feature table
feature_table <- read_delim(input_file, delim="\t", comment="#") 
colnames(feature_table)[1] <- "taxonomy"

feature_table <- feature_table %>%
  filter(str_detect(taxonomy, "UNCLASSIFIED|s__") & 
         str_detect(taxonomy, "t__", negate = TRUE)) %>%
  mutate(Species=str_replace_all(taxonomy, '\\w__', "")) %>%
  separate(Species, into=taxon_levels, sep="\\|") %>%
  mutate(across(where(is.character), function(x) replace_na(x, "UNCLASSIFIED"))) %>%
  mutate(Species=str_replace_all(Species, "_", " ")) %>%
  select(-taxonomy, -Kingdom, -Phylum, -Class, -Order, -Family, -Genus) %>%
  select(Species, everything()) %>%
  as.data.frame

rownames(feature_table) <- feature_table$Species
feature_table <- feature_table[,-match("Species", colnames(feature_table))]

# Set max abundance equal to 1
tab2 <- (feature_table %>% t) / 100

# read in sample read counts
counts <- read_delim(read_count_file, delim = "\t", 
                     col_names = c("Sample_ID", "Reads")) %>%
  as.data.frame

# Set rownames as sample names
rownames(counts) <- counts$Sample_ID
# Drop the Sample_ID column
counts <- counts[, -1, drop = FALSE]

tab2 <- tab2[rownames(counts),]

# Convert relative abundance to raw count
species_table <- map2(tab2 %>% as.data.frame, 
                      colnames(tab2), function(col, specie) {
                        df <- col * counts
                        colnames(df) <- specie
                        return(df) 
                      }) %>% list_cbind() %>% t

table2write <- species_table  %>%
  as.data.frame() %>%
  rownames_to_column("Species")

write_tsv(x = table2write, file = "metaphlan_species_table_GLmetagenomics.tsv")
```

**Input Data:**

- metaphlan-taxonomy_GLmetagenomics.tsv (MetaPhlan taxonomy table from [Step 20g](#20g-combine-metaphlan-taxonomy-tables))
- reads_per_sample.tsv (a 2-column tab delimited file with sample names and read counts as columns 1 and 2, respectively from [Step 20hi](#20hi-get-sample-read-counts))

**Output Data:**

- **metaphlan_species_table_GLmetagenomics.tsv** (a file containing the MetaPhlan species table)

#### 20i. Filter MetaPhlan Species Count Table

```R
feature_table_file <- "metaphlan_species_table_GLmetagenomics.tsv"
output_file <- "metaphlan_filtered_species_table_GLmetagenomics.tsv"
threshold <- 0.5

# string used to define non-microbial taxa
non_microbial <- "UNCLASSIFIED"

# read in feature table
feature_table <- read_delim(feature_table_file) %>%
                 mutate(across(where(is.numeric), function(col) replace_na(col, 0))) %>%
                 as.data.frame()
feature_name <- colnames(feature_table)[1]
rownames(feature_table) <- feature_table[,1]
feature_table <- feature_table[, -1]

# read-based count table
table2write <- filter_rare(feature_table, non_microbial, threshold = threshold) %>%
  as.data.frame %>%
  rownames_to_column(feature_name)

write_tsv(x = table2write, file = output_file)
```

**Custom Functions Used:**
- [group_low_abund_taxa()](#group_low_abund_taxa)

**Parameter Definitions:**

- `threshold` - threshold for filtering out rare taxa, a percentage between 0 and 100.
- `output_file` - output filename
- `input_file` - input filename

**Input Data:**

- metaphlan_species_table_GLmetagenomics.tsv (path to MetaPhlan species count table from [Step 20hii](#20hii-process-metaphlan-taxonomy-table))

**Output Data:**

- **metaphlan_filtered_species_table_GLmetagenomics.tsv** (a file containing the filtered MetaPhlan species table)

#### 20j. MetaPhlan Taxonomy Barplots

```R
species_table_file <- "metaphlan_species_table_GLmetagenomics.tsv"
filtered_species_table_file <- "metaphlan_filtered_species_table_GLmetagenomics.tsv"
metadata_file <- "/path/to/sample/metadata"

make_barplot(metadata_file = metadata_file, feature_table_file = species_table_file, 
             feature_column = "Species", samples_column = "sample_id", group_column = "group",
             output_prefix = "metaphlan_unfiltered_species", assay_suffix = "_GLmetagenomics",
             publication_format = publication_format, custom_palette = custom_palette)

# Save static unfiltered plot
make_barplot(metadata_file = metadata_file, feature_table_file = filtered_species_table_file, 
             feature_column = "Species", samples_column = "sample_id", group_column = "group",
             output_prefix = "metaphlan_filtered_species", assay_suffix = "_GLmetagenomics",
             publication_format = publication_format, custom_palette = custom_palette)
```

**Custom Functions Used:**
- [make_barplot()](#make_barplot)

**Parameter Definitions:**

- `species_table_file` - a file containing the species count table
- `filtered_species_table_file` - a file containing the filtered species count table
- `metadata_file` - a file containing group information for each sample in the species count files

**Input Data:**

- `metaphlan_species_table_GLmetagenomics.tsv` (path to MetaPhlan species table from [Step 20h](#20h-create-metaphlan-species-count-table))
- `metaphlan_filtered_species_table_GLmetagenomics.tsv` (a file containing the filtered species count table, output from [Step 20i](#20i-filter-metaphlan-species-count-table))
- `/path/to/sample/metadata` (a file containing sample-wise metadata, mapping sample names to group metadata)

**Output Data:**

- metaphlan_unfiltered_species_barplot_GLmetagenomics.png (taxonomy barplot without filtering)
- **metaphlan_unfiltered_species_barplot_GLmetagenomics.html** (interactive taxonomy barplot without filtering)
- metaphlan_filtered_species_barplot_GLmetagenomics.png (taxonomy barplot after filtering rare and non-microbial taxa)
- **metaphlan_filtered_species_barplot_GLmetagenomics.html** (interactive taxonomy barplot after filtering rare and non-microbial taxa)

#### 20k. Filter Humann Output

```R
# read in humann tables
humann_uniref_table <- read_delim(file = "Gene-families-cpm_GLmetagenomics.tsv", delim = "\t")
humann_KO_table <- read_delim(file = "Gene-families-KO-cpm_GLmetagenomics.tsv", delim = "\t")
humann_pathway_table <- read_delim(file = "Pathway-abundances-cpm_GLmetagenomics.tsv", delim = "\t")

# rename headers
humann_uniref_table <-  humann_uniref_table  %>% 
  rename(Uniref90=`# Gene Family`) %>%
  mutate(Uniref90=str_replace_all(Uniref90, "UniRef90_", "")) %>%
  set_names(colnames(.) %>% str_replace_all("_Abundance-CPM", "")) %>%
  as.data.frame()
write_tsv(x = humann_uniref_table, file = "Gene-families-uniref_unfiltered_GLmetagenomics.tsv")

humann_KO_table <- humann_KO_table %>%
  rename(KO=`# Gene Family`) %>%
  set_names(colnames(.) %>% str_replace_all("_Abundance-CPM", "")) %>%
  as.data.frame()
write_tsv(x = humann_KO_table, file = "Gene-families-KO_unfiltered_GLmetagenomics.tsv")

humann_pathway_table <-  humann_pathway_table  %>% 
  rename(Pathway=`# Pathway`) %>%
  set_names(colnames(.) %>% str_replace_all("_Abundance-CPM", "")) %>%
  as.data.frame()
write_tsv(x = humann_pathway_table, file = "Pathway-abundances_unfiltered_GLmetagenomics.tsv")

# filter data
threshold <- 500

humann_uniref_table <- humann_uniref_table %>%
  mutate(across(where(is.numeric), function(col) replace_na(col, 0))) %>% column_to_rownames("Uniref90")
humann_uniref_filtered <- get_abundant_features(humann_uniref_table, cpm_threshold = threshold) %>%
  as.data.frame() %>% rownames_to_column("Uniref90")
write_tsv(x = table2write, file = "Gene-families-uniref_filtered_GLmetagenomics.tsv")

humann_KO_table <- humann_KO_table %>%
  mutate(across(where(is.numeric), function(col) replace_na(col, 0))) %>% column_to_rownames("KO")
humann_KO_filtered <- get_abundant_features(humann_KO_table, cpm_threshold = threshold) %>%
  as.data.frame() %>% rownames_to_column("KO")
write_tsv(x = table2write, file = "Gene-families-KO_filtered_GLmetagenomics.tsv")

humann_pathway_table <- humann_pathway_table %>%
  mutate(across(where(is.numeric), function(col) replace_na(col, 0))) %>% column_to_rownames("Pathway")
humann_pathway_filtered <- get_abundant_features(humann_pathway_table, cpm_threshold = threshold) %>%
  as.data.frame() %>% rownames_to_column("Pathway")
write_tsv(x = table2write, file = "Pathway-abundances_filtered_GLmetagenomics.tsv")

```

**Custom Functions Used:**
- [get_abundant_features()](#get_abundant_features)

**Parameter Definitions:**

- `threshold` - threshold for filtering out low abundance features, a value greater than 0

**Input Data:**

- Gene-families-cpm_GLmetagenomics.tsv (Humann taxonomy table from [Step 20e](#20e-normalize-gene-families-and-pathway-abundances-tables))
- Gene-families-KO-cpm_GLmetagenomics.tsv (Humann pathway table from [Step 20e](#20e-normalize-gene-families-and-pathway-abundances-tables))
- Pathway-abundances-cpm_GLmetagenomics.tsv (Humann KO function table from [Step 20f](#20f-generate-normalized-gene-family-table-grouped-by-kegg-orthologs-kos))

**Output Data:**

- Gene-families-KO_unfiltered_GLmetagenomics.tsv (KO term abundances normalized to copies-per-million, with cleaned headers)
- Gene-families-uniref_unfiltered_GLmetagenomics.tsv (gene-family abundances normalized to copies-per-million, with cleaned headers)
- Pathway-abundances_unfiltered_GLmetagenomics.tsv (pathway abundances normalized to copies-per-million, with cleaned headers)
- **Gene-families-KO_filtered_GLmetagenomics.tsv** (KO term abundances filtered for features with less than 500 CPM across samples) 
- **Gene-families-uniref_filtered_GLmetagenomics.tsv** (gene-family abundances filtered for features with less than 500 CPM across samples) 
- **Gene-families-KO_filtered_GLmetagenomics.tsv** (Pathway abundances filtered for features with less than 500 CPM across samples) 

#### 20l. Create Humann Function Heatmaps

```R
metadata_table < "/path/to/sample_metadata"

make_heatmap(metadata_table_file = metadata_table, 
             feature_table_file = "Gene-families-uniref_unfiltered_GLmetagenomics.tsv", 
             samples_column="sample_id", group_column = "group", 
             output_prefix = "Gene-families-uniref_unfiltered", 
             assay_suffix = "_GLmetagenomics", 
             custom_palette = custom_palette)

make_heatmap(metadata_table_file = metadata_table, 
             feature_table_file = "Gene-families-uniref_filtered_GLmetagenomics.tsv", 
             samples_column="sample_id", group_column = "group", 
             output_prefix = "Gene-families-uniref_filtered", 
             assay_suffix = "_GLmetagenomics", 
             custom_palette = custom_palette)

make_heatmap(metadata_table_file = metadata_table, 
             feature_table_file = "Gene-families-KO_unfiltered_GLmetagenomics.tsv", 
             samples_column="sample_id", group_column = "group", 
             output_prefix = "Gene-families-KO_unfiltered", 
             assay_suffix = "_GLmetagenomics", 
             custom_palette = custom_palette)

make_heatmap(metadata_table_file = metadata_table, 
             feature_table_file = "Gene-families-KO_filtered_GLmetagenomics.tsv", 
             samples_column="sample_id", group_column = "group", 
             output_prefix = "Gene-families-KO_filtered", 
             assay_suffix = "_GLmetagenomics", 
             custom_palette = custom_palette)

make_heatmap(metadata_table_file = metadata_table, 
             feature_table_file = "Pathway-abundances_unfiltered_GLmetagenomics.tsv", 
             samples_column="sample_id", group_column = "group", 
             output_prefix = "Pathway-abundances_unfiltered", 
             assay_suffix = "_GLmetagenomics", 
             custom_palette = custom_palette)

make_heatmap(metadata_table_file = metadata_table, 
             feature_table_file = "Pathway-abundances_filtered_GLmetagenomics.tsv", 
             samples_column="sample_id", group_column = "group", 
             output_prefix = "Pathway-abundances_filtered", 
             assay_suffix = "_GLmetagenomics", 
             custom_palette = custom_palette)
```

**Custom Functions Used:**
- [make_heatmap()](#make_heatmap)

**Parameter Definitions:**

- `metadata_file` - a file containing group information for each sample in the species count files

**Input Data:**

- `/path/to/sample/metadata` (a file containing sample-wise metadata, mapping sample names to group metadata)
- `Gene-families-uniref_unfiltered_GLmetagenomics.tsv` (gene-family abundances table, output from [Step 20k](#20k-filter-humann-output))
- `Gene-families-KO_unfiltered_GLmetagenomics.tsv` (KO term abundances table, output from [Step 20k](#20k-filter-humann-output))
- `Pathway-abundances_unfiltered_GLmetagenomics.tsv` (pathway abundances table, output from [Step 20k](#20k-filter-humann-output))
- `Gene-families-uniref_filtered_GLmetagenomics.tsv` (filtered gene-family abundances table, output from [Step 20k](#20k-filter-humann-output)) 
- `Gene-families-KO_filtered_GLmetagenomics.tsv` (filtered KO term abundances table, output from [Step 20k](#20k-filter-humann-output)) 
- `Pathway-abundances_filtered_GLmetagenomics.tsv` (filtered Pathway abundances table, output from [Step 20k](#20k-filter-humann-output)) 

**Output Data:**

- **Gene-families-uniref_unfiltered_heatmap_GLmetagenomics.png** (gene family abundances heatmap without filtering)
- **Gene-families-uniref_filtered_heatmap_GLmetagenomics.png** (gene family abundances heatmap after filtering rare and non-microbial taxa)
- **Gene-families-KO_unfiltered_heatmap_GLmetagenomics.png** (KO term abundances heatmap without filtering)
- **Gene-families-KO_filtered_heatmap_GLmetagenomics.png** (KO term abundances heatmap after filtering rare and non-microbial taxa)
- **Pathway-abundances_unfiltered_heatmap_GLmetagenomics.png** (pathway abundances heatmap without filtering)
- **Pathway-abundances_filtered_heatmap_GLmetagenomics.png** (pathway abundances heatmap after filtering rare and non-microbial taxa)
- **Gene-families-uniref_unfiltered_top_50_heatmap_GLmetagenomics.png** (gene family abundances heatmap without filtering)
- **Gene-families-uniref_filtered_top_50_heatmap_GLmetagenomics.png** (gene family abundances heatmap after filtering rare and non-microbial taxa)
- **Gene-families-KO_unfiltered_top_50_heatmap_GLmetagenomics.png** (KO term abundances heatmap without filtering)
- **Gene-families-KO_filtered_top_50_heatmap_GLmetagenomics.png** (KO term abundances heatmap after filtering rare and non-microbial taxa)
- **Pathway-abundances_unfiltered_top_50_heatmap_GLmetagenomics.png** (pathway abundances heatmap without filtering)
- **Pathway-abundances_filtered_top_50_heatmap_GLmetagenomics.png** (pathway abundances heatmap after filtering rare and non-microbial taxa)

---
