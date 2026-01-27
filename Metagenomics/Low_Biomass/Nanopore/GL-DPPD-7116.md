# Bioinformatics pipeline for Low biomass long-read metagenomics data

> **This document holds an overview and some example commands of how GeneLab processes low-biomass, long-read metagenomics datasets. Exact processing commands for specific datasets that have been released are provided with their processed data in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/).**  

---

**Date:** January MM, 2026  
**Revision:** -  
**Document Number:** GL-DPPD-7116  

**Submitted by:**  
Olabiyi A. Obayomi (GeneLab Analysis Team)  

**Approved by:**  
Samrawit Gebre (OSDR Project Manager)  
Jonathan Galazka (OSDR Project Scientist)  
Amanda Saravia-Butler (GeneLab Science Lead)  
Barbara Novak (GeneLab Data Processing Lead)  


---

# Table of contents

- [**Software used**](#software-used)
- [**General processing overview with example commands**](#general-processing-overview-with-example-commands)
  - [**Pre-processing**](#pre-processing)
    - [1. Basecalling](#1-basecalling)
    - [2. Demultiplexing](#2-demultiplexing)
      - [2a. Split fastq ](#2a-split-fastq)
      - [2b. Concatenate files for each sample](#2b-concatenate-files-for-each-sample)
    - [3. Raw Data QC](#3-raw-data-qc)
      - [3a. Raw Data QC](#3a-raw-data-qc)
      - [3b. Compile Raw Data QC](#3b-compile-raw-data-qc)
    - [4. Quality filtering](#4-quality-filtering)
      - [4a. Filter Raw Data](#4a-filter-raw-data)
      - [4a. Filtered Data QC](#4b-filtered-data-qc)
      - [4c. Compile Filtered Data QC](#4c-compile-filtered-data-qc)
    - [5. Trimming](#5-trimming)
      - [5a. Trim Filtered Data](#5a-trim-filtered-data)
      - [5b. Trimmed Data QC](#5b-trimmed-data-qc)
      - [5c. Compile Trimmed Data QC](#5c-compile-trimmed-data-qc)
    - [6. Human Read Removal](#6-human-read-removal)
      - [6a. Build Kraken2 Database](#6a-build-kraken2-database)
      - [6b. Remove Human Reads](#6b-remove-human-reads)
      - [6c. Compile Human Read Removal QC](#6c-compile-human-read-removal-qc)
    - [7. Contaminant Removal](#7-contaminant-removal)
      - [7a. Assemble Contaminants](#7a-assemble-contaminants)
      - [7b. Build Contaminant Index and Map Reads](#7b-build-contaminant-index-and-map-reads)
      - [7c. Sort and Index Contaminant Reads](#7c-sort-and-index-contaminant-alignments)
      - [7d. Gather Contaminant Mapping Metrics](#7d-gather-contaminant-mapping-metrics)
      - [7e. Generate Decontaminated Read Files](#7e-generate-decontaminated-read-files)
      - [7f. Contaminant Removal QC](#7f-contaminant-removal-qc)
      - [7g. Compile Contaminant Removal QC](#7g-compile-contaminant-removal-qc)
    - [8. Host Read Removal](#8-host-read-removal)
      - [8a. Build Kraken2 Database](#8a-build-kraken2-database)
      - [8b. Remove Host Reads](#8b-remove-host-reads)
      - [8c. Compile Host Read Removal QC](#8c-compile-host-read-removal-qc)
    - [9. R Environment Setup](#9-r-environment-setup)
      - [9a. Load Libraries](#9a-load-libraries)
      - [9b. Define Custom Functions](#9b-define-custom-functions)
      - [9c. Set global variables](#9c-set-global-variables)
  - [**Read-based processing**](#read-based-processing)
    - [10. Taxonomic profiling using kaiju](#10-taxonomic-profiling-using-kaiju)
      - [10a. Build Kaiju Database](#10a-build-kaiju-database)
      - [10b. Kaiju Taxonomic Classification](#10b-kaiju-taxonomic-classification)
      - [10c. Compile Kaiju Taxonomy Results](#10c-compile-kaiju-taxonomy-results)
      - [10d. Convert Kaiju Output To Krona Format](#10d-convert-kaiju-output-to-krona-format)
      - [10e. Compile Kaiju Krona Reports](#10e-compile-kaiju-krona-reports)
      - [10f. Create Kaiju Species Count Table](#10f-create-kaiju-species-count-table)
      - [10g. Read-in Tables](#10g-read-in-tables)
      - [10h. Taxonomy Barplots](#10h-taxonomy-barplots)
      - [10i. Feature Decontamination](#10i-feature-decontamination)
    - [11. Taxonomic Profiling Using Kraken2](#11-taxonomic-profiling-using-kraken2)
      - [11a. Download Kraken2 Database](#11a-download-kraken2-database)
      - [11b. Kraken2 Taxonomic Classification](#11b-kraken2-taxonomic-classification)
      - [11c. Compile Kraken2 Taxonomy Results](#11c-compile-kraken2-taxonomy-results)
        - [11ci. Create Merged Kraken2 Taxonomy Table](#11ci-create-merged-kraken2-taxonomy-table)
        - [11cii. Compile Kraken2 Taxonomy Reports](#11cii-compile-kraken2-taxonomy-reports)
      - [11d. Convert Kraken2 Output to Krona Format](#11d-convert-kraken2-output-to-krona-format)
      - [11e. Compile Kraken2 Krona Reports](#11e-compile-kraken2-krona-reports)
      - [11f. Create Kraken2 Species Count Table](#11f-create-kraken2-species-count-table)
      - [11g. Read-in Tables](#11g-read-in-tables)
      - [11h. Taxonomy Barplots](#11h-taxonomy-barplots)
      - [11i. Feature Decontamination](#11i-feature-decontamination)
  - [**Assembly-based processing**](#assembly-based-processing)
    - [12. Sample Assembly](#12-sample-assembly)
    - [13. Polish Assembly](#13-polish-assembly)
    - [14. Rename Contigs and Summarize Assemblies](#14-rename-contigs-and-summarize-assemblies)
      - [14a. Rename Contig Headers](#14a-rename-contig-headers)
      - [14b. Summarize Assemblies](#14b-summarize-assemblies)
    - [15. Gene Prediction](#15-gene-prediction)
      - [15a. Generate Gene Predictions](15a-generate-gene-predictions)
      - [15b. Remove Line Wraps In Gene Prediction Output](#15a-remove-line-wraps-in-gene-prediction-output)
    - [16. Functional Annotation](#16-functional-annotation)
      - [16a. Download Reference Database of HMM Models](#16a-download-reference-database-of-hmm-models)
      - [16b. Run KEGG Annotation](#16b-run-kegg-annotation)
      - [16c. Filter KO Outputs](#16c-filter-ko-outputs)
    - [17. Taxonomic Classification](#17-taxonomic-classification)
      - [17a. Pull and Unpack Pre-built Reference DB](#17a-pull-and-unpack-pre-built-reference-db)
      - [17b. Run Taxonomic Classification](#17b-run-taxonomic-classification)
      - [17c. Add Taxonomy Info From Taxids To Genes](#17c-add-taxonomy-info-from-taxids-to-genes)
      - [17d. Add Taxonomy Info From Taxids To Contigs](#17d-add-taxonomy-info-from-taxids-to-contigs)
      - [17e. Format Gene-level Output With awk and sed](#17e-format-gene-level-output-with-awk-and-sed)
      - [17f. Format Contig-level Output With awk and sed](#17f-format-contig-level-output-with-awk-and-sed)
    - [18. Read-Mapping](#17-read-mapping)
      - [18a. Align Reads to Sample Assembly](#18a-align-reads-to-sample-assembly)
      - [18b. Sort and Index Assembly Alignments](#18b-sort-and-index-assembly-alignments)
    - [19. Get Coverage Information and Filter Based On Detection](#19-get-coverage-information-and-filter-based-on-detection)
      - [19a. Filter Coverage Levels Based On Detection](#19a-filter-coverage-levels-based-on-detection)
      - [19b. Filter Gene and Contig Coverage Based On Detection](#19b-filter-gene-and-contig-coverage-based-on-detection)
    - [20. Combine Gene-level Coverage, Taxonomy, and Functional Annotations For Each Sample](#20-combine-gene-level-coverage-taxonomy-and-functional-annotations-for-each-sample)
    - [21. Combine Contig-level Coverage and Taxonomy For Each Sample](#21-combine-contig-level-coverage-and-taxonomy-for-each-sample)
    - [22. Generate Normalized, Gene- and Contig-level Coverage Summary Tables of KO-annotations and Taxonomy Across Samples](#22-generate-normalized-gene--and-contig-level-coverage-summary-tables-of-ko-annotations-and-taxonomy-across-samples)
      - [22a. Generate Gene-level Coverage Summary Tables](#22a-generate-gene-level-coverage-summary-tables)
      - [22b. Generate Contig-level Coverage Summary Tables](#22b-generate-contig-level-coverage-summary-tables)
    - [23. **M**etagenome-**A**ssembled **G**enome (MAG) recovery](#23-metagenome-assembled-genome-mag-recovery)
      - [23a. Bin Contigs](#23a-bin-contigs)
      - [23b. Bin Quality Assessment](#23b-bin-quality-assessment)
      - [23c. Filter MAGs](#23c-filter-mags)
      - [23d. MAG Taxonomic Classification](#23d-mag-taxonomic-classification)
      - [23e. Generate Overview Table Of All MAGs](#23e-generate-overview-table-of-all-mags)
    - [24. Generate MAG-level Functional Summary Overview](#24-generate-mag-level-functional-summary-overview)
      - [24a. Get KO Annotations Per MAG](#24a-get-ko-annotations-per-mag)
      - [24b. Summarize KO Annotations With KEGG-Decoder](#24b-summarize-ko-annotations-with-kegg-decoder)
    - [25. Decontamination and Visualization of Contig- and Gene-taxonomy and Gene-function Outputs](#25-decontamination-and-visualization-of-contig--and-gene-taxonomy-and-gene-function-outputs)
      - [25a. Gene-level taxonomy heatmaps](#25a-gene-level-taxonomy-heatmaps)
      - [25b. Gene-level taxonomy decontamination](#25b-gene-level-taxonomy-decontamination)
      - [25c. Gene-level KO functions heatmaps](#25c-gene-level-ko-functions-heatmaps)
      - [25d. Gene-level KO functions decontamination](#25d-gene-level-ko-functions-decontamination)
      - [25e. Contig-level heatmaps](#25e-contig-level-heatmaps)
      - [25f. Contig-level decontamination](#25f-contig-level-decontamination)


---

# Software used

|Program|Version|Relevant Links|
|:------|:-----:|------:|
|bbduk| 38.86 |[https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/)|
|bit| 1.8.53 |[https://github.com/AstrobioMike/bioinf_tools#bioinformatics-tools-bit](https://github.com/AstrobioMike/bioinf_tools#bioinformatics-tools-bit)|
|CAT| 5.2.3 |[https://github.com/dutilh/CAT#cat-and-bat](https://github.com/dutilh/CAT#cat-and-bat)|
|CheckM| 1.1.3 |[https://github.com/Ecogenomics/CheckM](https://github.com/Ecogenomics/CheckM)|
|Dorado| 1.1.1| [https://github.com/nanoporetech/dorado](https://github.com/nanoporetech/dorado)|
|filtlong| 0.2.1 |[https://github.com/rrwick/Filtlong](https://github.com/rrwick/Filtlong)|
|Flye| 2.9.5 | [https://github.com/mikolmogorov/Flye](https://github.com/mikolmogorov/Flye) |
|GTDB-Tk| 2.4.0 |[https://github.com/Ecogenomics/GTDBTk](https://github.com/Ecogenomics/GTDBTk)|
|HUMAnN| 3.9 |[https://github.com/biobakery/humann](https://github.com/biobakery/humann)|
|Kaiju| 1.10.1 | [https://bioinformatics-centre.github.io/kaiju/](https://bioinformatics-centre.github.io/kaiju/) |
|KEGG-Decoder| 1.2.2 |[https://github.com/bjtully/BioData/tree/master/KEGGDecoder#kegg-decoder](https://github.com/bjtully/BioData/tree/master/KEGGDecoder#kegg-decoder)
|KOFamScan| 1.3.0 |[https://github.com/takaram/kofam_scan](https://github.com/takaram/kofam_scan)|
|Kraken2| 2.1.6 | [https://github.com/DerrickWood/kraken2](https://github.com/DerrickWood/kraken2) |
|KrakenTools | 1.2 | [https://ccb.jhu.edu/software/krakentools/](https://ccb.jhu.edu/software/krakentools/) |
|Krona| 2.8.1 | [https://github.com/marbl/Krona/wiki](https://github.com/marbl/Krona/wiki)|
|MetaBAT| 2.15 |[https://bitbucket.org/berkeleylab/metabat/src/master/](https://bitbucket.org/berkeleylab/metabat/src/master/)|
|Minimap2| 2.28 | [https://github.com/lh3/minimap2](https://github.com/lh3/minimap2) |
|MultiQC| 1.27.1 |[https://multiqc.info/](https://multiqc.info/)|
|Medaka| 2.1.1 | [https://github.com/nanoporetech/medaka](https://github.com/nanoporetech/medaka) |
|MetaPhlAn| 4.1.0 |[https://github.com/biobakery/MetaPhlAn](https://github.com/biobakery/MetaPhlAn)|
|NanoPlot| 1.44.1 | [https://github.com/wdecoster/NanoPlot](https://github.com/wdecoster/NanoPlot)|
|Porechop| 0.2.4 | [https://github.com/rrwick/Porechop](https://github.com/rrwick/Porechop) |
|Prodigal| 2.6.3 |[https://github.com/hyattpd/Prodigal#prodigal](https://github.com/hyattpd/Prodigal#prodigal)|
|samtools| 1.22.1 |[https://github.com/samtools/samtools#samtools](https://github.com/samtools/samtools#samtools)|
| R | 4.5.1 | [https://www.r-project.org](https://www.r-project.org) |
|Bioconductor | 3.21 | [https://www.bioconductor.org](https://www.bioconductor.org) |
|decontam| 1.28.0 | [https://www.bioconductor.org/packages/release/bioc/html/decontam.html](https://www.bioconductor.org/packages/release/bioc/html/decontam.html) |
|optparse| 1.7.5 |[https://cran.r-project.org/web/packages/optparse/index.html](https://cran.r-project.org/web/packages/optparse/index.html) |
|pavian| 1.2.1 | [https://github.com/fbreitwieser/pavian](https://github.com/fbreitwieser/pavian) |
|pheatmap| 1.0.13 | [https://cran.r-project.org/package=pheatmap](https://cran.r-project.org/package=pheatmap) |
|phyloseq| 1.52.0 | [https://www.bioconductor.org/packages/release/bioc/html/phyloseq.html](https://www.bioconductor.org/packages/release/bioc/html/phyloseq.html) |
|tidyverse| 2.0.0 | [https://www.tidyverse.org](https://www.tidyverse.org) |

---

# General processing overview with example commands

> Exact processing commands and output files listed in **bold** below are included with each Low Biomass Metagenomics Seq processed dataset in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/).  

## Pre-processing

### 1. Basecalling

```bash
model="hac" # high accuracy model
input_directory=/path/to/pod5/or/fast5/data
kit_name=SQK-RPB004

dorado basecaller ${model} ${input_directory} \
  --no-trim \
  --device auto \
  --recursive \
  --kit-name ${kit_name} \
  --min-qscore 8 > basecalled.bam
```

**Parameter Definitions:**

- `model` - Positional argument specifying the basecalling model to use or a path to the model directory. `hac` chooses the high accuracy model.
- `input_directory` - Positional argument specifying the location of the raw data in POD5 or FAST5 format.
- `--no-trim` - Skips trimming of barcodes, adapters, and primers.
- `--device` - Specifies CPU or GPU device; specifying 'auto' chooses either 'cpu' or 'gpu' depending on detected presence of a GPU device.
- `--recursive` - Enables recursive scanning through input directory to load FAST5 and/or POD5 files.
- `--kit-name` - The nanopore barcoding kit used during sequencing preparation. Enables barcoding with the provided kit name; see [dorado documentation](https://software-docs.nanoporetech.com/dorado/1.1.1/barcoding/barcoding/) for a full list of accepted kit names.
- `--min-qscore` - Specifies the minimum Q-score, reads with a mean Q-score below this threshold are discarded (default to `8`).

**Input Data:**

- *pod5 and/or *fast5 (raw nanopore data)

**Output Data:**

- basecalled.bam (basecalled data in bam format)

<br>

---

### 2. Demultiplexing

#### 2a. Split Fastq

```bash
dorado demux \
  --output-dir /path/to/fastq/output \
  --emit-fastq \
  --emit-summary \
  --kit-name ${kit_name} \
  basecalled.bam
```

**Parameter Definitions:**

- `--output-dir` - Specifies the output folder that is the root of the nested output structure. 
- `--emit-fastq` - Specifies that output is fastq format.
- `--emit-summary` - Creates a summary listing each read and its classified barcode.
- `--kit-name` - The nanopore barcoding kit used during sequencing preparation. Enables barcoding with the provided kit name; see [dorado documentation](https://software-docs.nanoporetech.com/dorado/1.1.1/barcoding/barcoding/) for a full list of accepted kit names.
- `basecalled.bam` - Positional argument specifying the input bam file.

**Input Data:**

- basecalled.bam (basecalled nanopore data in bam format, output from [Step 1](#1-basecalling))

**Output Data:**

- /path/to/fastq/output/\*_barcode\*.fastq (demultiplexed reads in fastq format)
- /path/to/fastq/output/\*_unclassified.fastq (unclassified reads in fastq format)
- /path/to/fastq/output/barcoding_summary.txt (barcode summary file listing each read, the file it was assigned to, and its classified barcode)


#### 2b. Concatenate Files For Each Sample

```bash
# Change to directory containing split fastq files generated from step 2a. 
cd /path/to/fastq/output/ # output of step 2a

# Get unique barcode names from demultiplexed file names
BARCODES=($(ls -1 *fastq* | sed -E 's/.+_(barcode[0-9]+)_.+/\1/g' | sort -u))

# Concat separate barcode/sample fastq files into per sample fastq gzipped files
[ -d raw_data/ ] || mkdir raw_data/
for sample in ${BARCODES[*]}; do

  [ -d  ${sample}/ ] ||  mkdir ${sample}/  
  mv *_${sample}_*  ${sample}/ 

  cat ${sample}/* | gzip > raw_data/${sample}.fastq.gz

done
```

**Parameter Definitions:**

- `cat ${sample}/*` - Concatenates all fastq files with the same barcode into one fastq file.
- `| gzip` - Sends the concatenated fastq file output from the `cat` command to the `gzip` command to create a compressed fastq.gz file for each barcode.

**Input Data:**

- /path/to/fastq/output/ (directory containing spilt fastq files from [Step 2a](#2a-split-fastq))

**Output Data:**

-  raw_data/sample.fastq.gz (gzipped per sample/barcode fastq files)

<br>

---

### 3. Raw Data QC

#### 3a. Raw Data QC

```bash 
NanoPlot --only-report \
         --prefix sample_raw_ \
         --outdir /path/to/raw_nanoplot_output \
         --threads NumberOfThreads \
         --fastq \
         /path/to/raw_data/sample.fastq.gz
```

**Parameter Definitions:**

- `--only-report` - Output only the report files.
- `--prefix` - Adds a sample specific prefix to the name of each output file.
- `--outdir` – Specifies the output directory to store results.
- `--threads` - Number of parallel processing threads to use.
- `--fastq` - Specifies that the input data is in fastq format.
- `/path/to/raw_data/sample.fastq.gz` – The input reads, specified as a positional argument.

**Input Data:**

- /path/to/raw_data/sample.fastq.gz (concatenated raw reads, output from [Step 2b](#2b-concatenate-files-for-each-sample))

**Output Data:**

- **/path/to/raw_nanoplot_output/sample_raw_NanoPlot-report.html** (NanoPlot html summary)
- /path/to/raw_nanoplot_output/sample_raw_NanoPlot_\<date\>_\<time\>.log (NanoPlot log file)
- /path/to/raw_nanoplot_output/sample_raw_NanoStats.txt (text file containing basic statistics)

#### 3b. Compile Raw Data QC

```bash 
multiqc --zip-data-dir \
        --outdir raw_multiqc_report \
        --filename raw_multiqc_GLlblMetag \
        --interactive \
        /path/to/raw_nanoplot_output/
```

**Parameter Definitions:**

- `--zip-data-dir` - Compress the data directory.
- `--outdir` – Specifies the output directory to store results.
- `--filename` – Specifies the filename prefix of results.
- `--interactive` - Force multiqc to always create interactive javascript plots.
- `/path/to/raw_nanoplot_output/` – The directory holding the output data from the NanoPlot run, provided as a positional argument.

**Input Data:**

- /path/to/raw_nanoplot_output/*raw_NanoStats.txt (NanoPlot output data, from [Step 3a](#3a-raw-data-qc))

**Output Data:**

- **raw_multiqc_GLlblMetag.html** (multiqc output html summary)
- **raw_multiqc_GLlblMetag_data.zip** (zip archive containing multiqc output data)

<br>  

---

### 4. Quality Filtering

#### 4a. Filter Raw Data

```bash
filtlong --min_length 200 --min_mean_q 8 /path/to/raw_data/sample.fastq.gz > sample_filtered.fastq
```

**Parameter Definitions:**

- `--min_length` – Specifies the minimum read length to retain (default: `200`).
- `--min_mean_q` – Specifies the minimum mean read quality to retain (default: `8`).
- `/path/to/raw_data/sample.fastq.gz` - The path to the input fastq file, provided as a positional argument.
- `> sample_filtered.fastq` - Redirects the output to a sample_filtered.fastq file.

**Input Data:**

- /path/to/raw_data/sample.fastq.gz (concatenated raw reads, output from [Step 2b](#2b-concatenate-files-for-each-sample))

**Output Data:**

- *sample_filtered.fastq (quality filtered reads)


#### 4b. Filtered Data QC

```bash
NanoPlot --only-report \
         --prefix sample_filtered_ \
         --outdir /path/to/filtered_nanoplot_output \
         --threads NumberOfThreads \
         --fastq \
         sample_filtered.fastq
```

**Parameter Definitions:**

- `--only-report` - Output only the report files.
- `--prefix` - Adds a sample specific prefix to the name of each output file.
- `--outdir` – Specifies the output directory to store results.
- `--threads` - Number of parallel processing threads to use.
- `--fastq` - Specifies that the input data is in fastq format.
- `sample_filtered.fastq` – The input reads, specified as a positional argument.

**Input Data:**

- sample_filtered.fastq (filtered reads, output from [Step 4a](#4a-filter-raw-data))

**Output Data:**

- **/path/to/filtered_nanoplot_output/sample_filtered_NanoPlot-report_GLlblMetag.html** (NanoPlot html summary)
- /path/to/filtered_nanoplot_output/sample_filtered_NanoPlot_\<date\>_\<time\>.log (NanoPlot log file)
- /path/to/filtered_nanoplot_output/sample_filtered_NanoStats.txt (text file containing basic statistics)

#### 4c. Compile Filtered Data QC

```bash
multiqc  --zip-data-dir \ 
         --outdir filtered_multiqc_report \
         --filename filtered_multiqc_GLlblMetag \
         --interactive \
         /path/to/filtered_nanoplot_output/
```

**Parameter Definitions:**

- `--zip-data-dir` - Compress the data directory.
- `--outdir` – Specifies the output directory to store results.
- `--filename` – Specifies the filename prefix of results.
- `--interactive` - Force multiqc to always create interactive javascript plots.
- `/path/to/filtered_nanoplot_output/` – The directory holding the output data from the NanoPlot run, provided as a positional argument.

**Input Data:**

- /path/to/filtered_nanoplot_output/*filtered_NanoStats.txt (NanoPlot output data, from [Step 4b](#4b-filtered-data-qc))

**Output Data:**

- **filtered_multiqc_report/filtered_multiqc_GLlblMetag.html** (multiqc output html summary)
- **filtered_multiqc_report/filtered_multiqc_GLlblMetag_data.zip** (zip archive containing multiqc output data)

<br>

---

### 5. Trimming

#### 5a. Trim Filtered Data

```bash
porechop --input sample_filtered.fastq \
         --threads NumberOfThreads \
         --discard_middle \
         --output sample_trimmed.fastq.gz  > sample_porechop.log
```

**Parameter Definitions:**

- `--input` – Specifies the input sequence file in fastq format.
- `--threads` - Number of parallel processing threads to use.
- `--discard_middle` -  Reads with middle adapters will be discarded.
- `--output` - Specifies the trimmed reads output fastq filename.
- `> sample_porechop.log` - Redirects the standard output to a log file.

**Input Data:**

- sample_filtered.fastq (filtered reads output from [Step 4a](#4a-filter-raw-data))

**Output Data:**

- sample_trimmed.fastq.gz (filtered and trimmed reads)
- sample_porechop.log (porechop standard output containing trimming info)

#### 5b. Trimmed Data QC

```bash
NanoPlot --only-report \
         --prefix sample_trimmed_ \
         --outdir /path/to/trimmed_nanoplot_output \
         --threads NumberOfThreads \
         --fastq \
         sample_trimmed.fastq.gz
```

**Parameter Definitions:**

- `--only-report` - Output only the report files.
- `--prefix` - Adds a sample specific prefix to the name of each output file.
- `--outdir` – Specifies the output directory to store results.
- `--threads` - Number of parallel processing threads to use.
- `--fastq` - Specifies that the input data is in fastq format.
- `sample_trimmed.fastq.gz` – The input reads, specified as a positional argument.

**Input Data:**

- sample_trimmed.fastq.gz (filtered and trimmed reads, output from [Step 5a](#5a-trim-filtered-data))

**Output Data:**

- **/path/to/trimmed_nanoplot_output/sample_trimmed_NanoPlot-report_GLlblMetag.html** (NanoPlot html summary)
- /path/to/trimmed_nanoplot_output/sample_trimmed_NanoPlot_\<date\>_\<time\>.log (NanoPlot log file)
- /path/to/trimmed_nanoplot_output/sample_trimmed_NanoStats.txt (text file containing basic statistics)

#### 5c. Compile Trimmed Data QC

```bash
multiqc --zip-data-dir \ 
        --outdir trimmed_multiqc_report \
        --filename trimmed_multiqc_GLlblMetag \
        --interactive \
        /path/to/trimmed_nanoplot_output/
```

**Parameter Definitions:**

- `--zip-data-dir` - Compress the data directory.
- `--outdir` – Specifies the output directory to store results.
- `--filename` – Specifies the filename prefix of results.
- `--interactive` - Force multiqc to always create interactive javascript plots.
- `/path/to/trimmed_nanoplot_output/` – The directory holding the output data from the NanoPlot run, provided as a positional argument.

**Input Data:**

- /path/to/trimmed_nanoplot_output/*trimmed_NanoStats.txt (NanoPlot output data, output from [Step 5b](#5b-trimmed-data-qc))

**Output Data:**

- **trimmed_multiqc_GLlblMetag.html** (multiqc output html summary)
- **trimmed_multiqc_GLlblMetag_data.zip** (zip archive containing multiqc output data)

<br>

---

### 6. Human Read Removal

#### 6a. Build Kraken2 Database

> **Note:** It is recommended to use NCBI genome files with kraken2 because sequences not downloaded from 
NCBI may require explicit assignment of taxonomy information before they can be used to build the 
database, as mentioned in the [Kraken2 Documentation](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown).

```bash
# Download NCBI taxonomic information 
kraken2-build --download-taxonomy --db kraken2-human-db/

# Add genomic sequences to your database's genomic library
kraken2-build --add-to-library human.fasta --db kraken2-human-db/ \
              --no-masking --kmer-length 35 --minimizer-length 31

# Build the database
kraken2-build --build --db kraken2-human-db/

# Clean up intermediate files
kraken2-build --clean --db kraken2-human-db/
```
**Parameter Definitions:**
- `--download-taxonomy` - Instructs kraken2-build to download the NCBI taxonomic information.
- `--db` - Specifies the name of the directory for the kraken2 database
- `--add-to-library` - Instructs kraken2-build to add the contents of a file (`human.fasta`) to the kraken2 DB library
- `--no-masking` - Disables masking of low-complexity sequences. For additional 
                   information see the [kraken documentation for masking](https://github.com/DerrickWood/kraken2/wiki/Manual#masking-of-low-complexity-sequences).
- `--build` - Instructs kraken2-build to build the kraken2 DB from the library files
- `--clean` - Instructs kraken2-build to remove unneeded intermediate files.

**Input Data:**

- `human.fasta` (fasta file containing human genome, SPECIFY WHERE THIS GEONOME CAME FROM)

**Output Data:**

- kraken2_human_db/ - Kraken2 human database directory, containing hash.k2d, opts.k2d, and taxo.k2d files 

#### 6b. Remove Human Reads

```bash
kraken2 --db kraken2_human_db \
        --gzip-compressed \
        --threads NumberOfThreads \
        --use-names \
        --output sample-kraken2-output.txt \
        --report sample-kraken2-report.tsv \
        --unclassified-out sample_HRrm_GLlblMetag.fastq \
        sample_trimmed_fastq.gz

# gzip fastq output file
gzip sample_HRrm_GLlblMetag.fastq
```

**Parameter Definitions:**

- `--db` - Specifies the directory holding the kraken2 database.
- `--gzip-compressed` - Specifies that the input fastq files are gzip-compressed.
- `--threads` - Specifies the number of parallel processing threads to use.
- `--use-names` - Specifies adding taxa names in addition to taxon IDs.
- `--output` - Specifies the name of the kraken2 read-based output file (one line per read).
- `--report` - Specifies the name of the kraken2 report output file (one line per taxa, with number of reads assigned to it).
- `--unclassified-out` - Specifies the name of the output file containing reads that were not classified, i.e non-human reads.
- `sample_trimmed.fastq.gz` - Positional argument specifying the input read file.

**Input Data:**

- kraken2_human_db/ (kraken2 human database directory, output from [Step 6a](#6a-build-kraken2-database))
- sample_trimmed.fastq.gz (filtered and trimmed sample reads, output from [Step 5a](#5a-trim-filtered-data))

**Output Data:**

- sample-kraken2-output.txt (kraken2 read-based output file (one line per read))
- sample-kraken2-report.tsv (kraken2 report output file (one line per taxa, with number of reads assigned to it))
- **sample_HRrm_GLlblMetag.fastq.gz** (filtered and trimmed sample reads with human reads removed, gzipped fastq file)


#### 6c. Compile Human Read Removal QC

```bash
multiqc --zip-data-dir \ 
        --outdir HRrm_multiqc_report \
        --filename HRrm_multiqc_GLlblMetag \
        --interactive \
        /path/to/*kraken2-report.tsv
```

**Parameter Definitions:**

- `--zip-data-dir` - Compress the data directory.
- `--outdir` – Specifies the output directory to store results.
- `--filename` – Specifies the filename prefix of results.
- `--interactive` - Force multiqc to always create interactive javascript plots.
- `/path/to/*kraken2-report.tsv` – The kraken2 output report files, provided as a positional argument.

**Input Data:**

- /path/to/*kraken2-report.tsv (kraken2 report files, output from [Step 6b](#6b-remove-human-reads))

**Output Data:**

- **HRrm_multiqc_GLlblMetag.html** (multiqc output html summary)
- **HRrm_multiqc_GLlblMetag_data.zip** (zip archive containing multiqc output data)

<br>


---

### 7. Contaminant Removal

> A major issue with low biomass data is the high potential for contamination due to the low amount of DNA extracted from the samples. Because negative control/blank samples should by theory be contaminant free, any sequence detected in the negative control is a potential contaminant. To filter out contaminants found in negative control samples that may have been due to cross contamination in the lab, we use a read mapping approach. First negative/blank control sample reads are assembled then the filtered, trimmed, and human-removed reads from each low-biomass sample are mapped to the assembled contigs from the negative/blank control samples. Reads mapping to the assembled contigs are categorized as contaminants and are therefore filtered out and thus excluded from downstream analyses.

### 7a. Assemble Contaminants

```bash
flye --meta \
     --threads NumberOfThreads \
     --out-dir /path/to/contaminant_assembly \
     --nano-raw /path/to/blank_samples/\*_HRrm_GLlblMetag.fastq.gz

# rename output
mv assembly.fasta blank-assembly.fasta
mv flye.log blank-flye.log
```

**Parameter Definitions:**

- `--meta` – Use metagenome/uneven coverage mode.
- `--threads` - Number of parallel processing threads to use.
- `--out-dir` - Specifies the output directory.
- `--nano-raw` - Specifies that input is from Oxford Nanopore regular raw reads. This adds a polishing step for error correction after the assembly is generated.

**Input Data**

- *_HRrm_GLlblMetag.fastq.gz (one or more filtered, trimmed, and HRrm reads from blank (negative control) samples, output from [Step 6b](#6b-remove-human-reads))

**Output Data**

- /path/to/contaminant_assembly/blank-assembly.fasta (assembly built from reads in blank samples in fasta format)
- blank-flye.log (flye log file)

<br>

#### 7b. Build Contaminant Index and Map Reads

```bash
# Build contaminant index
minimap2 -t NumberOfThreads \
         -a \
         -x splice \
         -d blanks.mmi \
         /path/to/contaminant_assembly/blank-assembly.fasta

# Map reads to index
minimap2 -t NumberOfThreads \
         -a \
         -x splice \
         blanks.mmi \
         sample_HRrm_GLlblMetag.fastq.gz  > sample.sam 2> sample-mapping-info.txt
```

**Parameter Definitions:**

- `-t` - Number of parallel processing threads.
- `-a` – Output in SAM format.
- `-x splice` - Specifies preset for spliced alignment of long reads.
- `-d` - Specifies the output file for the index (specific to the build contaminant index command).
- `/path/to/contaminant_assembly/blank-assembly.fasta` - Specifies the input file in fasta format, provided as a positional argument (specific to the build contaminant index command).
- `blanks.mmi` - Specifies the index file in mmi format, provided as a positional argument (specific to the map reads command).
- `/path/to/trimmed_reads/sample_HRrm_GLlblMetag.fastq.gz` - Specifies the input file in fastq format, provided as a positional argument (specific to the map reads command).
- `> sample.sam` - Redirects the output of the map reads command to a separate SAM file (specific to the map reads command).

**Input Data**

- /path/to/contaminant_assembly/blank-assembly.fasta (contaminant assembly, output from [Step 7a](#7-assemble-contaminants))
- sample_HRrm_GLlblMetag.fastq.gz (filtered, trimmed, and HRrm reads, output from [Step 6b](#6b-remove-human-reads))

**Output Data**

- blanks.mmi (contaminant index in MMI format)
- sample.sam (reads aligned to contaminant assembly in SAM format)
- sample-mapping-info.txt (minimap2 mapping log file)

#### 7c. Sort and Index Contaminant Alignments
```bash
# Sort Sam, convert to bam and create index
samtools sort --threads NumberOfThreads \
              --output sample_sorted.bam \
              sample.sam

samtools index sample_sorted.bam sample_sorted.bam.bai
```

**Parameter Definitions:**

**samtools sort**
- `--threads` - Number of parallel processing threads to use.
- `--output` - Specifies the output file for the aligned and sorted reads.
- `sample.sam` - Specifies the input SAM file, provided as a positional argument.

**samtools index**
- `sample_sorted.bam` - The input BAM file, provided as a positional argument.
- `sample_sorted.bam.bai` - The output index file, provided as a positional argument.

**Input Data:**

- sample.sam (reads aligned to contaminant assembly, output from [Step 7b](#7b-build-contaminant-index-and-map-reads))

**Output Data:**

- sample_sorted.bam (sorted mapping to contaminant assembly file)
- sample_sorted.bam.bai (index of sorted mapping to contaminant assembly file)

#### 7d. Gather Contaminant Mapping Metrics

```bash

samtools flagstat sample_sorted.bam > sample_flagstats.txt  2> sample_flagstats.log
samtools stats --remove-dups sample_sorted.bam > sample_stats.txt   2> sample_stats.log
samtools idxstats sample_sorted.bam  > sample_idxstats.txt 2> sample_idxstats.log
```

**Parameter Definitions:**

- `flagstat` - Positional argument specifying the program for counting the number of alignments for each SAM FLAG type.
- `stats` - Positional argument specifying the program for producing comprehensive statistics from the alignment file.
- `idxstats` - Positional argument specifying the program for producing contig alignment summary statistics.
- `--remove-dups` - Excludes reads marked as duplicates from the comprehensive statistics.
- `sample_sorted.bam` - Positional argument specifying the input BAM file.
- `> sample_flagstats.txt` - Redirects the flagstat standard output to a text file.
- `2> sample_flagstats.log` - Redirects the flagstat standard error to a log file.
- `> sample_stats.txt` - Redirects the stats standard output to a text file.
- `2> sample_stats.log` - Redirects the stats standard error to a log file.
- `> sample_idxstats.txt` - Redirects the idxstats standard output to a text file.
- `2> sample_idxstats.log` - Redirects the idxstats standard error to a log file.

**Input Data:**

- sample_sorted.bam (sorted mapping to contaminant assembly file, output from [Step 7c](#7c-sort-and-index-contaminant-alignments))
- sample_sorted.bam.bai (index of sorted mapping to contaminant assembly file, output from [Step 7c](#7c-sort-and-index-contaminant-alignments))

**Output Data:**

- sample_flagstats.txt (SAM FLAG counts)
- sample_flagstats.log (log file containing the flagstat standard error)
- sample_stats.txt (comprehensive alignment statistics)
- sample_stats.log (log file containing the stats standard error)
- sample_idxstats.txt (contig alignment summary statistics)
- sample_idxstats.log (log file containing the idxstats standard error)

#### 7e. Generate Decontaminated Read Files
```bash
# Retain reads that do not map to contaminants
samtools fastq -t -f 4 -o sample_decontam_GLlblMetag.fastq.gz -0 sample_decontam_GLlblMetag.fastq.gz sample_sorted.bam 
```

**Parameter Definitions:**

- `fastq` - Positional argument specifying the program for generating fastq files from a SAM/BAM file.
- `-t` - Copy RG, BC, and QT tags to the FASTQ header line.
- `-f 4` - Only retain unmapped reads that have been marked with the SAM "segment unmapped" FLAG (4).
- `-o sample_decontam_GLlblMetag.fastq.gz` - Send reads flagged as either read1 or read2 to the named file (.gz ending ensures compressed output)
- `-0 sample_decontam_GLlblMetag.fastq.gz` - Send reads flagged as both read1 and read2 or neither to the same named file
- `sample_sorted.bam` - Positional argument specifying the input BAM file.

**Input Data:**

- sample_sorted.bam (sorted mapping to contaminant assembly file, output from [Step 7c](#7c-sort-and-index-contaminant-alignments))

**Output Data:**

- **sample_decontam_GLlblMetag.fastq.gz** (filtered, trimmed, and HRrm sample reads with contaminants removed in fastq format)

#### 7f. Contaminant Removal QC

```bash
NanoPlot --only-report \
         --prefix sample_noblank_ \
         --outdir /path/to/decontam_nanoplot_output \
         --threads NumberOfThreads \
         --fastq \
         sample_decontam_GLlblMetag.fastq.gz
```

**Parameter Definitions:**

- `--only-report` - Output only the report files.
- `--prefix` - Adds a sample specific prefix to the name of each output file.
- `--outdir` – Specifies the output directory to store results.
- `--threads` - Number of parallel processing threads to use.
- `--fastq` - Specifies that the input data is in fastq format.
- `sample_decontam_GLlblMetag.fastq.gz` – The input reads, specified as a positional argument.

**Input Data:**

- sample_decontam_GLlblMetag.fastq.gz (filtered, trimmed, and HRrm sample reads with all contaminants removed, output from [Step 7e](#7e-generate-decontaminated-read-files))

**Output Data:**

- **/path/to/decontam_nanoplot_output/sample_decontam_NanoPlot-report_GLlblMetag.html** (NanoPlot html summary)
- /path/to/decontam_nanoplot_output/sample_decontam_NanoPlot_\<date\>_\<time\>.log (NanoPlot log file)
- /path/to/decontam_nanoplot_output/sample_decontam_NanoStats.txt (text file containing basic statistics)


#### 7g. Compile Contaminant Removal QC

```bash
multiqc --zip-data-dir \ 
        --outdir decontam_multiqc_report \
        --filename decontam_multiqc_GLlblMetag \
        --interactive \
        /path/to/decontam_nanoplot_output/
```

**Parameter Definitions:**

- `--zip-data-dir` - Compress the data directory.
- `--outdir` – Specifies the output directory to store results.
- `--filename` – Specifies the filename prefix of results.
- `--interactive` - Force multiqc to always create interactive javascript plots.
- `/path/to/decontam_nanoplot_output/` – The directory holding the output data from the NanoPlot run, provided as a positional argument.

**Input Data:**

- /path/to/decontam_nanoplot_output/*decontam_NanoStats.txt (NanoPlot output data, output from [Step 7f](#7f-contaminant-removal-qc))

**Output Data:**

- **decontam_multiqc_GLlblMetag.html** (multiqc output html summary)
- **decontam_multiqc_GLlblMetag_data.zip** (zip archive containing multiqc output data)

<br>

---

### 8. Host Read Removal

If the samples were derived from a host organism other than human, potential host reads should be identified and removed. This step is optional.

#### 8a. Build Kraken2 Database

> **Note:** It is recommended to use NCBI genome files with kraken2 because sequences not downloaded from 
NCBI may require explicit assignment of taxonomy information before they can be used to build the 
database, as mentioned in the [Kraken2 Documentation](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown).

```bash
# Download NCBI taxonomic information 
kraken2-build --download-taxonomy --db kraken2-${hostname}-db/

# Add genomic sequences to your database's genomic library
kraken2-build --add-to-library ${hostname}.fasta --db kraken2-${hostname}-db/ \
              --no-masking --kmer-length 35 --minimizer-length 31

# Build the database
kraken2-build --build --db kraken2-${hostname}-db/

# Clean up intermediate files
kraken2-build --clean --db kraken2-${hostname}-db/
```
**Parameter Definitions:**
- `--download-taxonomy` - Instructs kraken2-build to download the NCBI taxonomic information.
- `--db` - Specifies the name of the directory for the kraken2 database
- `--add-to-library` - Instructs kraken2-build to add the contents of a file (`${hostname}.fasta`) to the kraken2 DB library
- `--no-masking` - Disables masking of low-complexity sequences. For additional 
                   information see the [kraken documentation for masking](https://github.com/DerrickWood/kraken2/wiki/Manual#masking-of-low-complexity-sequences).
- `--build` - Instructs kraken2-build to build the kraken2 DB from the library files
- `--clean` - Instructs kraken2-build to remove unneeded intermediate files.
- `{$hostname}` - Specifies the name of the host organism used to uniquely identify the kraken2 database

**Input Data:**

- `${hostname}.fasta` (fasta file containing host genome)

**Output Data:**

- kraken2_${hostname}_db/ - Kraken2 host database directory, containing hash.k2d, opts.k2d, and taxo.k2d files 


#### 8b. Remove Host Reads

```bash
kraken2 --db kraken2_host_db \
        --gzip-compressed \
        --threads NumberOfThreads \
        --use-names \
        --output sample-kraken2-output.txt \
        --report sample-kraken2-report.tsv \
        --unclassified-out sample_HostRm_GLlblMetag.fastq \
        sample_decontam_GLlblMetag.fastq.gz

# gzip fastq output file
gzip sample_HostRm_GLlblMetag.fastq
```

**Parameter Definitions:**

- `--db` - Specifies the directory holding the kraken2 database.
- `--gzip-compressed` - Specifies that the input fastq files are gzip-compressed.
- `--threads` - Number of parallel processing threads to use.
- `--use-names` - Specifies adding taxa names in addition to taxon IDs.
- `--output` - Specifies the name of the kraken2 read-based output file (one line per read).
- `--report` - Specifies the name of the kraken2 report output file (one line per taxa, with number of reads assigned to it).
- `--unclassified-out` - Specifies the name of the output file containing reads that were not classified, i.e non-human reads.
- `sample_decontam_GLlblMetag.fastq.gz` - Positional argument specifying the input read file.

**Input Data:**

- kraken2_host_db/ (kraken2 host database directory, output from [Step 8a](#8a-build-kraken2-database))
- sample_decontam_GLlblMetag.fastq.gz (filtered, trimmed, HRrm and contaminant-removed sample reads, output from [Step 7e](#7e-generate-decontaminated-read-files))

**Output Data:**

- sample-kraken2-output.txt (kraken2 read-based output file (one line per read))
- sample-kraken2-report.tsv (kraken2 report output file (one line per taxa, with number of reads assigned to it))
- **sample_HostRm_GLlblMetag.fastq.gz** (filtered, trimmed, HRrm and contaminant-removed sample reads with all host reads removed, gzipped fastq file)


#### 8c. Compile Host Read Removal QC

```bash
multiqc --zip-data-dir \ 
        --outdir HostRm_multiqc_report \
        --filename HostRm_multiqc_GLlblMetag \
        --interactive \
        /path/to/*kraken2-report.tsv
```

**Parameter Definitions:**

- `--zip-data-dir` - Compress the data directory.
- `--outdir` – Specifies the output directory to store results.
- `--filename` – Specifies the filename prefix of results.
- `--interactive` - Force multiqc to always create interactive javascript plots.
- `/path/to/*kraken2-report.tsv` – The kraken2 output report files, provided as a positional argument.

**Input Data:**

- /path/to/*kraken2-report.tsv (kraken2 report files, output from [Step 8b](#8b-remove-host-reads))

**Output Data:**

- **HostRm_multiqc_GLlblMetag.html** (multiqc output html summary)
- **HostRm_multiqc_GLlblMetag_data.zip** (zip archive containing multiqc output data)

<br>

---

### 9. R Environment Setup

> Taxonomy bar plots, heatmaps and feature decontamination with decontam are performed in R.

#### 9a. Load libraries

```R
library(decontam)
library(phyloseq)
library(tidyverse)
library(pheatmap)
library(pavian)
```

#### 9b. Define Custom Functions

##### get_last_assignment()
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

##### mutate_taxonomy()
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

##### process_kaiju_table()
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


##### merge_kraken_reports()
<details>
  <summary>merge and process multiple kraken outputs to one species table</summary>

  ```R
  library(pavian)

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
    species_table <- species_table[,-(which(colnames(species_table) == "species"))]
    species_table <- as.matrix(species_table)
    
    return(species_table)
  }
  ```
  **Custom Functions Used:**
  - [read_reports()]()


  **Function Parameter Definitions:**
  - `reports_dir` - path to a directory containing kraken2 reports 

  **Returns:** a kraken species count matrix, `species_table`, with samples and species as columns and rows, respectively.

</details>

##### get_abundant_features()
<details>
  <summary>Find abundant features based on the sum of feature values</summary>
  
  ```R
  get_abundant_features <- function(mat, cpm_threshold = 1000){
  
    features <- rowSums(mat) %>% sort()
    
    abund_features <- features[features > cpm_threshold] %>% names
    
    abund_features.m <- mat[abund_features, ]
    
    return(abund_features.m)
  }
  ```
  **Function Parameter Definitions:**
  - `mat` - a feature count matrix with features as rows and samples as columns
  - `cpm_threshold = 1000` - threshold to identify abundant features

  **Returns:** a species relative abundance matrix, `abund_features.m`, with samples and species as rows and columns, respectively.
</details>

##### count_to_rel_abundance()
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


##### filter_rare()
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
      mutate( across( where(is.numeric), function(x) (x/sum(x, na.rm = TRUE))*100 ) )
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

##### group_low_abund_taxa()
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
      message(glue::glue("Rare taxa were not grouped. please provide a higher 
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

##### make_plot()
<details>
  <summary>create bar plot of relative abundance</summary>

  ```R
  # Make bar plot
  make_plot <- function(abund_table, metadata, custom_palette, publication_format,
                        samples_column="Sample_ID", prefix_to_remove="barcode"){
  
    abund_table_wide <- abund_table %>%
        as.data.frame %>%
        rownames_to_column(samples_column) %>%
        inner_join(metadata) %>%
        select(!!!colnames(metadata), everything()) %>%
        mutate(!!samples_column := !!sym(samples_column) %>% str_remove(prefix_to_remove))
        
      
    abund_table_long <- abund_table_wide  %>%
        pivot_longer(-colnames(metadata), 
                     names_to = "Species",
                     values_to = "relative_abundance")
      
    p <- ggplot(abund_table_long, mapping = aes(x = !!sym(samples_column), 
                                                y = relative_abundance, fill = Species)) +
         geom_col() +
         scale_fill_manual(values = custom_palette) + 
         labs(x=NULL, y="Relative Abundance (%)") + 
         publication_format

    return(p)
  }
  ```

  **Function Parameter Definitions:**
  - `abund_table` - a relative bundance dataframe with rows summing to 100%
  - `metadata` - a metadata dataframe with samples as row and columns describing each sample
  - `custom_palette` - a vector of strings specifying a custom color palette for coloring plots
  - `publication_format` - a ggplot::theme object specifying a custom theme for plotting
  - `samples_column` - a character column specifying the column in `metadata` holding sample names, default is "Sample_ID"
  - `prefix_to_remove` - a string specifying a prefix or any character set to remove from sample names, default is "barcode"

  **Returns:** a relative abundance stacked bar plot, `p`

</details>

##### make_barplot()
<details>
  <summary>Creates barplots from a feature table file</summary>
  
  ```R
  make_barplot <- function(metadata_table_file, feature_table_file, 
                           feature_column = "species", samples_column = "sample_id", group_column = "group", 
                           output_prefix, assay_suffix = "_GLlblMetag",
                           publication_format, custom_palette) {
    # Prepare feature table
    feature_table <- read_csv(feature_table_file)
    rownames(feature_table) <- feature_table[[1]]
    feature_table <- feature_table[, -1]

    # Prepare metadata
    metadata <- read_delim(metdata_file, delim = ",") %>% as.data.frame
    row.names(metadata) <- metadata[, samples_column]

    # compute abundances from counts
    abund_table <- count_to_rel_abundance(feature_table)
    
    # create plot
    p <- make_plot(abund_table, metadata, custom_palette, publication_format, samples_column) +
         facet_wrap(~Description, nrow=1, scales = "free_x")

    number_of_species <- p$data$Species %>% unique() %>% length()
    # Don't save legend if the number of species to plot is gsreater than 30
    if(number_of_species > 30) {
      p <- p + theme(legend.position = "none")
    }

    return(p)

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
  - `assay_suffix` - a character string specifying the GeneLab assay suffix (default: "_GLlblMetag")
  - `publication_format` - a ggplot::theme object specifying a custom theme for plotting, from [Step 9c](#9c-set-global-variables)
  - `custom_palette` - a vector of strings specifying a custom color palette for coloring plots, from [Step 9c](#9c-set-global-variables)

  **Returns:** a relative abundance stacked bar plot, `p`

</details>

##### make_heatmap()
<details>
  <summary>Creates heatmaps from a feature table file</summary>
  
  ```R
  make_heatmap <- function(metadata, species_gene_table, 
                           samples_column = "sample_id", group_column = "group", 
                           output_prefix, assay_suffix = "_GLlblMetag",
                           custom_palette) {
    # Prepare feature table
    # feature_table <- read_csv(feature_table_file) %>% as.data.frame
    # rownames(feature_table) <- feature_table[[1]]
    # feature_table <- feature_table[, -1] %>% as.matrix()

    # # Prepare metadata
    # metadata <- read_delim(metadata_file, delim = ",") %>% as.data.frame
    # row.names(metadata) <- metadata[, samples_column]

    # # Get common samples and re-arrange feature table and metadata
    # common_samples <- intersect(colnames(feature_table), rownames(metadata))
    # feature_table <- feature_table[, common_samples]
    # metadata <- metadata[common_samples, ]
    # metadata <- metadata %>% arrange(!!sym(group_column))

    # Create column annotation
    col_annotation <- as.data.frame(metadata)[, group_column, drop = FALSE]

    # Calculate output plot width and height
    number_of_samples <- ncol(feature_table)
    width <- 1 * number_of_samples
    number_of_features <- nrow(feature_table)
    height <- 0.2 * number_of_features

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
  }
  ```
  **Function Parameter Definitions:**
  - `metadata_file` - path to a file with samples as rows and columns describing each sample
  - `feature_table_file` - path to a tab separated samples feature table i.e. species/functions 
                           table with species/functions as the first column and samples as other columns.
  - `samples_column` - a character string specifying the column in `metadata` holding sample names, default: "sample_id"
  - `group_column` - a character string specifying the column in `metadata` used to facet/group plots, default: "group"
  - `output_prefix` - a character string specifying the unique name to add to the output file names 
                      used to denote the data type/source, for example "unfiltered-kaiju_species"
  - `assay_suffix` - a character string specifying the GeneLab assay suffix (default: "_GLlblMetag")
  - `custom_palette` - a vector of strings specifying a custom color palette for coloring plots, from [Step 9c](#9c-set-global-variables)

</details>

##### run_decontam()
<details>
  <summary>Feature table decontamination with decontam</summary>

  ```R
  run_decontam <- function(feature_table, metadata, contam_threshold=0.1, 
                           prev_col = NULL, freq_col = NULL, ntc_name = "TRUE") {

    # retain metadata for only the samples present in the input feature table
    sub_metadata <- metadata[colnames(feature_table), ]
    # Modify NTC concentration
    # Often times the user may set the NTC concentration to zero because they think nothing 
    # should be in the negative control but decontam fails if the value is set to zero.
    # To prevent decontam from failing, we replace zero with a very small concentration value
    # 0.0000001
    if (!is.null(freq_col)) {

      sub_metadata <- sub_metadata %>%
        mutate(!!freq_col:=map_dbl(!!sym(freq_col), .f = function(conc) {
              if(conc == 0) return(0.0000001) else return(conc) 
            } 
          )
        )
      sub_metadata[, freq_col] <- as.numeric(sub_metadata[, freq_col])

    }

    # Create phyloseq object
    ps <- phyloseq(otu_table(feature_table, taxa_are_rows = TRUE), sample_data(sub_metadata))

    # In our phyloseq object, `prev_col` is the sample variable that holds the negative 
    # control information. We'll summarize the data as a logical variable, with TRUE for control 
    # samples, as that is the form required by isContaminant.
    # The line below assumes that control samples will always be named "Control_Sample"
    # in the `prev_col`.
    sd <- as.data.frame(sample_data(ps)) # Extract sample metadata
    sd[, "is.neg"] <- 0 # Initialize
    sd[, "is.neg"] <- sample_data(ps)[[prev_col]] == ntc_name # Assign boolean value
    sample_data(ps) <- sd

    # Run Decontam 
    if (!is.null(freq_col) && !is.null(prev_col)) {
      # Run decontam in both prevalence and frequency modes
      contamdf <- isContaminant(ps, neg="is.neg", conc=freq_col, threshold=contam_threshold) 
    } else if(!is.null(freq_col)) {
      # Run decontam in frequency mode
      contamdf <- isContaminant(ps, conc=freq_col, threshold=contam_threshold) 
    } else if(!is.null(prev_col)){
      # Run decontam in prevalence mode
      contamdf <- isContaminant(ps, neg="is.neg", threshold=contam_threshold)
    } else {
      cat("Both freq_col and prev_col cannot be set to NULL.\n")
      cat("Please supply either one or both column names in your metadata")
      cat("for frequency and prevalence based analysis, respectively\n")
      stop()
    }            
    return(contamdf)
  }
  ```

  **Function Parameter Definitions:**
  - `metadata` - a metadata dataframe with samples as row and columns describing each sample
  - `feature_table` -  feature [species, functions etc.] matrix to decontaminate with sample names as column and features as row
  - `prev_col` - a character column in metadata to be used for prevalence based analysis. Controls in this column should always be names "Control_Sample"
  - `freq_col` - a numeric column in metadata to be used for frequency based analysis
  - `contam_threshold` -  the probability threshold below which (strictly less than) the null-hypothesis 
                          (not a contaminant) should be rejected in favor of the alternate hypothesis (contaminant).

  **Returns:** dataframe, `contamdf`, containing detailed decontam results
</details>

##### feature_decontam() 
<details>
  <summary>decontaminate a feature table</summary>
  
  ```R
  library(tidyverse)
  library(glue)

  feature_decontam <- function(metadata_file, feature_table_file, 
                               feature_column = "Species", samples_column = "sample_id",
                               prevalence_column = "NTC", ntc_name = "TRUE", 
                               frequency_column = "concentration", 
                               threshold = 0.1, classification_method, 
                               output_prefix, assay_suffix = "_GLlblMetag") {
    # Prepare feature table
    feature_table <- read_csv(feature_table_file) %>%  as.data.frame
    rownames(feature_table) <- feature_table[[1]]
    feature_table <- feature_table[, -1]  %>% as.matrix()

    # Prepare metadata
    metadata <- read_csv(metadata_file) %>% as.data.frame
    row.names(metadata) <- metadata[, samples_column]

    # Run decontam
    contamdf <- run_decontam(feature_table, metadata, threshold, prev_col, freq_col, ntc_name) 

    contamdf <- as.data.frame(contamdf) %>% rownames_to_column(feature_column)

    # Write decontaminated feature table and decontam's primary results
    outfile <- glue("{output_prefix}{classification_method}_decontam_results{assay_suffix}.csv")
    write_csv(x = contamdf, file = outfile)

    # Get the list of contaminants identified by decontam
    contaminants <- contamdf %>%
                    filter(contaminant == TRUE) %>%
                    pull(!!sym(feature_column))

    # Drop contaminants(s) if detected
    if(length(contaminants) > 0){
      
      # Drop contaminant features identified by decontam
      decontaminated_table <- feature_table %>%
        as.data.frame %>%
        rownames_to_column(feature_column) %>%
        filter(str_detect(!!sym(feature_column),
                          pattern = str_c(contaminants,
                                          collapse = "|"),
                          negate = TRUE))

      rownames(decontaminated_table) <- decontaminated_table[[feature_column]]
      decontaminated_table <- decontaminated_table[,-1] %>% as.matrix

      outfile <- glue("{output_prefix}{classification_method}_decontam_species_table{assay_suffix}.csv")
      write_csv(x = decontaminated_table, file = outfile)

      return(decontaminated_table)

    } else {
      message("No contaminants were detected by Decontam")
      return(NULL)
    }
  }
  ```
  **Custom Functions Used:**
  - [run_decontam()](#run_decontam)

  **Function Parameter Definitions:**
  - `metadata_file` - path to a file with samples as rows and columns describing each sample
  - `feature_table_file` - path to a tab separated samples feature table i.e. species/functions 
                           table with species/functions as the first column and samples as other columns.
  - `feature_column` - a character string containing the feature column name in the feature table ['Species', 'species', 'KO_ID'].
  - `samples_column` - a character string specifying the column in `metadata` holding sample names, default: "sample_id"
  - `frequency_column` - a character string specifying the column in `metadata` to use for frequency based analysis, default: "concentration"
  - `prevalence_column` - a character string specifying the column in `metadata` to use for prevalence based analysis, default: "NTC"
  - `ntc_name` - a character string specifying the value in the prevalence column for all negative template control samples, default: "TRUE"
  - `threshold` - a number between 0 and 1 specfying the decontam threshold for both prevalence and frequency based analyses. default: 0.1
  - `output_prefix` - a character string specifying the unique name to add to the output file names 
                      used to denote the data type/source, for example "unfiltered-kaiju_species"
  - `classification_method` - a character string specifying the tool used to generate the classifications ['kaiju', 'kraken2', 'metaphlan', 'contig-taxonomy', 'gene-taxonomy', 'gene-function']
  - `assay_suffix` - a character string specifying the GeneLab assay suffix (default: "_GLlblMetag")

  **Output Data:**
  - {classification_method}_decontam_species_table_GLlblMetag.csv - decontaminated feature table file
  - {classification_method}_decontam_results_GLlblMetag.csv - Decontam results file

  **Returns:** dataframe, `decontaminated_table`, containing the decontaminated feature table

</details>

##### process_taxonomy()
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

##### fix_names()
<details>
  <summary>clean taxonomy names</summary>

  ```R
  fix_names<- function(taxonomy,stringToReplace="Othe",suffix=";Other"){
    
    for(index in seq_along(stringToReplace)){

      for (taxa_index in seq_along(taxonomy)) {    
        # Get the row indices of the current taxonomy columns
        # with rows matching the sting in `stringToReplace`
        indices <- grep(x = taxonomy[,taxa_index], pattern = stringToReplace)
        # Replace the value in that row with the value in the adjacent cell concated with `suffix`
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


##### read_assembly_coverage_table()
<details>
  <summary>Read Assembly-based coverage annotation table</summary>

  ```R
  read_assembly_coverage_table <- function(file_name, sample_names){
  
    df <- read_delim(file = file_name, delim = "\t", comment = "#")

    # Subset taxoxnomy portion (domain:species) of input table
    # and replace empty/Na domain assignments with "Unclassified"
    taxonomy_table <- df %>%
      select(domain:species) %>%
      mutate(domain=replace_na(domain, "Unclassified"))
    
    # Subset count table
    counts_table <- df %>% select(!!any_of(sample_names))

    # Mutate taxonomy mames
    taxonomy_table  <- process_taxonomy(taxonomy_table)
    taxonomy_table <- fix_names(taxonomy_table, "Other", ";_")

    # Column bind taxonomy dataframe with species count dataframe
    df <- bind_cols(taxonomy_table, counts_table)
    
    return(df)
  }
  ```

  **Custom Functions Used:**
  [process_taxonomy()](#process_taxonomy)  
  [fix_names()](#fix_names)   

  **Function Parameter Definitions:**

  - `file_name` - path to contig taxonomy assignment file to be read
  - `sample_names` - string of samples names to keep in the final dataframe

  **Returns:** dataframe, `df`, containing cleaned taxonomy names and sample species count

</details>


##### get_sample_names()
<details>
  <summary>retrieve sample names for which assemblies were generated</summary>

  ```R
  get_sample_names <- function (assembly_summary) {
    # Read in table and drop columns were all rows are NA
    overview_table <-  read_delim(file = assembly_summary, delim = "\t", comment = "#") %>%
                        select(where( ~all(!is.na(.)) )) 

    col_names <- names(overview_table) %>% str_remove_all("-assembly")
    sample_order <- col_names[-1] %>% sort()

    return(sample_order)
  }
  ```
  **Function Parameter Definitions:**

  - `assembly_summary` - path to assembly summary file

  **Returns:** a character vector, `sample_order`, of sorted sample names

</details>


#### 9c. Set global variables

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

## Read-based Processing

### 10. Taxonomic Profiling Using Kaiju

#### 10a. Build Kaiju Database

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


#### 10b. Kaiju Taxonomic Classification

```bash
kaiju -f kaiju-db/nr_euk/kaiju_db_nr_euk.fmi \
      -t kaiju-db/nodes.dmp \
      -z NumberOfThreads \
      -E 1e-05 \
      -i /path/to/sample_decontam_GLlblMetag.fastq.gz \
      -o sample_kaiju.out
```

**Parameter Definitions:**

- `-f` - Specifies the path to the kaiju database index file (.fmi).
- `-t` - Specifies the path to the kaiju taxonomy hierarchy file (nodes.dmp).
- `-z` - Number of parallel processing threads to use.
- `-E` - Specifies the minimum E-value to use for filter matches (an E-value of 1e-05 means that there's a 0.001% chance that the matches identified occurred randomly).
- `-i` - Specifies path to the input file.
- `-o` - Specifies the name of the output file.

**Input Data:**

- kaiju-db/nr_euk/kaiju_db_nr_euk.fmi (FM-index file containing the main Kaiju database index, output from [Step 10a](#10a-build-kaiju-database))
- kaiju-db/nodes.dmp (kaiju taxonomy hierarchy nodes file, output from [Step 10a](#10a-build-kaiju-database))
- sample_decontam_GLlblMetag.fastq.gz or sample_HostRm_GLlblMetag.fastq.gz (filtered and trimmed sample reads with both 
    contaminants and human reads (and optionally host reads) removed, gzipped fastq file, 
    output from [Step 7e](#7e-generate-decontaminated-read-files) or [Step 8b](#8b-remove-host-reads))

**Output Data:**

- sample_kaiju.out (kaiju output file)

#### 10c. Compile Kaiju Taxonomy Results

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

- kaiju-db/nodes.dmp (kaiju taxonomy hierarchy nodes file, output from [Step 10a](#10a-build-kaiju-database))
- kaiju-db/names.dmp (kaiju taxonomy names file, output from [Step 10a](#10a-build-kaiju-database))
- *kaiju.out (kaiju output files, output from [Step 10b](#10b-kaiju-taxonomic-classification))

**Output Data:**

- merged_kaiju_table.tsv (compiled kaiju summary table at the species level)

#### 9d. Convert Kaiju Output To Krona Format

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
- kaiju-db/names.dmp (kaiju taxonomy names file, output from [Step 10a](#10a-build-kaiju-database))
- kaiju-db/nodes.dmp (kaiju taxonomy hierarchy nodes file, output from [Step 10a](#9a-build-kaiju-database))
- sample_kaiju.out (kaiju output file, output from [Step 10b](#10b-kaiju-taxonomic-classification))

**Output Data:**

- sample.krona (krona formatted kaiju output)

#### 10e. Compile Kaiju Krona Reports

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

**find**

- `-type f` -  Specifies that the type of file to find is a regular file.
- `-name "*.krona"` - Specifies to find files ending with the .krona suffix.  

**sort**

- `-u` - Specifies to perform a unique sort.
- `-V` - Specifies to perform a mixed type of sorting with names containing numbers within text.
- `> krona_files.txt` - Redirects the sorted list to a separate text file.

**basename**

- `-a` - Support multiple arguments and treat each as a file name.
- `-s '.krona'` - Remove trailing '.krona' suffix.

**paste**

- `-d','` - Paste both krona and sample files together line by line delimited by comma ','.

**ktImportText**

- `-o` - Specifies the compiled output html file name.
- `${KTEXT_FILES[*]}` - An array positional argument with the following content: 
                        sample_1.krona,sample_1 sample_2.krona,sample_2 ... sample_n.krona,sample_n.

**Input Data:**
- *.krona (all sample .krona formatted files, output from [Step 10d](#10d-convert-kaiju-output-to-krona-format)) 

                      
**Output Data:**

- krona_files.txt (sorted list of all *.krona files)
- sample_names.txt (sorted list of all sample names)
- **kaiju-report_GLlblMetag.html** (compiled krona html report containing all samples)


#### 10f. Create Kaiju Species Count Table

```R
library(tidyverse)
feature_table <- process_kaiju_table(file_path="merged_kaiju_table_GLlblMetag.tsv")
table2write <- feature_table  %>%
                as.data.frame %>%
                rownames_to_column("Species")
write_csv(x = table2write, file = "kaiju_species_table_GLlblMetag.csv")
```

**Custom Functions Used:**

- [process_kaiju_table()](#process_kaiju_table)

**Parameter Definitions:**

- `file_path` - path to compiled kaiju table at the species taxon level
- `x`  - feature table dataframe to write to file
- `file` - path to where to write kaiju count table per sample

**Input Data:**

- merged_kaiju_table_GLlblMetag.tsv (compiled kaiju table at the species taxon level, from [Step 10c](#10c-compile-kaiju-taxonomy-results))

**Output Data:**

- **kaiju_species_table_GLlblMetag.csv** (kaiju species count table in csv format)


#### 10g. Filter Kaiju Species Count Table

```R
library(tidyverse)

input_file <- "kaiju_species_table_GLlblMetag.csv"
output_file <- "kaiju_filtered_species_table_GLlblMetag.csv"
threshold <- 0.5

# string used to define non-microbial taxa
non_microbial <- "UNCLASSIFIED|Unclassifed|unclassified|Homo sapien|cannot|uncultured|unidentified"

# read in feature table
feature_table <- read_csv(input_file) %>% as.data.frame
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

write_csv(x = table2write, file = output_file)
```

**Custom Functions Used:**

- [group_low_abund_taxa()](#group_low_abund_taxa)

**Parameter Definitions:**

- `threshold` - threshold for filtering out rare taxa, a percentage between 0 and 100.
- `output_file` - output filename
- `input_file` - input filename

**Input Data:**

- kaiju_species_table_GLlblMetag.csv (path to kaiju species table from [Step 9f](#9f-create-kaiju-species-count-table))

**Output Data:**

- **kaiju_filtered_species_table_GLlblMetag.csv** - a file containing the filtered species table

---

#### 10h. Taxonomy barplots

```R
library(tidyverse)

species_table_file <- "kaiju_species_table_GLlblMetag.csv"
filtered_species_table_file <- "kaiju_filtered_species_table_GLlblMetag.csv"
metadata_file <- "/path/to/sample/metadata"
number_samples <- 10 

# set width based on number of samples, with a cap at 50 inches
plot_width <- 2 * number_samples
if(plot_width > 50) { plot_width = 50 }

p <- make_barplot(metadata_file = metadata_file, feature_table_file = species_table_file, 
                  feature_column = "Species", samples_column = "sample_id", group_column = "group",
                  publication_format = publication_format, custom_palette = custom_palette)

ggsave(filename = "kaiju_unfiltered_species_barplot_GLlblMetag.png", plot = p,
       device = "png", width = plot_width, height = 10, units = "in", dpi = 300, limitsize = FALSE)

# Save static unfiltered plot
p <- make_barplot(metadata_file = metadata_file, feature_table_file = filtered_species_table_file, 
                  feature_column = "Species", samples_column = "sample_id", group_column = "group",
                  publication_format = publication_format, custom_palette = custom_palette)

# Save interactive unfilterted plot
htmlwidgets::saveWidget(ggplotly(p), glue("kaiju_unfiltered_species_barplot_GLlblMetag.html"), selfcontained = TRUE)

# Save static filtered plot
ggsave(filename = glue("kaiju_filtered_species_barplot_GLlblMetag.png"), plot = p,
      device = 'png', width = plot_width, height = 10, units = 'in', dpi = 300, limitsize = FALSE)

# Save interactive filtered plot
htmlwidgets::saveWidget(ggplotly(p), glue("kaiju_filtered_species_barplot_GLlblMetag.html"), selfcontained = TRUE)
```

**Custom Functions Used:**
- [make_barplot](#make_barplot)

**Parameter Definitions:**

- `species_table_file` - a file containing the species count table
- `filtered_species_table_file` - a file containing the filtered species count table
- `metadata_file` - a file containing group information for each sample in the species count files
- `number_samples` - the total number of samples in the species count files, adjust based on number of input samples.

**Input Data:**

- `kaiju_species_table_GLlblMetag.csv` (a file containing the species count table, output from [Step 9f](#9f-create-kaiju-species-count-table))
- `kaiju_filtered_species_table_GLlblMetag.csv` (a file containing the filtered species count table, output from [Step 9g](#9g-filter-kaiju-species-count-table))
- `/path/to/sample/metadata` (a file containing sample-wise metadata, mapping samplenames to group metadata)


**Output Data:**

- kaiju_unfiltered_species_barplot_GLlblMetag.png (taxonomy barplot without filtering)
- **kaiju_unfiltered_species_barplot_GLlblMetag.html** (interactive taxonomy barplot without filtering)
- kaiju_filtered_species_barplot_GLlblMetag.png (taxonomy barplot after filtering rare and non-microbial taxa)
- **kaiju_filtered_species_barplot_GLlblMetag.html** (interactive taxonomy barplot after filtering rare and non-microbial taxa)


#### 10i. Feature decontamination

> Feature (species) decontamination with decontam. Decontam is an R package that statistically identifies contaminating features in a feature table

```R
library(tidyverse)
library(decontam)
library(phyloseq)

feature_table_file <- "kaiju_filtered_species_table_GLlblMetag.csv"
metadata_table <- "/path/to/sample/metadata"
number_samples <- NumberOfSamples # integer indicating how many samples are in the file

# set width based on number of samples, with a cap at 50 inches
plot_width <- 2 * number_samples
if(plot_width > 50) { plot_width = 50 }

decontaminated_table <- feature_decontam(metadata_file = metadata_table, 
                                         feature_table_file = feature_table_file, 
                                         feature_column = "species", 
                                         samples_column = "sample_id",
                                         prevalence_column = "NTC", 
                                         ntc_name = "TRUE", 
                                         frequency_column = "concentration", 
                                         threshold = 0.1, 
                                         classification_method = "kaiju", 
                                         output_prefix = "", 
                                         assay_suffix = "_GLlblMetag")

# Convert count matrix to relative abundance matrix
decontaminated_species_table <- count_to_rel_abundance(decontaminated_table)

# Make plot after filtering out contaminants
p <- make_plot(decontaminated_species_table, metadata, custom_palette, publication_format)

ggsave(filename = "kaiju_decontam_species_barplot_GLlblMetag.png", plot = p,
         device = "png", width = plot_width, height = 10, units = "in", dpi = 300)

# Save interactive filtered plot
htmlwidgets::saveWidget(ggplotly(p), glue("kaiju_decontam_species_barplot_GLlblMetag.html"), selfcontained = TRUE)
```

**Custom Functions Used:**
- [feature_decontam()](#feature_decontam)
- [make_plot()](#make_plot)
- [count_to_rel_abundance()](#count_to_rel_abundance)

**Parameter Definitions:**
- `metadata_table` - path to a file with samples as rows and columns describing each sample
- `feature_table_file` - path to a tab separated samples feature table i.e. species/functions 
                         table with species/functions as the first column and samples as other columns.
- `number_samples` - the total number of samples in the species count files, adjust based on number of input samples in the feature_table_file

**Input Data:**

- `kaiju_filtered_species_table_GLlblMetag.csv`(path to filtered species count per sample, output from [Step 9g](#9g-filter-kaiju-species-count-table))
- `/path/to/sample/metadata` (a file containing sample-wise metadata, mapping samplenames to group metadata)

**Output Data:**

- **kaiju_decontam_results_GLlblMetag.csv** (decontam's result table, output from [feature_decontam() function](#feature_decontam))
- **kaiju_decontam_species_table_GLlblMetag.csv** (decontaminated species table, output from [feature_decontam() function](#feature_decontam))
- **kaiju_decontam_species_barplot_GLlblMetag.png** (barplot after filtering out contaminants)
- **kaiju_decontam_species_barplot_GLlblMetag.html** (barplot after filtering out contaminants)

<br>

---

### 11. Taxonomic Profiling Using Kraken2

#### 11a. Download Kraken2 Database

```bash 
## Download all microbial (including eukaryotes) - https://benlangmead.github.io/aws-indexes/k2

# Downloading and building kraken2's pluspfp database which contains the standard database (Refseq archaea, bacteria, viral, plasmid, human1, UniVec_Core) + plants + protists + fungi

mkdir kraken2-db/ && cd kraken2-db/

# Inspect file
INSPECT_URL=https://genome-idx.s3.amazonaws.com/kraken/pluspfp_20250714/inspect.txt
wget ${INSPECT_URL}

# Library report
LIRARY_REPORT_URL=https://genome-idx.s3.amazonaws.com/kraken/pluspfp_20250714/library_report.tsv
wget ${LIRARY_REPORT_URL}

# Md5sums
MD5_URL=https://genome-idx.s3.amazonaws.com/kraken/pluspfp_20250714/pluspfp.md5 
wget ${MD5_URL}

# Download and unzip the main database files
DB_URL=https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_20250714.tar.gz 
wget -O k2_pluspfp.tar.gz --timeout=3600 --tries=0 --continue ${DB_URL} && \
tar -xvzf k2_pluspfp.tar.gz
```

**Parameter Definitions:**

**wget**

- `O` - Name of file to download the url content to.
- `--timeout=3600` - Specifies the network timeout in seconds.
- `--tries=0` - Retry download infinitely.
- `--continue` -  Continue getting a partially-downloaded file.
- `*_URL` - Position arguement specifying the url to download a particular resource from.


**Input Data:**

- `INSPECT_URL=` - url specifying the location of kraken2 inspect file
- `LIRARY_REPORT_URL=` - url specifying the location of kraken2 library report file
- `MD5_URL=` - url specifying the location of the md5 file of the kraken database
- `DB_URL=` - url specifying the location of the main kraken database archive in .tar.gz format

**Output Data:**

- kraken2-db/  (a directory containing kraken2 database files)

#### 11b. Kraken2 Taxonomic Classification

```bash
kraken2 --db kraken2-db/ \
        --gzip-compressed \
        --threads NumberOfThreads \
        --use-names \
        --output sample-kraken2-output.txt \
        --report sample-kraken2-report.tsv \
        /path/to/sample_decontam_GLlblMetag.fastq.gz
```

**Parameter Definitions:**

- `--db` - Specifies the directory holding the kraken2 database files. 
- `--gzip-compressed` - Specifies the input files are gzip-compressed.
- `--threads` - Number of parallel processing threads to use.
- `--use-names` - Specifies to add taxa names in addition to taxids.
- `--output` - Specifies the name of the kraken2 read-based output file.
- `--report` - Specifies the name of the kraken2 report output file.
- `sample_decontam_GLlblMetag.fastq.gz` - Positional argument specifying the input file.

**Input Data:**

- kraken2-db/ (a directory containing kraken2 database files, output from [Step 11a](#11a-download-kraken2-database))
- sample_decontam_GLlblMetag.fastq.gz or sample_HostRm_GLlblMetag.fastq.gz (filtered and trimmed sample reads with both 
    contaminants and human reads (and, optionally, host reads) removed, gzipped fasta file, 
    output from [Step 7e](#7e-generate-decontaminated-read-files) or [Step 8b](#8b-remove-host-reads))

**Output Data:**

- sample-kraken2-output.txt (kraken2 read-based output file (one line per read))
- sample-kraken2-report.tsv (kraken2 report output file (one line per taxa, with number of reads assigned to it))


#### 11c. Compile Kraken2 Taxonomy Results

##### 11ci. Create Merged Kraken2 Taxonomy Table

```R
species_table <- merge_kraken_reports(reports-dir='/path/to/kraken2/reports')
write_csv(x = species_table, file = "kraken2_species_table_GLlblMetag.csv")
```

**Custom Functions Used:**

- [merge_kraken_reports()](#merge_kraken_reports)

**Parameter Definitions:**

- `file_path` - path to compiled kaiju table at the species taxon level
- `x`  - feature table dataframe to write to file
- `file` - path to where to write kaiju count table per sample

**Parameter Definitions:**

- `--output` - Specifies the name of the kraken2 compiled results output file.
- `--report-files` - Specifies the name of each input kraken2 report file to compile.
- `--sample-names` - Specifies the name of each sample. 

**Input Data:**

- \*-kraken2-report.tsv (kraken report from each sample to compile, outputs from [Step 11b](#11b-taxonomic-classification))

**Output Data:**

- **kraken2_species_table_GLlblMetag.csv** (kraken species count table in csv format)

##### 11cii. Compile Kraken2 Taxonomy Reports

```bash
multiqc --zip-data-dir \ 
        --outdir kraken2_multiqc_report \
        --filename kraken2_multiqc_GLlblMetag \
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

- \*-kraken2-report.tsv (kraken report from each sample to compile, outputs from [Step 11b](#11b-taxonomic-classification))

**Output Data:**

- **kraken2_multiqc_GLlblMetag.html** (multiqc output html summary)
- **kraken2_multiqc_GLlblMetag_data.zip** (zip archive containing multiqc output data)


#### 11d. Convert Kraken2 Output to Krona Format

```bash
kreport2krona.py --report-file sample-kraken2-report.tsv  \
                 --output sample.krona
```

**Parameter Definitions:**

- `--report-file` - Specifies the name of the input kraken2 report file.
- `--output` - Specifies the name of the krona output file.

**Input Data:**

- sample-kraken2-report.tsv (kraken report, output from [Step 11b](#11b-taxonomic-classification))

**Output Data:**

- sample.krona (krona formatted kraken2 output)


#### 11e. Compile Kraken2 Krona Reports

```bash
# Find, list and write all .krona files to file 
find . -type f -name "*.krona" | sort -uV > krona_files.txt

FILES=($(find . -type f -name "*.krona"))
basename --multiple --suffix='.krona' ${FILES[*]} | sort -uV  > sample_names.txt

# Create ktImportText input format files
KTEXT_FILES=($(paste -d',' "krona_files.txt" "sample_names.txt"))

# Create html   
ktImportText -o kraken2-report_GLlblMetag.html ${KTEXT_FILES[*]}
```

**Parameter Definitions:**

**find**

- `-type f` -  Specifies that the type of file to find is a regular file.
- `-name "*.krona"` - Specifies to find files ending with the .krona suffix.  

**sort**

- `-u` - Specifies to perform a unique sort.
- `-V` - Specifies to perform a mixed type of sorting with names containing numbers within text.
- `> {}.txt` - Redirects the sorted list to a separate text file.

**basename**

- `--multiple` - Support multiple arguments and treat each as a file name.
- `--suffix='.krona'` - Remove a trailing '.krona' suffix.

**paste**

- `-d','` - Paste both krona and sample files together line by line delimited by comma ','.

**ktImportText**

- `-o` - Specifies the compiled output html file name.
- `${KTEXT_FILES[*]}` - An array positional arguement with the following content: sample_1.krona,sample_1 sample_2.krona,sample_2 .. sample_n.krona,sample_n.

**Input Data:**

- *.krona (all sample .krona formatted files, output from [Step 11d](#11d-convert-kraken2-output-to-krona-format)) 

                      
**Output Data:**

- krona_files.txt (sorted list of all *.krona files)
- sample_names.txt (sorted list of all sample names)
- **kraken2-report_GLlblMetag.html** (compiled krona html report containing all samples)

---

#### 11f. Filter Kraken2 Species Count Table

```R
library(tidyverse)

input_file <- "kraken2_species_table_GLlblMetag.csv"
output_file <- "kraken2_filtered_species_table_GLlblMetag.csv"
threshold <- 0.5

# string used to define non-microbial taxa
non_microbial <- "UNCLASSIFIED|Unclassifed|unclassified|Homo sapien|cannot|uncultured|unidentified"

# read in feature table
feature_table <- read_csv(input_file) %>% as.data.frame
feature_name <- colnames(feature_table)[1]
rownames(feature_table) <- feature_table[,1]
feature_table <- feature_table[, -1]

# read-based count table
table2write <- filter_rare(feature_table, non_microbial, threshold = threshold) %>%
  as.data.frame %>%
  rownames_to_column(feature_name)

write_csv(x = table2write, file = output_file)
```

**Custom Functions Used:**

- [group_low_abund_taxa()](#group_low_abund_taxa)

**Parameter Definitions:**

- `threshold` - threshold for filtering out rare taxa, a percentage between 0 and 100.
- `output_file` - output filename
- `input_file` - input filename

**Input Data:**

- kraken2_species_table_GLlblMetag.csv (path to kaiju species table from [Step 11ci.](#11ci-create-merged-kraken2-taxonomy-table))

**Output Data:**

- **kraken2_filtered_species_table_GLlblMetag.csv** - a file containing the filtered species table

---


#### 11g. Taxonomy barplots

```R
library(tidyverse)

species_table_file <- "kraken2_species_table_GLlblMetag.csv"
filtered_species_table_file <- "kraken2_filtered_species_table_GLlblMetag.csv"
metadata_file <- "/path/to/sample/metadata"
number_samples <- 10 

# set width based on number of samples, with a cap at 50 inches
plot_width <- 2 * number_samples
if(plot_width > 50) { plot_width = 50 }

p <- make_barplot(metadata_file = metadata_file, feature_table_file = species_table_file, 
                  feature_column = "species", samples_column = "sample_id", group_column = "group",
                  publication_format = publication_format, custom_palette = custom_palette)

ggsave(filename = "kraken2_unfiltered_species_barplot_GLlblMetag.png", plot = p,
       device = "png", width = plot_width, height = 10, units = "in", dpi = 300, limitsize = FALSE)

# Save static unfiltered plot
p <- make_barplot(metadata_file = metadata_file, feature_table_file = filtered_species_table_file, 
                  feature_column = "Species", samples_column = "sample_id", group_column = "group",
                  publication_format = publication_format, custom_palette = custom_palette)

# Save interactive unfilterted plot
htmlwidgets::saveWidget(ggplotly(p), glue("kraken2_unfiltered_species_barplot_GLlblMetag.html"), selfcontained = TRUE)

# Save static filtered plot
ggsave(filename = glue("kraken2_filtered_species_barplot_GLlblMetag.png"), plot = p,
      device = 'png', width = plot_width, height = 10, units = 'in', dpi = 300, limitsize = FALSE)

# Save interactive filtered plot
htmlwidgets::saveWidget(ggplotly(p), glue("kraken2_filtered_species_barplot_GLlblMetag.html"), selfcontained = TRUE)
```
**Custom Functions Used:**
- [make_barplot()](#make_plot)

**Parameter Definitions:**

- `species_table_file` - a file containing the species count table
- `filtered_species_table_file` - a file containing the filtered species count table
- `metadata_file` - a file containing group information for each sample in the species count files
- `number_samples` - the total number of samples in the species count files, adjust based on number of input samples.

**Input Data:**

- `kraken2_species_table_GLlblMetag.csv` (path to kaiju species table from [Step 11ci.](#11ci-create-merged-kraken2-taxonomy-table))
- `kraken2_filtered_species_table_GLlblMetag.csv` (a file containing the filtered species count table, output from [Step 11f](#11f-filter-kraken2-species-count-table))
- `/path/to/sample/metadata` (a file containing sample-wise metadata, mapping samplenames to group metadata)

**Output Data:**

- kraken2_unfiltered_species_barplot_GLlblMetag.png (taxonomy barplot without filtering)
- **kraken2_unfiltered_species_barplot_GLlblMetag.html** (interactive taxonomy barplot without filtering)
- kraken2_filtered_species_barplot_GLlblMetag.png (taxonomy barplot after filtering rare and non-microbial taxa)
- **kraken2_filtered_species_barplot_GLlblMetag.html** (interactive taxonomy barplot after filtering rare and non-microbial taxa)


#### 11h. Feature decontamination

> Feature (species) decontamination with decontam. Decontam is an R package that statistically 
  identifies contaminating features in a feature table

```R
library(tidyverse)
library(decontam)
library(phyloseq)

feature_table_file <- "kraken2_filtered_species_table_GLlblMetag.csv"
metadata_table <- "/path/to/sample/metadata"
number_samples <- NumberOfSamples # integer indicating how many samples are in the file

# set width based on number of samples, with a cap at 50 inches
plot_width <- 2 * number_samples
if(plot_width > 50) { plot_width = 50 }

decontaminated_table <- feature_decontam(metadata_file = metadata_table, 
                                         feature_table_file = feature_table_file, 
                                         feature_column = "species", 
                                         samples_column = "sample_id",
                                         prevalence_column = "NTC", 
                                         ntc_name = "TRUE", 
                                         frequency_column = "concentration", 
                                         threshold = 0.1, 
                                         classification_method = "kraken2", 
                                         output_prefix = "", 
                                         assay_suffix = "_GLlblMetag")

# Convert count matrix to relative abundance matrix
decontaminated_species_table <- count_to_rel_abundance(decontaminated_table)

# Make plot after filtering out contaminants
p <- make_plot(decontaminated_species_table, metadata, custom_palette, publication_format)

ggsave(filename = "kraken2_decontam_species_barplot_GLlblMetag.png", plot = p,
         device = "png", width = plot_width, height = 10, units = "in", dpi = 300)

# Save interactive filtered plot
htmlwidgets::saveWidget(ggplotly(p), glue("kraken2_decontam_species_barplot_GLlblMetag.html"), selfcontained = TRUE)
```

**Custom Functions Used:**
- [feature_decontam()](#feature_decontam)
- [make_plot()](#make_plot)
- [count_to_rel_abundance()](#count_to_rel_abundance)

**Parameter Definitions:**
- `metadata_table` - path to a file with samples as rows and columns describing each sample
- `feature_table_file` - path to a tab separated samples feature table i.e. species/functions 
                          table with species/functions as the first column and samples as other columns.
- `number_samples` - the total number of samples in the species count files, adjust based on number of input samples.

**Input Data:**

- `kraken2_filtered_species_table_GLlblMetag.csv`(path to filtered species count per sample, output from [Step 11f](#11f-filter-kraken2-species-count-table))
- `/path/to/sample/metadata` (a file containing sample-wise metadata, mapping samplenames to group metadata)

**Output Data:**

- **kraken2_decontam_results_GLlblMetag.csv** (decontam's result table, output from [feature_decontam() function](#feature_decontam))
- **kraken2_decontam_species_table_GLlblMetag.csv** (decontaminated species table, output from [feature_decontam() function](#feature_decontam))
- **kraken2_decontam_species_barplot_GLlblMetag.png** (barplot after filtering out contaminants)
- **kraken2_decontam_species_barplot_GLlblMetag.html** (barplot after filtering out contaminants)

<br>

---

## Assembly-based Processing

### 12. Sample Assembly

```bash
flye --meta \
     --threads NumberOfThreads \
     --out-dir sample/ \
     --nano-hq \
     /path/to/sample_decontam_GLlblMetag.fastq.gz

# rename output files            
mv sample/assembly.fasta sample_assembly.fasta
mv sample/flye.log sample_assembly.log
```

**Parameter Definitions:**

- `--meta` – Use metagenome/uneven coverage mode.
- `--threads` - Number of parallel processing threads to use.
- `--out-dir` - Specifies the name of the output directory.
- `--nano-hq` - Specifies that input is from Oxford Nanopore high-quality reads (Guppy5+ SUP or Q20, <5% error). This skips a genome polishing step since the assembly will be polished with medaka in the next step.
- `/path/to/sample_decontam_GLlblMetag.fastq.gz` - Path to the input file, specified as a positional argument.

**Input Data**

- sample_decontam_GLlblMetag.fastq.gz or sample_HostRm_GLlblMetag.fastq.gz (filtered and trimmed sample reads with both 
    contaminants and human reads (and optionally host reads) removed, gzipped fasta file, 
    output from [Step 7e](#7e-generate-decontaminated-read-files) or [Step 8b](#8b-remove-host-reads))

**Output Data**

- sample_assembly.fasta (sample assembly fasta)
- sample_assembly.log (flye log file)

<br>

---

### 13. Polish Assembly

```bash
medaka_consensus -t NumberOfThreads \
                 -i /path/to/sample_decontam_GLlblMetag.fastq.gz \
                 -d /path/to/assemblies/sample_assembly.fasta \
                 -o sample/
  
mv sample/consensus.fasta sample_polished.fasta
```

**Parameter Definitions:**

- `-t` - Number of parallel processing threads to use.
- `-i` - Specifies path to input read files used in creating the assembly.
- `-d` - Specifies path to the assembly fasta file.
- `-o` - Specifies the output directory.

**Input Data:**

- sample_decontam_GLlblMetag.fastq.gz or sample_HostRm_GLlblMetag.fastq.gz (filtered and trimmed sample reads with both 
    contaminants and human reads (and optionally host reads) removed, gzipped fasta file, 
    output from [Step 7e](#7e-generate-decontaminated-read-files) or [Step 8b](#8b-remove-host-reads))
- /path/to/assemblies/sample_assembly.fasta (sample assembly, output from [Step 12](#12-sample-assembly))

**Output Data:**

- sample_polished.fasta (polished sample assembly)

---

### 14. Rename Contigs and Summarize Assemblies

#### 14a. Rename Contig Headers

```bash
bit-rename-fasta-headers -i sample_polished.fasta \
                         -w c_sample \
                         -o sample-assembly_GLlblMetag.fasta
```

**Parameter Definitions:**  

- `-i` – Specifies the input fasta file.
- `-w` – Specifies the wanted header prefix (a number will be appended for each contig), starts with a "c" to ensure they won't start with a number which can be problematic.
- `-o` – Specifies the output fasta file.


**Input Data:**

- sample_polished.fasta (polished assembly file from [Step 13](#13-polish-assembly))

**Output files:**

- **sample-assembly_GLlblMetag.fasta** (contig-renamed assembly file)


#### 14b. Summarize Assemblies

```bash
bit-summarize-assembly -o assembly-summaries_GLlblMetag.tsv \
                       *-assembly_GLlblMetag.fasta
```

**Parameter Definitions:**  

- `-o` – Specifies the output summary table.
- `*-assembly.fasta` - Specifies the input assemblies to summarize, provided as positional arguments.

**Input Data:**

- *-assembly.fasta (contig-renamed assembly files from [Step 14a](#14a-renaming-contig-headers))

**Output files:**

- **assembly-summaries_GLlblMetag.tsv** (table of assembly summary statistics)

<br>

---

### 15. Gene Prediction

#### 15a. Generate Gene Predictions

```bash
prodigal -a sample-genes.faa \
         -d sample-genes.fasta \
         -f gff \
         -p meta \
         -c \
         -q \
         -o sample-genes_GLlblMetag.gff \
         -i sample-assembly_GLlblMetag.fasta
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

**Input Data:**

- sample-assembly_GLlblMetag.fasta (contig-renamed assembly file from [Step 14a](#14a-renaming-contig-headers))

**Output Data:**

- sample-genes.faa (gene-calls amino-acid fasta file)
- sample-genes.fasta (gene-calls nucleotide fasta file)
- **sample-genes_GLlblMetag.gff** (gene-calls in general feature format)

<br>

#### 15b. Remove Line Wraps In Gene Prediction Output

```bash
bit-remove-wraps sample-genes.faa > sample-genes.faa.tmp 2> /dev/null
mv sample-genes.faa.tmp sample-genes_GLlblMetag.faa

bit-remove-wraps sample-genes.fasta > sample-genes.fasta.tmp 2> /dev/null
mv sample-genes.fasta.tmp sample-genes_GLlblMetag.fasta
```

**Input Data:**

- sample-genes.faa (gene-calls amino-acid fasta file, output from [Step 15a](#15a-gene-prediction))
- sample-genes.fasta (gene-calls nucleotide fasta file, output from [Step 15a](#15a-gene-prediction))

**Output Data:**

- **sample-genes_GLlblMetag.faa** (gene-calls amino-acid fasta file with line wraps removed)
- **sample-genes_GLlblMetag.fasta** (gene-calls nucleotide fasta file with line wraps removed)

<br>

---

### 16. Functional Annotation

> **Note:**  
> The annotation process overwrites the same temporary directory by default. When running multiple 
processses at a time, it is necessary to specify a specific temporary directory with the 
`--tmp-dir` argument as shown below.


#### 16a. Download Reference Database of HMM Models

> **Note:** This step only needs to be done once.

```bash
curl -LO ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz
curl -LO ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz
tar -xzvf profiles.tar.gz
gunzip ko_list.gz 
```

#### 16b. Run KEGG Annotation

```bash
exec_annotation -p profiles/ \
                -k ko_list \
                --cpu NumberOfThreads \
                -f detail-tsv \
                -o sample-KO-tab.tmp \
                --tmp-dir sample-tmp-KO \
                --report-unannotated \
                sample-genes_GLlblMetag.faa 
```

**Parameter Definitions:**

- `-p` – Specifies the directory holding the downloaded reference HMMs.
- `-k` – Specifies the downloaded reference KO  (Kegg Orthology) terms. 
- `--cpu` – Specifies the number of searches to run in parallel.
- `-f` – Specifies the output format.
- `-o` – Specifies the output file name.
- `--tmp-dir` – Specifies the temporary directory to write to (needed if running more than one process concurrently, see Note above).
- `--report-unannotated` – Specifies to generate an output for each entry, event when no KO is assigned.
- `sample-genes.faa` – Specifies the input file, provided as a positional argument. 


**Input Data:**

- sample-genes_GLlblMetag.faa (amino-acid fasta file, output from [Step 15b](#15b-remove-line-wraps-in-gene-prediction-output))
- profiles/ (reference directory holding the KO HMMs, downloaded in [Step 16a](16a-download-reference-database-of-hmm-models))
- ko_list (reference list of KOs to scan for, downloaded in [Step 16a](16a-download-reference-database-of-hmm-models))

**Output Data:**

- sample-KO-tab.tmp (table of KO annotations assigned to gene IDs)


#### 16c. Filter KO Outputs
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

**Input Data:**

- sample-KO-tab.tmp (table of KO annotations assigned to gene IDs, output from [Step 16b](#16b-run-kegg-annotation))

**Output Data:**

- sample-annotations.tsv (table of KO annotations assigned to gene IDs)

<br>

---

### 17. Taxonomic Classification 

#### 17a. Pull and Unpack Pre-built Reference DB 

> **Note:** This step only needs to be done once.

```bash
wget tbb.bio.uu.nl/bastiaan/CAT_prepare/CAT_prepare_20200618.tar.gz
tar -xvzf CAT_prepare_20200618.tar.gz
```

#### 17b. Run Taxonomic Classification

```bash
CAT contigs -c sample-assembly.fasta \
            -d CAT_prepare_20200618/2020-06-18_database/ \
            -t CAT_prepare_20200618/2020-06-18_taxonomy/ \
            -p sample-genes.faa \
            -o sample-tax-out.tmp \
            -n NumberOfThreads \
            -r 3 \
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

**Input Data:**

- CAT_prepare_20200618/2020-06-18_database/ (directory holding the CAT reference sequence database, output from [Step 17a](17a-pull-and-unpack-pre-built-reference-db))
- CAT_prepare_20200618/2020-06-18_taxonomy/ (directory holding the CAT reference taxonomy database, output from [Step 17a](17a-pull-and-unpack-pre-built-reference-db))
- sample-assembly.fasta (contig-renamed assembly file from [Step 14a](#14a-rename-contig-headers))
- sample-genes.faa (amino-acid fasta file, output from [Step 15b](#15b-remove-line-wraps-in-gene-prediction-output))

**Output Data:**

- sample-tax-out.tmp.ORF2LCA.txt (gene-calls taxonomy file)
- sample-tax-out.tmp.contig2classification.txt (contig taxonomy file)


#### 17c. Add Taxonomy Info From Taxids To Genes

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

**Input Data:**

- sample-tax-out.tmp.ORF2LCA.txt (gene-calls taxonomy file, output from [Step 17b](#17b-run-taxonomic-classification))
- CAT_prepare_20200618/2020-06-18_taxonomy/ (directory holding the CAT reference taxonomy database, output from [Step 17a](#17a-pull-and-unpack-pre-built-reference-db))

**Output Data:**

- sample-gene-tax-out.tmp (gene-calls taxonomy file with lineage info added)


#### 17d. Add Taxonomy Info From Taxids To Contigs

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

**Input Data:**

- sample-tax-out.tmp.contig2classification.txt (contig taxonomy file, output from [Step 17b](#17b-run-taxonomic-classification))
- CAT_prepare_20200618/2020-06-18_taxonomy/ (directory holding the CAT reference taxonomy database, output from [Step 17a](#17a-pull-and-unpack-pre-built-reference-db))

**Output Data:**

- sample-contig-tax-out.tmp (contig taxonomy file with lineage info added)


#### 17e. Format Gene-level Output With awk and sed

```bash
awk -F $'\t' ' BEGIN { OFS=FS } { if ( $3 == "lineage" ) { print $1,$3,$5,$6,$7,$8,$9,$10,$11 } \
    else if ( $2 == "ORF has no hit to database" || $2 ~ /^no taxid found/ ) \
    { print $1,"NA","NA","NA","NA","NA","NA","NA","NA" } else { n=split($3,lineage,";"); \
    print $1,lineage[n],$5,$6,$7,$8,$9,$10,$11 } } ' sample-gene-tax-out.tmp | \
    sed 's/no support/NA/g' | sed 's/superkingdom/domain/' | sed 's/# ORF/gene_ID/' | \
    sed 's/lineage/taxid/'  > sample-gene-tax-out.tsv
```

**Input Data:**

- sample-gene-tax-out.tmp (gene-calls taxonomy file with lineage info added, output from [Step 17c](#17c-add-taxonomy-info-from-taxids-to-genes))

**Output Data:**

- sample-gene-tax-out.tsv (reformatted gene-calls taxonomy file with lineage info)


#### 17f. Format Contig-level Output With awk and sed

```bash
awk -F $'\t' ' BEGIN { OFS=FS } { if ( $2 == "classification" ) { print $1,$4,$6,$7,$8,$9,$10,$11,$12 } \
    else if ( $2 == "no taxid assigned" ) { print $1,"NA","NA","NA","NA","NA","NA","NA","NA" } \
    else { n=split($4,lineage,";"); print $1,lineage[n],$6,$7,$8,$9,$10,$11,$12 } } ' sample-contig-tax-out.tmp | \
    sed 's/no support/NA/g' | sed 's/superkingdom/domain/' | sed 's/^# contig/contig_ID/' | \
    sed 's/lineage/taxid/' > sample-contig-tax-out.tsv

  # clearing intermediate files
rm sample*.tmp*
```

**Input Data:**

- sample-contig-tax-out.tmp (contig taxonomy file with lineage info added, output from [Step 17d](#17d-add-taxonomy-info-from-taxids-to-contigs))

**Output Data:**

- sample-contig-tax-out.tsv (reformatted contig taxonomy file with lineage info)

<br>

---

### 18. Read-Mapping

#### 18a. Align Reads to Sample Assembly

```bash
minimap2 -a \
         -x map-ont \
         -t NumberOfThreads \
         sample_assembly.fasta \
         sample_decontam_GLlblMetag.fastq.gz \
         > sample.sam  2> sample-mapping-info.txt
```

**Parameter Definitions:**

- `-a` – Output in SAM format.
- `-x map-ont` - Specifies preset for mapping Nanopore reads to a reference.
- `-t` - Number of parallel processing threads to use
- `sample_assembly.fasta` – Assembly fasta file, provided as a positional argument.
- `sample_decontam_GLlblMetag.fastq.gz` - Input sequence data file, provided as a positional argument.
- `> sample.sam` - Redirects the output to a separate file.
- `2> sample-mapping-info.txt` - Redirects the standar error to a separate file.

**Input Data**

- sample-assembly.fasta (contig-renamed assembly file, output from [Step 14a](#14a-rename-contig-headers))
- sample_decontam_GLlblMetag.fastq.gz or sample_HostRm_GLlblMetag.fastq.gz (filtered and trimmed sample reads with both 
    contaminants and human reads (and optionally host reads) removed, gzipped fasta file, 
    output from [Step 7e](#7e-generate-decontaminated-read-files) or [Step 8b](#8b-remove-host-reads))

**Output Data**

- sample.sam (reads aligned to sample assembly in SAM format)
- **sample-mapping-info_GLlblMetag.txt** (read mapping information)


#### 18b. Sort and Index Assembly Alignments

```bash
# Sort Sam, convert to bam and create index
samtools sort --threads NumberOfThreads \
              -o sample_sorted_GLlblMetag.bam \
              sample.sam > sample_sort.log 2>&1

samtools index sample_sorted_GLlblMetag.bam sample_sorted_GLlblMetag.bam.bai
```

**Parameter Definitions:**

**samtools sort**
- `--threads` - Number of parallel processing threads to use.
- `-o` - Specifies the output file for the sorted aligned reads.
- `sample.sam` - Positional argument specifying the input SAM file.
- `> sample_sort.log 2>&1` - Redirects the standard output and standard error to a separate file.

**samtools index**
- `sample_sorted.bam` - Positional argument specifying the input BAM file to be sorted.
- `sample_sorted.bam.bai` - Positional argument specifying the name of the index file.

**Input Data:**

- sample.sam (reads aligned to sample assembly, output from [Step 18a](#18a-align-reads-to-sample-assembly))

**Output Data:**

- **sample_sorted_GLlblMetag.bam** (sorted mapping to sample assembly, in BAM format)
- **sample_sorted_GLlblMetag.bam.bai** (index of sorted mapping to sample assembly)

<br>

---

### 19. Get Coverage Information and Filter Based On Detection
> **Note:**  
> “Detection” is a measure of what proportion of a reference sequence recruited reads 
(see the discussion of detection [here](http://merenlab.org/2017/05/08/anvio-views/#detection)). 
Filtering based on detection is one way of helping to mitigate non-specific read-recruitment.

#### 19a. Filter Coverage Levels Based On Detection

```bash
# pileup.sh comes from the bbduk.sh package
pileup.sh -in sample.bam \
          fastaorf=sample-genes_GLlblMetag.fasta \
          outorf=sample-gene-cov-and-det.tmp \
          out=sample-contig-cov-and-det.tmp
```

**Parameter Definitions:**  

- `-in` – Specifies the input BAM file.
- `fastaorf=` – Specifies the input gene-calls nucleotide fasta file.
- `outorf=` – Specifies the output gene-coverage tsv file name.
- `out=` – Specifies the output contig-coverage tsv file name.

**Input Data:**

- sample.bam (sorted mapping to sample assembly BAM file, output from [Step 18b](#18b-sort-and-index-assembly-alignments))
- sample-genes.fasta (gene-calls nucleotide fasta file, output from [Step 15a](#15a-gene-prediction))


**Output Data:**

- sample-gene-cov-and-det.tmp (gene-coverage tsv file)
- sample-contig-cov-and-det.tmp (contig-coverage tsv file)


#### 19b. Filter Gene and Contig Coverage Based On Detection

> *The following commands filter gene and contig coverage tsv files to only keep genes and contigs with at least 50% detection (as defined above) then parse the tables to retain only gene IDs and respective coverage.*

```bash
# Filtering gene coverage
grep -v "#" sample-gene-cov-and-det.tmp | \
awk -F $'\t' ' BEGIN { OFS=FS } { if ( $10 <= 0.5 ) $4 = 0 } \
     { print $1,$4 } ' > sample-gene-cov.tmp

cat <( printf "gene_ID\tcoverage\n" ) sample-gene-cov.tmp > sample-gene-coverages_GLlblMetag.tsv

# Filtering contig coverage
grep -v "#" sample-contig-cov-and-det.tmp | \
awk -F $'\t' ' BEGIN { OFS=FS } { if ( $5 <= 50 ) $2 = 0 } \
     { print $1,$2 } ' > sample-contig-cov.tmp

cat <( printf "contig_ID\tcoverage\n" ) sample-contig-cov.tmp > sample-contig-coverages_GLlblMetag.tsv

# removing intermediate files
rm sample-*.tmp
```

**Input Data:**

- sample-gene-cov-and-det.tmp (temporary gene-coverage tsv file, output from [Step 19a](#19a-filter-coverage-levels-based-on-detection))
- sample-contig-cov-and-det.tmp (temporary contig-coverage tsv file, output from [Step 19a](#19a-filter-coverage-levels-based-on-detection))

**Output Data:**

- sample-gene-coverages_GLlblMetag.tsv (table with gene-level coverages)
- sample-contig-coverages_GLlblMetag.tsv (table with contig-level coverages)

<br>

---

### 20. Combine Gene-level Coverage, Taxonomy, and Functional Annotations For Each Sample
> **Note:**  
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

cat sample-header.tmp sample-gene-tab.tmp > sample-gene-coverage-annotation-and-tax_GLlblMetag.tsv

# removing intermediate files
rm sample*tmp sample-gene-coverages.tsv sample-annotations.tsv sample-gene-tax-out.tsv
```

**Input Data:**

- sample-gene-coverages.tsv (table with gene-level coverages, output from [Step 19b](#19b-filter-gene-and-contig-coverage-based-on-detection))
- sample-annotations.tsv (table of KO annotations assigned to gene IDs, output from [Step 16c](#16c-filter-ko-outputs))
- sample-gene-tax-out.tsv (reformatted gene-calls taxonomy file with lineage info, output from [Step 17e](#17e-format-gene-level-output-with-awk-and-sed))


**Output Data:**

- **sample-gene-coverage-annotation-and-tax_GLlblMetag.tsv** (table with combined gene coverage, annotation, and taxonomy info)

<br>

---

### 21. Combine Contig-level Coverage and Taxonomy For Each Sample
> **Note:**  
> Just uses `paste`, `sed`, and `awk` standard Unix commands to combine contig-level coverage and taxonomy into one table for each sample.

```bash
paste <( tail -n +2 sample-contig-coverages.tsv | sort -V -k 1 ) \
      <( tail -n +2 sample-contig-tax-out.tsv | sort -V -k 1 | cut -f 2- ) \
      > sample-contig.tmp

paste <( head -n 1 sample-contig-coverages.tsv ) \
      <( head -n 1 sample-contig-tax-out.tsv | cut -f 2- ) \
      > sample-contig-header.tmp
      
cat sample-contig-header.tmp sample-contig.tmp > sample-contig-coverage-and-tax_GLlblMetag.tsv

# removing intermediate files
rm sample*tmp sample-contig-coverages.tsv sample-contig-tax-out.tsv
```

**Input Data:**

- sample-contig-coverages.tsv (table with contig-level coverages, output from [Step 19b](#19b-filter-gene-and-contig-coverage-based-on-detection))
- sample-contig-tax-out.tsv (reformatted contig taxonomy file with lineage info, output from [Step 17f](#17f-format-contig-level-output-with-awk-and-sed))


**Output Data:**

- **sample-contig-coverage-and-tax_GLlblMetag.tsv** (table with combined contig coverage and taxonomy info)

<br>

---

### 22. Generate Normalized, Gene- and Contig-level Coverage Summary Tables of KO-annotations and Taxonomy Across Samples

> **Note:**  
> * To combine across samples to generate these summary tables, we need the same "units". This is done for annotations 
based on the assigned KO terms, and all non-annotated functions are included together as "Not annotated". It is done for 
taxonomic classifications based on taxids (full lineages included in the table), and any genes not classified are included 
together as "Not classified". 
> * The values we are working with are coverage per gene (so they are number of bases recruited to the gene normalized 
by the length of the gene). These have been normalized by making the total coverage of a sample 1,000,000 and setting 
each individual gene-level coverage its proportion of that 1,000,000 total. So basically percent, but out of 1,000,000 
instead of 100 to make the numbers more friendly. 

#### 22a. Generate Gene-level Coverage Summary Tables

```bash
bit-GL-combine-KO-and-tax-tables *-gene-coverage-annotation-and-tax_GLlblMetag.tsv \
                                 -o Combined

# add assay specific suffix
mv "Combined-gene-level-KO-function-coverages-CPM.tsv Combined-gene-level-KO-function-coverages-CPM_GLlblMetag.tsv"
mv "Combined-gene-level-KO-function-coverages-CPM.tsv Combined-gene-level-KO-function-coverages-CPM_GLlblMetag.tsv"
mv "Combined-gene-level-KO-function-coverages.tsv Combined-gene-level-KO-function-coverages_GLlblMetag.tsv"
mv "Combined-gene-level-taxonomy-coverages.tsv Combined-gene-level-taxonomy-coverages_GLlblMetag.tsv"
```

**Parameter Definitions:**  

- `*-gene-coverage-annotation-and-tax_GLlbsMetag.tsv` - Positional arguments specifying the input tsv files, can be provided as a space-delimited list of files, or with wildcards like above.

- `-o` – Specifies the output file prefix.


**Input Data:**

- *-gene-coverage-annotation-and-tax_GLlblMetag.tsv (tables with combined gene coverage, annotation, and taxonomy info generated for individual samples, output from [Step 20](#20-combine-gene-level-coverage-taxonomy-and-functional-annotations-for-each-sample))

**Output Data:**

- **Combined-gene-level-KO-function-coverages-CPM_GLlblMetag.tsv** (table with all samples combined based on KO annotations; normalized to coverage per million genes covered)
- **Combined-gene-level-taxonomy-coverages-CPM_GLlblMetag.tsv** (table with all samples combined based on gene-level taxonomic classifications; normalized to coverage per million genes covered)
- **Combined-gene-level-KO-function-coverages_GLlblMetag.tsv** (table with all samples combined based on KO annotations)
- **Combined-gene-level-taxonomy-coverages_GLlblMetag.tsv** (table with all samples combined based on gene-level taxonomic classifications)

#### 22b. Generate Contig-level Coverage Summary Tables

```bash
bit-GL-combine-contig-tax-tables *-contig-coverage-and-tax.tsv -o Combined
```

**Parameter Definitions:**  

- `*-contig-coverage-and-tax.tsv` - Positional arguments specifying the input tsv files, can be provided as a space-delimited list of files, or with wildcards like above.
- `-o` – Specifies the output file prefix.


**Input Data:**

- *-contig-coverage-and-tax.tsv (tables with combined contig coverage and taxonomy info generated for individual samples, output from [Step 21](#21-combine-contig-level-coverage-and-taxonomy-for-each-sample))

**Output Data:**

- **Combined-contig-level-taxonomy-coverages-CPM_GLlblMetag.tsv** (table with all samples combined based on contig-level taxonomic classifications; normalized to coverage per million contigs covered)
- **Combined-contig-level-taxonomy-coverages_GLlblMetag.tsv** (table with all samples combined based on contig-level taxonomic classifications)

<br>

---

### 23. **M**etagenome-**A**ssembled **G**enome (MAG) Recovery

#### 23a. Bin Contigs

```bash
jgi_summarize_bam_contig_depths --outputDepth sample-metabat-assembly-depth.tsv \
                                --percentIdentity 97 \
                                --minContigLength 1000 \
                                --minContigDepth 1.0  \
                                --referenceFasta sample-assembly.fasta \
                                sample.bam

metabat2  --inFile sample-assembly.fasta \
          --outFile sample \
          --abdFile sample-metabat-assembly-depth.tsv \
          -t NumberOfThreads

mkdir sample-bins
mv sample*bin*.fasta sample-bins
zip -r sample-bins.zip sample-bins
```

**Parameter Definitions:**  

**jgi_summarize_bam_contig_depths**

-  `--outputDepth` – Specifies the output depth file name.
-  `--percentIdentity` – Minimum end-to-end percent identity of a mapped read to be included.
-  `--minContigLength` – Minimum contig length to include.
-  `--minContigDepth` – Minimum contig depth to include.
-  `--referenceFasta` – Specifies the input assembly fasta file.
-  `sample.bam` – Input alignment BAM file, specified as a positional argument.

**metabat2**

-  `--inFile` - Specifies the input assembly fasta file.
-  `--outFile` - Specifies the prefix of the identified bins output files.
-  `--abdFile` - The depth file generated by the previous `jgi_summarize_bam_contig_depths` command.
-  `-t` - Number of parallel processing threads to use.


**Input Data:**

- sample-assembly.fasta (contig-renamed assembly file from [Step 14a](#14a-renaming-contig-headers))
- sample.bam (sorted mapping to sample assembly BAM file, output from [Step 18b](#18b-sort-and-index-assembly-alignments))

**Output Data:**

- **sample-metabat-assembly-depth.tsv** (tab-delimited summary of coverages)
- sample-bins/sample-bin\*.fasta (fasta files of recovered bins)
- **sample-bins.zip** (zip file containing fasta files of recovered bins)

#### 23b. Bin quality assessment 
> Utilizes the default `checkm` database available [here](https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz), `checkm_data_2015_01_16.tar.gz`.

```bash
checkm lineage_wf -f bins-overview_GLlblMetag.tsv \
                  --tab_table \
                  -x fasta \
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

**Input Data:**

- sample-bins/sample-bin\*.fasta (fasta files of recovered bins, output from [Step 23a](#23a-bin-contigs))

**Output Data:**

- **bins-overview_GLlblMetag.tsv** (tab-delimited file with quality estimates per bin)
- checkm-output-dir/ (directory holding detailed checkm outputs)

#### 23c. Filter MAGs

```bash
cat <( head -n 1 bins-overview_GLlblMetag.tsv ) \
    <( awk -F $'\t' ' $12 >= 90 && $13 <= 10 && $14 == 0 ' bins-overview_GLlblMetag.tsv | sed 's/bin./MAG-/' ) \
    > checkm-MAGs-overview.tsv
    
# copying bins into a MAGs directory in order to run tax classification
awk -F $'\t' ' $12 >= 90 && $13 <= 10 && $14 == 0 ' bins-overview_GLlblMetag.tsv | cut -f 1 > MAG-bin-IDs.tmp

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
  zip -r ${SAMPLE}-MAGs.zip ${SAMPLE}-MAGs
done
```

**Input Data:**

- bins-overview_GLlblMetag.tsv (tab-delimited file with quality estimates per bin from [Step 23b](#23b-bin-quality-assessment))

**Output Data:**

- checkm-MAGs-overview.tsv (tab-delimited file with quality estimates per MAG)
- MAGs/\*.fasta (directory holding high-quality MAGs)
- **\*-MAGs.zip** (zip files containing directories of high-quality MAGs)


#### 23d. MAG Taxonomic Classification
> Uses default `gtdbtk` database setup with program's `download.sh` command.

```bash
gtdbtk classify_wf --genome_dir MAGs/ \
                   -x fasta \
                   --out_dir gtdbtk-output-dir \
                   --skip_ani_screen
```

**Parameter Definitions:**  

-  `classify_wf` – Specifies the workflow being utilized.
-  `--genome_dir` – Specifies the directory holding the MAGs to classify.
-  `-x` – Specifies the extension that is on the MAG fasta files that are being taxonomically classified.
-  `--out_dir` – Specifies the output directory name.
-  `--skip_ani_screen`  - Specifies to skip ani_screening step to classify genomes using mash and skani.

**Input Data:**

- MAGs/\*.fasta (directory holding high-quality MAGs, output from [Step 23c](#23c-filter-mags))

**Output Data:**

- gtdbtk-output-dir/gtdbtk.\*.summary.tsv (files with assigned taxonomy and info)

#### 23e. Generate Overview Table Of All MAGs

```bash
# combine summaries
for MAG in $(cut -f 1 assembly-summaries_GLlblMetag.tsv | tail -n +2); do

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

paste assembly-summaries_GLlblMetag.tsv \
checkm-estimates-with-headers.tmp \
gtdb-taxonomies-with-headers.tmp \
    > MAGs-overview.tmp

# Ordering by taxonomy
head -n 1 MAGs-overview.tmp > MAGs-overview-header.tmp

tail -n +2 MAGs-overview.tmp | sort -t \$'\t' -k 14,20 > MAGs-overview-sorted.tmp

cat MAGs-overview-header.tmp MAGs-overview-sorted.tmp \
    > MAGs-overview_GLlblMetag.tsv
```

**Input Data:**

- assembly-summaries_GLlblMetag.tsv (table of assembly summary statistics, output from [Step 14b](#14b-summarize-assemblies))
- MAGs/\*.fasta (directory holding high-quality MAGs, output from [Step 23c](#23c-filter-mags))
- checkm-MAGs-overview.tsv (tab-delimited file with quality estimates per MAG, output from [Step 23c](#23c-filter-mags))
- gtdbtk-output-dir/gtdbtk.\*.summary.tsv (directory of files with assigned taxonomy and info, output from [Step 23d](#23d-mag-taxonomic-classification))

**Output Data:**

- **MAGs-overview_GLlblMetag.tsv** (a tab-delimited overview of all recovered MAGs)


<br>

---

### 24. Generate MAG-level Functional Summary Overview

#### 24a. Get KO Annotations Per MAG
> This utilizes the helper script [`parse-MAG-annots.py`](../Workflow_Documentation/NF_MGIllumina/workflow_code/bin/parse-MAG-annots.py) 

```bash
for file in $( ls MAGs/*.fasta )
do

    MAG_ID=$( echo ${file} | cut -f 2 -d "/" | sed 's/.fasta//' )
    sample_ID=$( echo ${MAG_ID} | sed 's/-MAG-[0-9]*$//' )

    grep "^>" ${file} | tr -d ">" > ${MAG_ID}-contigs.tmp

    python parse-MAG-annots.py -i annotations-and-taxonomy/${sample_ID}-gene-coverage-annotation-and-tax.tsv \
                               -w ${MAG_ID}-contigs.tmp \
                               -M ${MAG_ID} \
                               -o MAG-level-KO-annotations_GLlblMetag.tsv

    rm ${MAG_ID}-contigs.tmp

done
```

**Parameter Definitions:**  

- `-i` – Specifies the input sample TSV file containing sample coverage, annotation, and taxonomy info.
- `-w` – Specifies the appropriate temporary file holding all the contigs in the current MAG.
- `-M` – Specifies the current MAG unique identifier.
- `-o` – Specifies the output file name.

**Input Data:**

- \*-gene-coverage-annotation-and-tax.tsv (tables with combined gene coverage, annotation, and taxonomy info generated for individual samples, output from [Step 19](#19-combine-gene-level-coverage-taxonomy-and-functional-annotations-for-each-sample))
- MAGs/\*.fasta (directory holding high-quality MAGs, output from [Step 23c](#23c-filter-mags))

**Output Data:**

- **MAG-level-KO-annotations_GLlblMetag.tsv** (tab-delimited table holding MAGs and their KO annotations)


#### 24b. Summarize KO Annotations With KEGG-Decoder

```bash
KEGG-decoder -v interactive \
             -i MAG-level-KO-annotations_GLlblMetag.tsv \
             -o MAG-KEGG-Decoder-out_GLlblMetag.tsv
```

**Parameter Definitions:**  

- `-v interactive` – Specifies to create an interactive html output.
- `-i` – Specifies the input tab-delimited table holding MAGs and their KO annotations.
- `-o` – Specifies the output table.

**Input Data:**

- MAG-level-KO-annotations_GLlblMetag.tsv (tab-delimited table holding MAGs and their KO annotations, output from [Step 24a](#24a-getting-ko-annotations-per-mag))

**Output Data:**

- **MAG-KEGG-Decoder-out_GLlblMetag.tsv** (tab-delimited table holding MAGs and their proportions of 
                                           genes held known to be required for specific pathways/metabolisms)
- **MAG-KEGG-Decoder-out_GLlbnMetag.html** (interactive heatmap html file of the above output table)

<br>

---

### 25. Decontamination and Visualization of Contig- and Gene-taxonomy and Gene-function Outputs

#### 25a. Gene-level taxonomy heatmaps

```R
library(tidyverse)

metadata_file <- "/path/to/sample/metadata"
feature_data_file <- "Combined-gene-level-taxonomy-coverages-CPM_GLlblMetag.tsv"

# Prepare metadata
metadata <- read_delim(metadata_file, delim = ",") %>% as.data.frame
sample_names = metadata[, samples_column]
row.names(metadata) <- sample_names

# Prepare feature table
gene_taxonomy_table <- read_assembly_coverage_table(feature_table_file, sample_names) %>% as.data.frame

# Summarize gene table
species_gene_table <- gene_taxonomy_table %>%
  select(species, !!any_of(sample_names)) %>% 
  group_by(species) %>% 
  summarise(across(everything(), sum)) %>% 
  filter(species != "Unclassified;_;_;_;_;_;_") %>% # Drop unclassifed
  as.data.frame

rownames(species_gene_table) <- species_gene_table[[1]]
species_gene_table <- species_gene_table[, -1] %>% as.matrix()

# Get common samples and re-arrange feature table and metadata
common_samples <- intersect(colnames(species_gene_table), rownames(metadata))
species_gene_table <- species_gene_table[, common_samples]
metadata <- metadata[common_samples, ]
metadata <- metadata %>% arrange(!!sym(group_column))

table2write = species_gene_table %>% as.data.frame %>% rownames_to_column("species")
# Write out gene taxonomy table
write_csv(x = table2write, file = "gene_taxonomy_table.csv")

make_heatmap(metadata, species_gene_table, 
             samples_column="sample_id", group_column = "group", 
             output_prefix = "Combined-gene-level-taxonomy", 
             assay_suffix = "_GLlblMetag", 
             custom_palette = custom_palette)

```

**Custom Functions Used:**
- [read_assembly_coverage_table()](#read_assembly_coverage_table)
- [make_heatmap()](#make_heatmap)

**Input data:**
- /path/to/sample/metadata (a file with samples as rows and columns describing each sample)
- Combined-gene-level-taxonomy-coverages-CPM_GLlblMetag.tsv (table with all samples 
    combined based on gene-level taxonomic classifications, output from 
    [Step 22a](#22a-generating-gene-level-coverage-summary-tables)) 

**Output data:**
- gene_taxonomy_table.csv (aggregated gene taxonomy table with samples in columns and species in rows)
- **Combined-gene-level-taxonomy_heatmap_GLlblMetag.png** (heatmap of all gene taxonomy assignments)

#### 25b. Gene-level taxonomy decontamination

```R
library(tidyverse)
library(decontam)
library(phyloseq)

feature_table_file <- "gene_taxonomy_table.csv"
metadata_table <- "/path/to/sample/metadata"
number_samples <- NumberOfSamples # integer indicating how many samples are in the file

# set width based on number of samples, with a cap at 50 inches
plot_width <- 2 * number_samples
if(plot_width > 50) { plot_width = 50 }

# Prepare metadata
metadata <- read_delim(metadata_file, delim = ",") %>% as.data.frame
sample_names = metadata[, samples_column]
row.names(metadata) <- sample_names

decontaminated_table <- feature_decontam(metadata_file = metadata_table, 
                                         feature_table_file = feature_table_file, 
                                         feature_column = "species", 
                                         samples_column = "sample_id",
                                         prevalence_column = "NTC", 
                                         ntc_name = "TRUE", 
                                         frequency_column = "concentration", 
                                         threshold = 0.1, 
                                         classification_method = "Combined-gene-level-taxonomy", 
                                         output_prefix = "", 
                                         assay_suffix = "_GLlblMetag")

# Get common samples and re-arrange feature table and metadata
common_samples <- intersect(colnames(decontaminated_table), rownames(metadata))
decontaminated_table <- decontaminated_table[, common_samples]
metadata <- metadata[common_samples, ]
metadata <- metadata %>% arrange(!!sym(group_column))

make_heatmap(metadata, decontaminated_table, 
             samples_column = "sample_id", group_column = "group", 
             output_prefix = "Combined-gene-level-taxonomy_decontam", 
             assay_suffix = "_GLlblMetag",
             custom_palette)

```

**Custom Functions Used:**
- [feature_decontam()](#feature_decontam)
- [make_heatmap()](#make_plot)

**Parameter Definitions:**
- `metadata_table` - path to a file with samples as rows and columns describing each sample
- `feature_table_file` - path to a tab separated samples feature table containing gene-level coverage data 
                         species/functions as the first column and samples as other columns.
- `number_samples` - the total number of samples in the feature_table_file, adjust based on number of input samples

**Input Data:**

- `gene_taxonomy_table.csv`(aggregated gene taxonomy table with samples in columns and species in rows, from [Step 25a](#25a-gene-level-taxonomy-heatmaps))
- `/path/to/sample/metadata` (a file containing sample-wise metadata, mapping samplenames to group metadata)

**Output Data:**

- **Combined-gene-level-taxonomy_decontam_results_GLlblMetag.csv** (decontam's results table)
- **Combined-gene-level-taxonomy_decontam_species_table_GLlblMetag.csv** (decontaminated species table)
- **Combined-gene-level-taxonomy_decontam_heatmap_GLlblMetag.png** (gene-level taxonomy heatmap after filtering out contaminants)

#### 25c. Gene-level KO functions heatmaps

```R
library(tidyverse)
library(pheatmap)

metadata_file <- "/path/to/sample/metadata"
feature_data_file <- "Combined-gene-level-KO-function-coverages-CPM_GLlblMetag.ts"

# Abundant functions with CPM > 2000
abundance_threshold <- 2000

# Prepare metadata
metadata <- read_delim(metadata_file, delim = ",") %>% as.data.frame
sample_names = metadata[, samples_column]
row.names(metadata) <- sample_names

# Read-in KO functions table and drop unannotated
functions_table <- read_delim(file = feature_table_file, delim = "\t", comment = "#") %>%
                   select(KO_ID, KO_function, !!any_of(sample_names)) %>%
                   filter(KO_ID != "Not annotated")

# Convert the sample level data into a matrix
functions.m <- functions_table %>% select(any_of(sample_names)) %>% as.matrix()
rownames(functions.m) <- functions_table$KO_ID

# convert to dataframe without unannotated/unclassified species for output
table2write <- functions.m %>% as.data.frame %>%
               rownames_to_column("KO_ID")
# Write out  taxonomy table
write_csv(x = table2write  , file = "genes-KO-functions_table.csv")

# Get common samples and re-arrange feature table and metadata
common_samples <- intersect(colnames(functions_table), rownames(metadata))
functions_table <- functions_table[, common_samples]
metadata <- metadata[common_samples, ]
metadata <- metadata %>% arrange(!!sym(group_column))

make_heatmap(metadata, table2write,
             samples_column="sample_id", group_column = "group", 
             output_prefix = "Combined-gene-level-KO-function", 
             assay_suffix = "_GLlblMetag", 
             custom_palette = custom_palette)

```

**Custom Functions Used:**
- [make_heatmap()](#make_heatmap)

**Input data:**
- /path/to/sample/metadata (a file with samples as rows and columns describing each sample)
- Combined-gene-level-KO-function-coverages-CPM_GLlblMetag.tsv (table with all samples combined 
    based on KO annotations; normalized to coverage per million genes covered, output from 
    [Step 22a](#22a-generate-gene-level-coverage-summary-tables)

**Output data:**
- genes-KO-functions_table.csv (aggregated and subsetted gene KO function table)
- **Combined-gene-level-KO-function_heatmap_GLlblMetag.png** (heatmap of all gene-level KO function assignments)

#### 25d. Gene-level KO functions decontamination

```R
library(tidyverse)
library(decontam)
library(phyloseq)

feature_table_file <- "genes-KO-functions_table.csv"
metadata_table <- "/path/to/sample/metadata"
number_samples <- NumberOfSamples # integer indicating how many samples are in the file

# set width based on number of samples, with a cap at 50 inches
plot_width <- 2 * number_samples
if(plot_width > 50) { plot_width = 50 }

# Prepare metadata
metadata <- read_delim(metadata_file, delim = ",") %>% as.data.frame
sample_names = metadata[, samples_column]
row.names(metadata) <- sample_names

decontaminated_table <- feature_decontam(metadata_file = metadata_table, 
                                         feature_table_file = feature_table_file, 
                                         feature_column = "KO_ID", 
                                         samples_column = "sample_id",
                                         prevalence_column = "NTC", 
                                         ntc_name = "TRUE", 
                                         frequency_column = "concentration", 
                                         threshold = 0.1, 
                                         classification_method = "Combined-gene-level-KO-function", 
                                         output_prefix = "", 
                                         assay_suffix = "_GLlblMetag")

# Get common samples and re-arrange feature table and metadata
common_samples <- intersect(colnames(decontaminated_table), rownames(metadata))
decontaminated_table <- decontaminated_table[, common_samples]
metadata <- metadata[common_samples, ]
metadata <- metadata %>% arrange(!!sym(group_column))

make_heatmap(metadata, decontaminated_table, 
             samples_column = "sample_id", group_column = "group", 
             output_prefix = "Combined-gene-level-KO-function_decontam", 
             assay_suffix = "_GLlblMetag",
             custom_palette)

```

**Custom Functions Used:**
- [feature_decontam()](#feature_decontam)
- [make_heatmap()](#make_plot)

**Parameter Definitions:**
- `metadata_table` - path to a file with samples as rows and columns describing each sample
- `feature_table_file` - path to a tab separated samples feature table containing gene-level KO functions coverage data 
                         with KO_ID as the first column and samples as other columns.
- `number_samples` - the total number of samples in the feature_table_file, adjust based on number of input samples

**Input Data:**

- `genes-KO-functions_table.csv`(aggregated gene KO functions table table with samples in columns and KO_ID in rows, from [Step 25c](#25c-gene-level-ko-functions-heatmaps))
- `/path/to/sample/metadata` (a file containing sample-wise metadata, mapping samplenames to group metadata)

**Output Data:**

- **Combined-gene-level-KO-function_decontam_results_GLlblMetag.csv** (decontam's results table)
- **Combined-gene-level-KO-function_decontam_species_table_GLlblMetag.csv** (decontaminated gene-level KO functions table)
- **Combined-gene-level-KO-function_decontam_heatmap_GLlblMetag.png** (gene-level KO functions heatmap after filtering out contaminants)


#### 25f. Contig-level Heatmaps

```R
library(tidyverse)

metadata_file <- "/path/to/sample/metadata"
feature_data_file <- "Combined-contig-level-taxonomy-coverages-CPM_GLlblMetag.tsv"

# Prepare metadata
metadata <- read_delim(metadata_file, delim = ",") %>% as.data.frame
sample_names = metadata[, samples_column]
row.names(metadata) <- sample_names

# Prepare feature table
contig_taxonomy_table <- read_assembly_coverage_table(feature_table_file, sample_names) %>% as.data.frame

# Summarize contig table
species_contig_table <- contig_taxonomy_table %>%
  select(species, !!any_of(sample_names)) %>%
  group_by(species) %>%
  summarise(across(everything(), sum)) %>% 
  filter(species != "Unclassified;_;_;_;_;_;_") %>% # Drop unclassifed
  as.data.frame

rownames(species_contig_table) <- species_contig_table[[1]]
species_contig_table <- species_contig_table[, -1] %>% as.matrix()

# Get common samples and re-arrange feature table and metadata
common_samples <- intersect(colnames(species_contig_table), rownames(metadata))
species_contig_table <- species_contig_table[, common_samples]
metadata <- metadata[common_samples, ]
metadata <- metadata %>% arrange(!!sym(group_column))

table2write = species_contig_table %>% as.data.frame %>% rownames_to_column("species")
# Write out contig taxonomy table
write_csv(x = table2write, file = "contig_taxonomy_table.csv")

make_heatmap(metadata, species_contig_table, 
             samples_column="sample_id", group_column = "group", 
             output_prefix = "Combined-contig-level-taxonomy", 
             assay_suffix = "_GLlblMetag", 
             custom_palette = custom_palette)
```

**Custom Functions Used:**
- [read_assembly_coverage_table()](#read_assembly_coverage_table)
- [make_heatmap()](#make_heatmap)

**Input data:**
- /path/to/sample/metadata (a file with samples as rows and columns describing each sample)
- Combined-contig-level-taxonomy-coverages-CPM_GLlblMetag.tsv (table with all samples 
    combined based on contig-level taxonomic classifications, output from 
    [Step 21b](#21b-generate-contig-level-coverage-summary-tables)) 

**Output data:**
- contig_taxonomy_table.csv (aggregated contig taxonomy table with samples in columns and species in rows)
- **Combined-contig-level-taxonomy_heatmap_GLlblMetag.png** (heatmap of all contig taxonomy assignments)

#### 25g. Contig-level decontamination

```R
library(tidyverse)
library(decontam)
library(phyloseq)

feature_table_file <- "contig_taxonomy_table.csv"
metadata_table <- "/path/to/sample/metadata"
number_samples <- NumberOfSamples # integer indicating how many samples are in the file

# set width based on number of samples, with a cap at 50 inches
plot_width <- 2 * number_samples
if(plot_width > 50) { plot_width = 50 }

# Prepare metadata
metadata <- read_delim(metadata_file, delim = ",") %>% as.data.frame
sample_names = metadata[, samples_column]
row.names(metadata) <- sample_names

decontaminated_table <- feature_decontam(metadata_file = metadata_table, 
                                         feature_table_file = feature_table_file, 
                                         feature_column = "species", 
                                         samples_column = "sample_id",
                                         prevalence_column = "NTC", 
                                         ntc_name = "TRUE", 
                                         frequency_column = "concentration", 
                                         threshold = 0.1, 
                                         classification_method = "Combined-contig-level-taxonomy", 
                                         output_prefix = "", 
                                         assay_suffix = "_GLlblMetag")

# Get common samples and re-arrange feature table and metadata
common_samples <- intersect(colnames(decontaminated_table), rownames(metadata))
decontaminated_table <- decontaminated_table[, common_samples]
metadata <- metadata[common_samples, ]
metadata <- metadata %>% arrange(!!sym(group_column))

make_heatmap(metadata, decontaminated_table, 
             samples_column = "sample_id", group_column = "group", 
             output_prefix = "Combined-contig-level-taxonomy_decontam", 
             assay_suffix = "_GLlblMetag",
             custom_palette)

```

**Custom Functions Used:**
- [feature_decontam()](#feature_decontam)
- [make_heatmap()](#make_plot)

**Parameter Definitions:**
- `metadata_table` - path to a file with samples as rows and columns describing each sample
- `feature_table_file` - path to a tab separated samples feature table containing contig-level coverage data 
                         species/functions as the first column and samples as other columns.
- `number_samples` - the total number of samples in the feature_table_file, adjust based on number of input samples

**Input Data:**

- `contig_taxonomy_table.csv`(aggregated contig taxonomy table with samples in columns and species in rows, from [Step 25f](#25f-contig-level-heatmaps))
- `/path/to/sample/metadata` (a file containing sample-wise metadata, mapping samplenames to group metadata)

**Output Data:**

- **Combined-contig-level-taxonomy_decontam_results_GLlblMetag.csv** (decontam's results table)
- **Combined-contig-level-taxonomy_decontam_species_table_GLlblMetag.csv** (decontaminated contig-level species table)
- **Combined-contig-level-taxonomy_decontam_heatmap_GLlblMetag.png** (contig-level heatmap after filtering out contaminants)

