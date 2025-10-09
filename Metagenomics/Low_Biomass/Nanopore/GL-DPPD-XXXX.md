# Bioinformatics pipeline for Low biomass long-read metagenomics data

> **This document holds an overview and some example commands of how GeneLab processes low-biomass, long-read metagenomics datasets. Exact processing commands for specific datasets that have been released are provided with their processed data in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/).**  

---

**Date:** XXX NN, 2025  
**Revision:** -  
**Document Number:** GL-DPPD-XXXX  

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
      - [5c. Compile Trimmed Data QC](#5c-compile-filtered-data-qc)
    - [6. Contaminant Removal](#6-contaminant-removal)
      - [6a. Assemble Contaminants](#6a-assemble-contaminants)
      - [6b. Build Contaminant Index and Map Reads](#6b-build-contaminant-index-and-map-reads)
      - [6c. Sort and Index Contaminant Reads](#6c-sort-and-index-contaminant-alignments)
      - [6d. Gather Contaminant Mapping Metrics](#6d-gather-contaminant-mapping-metrics)
      - [6e. Generate Decontaminated Read Files](#6e-generate-decontaminated-read-files)
      - [6f. Contaminant Removal QC](#6f-contaminant-removal-qc)
      - [6g. Compile Contaminant Removal QC](#6g-compile-contaminant-removal-qc)
    - [7. Host Removal](#7-host-removal)
      - [7a. Build or download host database](#7a-build-or-download-host-database)
        - [7a.i. Download from URL](#7ai-download-from-url)
        - [7a.ii. Build from custom reference](#7aii-build-from-custom-reference)
        - [7a.iii. Build from host name](#7aiii-build-from-host-name)
      - [7b. Remove Host Reads](#7b-remove-host-reads)
    - [8. R Environment Setup](#8-r-environment-setup)
      - [8a. Load libraries](#8a-load-libraries)
      - [8b. Define Custom Functions](#8b-define-custom-functions)
      - [8c. Set global variables](#8c-set-global-variables)
  - [**Read-based processing**](#read-based-processing)
    - [9. Taxonomic profiling using kaiju](#9-taxonomic-profiling-using-kaiju)
      - [9a. Build kaiju database](#9a-build-kaiju-database)
      - [9b. Kaiju Taxonomic Classification](#9b-kaiju-taxonomic-classification)
      - [9c. Compile kaiju taxonomy results](#9c-compile-kaiju-taxonomy-results)
      - [9d. Convert kaiju output to krona format](#9d-convert-kaiju-output-to-krona-format)
      - [9e. Compile kaiju krona report](#9e-compile-kaiju-krona-report)
      - [9f. Create kaiju species count table](#9f-create-kaiju-species-count-table)
      - [9g. Read-in tables](#9g-read-in-tables)
      - [9h. Taxonomy barplots](#9h-taxonomy-barplots)
      - [9i. Feature decontamination](#9i-feature-decontamination)
    - [10. Taxonomic Profiling using Kraken2](#10-taxonomic-profiling-using-kraken2)
      - [10a. Download kraken2 database](#10a-download-kraken2-database)
      - [10b. Taxonomic Classification](#10b-taxonomic-classification)
      - [10c. Convert Kraken2 output to Krona format](#10c-convert-kraken2-output-to-krona-format)
      - [10d. Compile kraken2 krona report](#10d-compile-kraken2-krona-report)
      - [10e. Create kraken species count table](#10e-create-kraken-species-count-table)
      - [10f. Read-in tables](#10f-read-in-tables)
      - [10g. Taxonomy barplots](#10g-taxonomy-barplots)
      - [10h. Feature decontamination](#10h-feature-decontamination)
  - [**Assembly-based processing**](#assembly-based-processing)
    - [11. Sample assembly](#11-sample-assembly)
    - [12. Polish assembly](#12-polish-assembly)
    - [13. Renaming contigs and summarizing assemblies](#13-renaming-contigs-and-summarizing-assemblies)
      - [13a. Renaming contig headers](#13a-renaming-contig-headers)
      - [13b. Summarizing assemblies](#13b-summarizing-assemblies)
    - [14. Gene prediction](#14-gene-prediction)
      - [14a. Remove line wraps in gene prediction output](#14a-remove-line-wraps-in-gene-prediction-output)
    - [15. Functional annotation](#15-functional-annotation)
      - [15a. Downloading reference database of HMM models (only needs to be done once)](#15a-downloading-reference-database-of-hmm-models-only-needs-to-be-done-once)
      - [15b. Running KEGG annotation](#15b-running-kegg-annotation)
      - [15c. Filtering output to retain only those passing the KO-specific score and top hits](#15c-filtering-output-to-retain-only-those-passing-the-ko-specific-score-and-top-hits)
    - [16. Taxonomic classification](#16-taxonomic-classification)
      - [16a. Pulling and un-packing pre-built reference db (only needs to be done once)](#16a-pulling-and-un-packing-pre-built-reference-db-only-needs-to-be-done-once)
      - [16b. Running taxonomic classification](#16b-running-taxonomic-classification)
      - [16c. Adding taxonomy info from taxids to genes](#16c-adding-taxonomy-info-from-taxids-to-genes)
      - [16d. Adding taxonomy info from taxids to contigs](#16d-adding-taxonomy-info-from-taxids-to-contigs)
      - [16e. Formatting gene-level output with awk and sed](#16e-formatting-gene-level-output-with-awk-and-sed)
      - [16f. Formatting contig-level output with awk and sed](#16f-formatting-contig-level-output-with-awk-and-sed)
    - [17. Read-mapping](#17-read-mapping)
      - [17a. Align Reads to Sample Assembly](#17a-align-reads-to-sample-assembly)
      - [17b. Sort and Index Assembly Alignments](#17b-sort-and-index-assembly-alignments)
    - [18. Getting coverage information and filtering based on detection](#18-getting-coverage-information-and-filtering-based-on-detection)
      - [18a. Filtering coverage levels based on detection](#18a-filtering-coverage-levels-based-on-detection)
      - [18b. Filtering gene and contig coverage based on requiring 50% detection and parsing down to just gene ID and coverage](#18b-filtering-gene-and-contig-coverage-based-on-requiring-50-detection-and-parsing-down-to-just-gene-id-and-coverage)
    - [19. Combining gene-level coverage, taxonomy, and functional annotations into one table for each sample](#19-combining-gene-level-coverage-taxonomy-and-functional-annotations-into-one-table-for-each-sample)
    - [20. Combining contig-level coverage and taxonomy into one table for each sample](#20-combining-contig-level-coverage-and-taxonomy-into-one-table-for-each-sample)
    - [21. Generating normalized, gene- and contig-level coverage summary tables of KO-annotations and taxonomy across samples](#21-generating-normalized-gene--and-contig-level-coverage-summary-tables-of-ko-annotations-and-taxonomy-across-samples)
      - [21a. Generating gene-level coverage summary tables](#21a-generating-gene-level-coverage-summary-tables)
      - [21b. Gene-level taxonomy heatmaps](#21b-gene-level-taxonomy-heatmaps)
      - [21c. Gene-level taxonomy decontamination](#21c-gene-level-taxonomy-decontamination)
      - [21d. Gene-level KO functions heatmaps](#21d-gene-level-ko-functions-heatmaps)
      - [21e. Gene-level KO functions decontamination](#21e-gene-level-ko-functions-decontamination)
      - [21f. Generating contig-level coverage summary tables](#21f-generating-contig-level-coverage-summary-tables)
      - [21g. Contig-level Heatmaps](#21g-contig-level-heatmaps)
      - [21h. Contig-level decontamination](#21h-contig-level-decontamination)
    - [22. **M**etagenome-**A**ssembled **G**enome (MAG) recovery](#22-metagenome-assembled-genome-mag-recovery)
      - [22a. Binning contigs](#22a-binning-contigs)
      - [22b. Bin quality assessment](#22b-bin-quality-assessment)
      - [22c. Filtering MAGs](#22c-filtering-mags)
      - [22d. MAG taxonomic classification](#22d-mag-taxonomic-classification)
      - [22e. Generating overview table of all MAGs](#22e-generating-overview-table-of-all-mags)
    - [23. Generating MAG-level functional summary overview](#23-generating-mag-level-functional-summary-overview)
      - [23a. Getting KO annotations per MAG](#23a-getting-ko-annotations-per-mag)
      - [23b. Summarizing KO annotations with KEGG-Decoder](#23b-summarizing-ko-annotations-with-kegg-decoder)



---

# Software used

|Program|Version|Relevant Links|
|:------|:-----:|------:|
|bbduk| 38.86 |[https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/)|
|CAT| 5.2.3 |[https://github.com/dutilh/CAT#cat-and-bat](https://github.com/dutilh/CAT#cat-and-bat)|
|CheckM| 1.1.3 |[https://github.com/Ecogenomics/CheckM](https://github.com/Ecogenomics/CheckM)|
|Dorado| 1.1.1| [https://github.com/nanoporetech/dorado](https://github.com/nanoporetech/dorado)|
|filtlong| 0.2.1 |[https://github.com/rrwick/Filtlong](https://github.com/rrwick/Filtlong)|
|Flye| 2.9.5 | [https://github.com/mikolmogorov/Flye](https://github.com/mikolmogorov/Flye) |
|GTDB-Tk| 2.4.0 |[https://github.com/Ecogenomics/GTDBTk](https://github.com/Ecogenomics/GTDBTk)|
|Kaiju| 1.10.1 | [https://bioinformatics-centre.github.io/kaiju/](https://bioinformatics-centre.github.io/kaiju/) |
|KEGG-Decoder| 1.2.2 |[https://github.com/bjtully/BioData/tree/master/KEGGDecoder#kegg-decoder](https://github.com/bjtully/BioData/tree/master/KEGGDecoder#kegg-decoder)
|KOFamScan| 1.3.0 |[https://github.com/takaram/kofam_scan](https://github.com/takaram/kofam_scan)|
|Kraken2| 2.1.6 | [https://github.com/DerrickWood/kraken2](https://github.com/DerrickWood/kraken2) |
|KrakenTools | 1.2 | [https://ccb.jhu.edu/software/krakentools/](https://ccb.jhu.edu/software/krakentools/) |
|Krona| 2.8.1 | [https://github.com/marbl/Krona/wiki](https://github.com/marbl/Krona/wiki)|
|MetaBAT| 2.15 |[https://bitbucket.org/berkeleylab/metabat/src/master/](https://bitbucket.org/berkeleylab/metabat/src/master/)|
|Minimap2| 2.2.8 | [https://github.com/lh3/minimap2](https://github.com/lh3/minimap2) |
|MultiQC| 1.27.1 |[https://multiqc.info/](https://multiqc.info/)|
|Medaka| 2.0.1 | [https://github.com/nanoporetech/medaka](https://github.com/nanoporetech/medaka) |
|NanoPlot| 1.44.1 | [https://github.com/wdecoster/NanoPlot](https://github.com/wdecoster/NanoPlot)|
|Porechop| 0.2.4 | [https://github.com/rrwick/Porechop](https://github.com/rrwick/Porechop) |
|Prodigal| 2.6.3 |[https://github.com/hyattpd/Prodigal#prodigal](https://github.com/hyattpd/Prodigal#prodigal)|
|samtools| 1.20 |[https://github.com/samtools/samtools#samtools](https://github.com/samtools/samtools#samtools)|
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
  --min-qscore 7 > basecalled.bam
```

**Parameter Definitions:**

- `--no-trim` - skips trimming of barcodes, adapters, and primers
- `--device` - specifies CPU or GPU device; specifying 'auto' chooses either 'cpu' or 'gpu' depending on detected presence of a GPU device
- `--recursive` - enables recursive scanning through input directory to load FAST5 and/or POD5 files
- `--kit-name` - The nanopore barcoding kit used during sequencing preparation. Enables barcoding with the provided kit name; see [dorado documentation](https://software-docs.nanoporetech.com/dorado/1.1.1/barcoding/barcoding/) for a full list of accepted kit names
- `--min-qscore` - specifies the minimum Q-score, reads with a mean Q-score below this threshold are discarded (default to `7` for this pipeline)
- `model` - positional argument specifying the basecalling model to use or a path to the model directory. `hac` chooses the high accuracy model.
- `input_directory` - positional argument specifying the location of the raw data in POD5 or FAST5 format

**Input Data:**

- *pod5 and/or *fast5 (raw nanopore data)

**Output Data:**

- basecalled.bam (basecalled data in bam format)

<br>

---

### 2. Demultiplexing

#### 2a. Split fastq

```bash
dorado demux \
  --output-dir /path/to/fastq/output \
  --emit-fastq \
  --emit-summary \
  --kit-name ${kit_name} \
  basecalled.bam
```

**Parameter Definitions:**

- `--output-dir` - specifies the output folder that is the root of the nested output structure. 
- `--emit-fastq` - specifies that output is fastq format
- `--emit-summary` - creates a summary listing each read and its classified barcode.
- `--kit-name` - The nanopore barcoding kit used during sequencing preparation. Enables barcoding with the provided kit name; see [dorado documentation](https://software-docs.nanoporetech.com/dorado/1.1.1/barcoding/barcoding/) for a full list of accepted kit names

**Input Data:**

- basecalled.bam (basecalled nanopore data in bam format, output from [Step 1](#1-basecalling))

**Output Data:**

- /path/to/fastq/output/\*_barcode\*.fastq (demultiplexed reads in fastq format)
- /path/to/fastq/output/\*_unclassified.fastq (unclassified reads in fastq format)
- /path/to/fastq/output/barcoding_summary.txt (barcode summary file listing each read, the file it was assigned to, and its classified barcode )


#### 2b. Concatenate files for each sample

```bash
# Change to directory containing split fastq files generated from step 2a. split fastq above
cd /path/to/fastq/output/ # output of step 2a
# Get unique barcode names from demultiplexed file names
BARCODES=($(ls -1 *fastq* |sed -E 's/.+_(barcode[0-9]+)_.+/\1/g' | sort -u))

# Concat separate barcode/sample fastq files into per sample fastq gzippped files
[ -d raw_data/ ] || mkdir raw_data/
for sample in ${BARCODES[*]}; do

  [ -d  ${sample}/ ] ||  mkdir ${sample}/  
  mv *_${sample}_*  ${sample}/ 

  cat ${sample}/* | gzip --to-stdout raw_data/${sample}.fastq.gz

done
```

**Parameter Definitions:**

- `| gzip --to-stdout` - sends output from `cat` to `gzip` to create compressed fastq.gz file

**Input Data:**

- /path/to/fastq/output/ (directory containing spilt fastq files from [Step 2a](#2a-split-fastq))

**Output Data:**

-  raw_data/sample.fastq.gz (gzipped per sample/barcode fastq files)

<br>

---

### 3.  Raw Data QC

#### 3a. Raw Data QC

```bash 
NanoPlot --only-report \
         --prefix sample_ \
         --outdir /path/to/raw_nanoplot_output \
         --threads NumberOfThreads \
         --fastq /path/to/raw_data/sample.fastq.gz
```

**Parameter Definitions:**

- `--outdir` – specifies the output directory to store results
- `--only-report` - output only the report files
- `--prefix` - adds a sample specific prefix to the name of each output file
- `--threads` - number of parallel processing threads to use
- `--fastq` - specifies that the input data is in a fastq format
- `/path/to/raw_data/sample.fastq.gz` – the input reads are specified as a positional argument

**Input Data:**

- /path/to/raw_data/sample.fastq.gz (raw reads, output from [Step 2b](#2b-concatenate-files-for-each-sample))

**Output Data:**

- **/path/to/raw_nanoplot_output/sample_NanoPlot-report.html** (NanoPlot html summary)
- /path/to/raw_nanoplot_output/sample_NanoPlot_<date>_<time>.log (NanoPlot log file)
- /path/to/raw_nanoplot_output/sample_NanoStats.txt (text file containing basic statistics)

#### 3b. Compile Raw Data QC

```bash 
multiqc --zip-data-dir \
        --outdir raw_multiqc_report \
        --filename raw_multiqc \
        --interactive /path/to/raw_nanoplot_output/
```

**Parameter Definitions:**

- `--zip-data-dir` - compress the data directory
- `--outdir` – the output directory to store results
- `--filename` – the filename prefix of results
- `--interactive` - force multiqc to always create interactive javascript plots
- `/path/to/raw_nanoplot_output/` – the directory holding the output data from the NanoPlot run, provided as a positional argument

**Input Data:**

- /path/to/raw_nanoplot_output/*NanoStats.txt (NanoPlot output data, from [Step 3a](#3a-raw-data-qc))

**Output Data:**

- **raw_multiqc.html** (multiqc output html summary)
- **raw_multiqc_data.zip** (zip archive containing multiqc output data)

<br>  

---

### 4. Quality filtering

#### 4a. Filter Raw Data

```bash
filtlong --min_length 200 --min_mean_q 8 /path/to/raw_data/sample.fastq.gz > sample_filtered.fastq
```

**Parameter Definitions:**

- `--min_length` – specifies the minimum read length to retain (default to `200` for this pipeline)
- `--min_mean_q` – specifies the minimum mean read quality (default to `8` for this pipeline)

**Input Data:**

- /path/to/raw_data/sample.fastq.gz (raw reads, output from [Step 2b](#2b-concatenate-files-for-each-sample))

**Output Data:**

- *sample_filtered.fastq (quality filtered reads)


#### 4b. Filtered Data QC

```bash
NanoPlot --only-report \
         --prefix sample_ \
         --outdir /path/to/filtered_nanoplot_output \
         --threads NumberOfThreads \
         --fastq sample_filtered.fastq
```

**Parameter Definitions:**

- `--outdir` – specifies the output directory to store results
- `--only-report` - output only the report files
- `--prefix` - adds a sample specific prefix to the name of each output file
- `--threads` - number of parallel processing threads to use
- `sample_filtered.fastq` – the input reads are specified as a positional argument

**Input Data:**

- sample_filtered.fastq (filtered reads, output from [Step 4a](#4a-filter-raw-data))

**Output Data:**

- **/path/to/filtered_nanoplot_output/sample_NanoPlot-report.html** (NanoPlot html summary)
- /path/to/filtered_nanoplot_output/sample_NanoPlot_<date>_<time>.log (NanoPlot log file)
- /path/to/filtered_nanoplot_output/sample_NanoStats.txt (text file containing basic statistics)

#### 4c. Compile Filtered Data QC

```bash
multiqc  --zip-data-dir \ 
         --outdir filtered_multiqc_report \
         --filename filtered_multiqc \
         --interactive /path/to/filtered_nanoplot_output/
```

**Parameter Definitions:**

- `--zip-data-dir` - compress the data directory
- `--outdir` – the output directory to store results
- `--filename` – the filename prefix of results
- `--interactive` - force multiqc to always create interactive javascript plots
- `/path/to/filtered_nanoplot_output/` – the directory holding the output data from the NanoPlot run, provided as a positional argument

**Input Data:**

- /path/to/filtered_nanoplot_output/*NanoStats.txt (NanoPlot output data, from [Step 4b](#4b-filtered-data-qc))

**Output Data:**

- **filtered_multiqc_report/filtered_multiqc.html** (multiqc output html summary)
- **filtered_multiqc_report/filtered_multiqc_data.zip** (zip archive containing multiqc output data)

<br>

---

### 5. Trimming

#### 5a. Trim Filtered Data

```bash
porechop --input sample_filtered.fastq --threads NumberOfThreads \
         --discard_middle --output sample_trimmed.fastq  > sample_porechop.log
```

**Parameter Definitions:**

- `--input` – the input read file in fastq format
- `--threads` - number of parallel processing threads to use
- `--discard_middle` -  reads with middle adapters will be discarded
- `--output` - trimmed reads output fastq filename
- `> sample_porechop.log` - capture stdout in a log file

**Input Data:**

- sample_filtered.fastq (filtered reads output from [Step 4a](#4a-filter-raw-data))

**Output Data:**

- **sample_trimmed.fastq** (filtered and trimmed reads)

#### 5b. Trimmed Data QC

```bash
NanoPlot --only-report \
         --prefix sample_ \
         --outdir /path/to/trimmed_nanoplot_output \
         --threads NumberOfThreads \
         --fastq sample_trimmed.fastq
```

**Parameter Definitions:**

- `--outdir` – specifies the output directory to store results
- `--only-report` - output only the report files
- `--prefix` - adds a sample specific prefix to the name of each output file
- `--threads` - number of parallel processing threads to use
- `sample_trimmed.fastq` – the input reads are specified as a positional argument

**Input Data:**

- sample_trimmed.fastq (filtered and trimmed reads, output from [Step 5a](#5a-trim-filtered-data))

**Output Data:**

- **/path/to/trimmed_nanoplot_output/sample_NanoPlot-report.html** (NanoPlot html summary)
- /path/to/trimmed_nanoplot_output/sample_NanoPlot_<date>_<time>.log (NanoPlot log file)
- /path/to/trimmed_nanoplot_output/sample_NanoStats.txt (text file containing basic statistics)

#### 5c. Compile Trimmed Data QC

```bash
multiqc --zip-data-dir \ 
        --outdir trimmed_multiqc_report \
        --filename trimmed_multiqc \
        --interactive /path/to/trimmed_nanoplot_output/
```

**Parameter Definitions:**

- `--zip-data-dir` - compress the data directory
- `--outdir` – the output directory to store results
- `--filename` – the filename prefix of results
- `--interactive` - force multiqc to always create interactive javascript plots
- `/path/to/trimmed_nanoplot_output/` – the directory holding the output data from the NanoPlot run, provided as a positional argument

**Input Data:**

- /path/to/trimmed_nanoplot_output/*NanoStats.txt (NanoPlot output data, output from [Step 5b](#5b-trimmed-data-qc))

**Output Data:**

- **trimmed_multiqc.html** (multiqc output html summary)
- **trimmed_multiqc_data.zip** (zip archive containing multiqc output data)

<br>

---

### 6. Contaminant Removal

> A major issue with low biomass data is the high potential for contamination due to the low amount of DNA extracted in the sample for sequencing.  Because negative control/blank samples should by theory be contaminant free, any sequence detected in the negative control is a potential contaminant. To filter out contaminants found in negative control samples that may have been due to cross contamination in the lab, we use a read mapping approach. First negative/blank control samples reads are assembled then filtered and trimmed reads mapped to the assembled contigs. Reads mapping to the assembled contigs are categorized as contaminants and are therefore filtered out and thus excluded from further analyses.

### 6a. Assemble Contaminants

```bash
flye --meta --threads NumberOfThreads \
     --out-dir /path/to/contaminant_assembly \
     --nano-raw /path/to/blank_samples/\*_trimmed.fastq
```

**Parameter Definitions:**

- `--meta` – use metagenome/uneven coverage mode
- `--threads` - number of parallel processing threads to use
- `--out-dir` - output directory
- `--nano-raw` - specifies that input is from Oxford Nanopore regular raw reads. This adds a polishing step for error correction after the assembly is generated.

**Input Data**

- *_trimmed.fastq (one or more trimmed reads from blank samples, output from [Step 5a](#5a-trim-filtered-data))

**Output Data**

- /path/to/contaminant_assembly/assembly.fasta (Assembly built from reads in blank samples in fasta format)

<br>

---

#### 6b. Build Contaminant Index and Map Reads

```bash
# Build contaminant index
minimap2 -t NumberOfThreads -a -x splice -d blanks.mmi /path/to/contaminant_assembly/assembly.fasta

# Map reads to index
minimap2 -t NumberOfThreads -a -x splice blanks.mmi /path/to/trimmed_reads/sample_trimmed.fastq  > sample.sam
```

**Parameter Definitions:**

- `-t` - number of parallel processing threads
- `-a` – output in SAM format
- `-x splice` - specifies preset for spliced alignment of long reads
- `-d` - specifies the output file for the index

**Input Data**

- /path/to/contaminant_assembly/assembly.fasta (Contaminant assembly, output from [Step 6a](#6-assemble-contaminants))
- /path/to/trimmed_reads/sample_trimmed.fastq (Filtered and trimmed reads, output from [Step 5a](#5a-trim-filtered-data))

**Output Data**

- sample.sam (Reads aligned to contaminant assembly)

#### 6c. Sort and Index Contaminant Alignments
```bash
# Sort Sam, convert to bam and create index
samtools sort --threads NumberOfThreads -o sample_sorted.bam sample.sam > sample_sort.log 2>&1

samtools index sample_sorted.bam sample_sorted.bam.bai
```

**Parameter Definitions:**

**samtools sort**
- `--threads` - number of parallel processing threads to use
- `-o` - specifies the output file for the sorted reads
- `sample.sam` - positional argument specifying the input SAM file

**samtools index**
- `sample_sorted.bam` - positional argument specifying the input BAM file to be indexed
- `sample_sorted.bam.bai` - positional argument specifying the name of the index file

**Input Data:**

- sample.sam (Reads aligned to contaminant assembly, output from [Step 6b](#6b-build-contaminant-index-and-map-reads))

**Output Data:**

- sample_sorted.bam (sorted mapping to contaminant assembly)
- sample_sorted.bam.bai (index of sorted mapping to contaminant assembly)

#### 6d. Gather Contaminant Mapping Metrics

```bash

samtools flagstat sample_sorted.bam > sample_flagstats.txt  2> sample_flagstats.log
samtools stats --remove-dups sample_sorted.bam > sample_stats.txt   2> sample_stats.log
samtools idxstats sample_sorted.bam  > sample_idxstats.txt 2> sample_idxstats.log
```

**Parameter Definitions:**

- `flagstat` - positional argument specifying the program for counting the number of alignments for each SAM FLAG type
- `stats` - positional argument specifying the program for producing comprehensive statistics from the alignment file
- `idxstats` - positional argument specifying the program for producing contig alignment summary statistics
- `--remove-dups` - excludes reads marked as duplicates from comprehensive statistics
- `sample_sorted.bam` - positional argument specifying the input BAM file

**Input Data:**

- sample_sorted.bam (sorted mapping to contaminant assembly, output from [Step 6c](#6c-sort-and-index-contaminant-alignments))
- sample_sorted.bam.bai (index of sorted mapping to contaminant assembly, output from [Step 6c](#6c-sort-and-index-contaminant-alignments))

**Output Data:**

- sample_flagstats.txt (SAM FLAG counts)
- sample_stats.txt (comprehensive alignment statistics)
- sample_idxstats.txt (contig alignment summary statistics)

#### 6e. Generate Decontaminated Read Files
```bash
# Retain reads that do not map to contaminants
samtools fastq -t -f 4 sample_sorted.bam | gzip --to-stdout > sample_blank_removed.fastq.gz
```

**Parameter Definitions:**

- `fastq` - positional argument specifying the program for generating fastq files from a SAM/BAM file
- `-t` - copy RG, BC, and QT tags to the FASTQ header line
- `-f 4` - only retain unmapped reads that have been marked with the SAM "segment unmapped" FLAG (4)
- `sample_sorted.bam` - positional argument specifying the input BAM file
- `| gzip --to-stdout` - sends output from `samtools fastq` to `gzip` to create compressed fastq.gz file
- `> sample_blank_removed.fastq.gz` - specifies the name of the file used to store the fastq.gz output

**Input Data:**

- sample_sorted.bam (sorted mapping to contaminant assembly, output from [Step 6c](#6c-sort-and-index-contaminant-alignments))

**Output Data:**

- sample_blank_removed.fastq.gz (blank removed reads in fastq format)

#### 6f. Contaminant Removal QC

```bash
NanoPlot --only-report \
         --prefix sample_ \
         --outdir /path/to/noblank_nanoplot_output \
         --threads NumberOfThreads \
         --fastq sample_blank_removed.fastq.gz
```

**Parameter Definitions:**

- `--outdir` – specifies the output directory to store results
- `--only-report` - output only the report files
- `--prefix` - adds a sample specific prefix to the name of each output file
- `--threads` - number of parallel processing threads to use
- `--fastq` - specifies that the input data is in a fastq format
- `sample_blank_removed.fastq.gz` – the input reads are specified as a positional argument

**Input Data:**

- sample_blank_removed.fastq.gz (blank removed reads, output from [Step 6e](#6e-generate-decontaminated-read-files))

**Output Data:**

- **/path/to/noblank_nanoplot_output/sample_NanoPlot-report.html** (NanoPlot html summary)
- /path/to/noblank_nanoplot_output/sample_NanoPlot_<date>_<time>.log (NanoPlot log file)
- /path/to/noblank_nanoplot_output/sample_NanoStats.txt (text file containing basic statistics)


#### 6g. Compile Contaminant Removal QC

```bash
multiqc --zip-data-dir \ 
        --outdir noblank_multiqc_report \
        --filename noblank_multiqc \
        --interactive /path/to/noblank_nanoplot_output/
```

**Parameter Definitions:**

- `--zip-data-dir` - compress the data directory
- `--outdir` – the output directory to store results
- `--filename` – the filename prefix of results
- `--interactive` - force multiqc to always create interactive javascript plots
- `/path/to/noblank_nanoplot_output/` – the directory holding the output data from the NanoPlot run, provided as a positional argument

**Input Data:**

- /path/to/noblank_nanoplot_output/*NanoStats.txt (NanoPlot output data, output from [Step 6f](#6f-contaminant-removal-qc))

**Output Data:**

- **noblank_multiqc_report/noblank_multiqc.html** (multiqc output html summary)
- **noblank_multiqc_report/noblank_multiqc_data.zip** (zip archive containing multiqc output data)

<br>

---

### 7. Host Removal

#### 7a. Build or download host database

##### 7a.i. Download from URL

```bash
  # Downloading and unpacking database from ${host_url}
  wget -O host.tar.gz --timeout=3600 --tries=0 --continue  host_url

  mkdir kraken2_host_db/ && \
  tar -zxvf host.tar.gz -C kraken2_host_db/ && \
  rm -rf  host.tar.gz # Cleaning up
```

**Parameter Definitions:**

- `--timeout` - network timeout in seconds
- `--tries` - number of times to retry the download
- `--continue` - continue getting a partially downloaded file (if it exists)
- `host_url` - positional argument specifying the URl for the host database

**Output Data:**

- kraken2_host_db/ - Kraken2 database directory


##### 7a.ii. Build from custom reference

```bash
# Install taxonomy       
kraken2-build --download-taxonomy --db kraken2_host_db/
# Add sequence to your database's genomic library
kraken2-build --add-to-library host_assembly.fasta --db kraken2_host_db/ --no-masking
# Once your library is finalized, build the database
kraken2-build --build --db kraken2_host_db/
```

**Parameter Definitions:**

- `--download-taxonomy` - downloads taxonomic mapping information
- `--add-to-library host_assembly.fasta` - specifies to add assembly fasta to library
- `--db` - specifies the output directory for the kraken database
- `--build` - specifies to construct kraken2-formatted database

**Input Data:**

- `host_assembly.fasta` - host genome assembly in fasta format 

**Output Data:**

- kraken2_host_db/ - Kraken2 database directory


##### 7a.iii. Build from host name

```bash
# Build kraken reference from host_name
kraken2-build --download-library host_name  -db kraken2_host_db/ \
              --threads numberOfThreads  --no-masking
kraken2-build --download-taxonomy --db kraken2_host_db/
kraken2-build --build --db kraken2_host_db/ --threads numberOfThreads 
kraken2-build --clean --db kraken2_host_db/
```

**Parameter Definitions:**

- `--download-library` - specifies the reference name/type to download, host_name must 
                         be one of: "archaea", "bacteria", "plasmid", "viral", "human", 
                         "fungi", "plant", "protozoa", "nr", "nt", "UniVec", "UniVec_Core"
- `--db` - specifies the directory we are putting the database in
- `--threads` - number of parallel processing threads to use
- `--no-masking` - prevents masking of low-complexity sequences. For additional 
                   information see the [kraken documentation for masking](https://github.com/DerrickWood/kraken2/wiki/Manual#masking-of-low-complexity-sequences)
- `--download-taxonomy` - downloads taxonomic mapping information
- `--build` - specifies to construct kraken2-formatted database
- `--clean` - specifies to remove unnecessarily intermediate files

**Input Data:**

- `host_name` - host database name (one of those specified in `--download-library` above)

**Output Data:**

- kraken2_host_db/ - Kraken2 database directory

#### 7b. Remove host reads

```bash
kraken2 --db kraken2_host_db/ --gzip-compressed --threads NumberOfThreads --use-names \
        --output sample-kraken2-output.txt --report sample-kraken2-report.tsv \
        --unclassified-out sample_host_removed.fastq sample_blank_removed.fastq.gz
gzip sample_host_removed.fastq
```

**Parameter Definitions:**

- `--db` - specifies the directory holding the kraken2 database files created in [Step 7a](#7a-build-or-download-host-database)
- `--gzip-compressed` - specifies the input fastq files are gzip-compressed
- `--threads` - number of parallel processing threads to use
- `--use-names` - specifies adding taxa names in addition to taxon IDs
- `--output` - specifies the name of the kraken2 read-based output file (one line per read)
- `--report` - specifies the name of the kraken2 report output file (one line per taxa, with number of reads assigned to it)
- `--unclassified-out` - name of output file of reads that were not classified i.e non-host reads.
- `sample_blank_removed.fastq.gz` - positional argument specifying the input read file

**Input Data:**

- sample_blank_removed.fastq.gz (gzipped blank removed fastq file, output from [Step 6d](#6d-generate-decontaminated-read-files))

**Output Data:**

- sample-kraken2-output.txt (kraken2 read-based output file (one line per read))
- sample-kraken2-report.tsv (kraken2 report output file (one line per taxa, with number of reads assigned to it))
- **sample_host_removed.fastq.gz** (host-read removed, gzipped fastq file)

<br>

---


### 8. R Environment Setup

> Taxonomy bar plots, heatmaps and feature decontamination with decontam are performed in R.

#### 8a. Load libraries

```R
library(decontam)
library(phyloseq)
library(tidyverse)
library(pheatmap)
library(pavian)
```

#### 8b. Define Custom Functions

##### get_last_assignment()
<details>
  <summary>retrieves the last taxonomy assignment from a taxonomy string</summary>

  ```R
  get_last_assignment <- function(taxonomy_string, split_by=';', remove_prefix=NULL) {

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
      level_name <- gsub(pattern = remove_prefix, replacement = '', x = level_name)
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
    colnames(df)[col_index] <- 'taxonomy'
    df <- df %>% dplyr::mutate(across( where(is.numeric), function(x) tidyr::replace_na(x,0)  ) )%>% 
      dplyr::mutate(taxonomy=map_chr(taxonomy,.f = function(taxon_name=.x){
        last_assignment <- get_last_assignment(taxon_name) 
        last_assignment  <- gsub(pattern = "\\[|\\]|'", replacement = '',x = last_assignment)
        trimws(last_assignment, which = "both")
      })) %>% 
      as.data.frame(check.names=FALSE, StringAsFactor=FASLE)
    # Ensure the taxonomy names are unique by aggregating duplicates
    df <- aggregate(.~taxonomy,data = df, FUN = sum)
    return(df)
  }
  ```

  **Function Parameter Definitions:**
  - `df` - a dataframe containing the taxonomy assignments
  - `taxonomy_column=` - name of the column in `df` containing the taxonomy assignments, default="taxonomy"

  **Returns:** a dataframe with unique last taxonomy names stored in a column named "taxonomy"

</details>

##### process_kaiju_table()
<details>
  <summary>reformat kaiju output table</summary>

  ```R
  process_kaiju_table <- function(file_path, taxon_col="taxon_name") {
  
    abs_abun_df <-  read_delim(file = file_path,
                               delim = "\t",
                               col_names = TRUE) %>% # read input table
             select(sample, reads, taxonomy=!!sym(taxon_col)) %>%
             pivot_wider(names_from = "sample", values_from = "reads", 
                             names_sort = TRUE) %>% # convert long dataframe to wide dataframe
             mutate_taxonomy # mutate the taxonomy coxlumn such that it contains only lowest taxonomy assignment
  
    # Set the taxon names as row names, drop the taxonomy column and convert to a matrix
    rownames(abs_abun_df) <- abs_abun_df[,"taxonomy"]
    abs_abun_df <- abs_abun_df[,-(which(colnames(abs_abun_df) == "taxonomy"))]
    abs_abun_matrix <- as.matrix(abs_abun_df)
    
    return(abs_abun_matrix)
  }
  ```

  **Function Parameter Definitions:**
  - `file_path` - file path to the tab-delimited kaiju output table file
  - `taxon_col=`- name of the taxon column in the input data file, default="taxon_path"

  **Returns:** a dataframe with reformated kaiju output

</details>


##### process_kraken_table()
<details>
  <summary>merge and process multiple kraken outputs to one species table</summary>

  ```R
  process_kraken_table <- function(reports_dir) {

    reports <- read_reports(reports_dir)
    # Retrieve sample names from file names
    samples <- names(reports) %>%
                  str_split("-") %>%
                  map_chr(function(x) pluck(x, 1))
    merged_reports  <- merge_reports2(reports, col_names = samples)
    taxonReads <- merged_reports$taxonReads
    cladeReads <- merged_reports$cladeReads
    tax_data <- merged_reports[["tax_data"]]

    species_table <- tax_data %>% 
      bind_cols(cladeReads) %>%
      filter(taxRank %in% c("U","S")) %>% # select unclassified and species rows 
      select(-contains("tax")) %>%
      zero_if_na() %>% 
      filter(name != 0) %>%  # drop unknown taxonomies
      group_by(name) %>% 
      summarise(across(everything(), sum)) %>% 
      ungroup() %>% 
      as.data.frame() %>% 
      rename(species=name)

    # Set rownames as species name, drop species column
    # and convert table from dataframe to matrix
    species_names <- species_table[,"species"]
    rownames(species_table) <- species_names
    species_table <- species_table[,-(which(colnames(species_table) == "species"))]
    species_table <- as.matrix(species_table)
    
    return(species_table)
  }
  ```

  **Function Parameter Definitions:**
  - `reports_dir` - path to a directory containing kraken2 reports 

  **Returns:** a kraken species count matrix with samples and species as columns and rows, respectively.

</details>


##### count_to_rel_abundance()
<details>
  <summary>Convert species count matrix to relative abundance matrix</summary>

  ```R
  count_to_rel_abundance <- function(species_table) {

    abund_table <- species_table %>% 
                        as.data.frame %>% 
                        mutate( across(everything(), function(x) (x/sum(x, na.rm = TRUE))*100 ) )  %>% # calculate species relative abundance per sample
        select(
                where( ~all(!is.na(.)) )
              )  %>% # drop columns where none of the reads were classified or were non-microbial (NA)
              rownames_to_column("Species") 

    # Set rownames as species name and drop species column  
    rownames(abund_table) <- abund_table$Species
    abund_table <- abund_table[, -match(x = "Species", colnames(abund_table))] %>% t

    return(abund_table)
  }
  ```

  **Function Parameter Definitions:**
  - `species_table` - a species count matrix with samples and species as columns and rows, respectively.

  **Returns:** a species relative abundance matrix with samples and species as rows and column, respectively.

</details>


##### filter_rare()
<details>
  <summary>filter out rare and non_microbial taxonomy assignments</summary>

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
    clean_tab  <- clean_tab[,-1] 
    
    
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
  - `non_microbial` - a regex denoting the string used to identify a species as non-microbial or unwanted
  - `threshold=` - abundance threshold used to determine if the relative abundance is rare, value denotes a percentage between 0 and 100.

  **Returns:** a dataframe with rare and non_microbial/unwanted species removed
</details>


##### make_plot()
<details>
  <summary>create bar plot of relative abundance</summary>

  ```R
  # Make bar plot
make_plot <- function(abund_table, metadata, colors2use, publication_format,
                      samples_column="Sample_ID", prefix_to_remove="barcode"){
  
abund_table_wide <- abund_table %>% 
    as.data.frame() %>% 
    rownames_to_column(samples_column) %>% 
    inner_join(metadata) %>% 
    select(!!!colnames(metadata), everything()) %>% 
    mutate(!!samples_column := !!sym(samples_column) %>% str_remove(prefix_to_remove))
    
  
abund_table_long <- abund_table_wide  %>%
    pivot_longer(-colnames(metadata), 
                 names_to = "Species",
                 values_to = "relative_abundance")
  
p <- ggplot(abund_table_long, mapping = aes(x=!!sym(samples_column), 
                                              y=relative_abundance, fill=Species)) +
    geom_col() +
    scale_fill_manual(values = colors2use) + 
    labs(x=NULL, y="Relative Abundance (%)") + 
    publication_format

return(p)
}
  ```

  **Function Parameter Definitions:**
  - `abund_table` - a relative bundance dataframe with rows summing to 100%
  - `metadata` - a metadata dataframe with samples as row and columns describing each sample
  - `colors2use` - a vector of strings specifying a custom color palette for coloring plots
  - `publication_format` - a ggplot::theme object specifying a custom theme for plotting
  - `samples_column` - a character column specifying the column in `metadata` holding sample names, default is "Sample_ID"
  - `prefix_to_remove` - a string specifying a prefix or any character set to remove from sample names, default is "barcode"

  **Returns:** a relative abundance stacked bar plot

</details>


##### run_decontam()
<details>
  <summary>Feature table decontamination with decontam</summary>

  ```R
  run_decontam <- function(feature_table, metadata, contam_threshold=0.1, prev_col=NULL, freq_col=NULL) {

    sub_metadata <- metadata[colnames(feature_table),]
    # Modify NTC concentration
    # Often times the user may set the NTC concentration to zero because they think nothing 
    # should be in the negative control but decontam fails if the value is set to zero.
    # To prevent decontam from failing, we replace zero with a very small concentration value
    # 0.0000001
    if (!is.null(freq_col)) {

      sub_metadata <- sub_metadata %>% 
        mutate(!!freq_col:=map_dbl(!!sym(freq_col), .f= function(conc) { 
                                      if(conc == 0) return(0.0000001) else return(conc) 
                                    } 
                                  )
              )
      sub_metadata[, freq_col] <- as.numeric(sub_metadata[,freq_col])

    }

    # Create phyloseq object
    ps <- phyloseq(otu_table(feature_table, taxa_are_rows = TRUE), sample_data(sub_metadata))

    # In our phyloseq object, `prev_col` is the sample variable that holds the negative 
    # control information. We'll summarize the data as a logical variable, with TRUE for control 
    # samples, as that is the form required by isContaminant.
    # The line below assumes that control samples will always be named "Control_Sample"
    # in the `prev_col`.
    sample_data(ps)$is.neg <- sample_data(ps)[[prev_col]] == "Control_Sample"

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

  **Returns:** a dataframe of detailed decontam results
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

```
**Function Parameter Definitions:**

- `taxonomy` - is a taxonomy assignment dataframe with ranks [Phylum, Class .. Species] as columns and taxonomy assignments as rows
- `prefix`  - is a regular expression specifying a character sequence to remove
              from taxon names

**Returns:** a dataframe of reformated taxonomy names

</details>


##### format_taxonomy_table()
<details>
  <summary>format a taxonomy assignment table by appending a suffix to a known name</summary>

```R
format_taxonomy_table <- function(taxonomy,stringToReplace="Other",
                                  suffix=";Other") {
  
  for (taxa_index in seq_along(taxonomy)) {
    
    # Get the row indices of the current taxonomy columns 
    # with rows matching the sting in `stringToReplace`
    indices <- grep(x = taxonomy[,taxa_index], pattern = stringToReplace)
    # Replace the value in that row with the value in the adjacent cell concated with `suffix` 
    taxonomy[indices,taxa_index] <- 
      paste0(taxonomy[indices,taxa_index-1],
             rep(x = suffix, times=length(indices)))
    
  }
  return(taxonomy)
}

```
**Function Parameter Definitions:**
- `taxonomy` -  taxonomy dataframe with taxonomy ranks as column names
- `stringToReplace` - a regex string specifying what to replace
- `suffix` - string specifying the replacement value

**Returns:** a dataframe of reformated taxonomy names

</details>


##### fix_names()
<details>
  <summary>clean taxonomy names</summary>

```R
fix_names<- function(taxonomy,stringToReplace,suffix){
  
  for(index in seq_along(stringToReplace)){
    taxonomy <- format_taxonomy_table(taxonomy = taxonomy,
                                      stringToReplace=stringToReplace[index], 
                                      suffix=suffix[index])
  }
  return(taxonomy)
}

```
**Function Parameter Definitions:**
- `taxonomy` -  taxonomy dataframe with taxonomy ranks as column names
- `stringToReplace` - a regex string specifying what to replace
- `suffix` - string specifying the replacement value

**Returns:** a dataframe of reformated/cleaned taxonomy names

</details>


##### read_input_table()
<details>
  <summary>read an input table into a tibble</summary>

```R
read_input_table <- function(file_name){
  
   df <- read_delim(file = file_name, delim = "\t", comment = "#")
   return(df)
   
}
```
**Function Parameter Definitions:**

- `file_name` - path to file to be read
**Returns:** a tibble generated from the input file

</details>



##### read_contig_table()
<details>
  <summary>Read Assembly-based contig annotation table</summary>

  ```R
read_contig_table <- function(file_name, sample_names){
  
  df <- read_input_table(file_name)

  # Subset taxoxnomy portion (domain:species) of input table
  # and replace empty/Na domain assignments with "Unclassified"
  taxonomy_table <- df %>%
    select(domain:species) %>%
    mutate(domain=replace_na(domain, "Unclassified"))
  
  # Subset count table
  counts_table <- df %>% select(!!sample_names)

  # Mutate taxonomy mames
  taxonomy_table  <- process_taxonomy(taxonomy_table)
  taxonomy_table <- fix_names(taxonomy_table, "Other", ";_")

  # Column bind taxonomy dataframe with species count dataframe
  df <- bind_cols(taxonomy_table, counts_table)
  
  return(df)
}

```

**Function Parameter Definitions:**

- `file_name` - path to contig taxonomy assignment file to be read
- `sample_names` - string of samples names to keep in the final dataframe

**Returns:** a dataframe with cleaned taxonomy names and sample species count

</details>



##### get_sample_names()
<details>
  <summary>retrieve sample names for which assemblies were generated</summary>

  ```R
get_sample_names <- function (assembly_summary) {


  overview_table <-  read_input_table(assembly_summary) %>%
                       select(
                         where( ~all(!is.na(.)) )
                         ) # Drop columns were all its rows are NAs

col_names <- names(overview_table) %>% str_remove_all("-assembly")
sample_order <- col_names[-1] %>% sort()

return(sample_order)

}
```
**Function Parameter Definitions:**

- `assembly_summary` - path to assembly summary file

**Returns:** a character vector of sorted sample names

</details>


#### 8c. Set global variables

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
# Heatmap color gradient - here from white to red
colours <- colorRampPalette(c('white','red'))(255)
```

**Input Data:** 

*No input data required*

**Output Data:**

- `publication_format` (a ggplot::theme object specifying a custom theme for plotting)
- `custom_palette` (a vector of strings specifying a custom color palette for coloring plots)

<br>

---

## Read-based Processing

### 9. Taxonomic profiling using kaiju

#### 9a. Build kaiju database

```bash
# Make a directory that will hold the downloaded kaiju database
mkdir kaiju-db/ && cd kaiju-db/
# Download kaiju's reference database
kaiju-makedb -s nr_euk -t NumberOfThreads
# Cleaning up
rm nr_euk/kaiju_db_nr_euk.bwt nr_euk/kaiju_db_nr_euk.sa
```

**Parameter Definitions:**

- `-t` - number of parallel processing threads to use
- `-s nr_euk` - specifies to download NCBI's nr and additionally including fungi and microbial eukaryotes databases

**Input Data:**

*No input data required*

**Output Data:**

- kaiju-db/nr_euk/kaiju_db_nr_euk.fmi (fmi file)
- kaiju-db/nodes.dmp (nodes file)
- kaiju-db/names.dmp (names file)


#### 9b. Kaiju Taxonomic Classification

```bash
kaiju -f kaiju-db/nr_euk/kaiju_db_nr_euk.fmi -t kaiju-db/nodes.dmp \
    -z NumberOfThreads \
    -E 1e-05 \
    -i /path/to/decontaminated_reads/sample_host_removed.fastq.gz \
    -o sample_kaiju.out
```

**Parameter Definitions:**

- `-f` - specifies path to the kaiju database (.fmi) file
- `-t` - specifies path to the kaiju nodes.dmp file
- `-z` - number of parallel processing threads to use
- `-E` - specifies the minimum E-value in Greedy mode (default: 1e-05)
- `-i` - specifies path to the input file
- `-o` - specifies the name of output file

**Input Data:**

- kaiju-db/nr_euk/kaiju_db_nr_euk.fmi (fmi file, output from [Step 9a](#9a-build-kaiju-database))
- kaiju-db/nodes.dmp (nodes file, output from [Step 9a](#9a-build-kaiju-database))
- sample_host_removed.fastq.gz (gzipped decontaminated reads fastq file, output from [Step 7b](#7b-remove-host-reads))

**Output Data:**

- sample_kaiju.out (kaiju output file)

#### 9c. Compile kaiju taxonomy results

```bash
# Merge kaiju reports to one table at the species level
  kaiju2table -t nodes.dmp -n names.dmp -p -r species \
              -o merged_kaiju_table.tsv *_kaiju.out

# Convert file names to sample names
sed -i -E 's/.+\/(.+)_kaiju\.out/\1/g' merged_kaiju_table.tsv && \
sed -i -E 's/file/sample/' merged_kaiju_table.tsv
```

**Parameter Definitions:**

- `-n` - specifies path to the kaiju names.dmp file
- `-t` - specifies path to the kaiju nodes.dmp file
- `-r` - specifies taxonomic rank, must be one of: phylum, class, order, family, genus, species
- `-o` - specifies the name of krona formatted kaiju output file
- `*_kaiju.out` - positional argument specifying the path to the kaiju output file (output from [Step 9ai](#9ai-read-taxonomic-classification-using-kaiju))

**Input Data:**

- kaiju-db/nodes.dmp (nodes file, output from [Step 9a](#9a-build-kaiju-database))
- kaiju-db/names.dmp (names file, output from [Step 9a](#9a-build-kaiju-database))
- *kaiju.out (kaiju report files, output from [Step 9b](#9b-kaiju-taxonomic-classification))

**Output Data:**

- **merged_kaiju_table.tsv** (Compiled kaiju table at the species taxon level)

#### 9d. Convert kaiju output to krona format

```bash
kaiju2krona -u -n kaiju-db/names.dmp -t kaiju-db/nodes.dmp \
            -i sample_kaiju.out -o sample.krona
```

**Parameter Definitions:**

- `-u` - include count for unclassified reads in output
- `-n` - specifies path to the kaiju names.dmp file
- `-t` - specifies path to the kaiju nodes.dmp file
- `-i` - specifies path to the kaiju output file (output from [Step 9b](#9b-kaiju-taxonomic-classification))
- `-o` - specifies the name of krona formatted kaiju output file

**Input Data:**
- kaiju-db/nodes.dmp (nodes file, output from [Step 9a](#9a-build-kaiju-database))
- kaiju-db/names.dmp (names file, output from [Step 9a](#9a-build-kaiju-database))
- sample_kaiju.out (kaiju output file, output from [Step 9b](#9b-kaiju-taxonomic-classification))

**Output Data:**

- sample.krona (krona formatted kaiju output)

#### 9e. Compile kaiju krona report

```bash
# Find, list and write all .krona files to file 
find . -type f -name "*.krona" |sort -uV > krona_files.txt

FILES=($(find . -type f -name "*.krona"))
basename --multiple --suffix='.krona' ${FILES[*]} | sort -uV  > sample_names.txt

# Create ktImportText input format files
KTEXT_FILES=($(paste -d',' "krona_files.txt" "sample_names.txt"))

# Create html   
ktImportText  -o kaiju-report.html ${KTEXT_FILES[*]}
```

**Parameter Definitions:**

**find**

- `-type f` -  specifies that the type of file to find is a regular file
- `-name "*.krona"` - specifies to find files ending with the .krona suffix  

**sort**

- `-u` - specifies to perform a unique sort
- `-V` - specifies to perform a mixed type of sorting

**basename**

- `--multiple` - support multiple arguments and treat each as a file name
- `--suffix='.krona'` - remove a trailing '.krona' suffix

**paste**

- `-d','` - paste both krona and sample files together line by line delimited by comma ','

**ktImportText**

- `-o` - specifies the compiled output html file name
- `${KTEXT_FILES[*]}` - an array positional arguement with the following content: 
                     sample_1.krona,sample_1 sample_2.krona,sample_2 .. sample_n.krona,sample_n.

**Input Data:**
*.krona (all sample .krona formatted files, output from [Step 9d](#9d-convert-kaiju-output-to-krona-format)) 

                      
**Output Data:**

- **kaiju-report.html** (compiled krona html report output)


#### 9f. Create kaiju species count table

```R
library(tidyverse)
feature_table <- process_kaiju_table (file_path="merged_kaiju_table.tsv")
write_csv(x = feature_table, file = "kaiju_species_table.csv")
```

**Parameter Definitions:**

- `file_path` - path to compiled kaiju table at the species taxon level
- `x`  - feature table dataframe to write to file
- `file` - path to where to write kaiju count table per sample

**Input Data:**

- merged_kaiju_table.tsv (compiled kaiju table at the species taxon level, from [Step 9c](#10c-compile-kaiju-taxonomy-results))

**Output Data:**

- kaiju_species_table.csv (kaiju species count table in csv format)


#### 9g. Read-in tables

```R
library(tidyverse)

# Read-in metadata
metdata_file <- "/path/to/sample/metadata"
samples_column <- "Sample_ID"
metadata <- read_delim(file=metdata_file , delim = "\t") %>% as.data.frame()
row.names(metadata) <- metadata[,samples_column]

# Read-in feature table
species_table <- read_csv(file="kaiju_species_table.csv") %>%  as.data.frame()
```

**Parameter Definitions:**

- `file` - path to input tables
- `delim` - file delimiter 

**Input Data:**

- metadata_file  (path to sample-wise metadata file)
- kaiju_species_table.csv (path to kaiju species table from [Step 9f](#9f-create-kaiju-species-count-table))

**Output Data:**

- `metadata` - a dataframe of sample-wise metadata
- `species_table` - a dataframe of species count per sample
---

#### 9h. Taxonomy barplots

```R
library(tidyverse)

# Threshold to filter out potential false positive
# taxonomy assignments
filter_threshold <- 0.5
# Filter out Rare and non-microbial assignments.
# You can add as many species that you'd like to filter out
# using the following syntax "|species_name1|species_name2"
non_microbial <- "Unclassifed|unclassified|Homo sapien|cannot|uncultured|unidentified"

plot_width <- 18
plot_height <- 8

# Convert count matrix to relative abundance matrix
abund_table <- count_to_rel_abundance(species_table)

# Make plot without filtering
p <- make_plot(abund_table, metadata, custom_palette, publication_format)

ggsave(filename =  "unfiltered-kaiju_species_plot.png", plot = p,
       device = "png", width = plot_width, height = plot_height, units = "in", dpi = 300)


# Get species with relative abundance greater than `filter_threshold` in all samples
# Drop rare and non-microbial assignments
filtered_species_table  <- filter_rare(species_table, non_microbial, threshold=filter_threshold)


# Convert count matrix to relative abundance matrix
filtered_species_table <- count_to_rel_abundance(filtered_species_table)

# Write filtered table to file
table2write <- filtered_species_table %>%
                 t %>%
                as.data.frame() %>%
                rownames_to_column("Species")

write_csv(x = table2write, file = "filtered-kaiju_species_table.csv")

# Make plot after filtering
p <- make_plot(filtered_species_table , metadata, custom_palette, publication_format)

ggsave(filename = "filtered-kaiju_species_plot.png", plot = p,
         device = "png", width = plot_width, height = plot_height, units = "in", dpi = 300)
```

**Parameter Definitions:**

- `filter_threshold` - a decimal threshold from 0-1 for filtering out rare species i.e potential false epositives.
- `non_microbial` - a regex string  listing out assignmnets to drop before filtering based on the `filter_threshold` above. 

**Input Data:**

- `species_table` (a dataframe of species count per sample, output from [Step 9g](#9g-read-in-tables))
- `metadata` - (a dataframe of sample-wise metadata, output from [Step 9g](#9g-read-in-tables))

**Output Data:**

- **unfiltered-kaiju_species_plot.png** (barplot plot without filtering)
- **filtered-kaiju_species_table.csv** (filtered relative abundance table)
- **filtered-kaiju_species_plot.png** (barplot after filtering rare and non-microbial taxa)


#### 9i. Feature decontamination

> Feature (species) decontamination with decontam. Decontam is an R package that statistically identifies contaminating features in a feature table

```R
library(tidyverse)
library(decontam)
library(phyloseq)

feature_table <- read_csv("filtered-kaiju_species_table.csv") %>%
                  as.data.frame()

 rownames(feature_table) <- feature_table$Species
 feature_table <- feature_table[,-1]  %>% as.matrix()
# Set to 0.5 for a more aggressive approach where species more prevalent
# in the negative controls are considered contaminants
contam_threshold <- 0.1
# Control samples in this column should always be written as 
# "Control_Sample" and true samples as "True_Sample"
prev_col <- "Sample_or_Control"
freq_col <- "input_conc_ng"
plot_width <- 18
plot_height <- 8

contamdf <- run_decontam(feature_table, metadata, contam_threshold, prev_col, freq_col)

# Write decontam results table to file
write_csv(x = contamdf %>% rownames_to_column("Species"), file = "decontam-kaiju_results.csv")

# Get the list of contaminants identified by decontam
contaminants <- contamdf %>%
                as.data.frame %>%
                rownames_to_column("Species") %>%
                filter(contaminant == TRUE) %>% pull(Species)

# Drop contaminant features identified by decontam
decontaminated_table <- feature_table %>% 
                as.data.frame  %>% 
                rownames_to_column("Species") %>% 
                filter(str_detect(Species, 
                                  pattern = str_c(contaminants,
                                                  collapse = "|"),
                                  negate = TRUE)) %>%
                select(-Species) %>% as.matrix

# Convert count matrix to relative abundance matrix
decontaminated_species_table <- count_to_rel_abundance(decontaminated_table)

# Write decontaminated species table to file
write_csv(x = decontaminated_species_table, file = "decontaminated-kaiju_species_table.csv")

# Make plot after filtering out contaminants
p <- make_plot(decontaminated_species_table , metadata, custom_palette, publication_format)

ggsave(filename = "decontaminated-kaiju-species_plot.png", plot = p,
         device = "png", width = plot_width, height = plot_height, units = "in", dpi = 300)
```

**Input Data:**

- `filtered-kaiju_species_table.csv`(path to filtered species count per sample, output from [Step 9h](#9h-taxonomy-barplots))
- `metadata`(a dataframe of sample-wise metadata, output from [Step 9g](#9g-read-in-tables))

**Output Data:**

- **decontam-kaiju_results.csv** (decontam's result table)
- **decontaminated-kaiju_species_table.csv** (decontaminated species table)
- **decontaminated-kaiju-species_plot.png** (barplot after filtering out contaminants)

<br>

---

### 10. Taxonomic Profiling using Kraken2

#### 10a. Download kraken2 database

```bash 
## Download all microbial (including eukaryotes) - https://benlangmead.github.io/aws-indexes/k2

# Downloading and building kraken2's pluspfp database which contains that standard database + plants + protists + fungi..

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

- `O` - name of file to download the url content to
- `--timeout=3600` - specifies the network timeout in seconds
- `--tries=0` - retry downdload infinitely
- `--continue` -  continue getting a partially-downloaded file
- `*_URL` - position arguement specifying the url to download a particular resource from.


**Input Data:**

- `INSPECT_URL=` - url specifying the location of kraken2 inspect file
- `LIRARY_REPORT_URL=` -  url specifying the location of kraken2 library report file
- `MD5_URL=` -  url specifying the location of the md5 file of the kraken database
- `DB_URL=` - url specifying the location of the main kraken database archive in .tar.gz format

**Output Data:**

- kraken2-db/  (a directory containing kraken 2 database files)

#### 10b. Taxonomic Classification

```bash
kraken2 --db kraken2-db/ --gzip-compressed --threads NumberOfThreads --use-names \
        --output sample-kraken2-output.txt --report sample-kraken2-report.tsv \
        /path/to/decontaminated_reads/sample_host_removed.fastq.gz
```

**Parameter Definitions:**

- `--db` - specifies the directory holding the kraken2 database files 
- `--gzip-compressed` - specifies the input fastq files are gzip-compressed
- `--threads` - number of parallel processing threads to use
- `--use-names` - specifies adding taxa names in addition to taxids
- `--output` - specifies the name of the kraken2 read-based output file (one line per read)
- `--report` - specifies the name of the kraken2 report output file (one line per taxa, with number of reads assigned to it)
- `sample_host_removed.fastq.gz` - positional argument specifying the input read file

**Input Data:**

- kraken2-db/ (a directory containing kraken 2 database files, output from [Step 10a](#10a-download-kraken2-database))
- sample_host_removed.fastq.gz (gzipped reads fastq file, output from [Step 7b](#7b-remove-host-reads))

**Output Data:**

- sample-kraken2-output.txt (kraken2 read-based output file (one line per read))
- sample-kraken2-report.tsv (kraken2 report output file (one line per taxa, with number of reads assigned to it))

#### 10c. Convert Kraken2 output to Krona format

```bash
kreport2krona.py --report-file sample-kraken2-report.tsv  --output sample.krona
```

**Parameter Definitions:**

- `--output` - specifies the name of the krona output file
- `--report-file` - specifies the name of the input kraken2 report file

**Input Data:**

- sample-kraken2-report.tsv (kraken report, output from [Step 10b](#10b-taxonomic-classification))

**Output Data:**

- sample.krona (krona formatted kraken2 output)


#### 10d. Compile kraken2 krona report

```bash
# Find, list and write all .krona files to file 
find . -type f -name "*.krona" |sort -uV > krona_files.txt

FILES=($(find . -type f -name "*.krona"))
basename --multiple --suffix='.krona' ${FILES[*]} | sort -uV  > sample_names.txt

# Create ktImportText input format files
KTEXT_FILES=($(paste -d',' "krona_files.txt" "sample_names.txt"))

# Create html   
ktImportText  -o kraken-report.html ${KTEXT_FILES[*]}
```

**Parameter Definitions:**

**find**

- `-type f` -  specifies that the type of file to find is a regular file
- `-name "*.krona"` - specifies to find files ending with the .krona suffix  

**sort**

- `-u` - specifies to perform a unique sort
- `-V` - specifies to perform a mixed type of sorting

**basename**

- `--multiple` - support multiple arguments and treat each as a file name
- `--suffix='.krona'` - remove a trailing '.krona' suffix

**paste**

- `-d','` - paste both krona and sample files together line by line delimited by comma ','

**ktImportText**

- `-o` - specifies the compiled output html file name
- `${KTEXT_FILES[*]}` - an array positional arguement with the following content: 
                     sample_1.krona,sample_1 sample_2.krona,sample_2 .. sample_n.krona,sample_n.

**Input Data:**

- *.krona (all sample .krona formatted files, output from [Step 10c](#10c-convert-kraken2-output-to-krona-format)) 

                      
**Output Data:**

- **kraken-report.html** (compiled krona html report output)

#### 10e. Create kraken species count table

```R
library(tidyverse)
library(pavian)

reports_dir <- "/path/to/directory/with/*-kraken2-report.tsv"
species_table <- process_kraken_table(reports_dir)
write_csv(x = species_table, 
          file = "kraken_species_table.csv")
```

**Parameter Definitions:**

- `reports_dir` - a directory containing kraken2 default reports
- `x` - table to write
- `file` - file name to write table to.

**Input Data:**

- *-kraken2-report.tsv (kraken2 report output file, from [Step 10b](#10b-taxonomic-classification))

**Output Data:**

- **kraken_species_table.csv** (kraken species count table in csv format)

#### 10f. Read-in tables

```R
library(tidyverse)

# Read-in metadata

metdata_file <- "/path/to/sample/metadata"
samples_column <- "Sample_ID"
metadata <- read_delim(file=metdata_file , delim = "\t") %>% as.data.frame()
row.names(metadata) <- metadata[,samples_column]
# Read-in feature table
species_table <- read_csv(file="kraken_species_table.csv") %>%  as.data.frame()
rownames(species_table) <- species_table$species
# Drop the species column
species_table <- species_table[,-match("species", colnames(species_table))]
```

**Parameter Definitions:**

- `file` - path to input table
- `delim` - file delimiter 

**Input Data:**

- metadata_file  (path to sample-wise metadata file)
- kraken_species_table.csv (path to kraken species taable)

**Output Data:**

- metadata (a dataframe of sample-wise metadata)
- species_table (a dataframe of species count with rows and columns as species and sample names, respectively)


#### 10g. Taxonomy barplots

```R
library(tidyverse)

# Threshold to filter out potential false positive
# taxonomy assignments
filter_threshold <- 0.5
# Filter out Rare and non-microbial assignments.
# You can add as many species that you'd like to filter out
# using the following syntax "|species_name1|species_name2"
non_microbial <- "Unclassifed|unclassified|Homo sapien"

plot_width <- 18
plot_height <- 8

# Convert count matrix to relative abundance matrix
abund_table <- count_to_rel_abundance(species_table)

# Make plot without filtering
p <- make_plot(abund_table, metadata, custom_palette, publication_format)

ggsave(filename =  "unfiltered-kraken_species_plot.png", plot = p, device = "png", 
       width = plot_width, height = plot_height, units = "in", dpi = 300)


# Get species with relative abundance greater than `filter_threshold` in all samples
# Drop rare and non-microbial assignments
filtered_species_table  <- filter_rare(species_table, non_microbial, threshold=filter_threshold)


# Convert count matrix to relative abundance matrix
filtered_species_table <- count_to_rel_abundance(filtered_species_table)

# Write filtered table to file
table2write <- filtered_species_table %>%
                 t %>%
                 as.data.frame() %>%
                rownames_to_column("Species")

write_csv(x = table2write , file = "filtered-kraken_species_table.csv")

# Make plot after filtering
p <- make_plot(filtered_species_table , metadata, custom_palette, publication_format)

ggsave(filename = "filtered-kraken_species_plot.png", plot = p,
         device = "png", width = plot_width, height = plot_height, units = "in", dpi = 300)
```

**Parameter Definitions:**

- `filter_threshold` - a decimal threshold from 0-1 to filter out rare species i.e potential false positives
- `non_microbial` - a regex string listing out assignments to drop before filtering based on the `filter_threshold` above 

**Input Data:**

- `species_table` (a dataframe of species count per sample, output from [Step 10f](#10f-read-in-tables))
- `metadata` - (a dataframe of sample-wise metadata, output from [Step 10f](#10f-read-in-tables))

**Output Data:**

- **unfiltered-kraken_species_plot.png** (barplot plot without filtering)
- **filtered-kraken_species_table.csv** (filtered relative abundance table)
- **filtered-kraken_species_plot.png** (barplot after filtering rare and non-microbial taxa)


#### 10h. Feature decontamination

Feature decontamination with decontam. Decontam is an R package that statistically identifies contaminating features in a feature table.

```R
library(tidyverse)
library(decontam)
library(phyloseq)

feature_table <- read_csv("filtered-kraken_species_table.csv") %>%
                  as.data.frame()

 rownames(feature_table) <- feature_table$Species
 feature_table <- feature_table[,-1]  %>% as.matrix()

# Set to 0.5 for a more aggressive approach where species more prevalent
# in the negative controls are considered contaminants
contam_threshold <- 0.1
# Control samples in this column should always be written as
# "Control_Sample" and true samples as "True_Sample" for the function below to
# function properly.
prev_col <- "Sample_or_Control"
freq_col <- "input_conc_ng"
plot_width <- 18
plot_height <- 8

contamdf <- run_decontam(feature_table, metadata, contam_threshold, prev_col, freq_col)

# Write decontam result table to file
write_csv(x = contamdf %>% rownames_to_column("Species"), file = "decontam-kraken_results.csv")

# Get the list of contaminants identified by decontam
contaminants <- contamdf %>%
                as.data.frame %>%
                rownames_to_column("Species") %>%
                filter(contaminant == TRUE) %>% pull(Species)

# Drop contaminant features identified by decontam
decontaminated_table <- feature_table %>% 
                as.data.frame  %>% 
                rownames_to_column("Species") %>% 
                filter(str_detect(Species, 
                                  pattern = str_c(contaminants,
                                                  collapse = "|"),
                                  negate = TRUE)) %>%
                select(-Species) %>% as.matrix

# Convert count matrix to relative abundance matrix
decontaminated_species_table <- count_to_rel_abundance(decontaminated_table)

# Write decontaminated species table to file
write_csv(x = decontaminated_species_table, file = "decontaminated-kraken_species_table.csv")

# Make plot after filtering out contaminants
p <- make_plot(decontaminated_species_table , metadata, custom_palette, publication_format)

ggsave(filename = "decontaminated-kraken-species_plot.png", plot = p,
         device = "png", width = plot_width, height = plot_height, units = "in", dpi = 300)
```

**Input Data:**

- `filtered-kraken_species_table.csv`(path to species count per sample, output from [Step 10g](#10g-taxonomy-barplots))
- `metadata`(a dataframe of sample-wise metadata, output from step[Step 10f](#10f-read-in-tables))

**Output Data:**

- **decontam-kraken_results.csv** (decontam's result table)
- **decontaminated-kraken_species_table.csv** (decontaminated species table)
- **decontaminated-kraken-species_plot.png** (barplot after filtering out contaminants)

<br>

---

## Assembly-based processing

### 11. Sample assembly

```bash
flye --meta --threads NumberOfThreads --out-dir sample/ \
     --nano-hq /path/to/decontaminated_raw_data/sample_host_removed.fastq.gz

# rename output files            
mv sample/assembly.fasta sample_assembly.fasta
mv sample/flye.log sample_flye.log
```

**Parameter Definitions:**

- `--meta` – use metagenome/uneven coverage mode
- `--threads` - number of parallel processing threads to use
- `--out-dir` - Output directory
- `--nano-hq` - specifies that input is from Oxford Nanopore high-quality reads (Guppy5+ SUP or Q20, <5% error). This skips a genome polishing step since the assembly will be polished with medaka in the next step

**Input Data**

- sample_host_removed.fastq.gz (decontaminated raw data in fastq format, output from [Step 7b](#7b-remove-host-reads))

**Output Data**

- sample_assembly.fasta (sample assembly)
- sample_flye.log (log file)

<br>

---

### 12. Polish assembly

```bash
medaka_consensus -t NumberOfThreads -i /path/to/decontaminated_raw_data/sample_host_removed.fastq.gz \
  -d /path/to/assemblies/sample_assembly.fasta -o sample/
  
mv sample/consensus.fasta sample_polished.fasta
```

**Parameter Definitions:**

- `-t` - number of parallel processing threads to use
- `-i` - specifies path to input read files used in creating the assembly
- `-d` - specifies path to the assembly fasta file
- `-o` - specifies the output directory

**Input Data:**

- /path/to/decontaminated_raw_data/sample_host_removed.fastq.gz (decontaminated raw data in fastq format, output from [Step 7b](#8b-remove-host-reads))
- /path/to/assemblies/sample_assembly.fasta (sample assembly, output from [Step 11](#11-sample-assembly))

**Output Data:**

- sample_polished.fasta (polished sample assembly)

---

### 13. Renaming contigs and summarizing assemblies

#### 13a. Renaming contig headers

```bash
bit-rename-fasta-headers -i sample_polished.fasta -w c_sample -o sample_assembly.fasta
```

**Parameter Definitions:**  

- `-i` – input fasta file
- `-w` – wanted header prefix (a number will be appended for each contig), starts with a "c" to ensure they won't start with a number which can be problematic
- `-o` – output fasta file


**Input Data:**

- sample_polished.fasta (polished assembly file from [Step 12](#12-polish-assembly))

**Output files:**

- **sample-assembly.fasta** (contig-renamed assembly file)


#### 13b. Summarizing assemblies

```bash
bit-summarize-assembly -o assembly-summaries_GLmetagenomics.tsv *-assembly.fasta
```

**Parameter Definitions:**  

- `-o` – output summary table
- `*-assembly.fasta` - multiple input assemblies provided as positional arguments

**Input Data:**

- *-assembly.fasta (contig-renamed assembly files from [Step 13a](#13a-renaming-contig-headers))

**Output files:**

- **assembly-summaries_GLmetagenomics.tsv** (table of assembly summary statistics)

<br>

---

### 14. Gene prediction
```bash
prodigal -a sample-genes.faa -d sample-genes.fasta -f gff -p meta -c -q \
         -o sample-genes.gff -i sample-assembly.fasta
```

**Parameter Definitions:**

- `-a` – specifies the output amino acid sequences file
- `-d` – specifies the output nucleotide sequences file
- `-f` – specifies the output format gene-calls file
- `-p` – specifies which mode to run the gene-caller in 
- `-c` – no incomplete genes reported 
- `-q` – run in quiet mode (don’t output process on each contig) 
- `-o` – specifies the name of the output gene-calls file 
- `-i` – specifies the input assembly

**Input Data:**

- sample-assembly.fasta (contig-renamed assembly file from [Step 13a](#13a-renaming-contig-headers))

**Output Data:**

- sample-genes.faa (gene-calls amino-acid fasta file)
- sample-genes.fasta (gene-calls nucleotide fasta file)
- **sample-genes.gff** (gene-calls in general feature format)

<br>

#### 14a. Remove line wraps in gene prediction output
```bash
bit-remove-wraps sample-genes.faa > sample-genes.faa.tmp 2> /dev/null
mv sample-genes.faa.tmp sample-genes.faa

bit-remove-wraps sample-genes.fasta > sample-genes.fasta.tmp 2> /dev/null
mv sample-genes.fasta.tmp sample-genes.fasta
```

**Input Data:**

- sample-genes.faa (gene-calls amino-acid fasta file, output from [Step 14](#14-gene-prediction))
- sample-genes.fasta (gene-calls nucleotide fasta file, output from [Step 14](#14-gene-prediction))

**Output Data:**

- **sample-genes.faa** (gene-calls amino-acid fasta file with line wraps removed)
- **sample-genes.fasta** (gene-calls nucleotide fasta file with line wraps removed)

<br>

---

### 15. Functional annotation
> **Note:**  
> The annotation process overwrites the same temporary directory by default. When running multiple 
processses at a time, it is necessary to specify a specific temporary directory with the 
`--tmp-dir` argument as shown below.


#### 15a. Downloading reference database of HMM models (only needs to be done once)

```bash
curl -LO ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz
curl -LO ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz
tar -xzvf profiles.tar.gz
gunzip ko_list.gz 
```

#### 15b. Running KEGG annotation

```bash
exec_annotation -p profiles/ -k ko_list --cpu NumberOfThreads -f detail-tsv -o sample-KO-tab.tmp \
                --tmp-dir sample-tmp-KO --report-unannotated sample-genes.faa 
```

**Parameter Definitions:**

- `-p` – specifies the directory holding the downloaded reference HMMs
- `-k` – specifies the downloaded reference KO  (Kegg Orthology) terms 
- `--cpu` – specifies the number of searches to run in parallel
- `-f` – specifies the output format
- `-o` – specifies the output file name
- `--tmp-dir` – specifies the temporary directory to write to (needed if running more than one process concurrently, see Notes above)
- `--report-unannotated` – specifies to generate an output for each entry
- `sample-genes.faa` – the input file is specified as a positional argument 


**Input Data:**

- sample-genes.faa (amino-acid fasta file, from [Step 14](#14-gene-prediction))
- profiles/ (reference directory holding the KO HMMs)
- ko_list (reference list of KOs to scan for)

**Output Data:**

- sample-KO-tab.tmp (table of KO annotations assigned to gene IDs)


#### 15c. Filtering output to retain only those passing the KO-specific score and top hits

```bash
bit-filter-KOFamScan-results -i sample-KO-tab.tmp -o sample-annotations.tsv

# removing temporary files
rm -rf sample-tmp-KO/ sample-KO-annots.tmp
```

**Parameter Definitions:**  

- `-i` – specifies the input table
- `-o` – specifies the output table

**Input Data:**

- sample-KO-tab.tmp (table of KO annotations assigned to gene IDs from [Step 15b](#15b-running-kegg-annotation))

**Output Data:**

- sample-annotations.tsv (table of KO annotations assigned to gene IDs)

<br>

---

### 16. Taxonomic classification

#### 16a. Pulling and un-packing pre-built reference db (only needs to be done once)

```bash
wget tbb.bio.uu.nl/bastiaan/CAT_prepare/CAT_prepare_20200618.tar.gz
tar -xvzf CAT_prepare_20200618.tar.gz
```

#### 16b. Running taxonomic classification

```bash
CAT contigs -c sample-assembly.fasta -d CAT_prepare_20200618/2020-06-18_database/ \
            -t CAT_prepare_20200618/2020-06-18_taxonomy/ -p sample-genes.faa \
            -o sample-tax-out.tmp -n NumberOfThreads -r 3 --top 4 --I_know_what_Im_doing --no-stars
```

**Parameter Definitions:**  

- `-c` – specifies the input assembly fasta file
- `-d` – specifies the CAT reference sequence database
- `-t` – specifies the CAT reference taxonomy database
- `-p` – specifies the input protein fasta file
- `-o` – specifies the output prefix
- `-n` – specifies the number of CPU cores to use
- `-r` – specifies the number of top protein hits to consider in assigning tax
- `--top` – specifies the number of protein alignments to store
- `--I_know_what_Im_doing` – allows us to alter the `--top` parameter
- `--no-stars` - suppress marking of suggestive taxonomic assignments

**Input Data:**

- sample-assembly.fasta (assembly file from [Step 13a](#13a-renaming-contig-headers))
- sample-genes.faa (gene-calls amino-acid fasta file from [Step 14](#14-gene-prediction))

**Output Data:**

- sample-tax-out.tmp.ORF2LCA.txt (gene-calls taxonomy file)
- sample-tax-out.tmp.contig2classification.txt (contig taxonomy file)

#### 16c. Adding taxonomy info from taxids to genes

```bash
CAT add_names -i sample-tax-out.tmp.ORF2LCA.txt -o sample-gene-tax-out.tmp \
              -t CAT_prepare_20200618/2020-06-18_taxonomy/ --only_official --exclude-scores
```

**Parameter Definitions:**  

- `-i` – specifies the input taxonomy file
- `-o` – specifies the output file 
- `-t` – specifies the CAT reference taxonomy database
- `--only_official` – specifies to add only standard taxonomic ranks
- `--exclude-scores` - specifies to exclude bit-score support scores in the lineage

**Input Data:**

- sample-tax-out.tmp.ORF2LCA.txt (gene-calls taxonomy file from [Step 16b](#16b-running-taxonomic-classification))

**Output Data:**

- sample-gene-tax-out.tmp (gene-calls taxonomy file with lineage info added)

#### 16d. Adding taxonomy info from taxids to contigs

```bash
CAT add_names -i sample-tax-out.tmp.contig2classification.txt -o sample-contig-tax-out.tmp \
              -t CAT-ref/2020-06-18_taxonomy/ --only_official --exclude-scores
```

**Parameter Definitions:**  

- `-i` – specifies the input taxonomy file
- `-o` – specifies the output file 
- `-t` – specifies the CAT reference taxonomy database
- `--only_official` – specifies to add only standard taxonomic ranks
- `--exclude-scores` - specifies to exclude bit-score support scores in the lineage

**Input Data:**

- sample-tax-out.tmp.contig2classification.txt (contig taxonomy file from [Step 16b](#16b-running-taxonomic-classification))

**Output Data:**

- sample-contig-tax-out.tmp (contig taxonomy file with lineage info added)


#### 16e. Formatting gene-level output with awk and sed

```bash
awk -F $'\t' ' BEGIN { OFS=FS } { if ( $3 == "lineage" ) { print $1,$3,$5,$6,$7,$8,$9,$10,$11 } \
    else if ( $2 == "ORF has no hit to database" || $2 ~ /^no taxid found/ ) \
    { print $1,"NA","NA","NA","NA","NA","NA","NA","NA" } else { n=split($3,lineage,";"); \
    print $1,lineage[n],$5,$6,$7,$8,$9,$10,$11 } } ' sample-gene-tax-out.tmp | \
    sed 's/no support/NA/g' | sed 's/superkingdom/domain/' | sed 's/# ORF/gene_ID/' | \
    sed 's/lineage/taxid/'  > sample-gene-tax-out.tsv
```

**Input Data:**

- sample-gene-tax-out.tmp (gene-calls taxonomy file with lineage info added from [Step 16c](#16c-adding-taxonomy-info-from-taxids-to-genes))

**Output Data:**

- sample-gene-tax-out.tsv (gene-calls taxonomy file with lineage info added reformatted)

#### 16f. Formatting contig-level output with awk and sed

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

- sample-contig-tax-out.tmp (contig taxonomy file with lineage info added from [Step 16d](#16d-adding-taxonomy-info-from-taxids-to-contigs))

**Output Data:**

- sample-contig-tax-out.tsv (contig taxonomy file with lineage info added reformatted)

<br>

---

### 17. Read-Mapping

#### 17a. Align Reads to Sample Assembly

```bash
minimap2 -a -x map-ont \
        -t NumberOfThreads \
        sample_assembly.fasta sample_host_removed.fastq.gz \
        > sample.sam  2> sample-mapping-info.txt
```

**Parameter Definitions:**

- `-t` - number of parallel processing threads to use
- `-a` – output in SAM format
- `-x map-ont` - specifies preset for mapping Nanopore reads to a reference

**Input Data**

- /path/to/assemblies/sample_assembly.fasta (Sample assembly, output from [Step 13a](#13a-renaming-contig-headers))
- /path/to/trimmed_reads/sample_host_removed.fastq.gz (Filtered and trimmed reads, output from [Step 7b](#7b-remove-host-reads))

**Output Data**

- sample.sam (Reads aligned to sample assembly)

#### 17b. Sort and Index Assembly Alignments

```bash
# Sort Sam, convert to bam and create index
samtools sort --threads NumberOfThreads -o sample_sorted.bam sample.sam > sample_sort.log 2>&1

samtools index sample_sorted.bam sample_sorted.bam.bai
```

**Parameter Definitions:**

**samtools sort**
- `--threads` - number of parallel processing threads to use
- `-o` - specifies the output file for the sorted reads
- `sample.sam` - positional argument specifying the input SAM file

**samtools index**
- `sample_sorted.bam` - positional argument specifying the input BAM file to be sorted
- `sample_sorted.bam.bai` - positional argument specifying the name of the index file

**Input Data:**

- sample.sam (Reads aligned to sample assembly, output from [Step 17a](#17a-align-reads-to-sample-assembly))

**Output Data:**

- sample_sorted.bam (sorted mapping to sample assembly)
- sample_sorted.bam.bai (index of sorted mapping to sample assembly)

<br>

---

### 18. Getting coverage information and filtering based on detection
> **Note:**  
> “Detection” is a measure of what proportion of a reference sequence recruited reads 
(see the discussion of detection [here](http://merenlab.org/2017/05/08/anvio-views/#detection)). 
Filtering based on detection is one way of helping to mitigate non-specific read-recruitment.

#### 18a. Filtering coverage levels based on detection

```bash
  # pileup.sh comes from the bbduk.sh package
pileup.sh -in sample.bam fastaorf=sample-genes.fasta outorf=sample-gene-cov-and-det.tmp \
          out=sample-contig-cov-and-det.tmp
```

**Parameter Definitions:**  

- `-in` – the input bam file
- `fastaorf=` – input gene-calls nucleotide fasta file
- `outorf=` – the output gene-coverage tsv file
- `out=` – the output contig-coverage tsv file

**Input Data:**

- sample.bam (mapping file from [Step 17b](#17b-sort-and-index-assembly-alignments))
- sample-genes.fasta (gene-calls nucleotide fasta file from [Step 14](#14-gene-prediction))


**Output Data:**

- sample-gene-cov-and-det.tmp (gene-coverage tsv file)
- sample-contig-cov-and-det.tmp (contig-coverage tsv file)


#### 18b. Filtering gene and contig coverage based on requiring 50% detection and parsing down to just gene ID and coverage

```bash
# Filtering gene coverage
grep -v "#" sample-gene-cov-and-det.tmp | awk -F $'\t' ' BEGIN { OFS=FS } { if ( $10 <= 0.5 ) $4 = 0 } \
     { print $1,$4 } ' > sample-gene-cov.tmp

cat <( printf "gene_ID\tcoverage\n" ) sample-gene-cov.tmp > sample-gene-coverages.tsv

# Filtering contig coverage
grep -v "#" sample-contig-cov-and-det.tmp | awk -F $'\t' ' BEGIN { OFS=FS } { if ( $5 <= 50 ) $2 = 0 } \
     { print $1,$2 } ' > sample-contig-cov.tmp

cat <( printf "contig_ID\tcoverage\n" ) sample-contig-cov.tmp > sample-contig-coverages.tsv

# removing intermediate files
rm sample-*.tmp
```

**Input Data:**

- sample-gene-cov-and-det.tmp (temporary gene-coverage tsv file from [Step 18a](#18a-filtering-coverage-levels-based-on-detection))
- sample-contig-cov-and-det.tmp (temporary contig-coverage tsv file from [Step 18a](#18a-filtering-coverage-levels-based-on-detection))

**Output Data:**

- sample-gene-coverages.tsv (table with gene-level coverages)
- sample-contig-coverages.tsv (table with contig-level coverages)

<br>

---

### 19. Combining gene-level coverage, taxonomy, and functional annotations into one table for each sample
> **Note:**  
> Just uses `paste`, `sed`, and `awk`, which are all standard in any Unix-like environment.  

```bash
paste <( tail -n +2 sample-gene-coverages.tsv | sort -V -k 1 ) <( tail -n +2 sample-annotations.tsv | sort -V -k 1 | cut -f 2- ) \
      <( tail -n +2 sample-gene-tax-out.tsv | sort -V -k 1 | cut -f 2- ) > sample-gene-tab.tmp

paste <( head -n 1 sample-gene-coverages.tsv ) <( head -n 1 sample-annotations.tsv | cut -f 2- ) \
      <( head -n 1 sample-gene-tax-out.tsv | cut -f 2- ) > sample-header.tmp

cat sample-header.tmp sample-gene-tab.tmp > sample-gene-coverage-annotation-and-tax.tsv

  # removing intermediate files
rm sample*tmp sample-gene-coverages.tsv sample-annotations.tsv sample-gene-tax-out.tsv
```

**Input Data:**

- sample-gene-coverages.tsv (table with gene-level coverages from [Step 18b](#18b-filtering-gene-and-contig-coverage-based-on-requiring-50-detection-and-parsing-down-to-just-gene-id-and-coverage))
- sample-annotations.tsv (table of KO annotations assigned to gene IDs from [Step 15c](#15c-filtering-output-to-retain-only-those-passing-the-ko-specific-score-and-top-hits))
- sample-gene-tax-out.tsv (gene-level taxonomic classifications from [Step 16f](#16f-formatting-contig-level-output-with-awk-and-sed))


**Output Data:**

- **sample-gene-coverage-annotation-and-tax.tsv** (table with combined gene coverage, annotation, and taxonomy info)

<br>

---

### 20. Combining contig-level coverage and taxonomy into one table for each sample
> **Note:**  
> Just uses `paste`, `sed`, and `awk`, which are all standard in any Unix-like environment.  

```bash
paste <( tail -n +2 sample-contig-coverages.tsv | sort -V -k 1 ) \
      <( tail -n +2 sample-contig-tax-out.tsv | sort -V -k 1 | cut -f 2- ) > sample-contig.tmp

paste <( head -n 1 sample-contig-coverages.tsv ) <( head -n 1 sample-contig-tax-out.tsv | cut -f 2- ) \
      > sample-contig-header.tmp
      
cat sample-contig-header.tmp sample-contig.tmp > sample-contig-coverage-and-tax.tsv

  # removing intermediate files
rm sample*tmp sample-contig-coverages.tsv sample-contig-tax-out.tsv
```

**Input Data:**

- sample-contig-coverages.tsv (table with contig-level coverages from [Step 18b](#18b-filtering-gene-and-contig-coverage-based-on-requiring-50-detection-and-parsing-down-to-just-gene-id-and-coverage))
- sample-contig-tax-out.tsv (contig-level taxonomic classifications from [Step 16f](#16f-formatting-contig-level-output-with-awk-and-sed))


**Output Data:**

- **sample-contig-coverage-and-tax.tsv** (table with combined contig coverage and taxonomy info)

<br>

---

### 21. Generating normalized, gene- and contig-level coverage summary tables of KO-annotations and taxonomy across samples

> **Note:**  
> * To combine across samples to generate these summary tables, we need the same "units". This is done for annotations 
based on the assigned KO terms, and all non-annotated functions are included together as "Not annotated". It is done for 
taxonomic classifications based on taxids (full lineages included in the table), and any not classified are included 
together as "Not classified". 
> * The values we are working with are coverage per gene (so they are number of bases recruited to the gene normalized 
by the length of the gene). These have been normalized by making the total coverage of a sample 1,000,000 and setting 
each individual gene-level coverage its proportion of that 1,000,000 total. So basically percent, but out of 1,000,000 
instead of 100 to make the numbers more friendly. 

#### 21a. Generating gene-level coverage summary tables

```bash
bit-GL-combine-KO-and-tax-tables *-gene-coverage-annotation-and-tax.tsv -o Combined
```

**Parameter Definitions:**  

- `*-gene-coverage-annotation-and-tax.tsv` - positional arguments specifying the input tsv files, can be provided as a space-delimited list of files, or with wildcards like above

- `-o` – specifies the output prefix


**Input Data:**

- *-gene-coverage-annotation-and-tax.tsv (tables with combined gene coverage, annotation, and taxonomy info generated for individual samples from [Step 19](#19-combining-gene-level-coverage-taxonomy-and-functional-annotations-into-one-table-for-each-sample))

**Output Data:**

- **Combined-gene-level-KO-function-coverages-CPM_GLmetagenomics.tsv** (table with all samples combined based on KO annotations; normalized to coverage per million genes covered)
- **Combined-gene-level-taxonomy-coverages-CPM_GLmetagenomics.tsv** (table with all samples combined based on gene-level taxonomic classifications; normalized to coverage per million genes covered)
- **Combined-gene-level-KO-function-coverages_GLmetagenomics.tsv** (table with all samples combined based on KO annotations)
- **Combined-gene-level-taxonomy-coverages_GLmetagenomics.tsv** (table with all samples combined based on gene-level taxonomic classifications)


#### 21b. Gene-level taxonomy heatmaps

```R
library(tidyverse)
library(pheatmap)

# Abundant taxa with CPM > 1000
abundance_threshold <- 1000

sample_order <- get_sample_names("assembly-summaries_GLmetagenomics.tsv")
# Read-in gene table
gene_taxonomy_table <-  read_contig_table("Combined-gene-level-taxonomy-coverages-CPM_GLmetagenomics.tsv", sample_order)

# Summarize gene table
species_gene_table <- gene_taxonomy_table %>%
  select(species, !!sample_order) %>% 
  group_by(species) %>% 
  summarise(across(everything(), sum)) 

# Convert gene dataframe table to a matrix table
gene.m <- species_gene_table %>% as.data.frame()
# Write out gene taxonomy table
write_csv(x = gene.m, file = "gene_taxonomy_table.csv")

rownames(gene.m) <- gene.m[['species']]
gene.m <- gene.m[,-match("species", colnames(gene.m))] %>% as.matrix()


#------ All gene taxonomy assignments

# Drop unclassified assignments
mat2plot <- gene.m[-match("Unclassified;_;_;_;_;_;_", rownames(gene.m)),]

png(filename = "All-genes-taxonomy-heatmap_GLmetagenomics.png", 
    width = plot_width, height = plot_height, units = "in", res=300)
pheatmap(mat = mat2plot,
         cluster_cols = FALSE, 
         cluster_rows = FALSE, 
         col = colours, 
         angle_col = 0, 
         display_numbers = TRUE,
         fontsize = 12,
         number_format = "%.0f")
dev.off()


#------ Abundant gene taxonomy assignments

taxa <- rowSums(gene.m) %>% sort()
abund_taxa <- taxa[ taxa > abundance_threshold ] %>% names
abund_gene.m <- gene.m[abund_taxa,]


# Drop unclassified assignments
mat2plot <- abund_gene.m[-match("Unclassified;_;_;_;_;_;_", rownames(abund_gene.m)),]

png(filename = "Abundant-genes-taxonomy-heatmap_GLmetagenomics.png", 
    width = plot_width, height = plot_height, units = "in", res=300)
pheatmap(mat = mat2plot,
         cluster_cols = FALSE, 
         cluster_rows = FALSE, 
         col = colours, 
         angle_col = 0, 
         display_numbers = TRUE,
         fontsize = 12,
         number_format = "%.0f")
dev.off()
```

**Input data:**
- assembly-summaries_GLmetagenomics.tsv (table of assembly summary statistics, output from [Step 13b](#13b-summarizing-assemblies))
- Combined-gene-level-taxonomy-coverages-CPM_GLmetagenomics.tsv (table with all samples combined based on gene-level taxonomic classifications, output from [Step 21a](#21a-generating-gene-level-coverage-summary-tables)) 

**Output data:**
- gene_taxonomy_table.csv (aggregated gene taxonomy table)
- **All-genes-taxonomy-heatmap_GLmetagenomics.png** (heatmap of all genes taxonomy assignments)
- **Abundant-genes-taxonomy-heatmap_GLmetagenomics.png** (heatmap of abundant genes taxonomy assignments)

#### 21c. Gene-level taxonomy decontamination

```R
library(tidyverse)
library(decontam)
library(phyloseq)

# Set to 0.5 for a more aggressive approach where species more prevalent
# in the negative controls are considered contaminants
contam_threshold <- 0.1
# Control samples in this column should always be written as 
# "Control_Sample" and true samples as "True_Sample"
prev_col <- "Sample_or_Control"
freq_col <- "input_conc_ng"
plot_width <- 18
plot_height <- 8

# Read-in metadata
metdata_file <- "/path/to/sample/metadata"
samples_column <- "Sample_ID"
metadata <- read_delim(file=metdata_file , delim = "\t") %>% as.data.frame()
row.names(metadata) <- metadata[,samples_column]

# Read-in featusre table
gene.m <- read_csv("gene_taxonomy_table.csv")
rownames(gene.m) <- gene.m[['species']]
gene.m <- gene.m[,-match("species", colnames(gene.m))] %>% as.matrix()
feature_table <- gene.m


contamdf <- run_decontam(feature_table, metadata, contam_threshold, prev_col, freq_col)

# Write decontam results table to file
write_csv(x = contamdf %>% rownames_to_column("Species"), file = "decontam-gene-taxonomy_results.csv")

# Get the list of contaminats identified by decontam
contaminants <- contamdf %>%
                as.data.frame %>%
                rownames_to_column("Species") %>%
                filter(contaminant == TRUE) %>% pull(Species)

# Drop contaminant features identified by decontam
decontaminated_table <- feature_table %>% 
                as.data.frame  %>% 
                rownames_to_column("Species") %>% 
                filter(str_detect(Species, 
                                  pattern = str_c(contaminants,
                                                  collapse = "|"),
                                  negate = TRUE)) %>%
                select(-Species) %>% as.matrix


# Write decontaminated species table to file
write_csv(x = decontaminated_table, file = "decontaminated-gene-taxonomy_table.csv")

# Get the index of species (contaminants and unclassified) to drop
non_microbial <- "Unclassified;_;_;_;_;_;_"
species_to_drop_index <- grep(x = rownames(feature_table), 
                              str_c(c(contaminants,non_microbial), 
                                    collapse = "|"))

mat2plot <- feature_table[-species_to_drop_index,]
png(filename = "decontaminated-gene-taxonomy-heatmap_GLmetagenomics.png", 
    width = plot_width, height = plot_height, units = "in", res=300)
pheatmap(mat = mat2plot,
         cluster_cols = FALSE, 
         cluster_rows = FALSE, 
         col = colours, 
         angle_col = 0, 
         display_numbers = TRUE,
         fontsize = 14, 
         number_format = "%.0f")
dev.off()

```

**Input data:**

- metadata_file  (path to sample-wise metadata file)
- gene_taxonomy_table.csv (aggregated gene taxonomy table, output from [Step 21b](#21b-gene-level-taxonomy-heatmaps))

**Output data:**

- **decontam-gene-taxonomy_results.csv** (decontam's results table)
- **decontaminated-gene-taxonomy_table.csv** (decontaminated species table)
- **decontaminated-gene-taxonomy-heatmap_GLmetagenomics.png** (heatmap after filtering out contaminants)



#### 21d. Gene-level KO functions heatmaps

```R
library(tidyverse)
library(pheatmap)

# Abundant functions with CPM > 2000
abundance_threshold <- 2000

sample_order <- get_sample_names("assembly-summaries_GLmetagenomics.tsv")
# Read-in KO functions table
functions_table <- read_input_table("Combined-gene-level-KO-function-coverages-CPM_GLmetagenomics.tsv") %>%
                    select(KO_ID, KO_function, !!sample_order)

# Subset table and then convert from datafame to matrix
functions.m <- functions_table[,sample_order] %>% as.matrix()
rownames(functions.m) <- functions_table$KO_ID
table2write <-  functions.m %>% 
                      as.data.frame() %>% rownames_to_column("KO_ID") %>%
                      filter(KO_ID != "Not annotated") # Drop unannotated / unclassified
# Write out  taxonomy table
write_csv(x = table2write  , file = "genes-KO-functions_table.csv")


#------ All KO functions assignments

# Drop unclassified assignments
mat2plot <- functions.m[-match("Not annotated", rownames(functions.m),]

png(filename = "All-genes-KO-functions-heatmap_GLmetagenomics.png", 
    width = plot_width, height = plot_height, units = "in", res=300)
pheatmap(mat = mat2plot,
         cluster_cols = FALSE, 
         cluster_rows = FALSE, 
         col = colours, 
         angle_col = 0, 
         display_numbers = TRUE,
         fontsize = 12,
         number_format = "%.0f")
dev.off()


#------ Abundant KO functions assignments

functions <- rowSums(functions.m) %>% sort()
abund_functions <- functions[ functions > abundance_threshold ] %>% names
abund_functions.m <- functions.m[abund_functions,]


# Drop unannotated assignments
mat2plot <- abund_functions.m[-match("Not annotated", rownames(abund_functions.m)),]

png(filename = "Abundant-genes-KO-functions-heatmap_GLmetagenomics.png", 
    width = plot_width, height = plot_height, units = "in", res=300)
pheatmap(mat = mat2plot,
         cluster_cols = FALSE, 
         cluster_rows = FALSE, 
         col = colours, 
         angle_col = 0, 
         display_numbers = TRUE,
         fontsize = 12,
         number_format = "%.0f")
dev.off()
```

**Parameter Definitions:**  


**Input data:**
- assembly-summaries_GLmetagenomics.tsv (table of assembly summary statistics, output from [Step 13b](#13b-summarizing-assemblies))
- Combined-gene-level-KO-function-coverages-CPM_GLmetagenomics.tsv (table with all samples combined based on KO annotations; normalized to coverage per million genes covered, output from [Step 21a](#21a-generating-gene-level-coverage-summary-tables))

**Output data:**
- genes-KO-functions_table.csv (aggregated and subsetted gene KO functions table)
- **All-genes-KO-functions-heatmap_GLmetagenomics.png** (heatmap of gene-wise KO function assignments)
- **Abundant-genes-KO-functions-heatmap_GLmetagenomics.png** (heatmap of gene-wise abundant KO function assignments)

#### 21e. Gene-level KO functions decontamination

```R
library(tidyverse)
library(decontam)
library(phyloseq)

# Set to 0.5 for a more aggressive approach where species more prevalent
# in negative controls are considered contaminants
contam_threshold <- 0.1 
# Control samples in this column should always be written as "Control_Sample" and true samples as "True_Sample"
prev_col <- "Sample_or_Control"
freq_col <- "input_conc_ng"
plot_width <- 18
plot_height <- 8

# Read-in metadata
metdata_file <- "/path/to/sample/metadata"
samples_column <- "Sample_ID"
metadata <- read_delim(file=metdata_file , delim = "\t") %>% as.data.frame()
row.names(metadata) <- metadata[,samples_column]

# Read-in feature table
functions.m <- read_csv("genes-KO-functions_table.csv")
rownames(functions.m) <- functions.m[['KO_ID']]
gene.m <- functions.m[,-match("KO_ID", colnames(functions.m))] %>% as.matrix()
feature_table <- functions.m


contamdf <- run_decontam(feature_table, metadata, contam_threshold, prev_col, freq_col)

# Write decontam results table to file
write_csv(x = contamdf %>% rownames_to_column("KO_ID"), file = "decontam-gene-KO-functions_results.csv")

# Get the list of contaminants identified by decontam
contaminants <- contamdf %>%
                as.data.frame %>%
                rownames_to_column("KO_ID") %>%
                filter(contaminant == TRUE) %>% pull(KO_ID)

# Drop contaminant features identified by decontam
decontaminated_table <- feature_table %>% 
                as.data.frame  %>% 
                rownames_to_column("KO_ID") %>% 
                filter(str_detect(Species, 
                                  pattern = str_c(contaminants,
                                                  collapse = "|"),
                                  negate = TRUE)) %>%
                select(-KO_ID) %>% as.matrix


# Write decontaminated species table to file
write_csv(x = decontaminated_table, file = "decontaminated-gene-KO-functions_table.csv")

# Get the index of species (contaminants and unclassified) to drop
unclassified <- "Not annotated"
functions_to_drop_index <- grep(x = rownames(feature_table), 
                              str_c(c(contaminants,unclassified), 
                                    collapse = "|"))

mat2plot <- feature_table[-functions_to_drop_index,]
png(filename = "decontaminated-gene-KO-functions-heatmap_GLmetagenomics.png", 
    width = plot_width, height = plot_height, units = "in", res=300)
pheatmap(mat = mat2plot,
         cluster_cols = FALSE, 
         cluster_rows = FALSE, 
         col = colours, 
         angle_col = 0, 
         display_numbers = TRUE,
         fontsize = 14, 
         number_format = "%.0f")
dev.off()

```

**Input data:**

- metadata_file  (path to sample-wise metadata file)
- gene_taxonomy_table.csv (agggregated gene taxomy table, output from [Step 21b](#21b-gene-level-taxonomy-heatmaps))

**Output data:**

- **decontam-gene-KO-functions_results.csv** (decontam's results table)
- **decontaminated-gene-KO-functions_table.csv** (decontaminated functions table)
- **decontaminated-gene-KO-functions-heatmap_GLmetagenomics.png** (heatmap after filtering out contaminants)



#### 21f. Generating contig-level coverage summary tables

```bash
bit-GL-combine-contig-tax-tables *-contig-coverage-and-tax.tsv -o Combined
```

**Parameter Definitions:**  

- `*-contig-coverage-and-tax.tsv` - positional arguments specifying the input tsv files, can be provided as a space-delimited list of files, or with wildcards like above
- `-o` – specifies the output prefix


**Input Data:**

- *-contig-coverage-annotation-and-tax.tsv (tables with combined contig coverage, annotation, and taxonomy info generated for individual samples from [Step 20](#20-combining-contig-level-coverage-and-taxonomy-into-one-table-for-each-sample))

**Output Data:**

- **Combined-contig-level-taxonomy-coverages-CPM_GLmetagenomics.tsv** (table with all samples combined based on contig-level taxonomic classifications; normalized to coverage per million genes covered)
- **Combined-contig-level-taxonomy-coverages_GLmetagenomics.tsv** (table with all samples combined based on contig-level taxonomic classifications)

<br>


#### 21g. Contig-level Heatmaps

```R
plot_width <- 20
plot_height <- 30
sample_order <- get_sample_names("assembly-summaries_GLmetagenomics.tsv")

contig_table <-  read_contig_table("Combined-contig-level-taxonomy-coverages-CPM_GLmetagenomics.tsv", sample_order)
species_contig_table <- contig_table %>% select(species, !!sample_order)

contig.m <- species_contig_table %>%
  group_by(species) %>%
  summarise(across(everything(), sum)) %>%
  filter(species != "Unclassified;_;_;_;_;_;_") %>% # Drop unclassifed
  as.data.frame()

# Write out contig taxonomy table
write_csv(x = contig.m, file = "contig_taxonomy_table.csv")

rownames(contig.m) <- contig.m[['species']]
contig.m <- contig.m[,-match("species", colnames(contig.m))] %>% as.matrix()

#------ All contig taxonomy assignments

# Drop unclassified assignments
mat2plot <- contig.m[-match("Unclassified;_;_;_;_;_;_", rownames(contig.m)),]

png(filename = "All-contig-taxonomy-heatmap_GLmetagenomics.png", 
    width = plot_width, height = plot_height, units = "in", res=300)
pheatmap(mat = mat2plot,
         cluster_cols = FALSE, 
         cluster_rows = FALSE, 
         col = colours, 
         angle_col = 0, 
         display_numbers = TRUE,
         fontsize = 12,
         number_format = "%.0f")
dev.off()


#------ Abundant contig taxonomy assignments

taxa <- rowSums(contig.m) %>% sort()
abund_taxa <- taxa[ taxa > abundance_threshold ] %>% names
abund_contig.m <- contig.m[abund_taxa,]

mat2plot <- abund_contig.m[-match("Unclassified;_;_;_;_;_;_", rownames(abund_contig.m)),]

png(filename = "Abundant-contig-taxonomy-heatmap_GLmetagenomics.png", 
    width = plot_width, height = plot_height, units = "in", res=300)
pheatmap(mat = mat2plot,
         cluster_cols = FALSE, 
         cluster_rows = FALSE, 
         col = colours, 
         angle_col = 0, 
         display_numbers = TRUE,
         fontsize = 12,
         number_format = "%.0f")
dev.off()
```


**Parameter Definitions:**  


**Input data:**

- assembly-summaries_GLmetagenomics.tsv (table of assembly summary statistics, output from [Step 13b](#13b-summarizing-assemblies))
- Combined-contig-level-taxonomy-coverages-CPM_GLmetagenomics.tsv (table with all samples combined based on contig-level taxonomic classifications; normalized to coverage per million genes covered from [Step 21f](#21f-generating-contig-level-coverage-summary-tables))

**Output data:**

- contig_taxonomy_table.csv (aggregated contig taxonomy)
- **All-contig-taxonomy-heatmap_GLmetagenomics.png** (All contig level taxonomy heatmap)
- **Abundant-contig-taxonomy-heatmap_GLmetagenomics.png** (Abundant contig level taxonomy heatmap)


#### 21h. Contig-level decontamination

```R
library(tidyverse)
library(decontam)
library(phyloseq)

# Set to 0.5 for a more aggressive approach where species more prevalent
# in the negative controls are considered contaminants
contam_threshold <- 0.1
# Control samples in this column should always be written as
# "Control_Sample" and true samples as "True_Sample"
prev_col <- "Sample_or_Control"
freq_col <- "input_conc_ng"
plot_width <- 18
plot_height <- 8

# Read-in metadata
metdata_file <- "/path/to/sample/metadata"
samples_column <- "Sample_ID"
metadata <- read_delim(file=metdata_file , delim = "\t") %>% as.data.frame()
row.names(metadata) <- metadata[,samples_column]

# Read-in feature table
contig.m <- read_csv("contig_taxonomy_table.csv")
rownames(contig.m) <- contig.m[['species']]
contig.m <- contig.m[,-match("species", colnames(contig.m))] %>% as.matrix()
feature_table <- contig.m


contamdf <- run_decontam(feature_table, metadata, contam_threshold, prev_col, freq_col)

# Write decontam results table to file
write_csv(x = contamdf %>% rownames_to_column("Species"), file = "decontam-gene-taxonomy_results.csv")

# Get a list of contaminants identified by decontam
contaminants <- contamdf %>%
                as.data.frame %>%
                rownames_to_column("Species") %>%
                filter(contaminant == TRUE) %>% pull(Species)

# Drop contaminant features identified by decontam
decontaminated_table <- feature_table %>% 
                as.data.frame  %>% 
                rownames_to_column("Species") %>% 
                filter(str_detect(Species, 
                                  pattern = str_c(contaminants,
                                                  collapse = "|"),
                                  negate = TRUE)) %>%
                select(-Species) %>% as.matrix


# Write decontaminated species table to file
write_csv(x = decontaminated_table, file = "decontaminated-gene-taxonomy_table.csv")

# Get the index of species (contaminants and unclassified) to drop
non_microbial <- "Unclassified;_;_;_;_;_;_"
species_to_drop_index <- grep(x = rownames(feature_table), 
                              str_c(c(contaminants,non_microbial), 
                                    collapse = "|"))

mat2plot <- feature_table[-species_to_drop_index,]
png(filename = "decontaminated-contig-taxonomy-heatmap_GLmetagenomics.png", 
    width = plot_width, height = plot_height, units = "in", res=300)
pheatmap(mat = mat2plot,
         cluster_cols = FALSE, 
         cluster_rows = FALSE, 
         col = colours, 
         angle_col = 0, 
         display_numbers = TRUE,
         fontsize = 14, 
         number_format = "%.0f")
dev.off()

```

**Input data:**

- metadata_file  (path to sample-wise metadata file)
- contig_taxonomy_table.csv (aggregated contig taxonomy table, output from [Step 21g](#21g-contig-level-heatmaps))

**Output data:**

- **decontam-contig-taxonomy_results.csv** (decontam's results table)
- **decontaminated-contig-taxonomy_table.csv** (decontaminated species table)
- **decontaminated-contig-taxonomy-heatmap_GLmetagenomics.png** (heatmap after filtering out contaminants)


---

### 22. **M**etagenome-**A**ssembled **G**enome (MAG) recovery

#### 22a. Binning contigs

```bash
jgi_summarize_bam_contig_depths --outputDepth sample-metabat-assembly-depth.tsv --percentIdentity 97 --minContigLength 1000 --minContigDepth 1.0  --referenceFasta sample-assembly.fasta sample.bam

metabat2  --inFile sample-assembly.fasta --outFile sample --abdFile sample-metabat-assembly-depth.tsv -t NumberOfThreads

mkdir sample-bins
mv sample*bin*.fasta sample-bins
zip -r sample-bins.zip sample-bins
```

**Parameter Definitions:**  

-  `--outputDepth` – specifies the output depth file
-  `--percentIdentity` – minimum end-to-end percent identity of a mapped read to be included
-  `--minContigLength` – minimum contig length to include
-  `--minContigDepth` – minimum contig depth to include
-  `--referenceFasta` – the assembly fasta file generated in step 5a
-  `sample.bam` – final positional arguments are the bam files generated in step 9
-  `--inFile` - the assembly fasta file generated in step 5a
-  `--outFile` - the prefix of the identified bins output files
-  `--abdFile` - the depth file generated by the previous `jgi_summarize_bam_contig_depths` command
-  `-t` - number of parallel processing threads to use


**Input Data:**

- sample-assembly.fasta (assembly fasta file created in [Step 13a](#13a-renaming-contig-headers))
- sample.bam (bam file created in [Step 17b](#17b-sort-and-index-assembly-alignments))

**Output Data:**

- **sample-metabat-assembly-depth.tsv** (tab-delimited summary of coverages)
- sample-bins/sample-bin\*.fasta (fasta files of recovered bins)
- **sample-bins.zip** (zip file containing fasta files of recovered bins)

#### 22b. Bin quality assessment
Utilizes the default `checkm` database available [here](https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz), `checkm_data_2015_01_16.tar.gz`.

```bash
checkm lineage_wf -f bins-overview_GLmetagenomics.tsv --tab_table -x fa ./ checkm-output-dir
```

**Parameter Definitions:**  

-  `lineage_wf` – specifies the workflow being utilized
-  `-f` – specifies the output summary file
-  `--tab_table` – specifies the output summary file should be a tab-delimited table
-  `-x` – specifies the extension that is on the bin fasta files that are being assessed
-  `./` – first positional argument at end specifies the directory holding the bins generated in step 14a
-  `checkm-output-dir` – second positional argument at end specifies the primary checkm output directory with detailed information

**Input Data:**

- sample-bins/sample-bin\*.fasta (bin fasta files generated in [Step 22a](#22a-binning-contigs))

**Output Data:**

- **bins-overview_GLmetagenomics.tsv** (tab-delimited file with quality estimates per bin)
- checkm-output-dir (directory holding detailed checkm outputs)

#### 22c. Filtering MAGs

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
  zip -r ${SAMPLE}-MAGs.zip ${SAMPLE}-MAGs
done
```

**Input Data:**

- bins-overview_GLmetagenomics.tsv (tab-delimited file with quality estimates per bin from [Step 22b](#22b-bin-quality-assessment))

**Output Data:**

- checkm-MAGs-overview.tsv (tab-delimited file with quality estimates per MAG)
- MAGs/\*.fasta (directory holding high-quality MAGs)
- **\*-MAGs.zip** (zip files containing directories of high-quality MAGs)


#### 22d. MAG taxonomic classification
Uses default `gtdbtk` database setup with program's `download.sh` command.

```bash
gtdbtk classify_wf --genome_dir MAGs/ -x fa --out_dir gtdbtk-output-dir  --skip_ani_screen
```

**Parameter Definitions:**  

-  `classify_wf` – specifies the workflow being utilized
-  `--genome_dir` – specifies the directory holding the MAGs generated in step 14c
-  `-x` – specifies the extension that is on the MAG fasta files that are being taxonomically classified
-  `--out_dir` – specifies the output directory
-  `--skip_ani_screen`  - specifies to skip ani_screening step to classify genomes using mash and skani

**Input Data:**

- MAGs/\*.fasta (directory holding high-quality MAGs from [Step 22c](#22c-filtering-mags))

**Output Data:**

- gtdbtk-output-dir/gtdbtk.\*.summary.tsv (files with assigned taxonomy and info)

#### 22e. Generating overview table of all MAGs

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

**Input Data:**

- assembly-summaries_GLmetagenomics.tsv (table of assembly summary statistics from [Step 13b](#13b-summarizing-assemblies))
- MAGs/\*.fasta (directory holding high-quality MAGs from [Step 22c](#23c-filtering-mags))
- checkm-MAGs-overview.tsv (tab-delimited file with quality estimates per MAG from [Step 22c](#22c-filtering-mags))
- gtdbtk-output-dir/gtdbtk.\*.summary.tsv (directory of files with assigned taxonomy and info from [Step 22d](#22d-mag-taxonomic-classification))

**Output Data:**

- **MAGs-overview_GLmetagenomics.tsv** (a tab-delimited overview of all recovered MAGs)


<br>

---

### 23. Generating MAG-level functional summary overview

#### 23a. Getting KO annotations per MAG
This utilizes the helper script [`parse-MAG-annots.py`](../Workflow_Documentation/NF_MGIllumina/workflow_code/bin/parse-MAG-annots.py) 

```bash
for file in $( ls MAGs/*.fasta )
do

    MAG_ID=$( echo ${file} | cut -f 2 -d "/" | sed 's/.fasta//' )
    sample_ID=$( echo ${MAG_ID} | sed 's/-MAG-[0-9]*$//' )

    grep "^>" ${file} | tr -d ">" > ${MAG_ID}-contigs.tmp

    python parse-MAG-annots.py -i annotations-and-taxonomy/${sample_ID}-gene-coverage-annotation-and-tax.tsv \
                               -w ${MAG_ID}-contigs.tmp -M ${MAG_ID} \
                               -o MAG-level-KO-annotations_GLmetagenomics.tsv

    rm ${MAG_ID}-contigs.tmp

done
```

**Parameter Definitions:**  

- `-i` – specifies the input sample gene-coverage-annotation-and-tax.tsv file generated in step 11
-  `-w` – specifies the appropriate temporary file holding all the contigs in the current MAG
- `-M` – specifies the current MAG unique identifier
- `-o` – specifies the output file

**Input Data:**

- \*-gene-coverage-annotation-and-tax.tsv (sample gene-coverage-annotation-and-tax.tsv file generated in [Step 19](#19-combining-gene-level-coverage-taxonomy-and-functional-annotations-into-one-table-for-each-sample))
- MAGs/\*.fasta (directory holding high-quality MAGs from [Step 22c](#22c-filtering-mags))

**Output Data:**

- **MAG-level-KO-annotations_GLmetagenomics.tsv** (tab-delimited table holding MAGs and their KO annotations)


#### 23b. Summarizing KO annotations with KEGG-Decoder

```bash
KEGG-decoder -v interactive -i MAG-level-KO-annotations_GLmetagenomics.tsv -o MAG-KEGG-Decoder-out_GLmetagenomics.tsv
```

**Parameter Definitions:**  

- `-v interactive` – specifies to create an interactive html output
- `-i` – specifies the input MAG-level-KO-annotations_GLmetagenomics.tsv file generated in [Step 23a](#23a-getting-ko-annotations-per-mag)
- `-o` – specifies the output table

**Input Data:**

- MAG-level-KO-annotations_GLmetagenomics.tsv (tab-delimited table holding MAGs and their KO annotations, generated in [Step 23a](#23a-getting-ko-annotations-per-mag))

**Output Data:**

- **MAG-KEGG-Decoder-out_GLmetagenomics.tsv** (tab-delimited table holding MAGs and their proportions of genes held known to be required for specific pathways/metabolisms)

- **MAG-KEGG-Decoder-out_GLmetagenomics.html** (interactive heatmap html file of the above output table)

<br>

---

