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
    - [6. Assemble Contaminants](#6-assemble-contaminants)
    - [7. Contaminant Removal](#7-remove-contaminants)
      - [7a. Build Contaminant Index and Map Reads](#7a-build-contaminant-index-and-map-reads)
      - [7b. Sort and Index Contaminant Reads](#7b-sort-and-index-contaminant-alignments)
      - [7c. Gather Contaminant Mapping Metrics](#7c-gather-contaminant-mapping-metrics)
      - [7d. Generate Decontaminated Read Files](#7d-generate-decontaminated-read-files)
      - [7e. Contaminant Removal QC](#7e-contaminant-removal-qc)
      - [7f. Compile Contaminant Removal QC](#7f-compile-contaminant-removal-qc)
    - [8. Host Removal](#8-host-removal)
      - [8a. Build or download host database](#8a-build-or-download-host-database)
        - [8a.i. Download from URL](#8ai-download-from-url)
        - [8a.ii. Build from custom reference](#8aii-build-from-custom-reference)
        - [8a.iii. Build from host name](#8aiii-build-from-host-name)
      - [8b. Remove Host Reads](#8b-remove-host-reads)
    - [9. R Environment Setup](#9-r-environment-setup)
      - [9a. Load libraries](#9a-load-libraries)
      - [9b. Define Custom Functions](#9b-define-custom-functions)
      - [9c. Set global variables](#9c-set-global-variables)
  - [**Read-based processing**](#read-based-processing)
    - [10. Taxonomic profiling using kaiju](#10-taxonomic-profiling-using-kaiju)
      - [10a. Build kaiju database](#10a-build-kaiju-database)
      - [10b. Kaiju Taxonomic Classification](#10b-kaiju-taxonomic-classification)
      - [10c. Compile kaiju taxonomy results](#10c-compile-kaiju-taxonomy-results)
      - [10d. Convert kaiju output to krona format](#10d-convert-kaiju-output-to-krona-format)
      - [10e. Compile kaiju krona report](#10e-compile-kaiju-krona-report)
      - [10f. Create kaiju species count table](#10f-create-kaiju-species-count-table)
      - [10g. Read-in tables](#10g-read-in-tables)
      - [10h. Taxonomy barplots](#10h-taxonomy-barplots)
    - [11. Taxonomic Profiling using Kraken2](#11-taxonomic-profiling-using-kraken2)
      - [11a. Download kraken2 database](#11a-download-kraken2-database)
      - [11b. Taxonomic Classification](#11b-taxonomic-classification)
      - [11c. Convert Kraken2 output to Krona format](#11c-convert-kraken2-output-to-krona-format)
      - [11d. Compile kraken2 krona report](#11d-compile-kraken2-krona-report)
      - [11e. Create kraken species count table](#11e-create-kraken-species-count-table)
      - [11f. Read-in tables](#11f-read-in-tables)
      - [11g. Taxonomy barplots](#11g-taxonomy-barplots)
  - [**Assembly-based processing**](#assembly-based-processing)
    - [12. Sample assembly](#12-sample-assembly)
    - [13. Polish assembly](#13-polish-assembly)
    - [14. Renaming contigs and summarizing assemblies](#14-renaming-contigs-and-summarizing-assemblies)
    - [15. Gene prediction](#15-gene-prediction)
    - [16. Functional annotation](#16-functional-annotation)
    - [17. Taxonomic classification](#17-taxonomic-classification)
    - [18. Read-mapping](#18-read-mapping)
    - [19. Getting coverage information and filtering based on detection](#19-getting-coverage-information-and-filtering-based-on-detection)
    - [20. Combining gene-level coverage, taxonomy, and functional annotations into one table for each sample](#20-combining-gene-level-coverage-taxonomy-and-functional-annotations-into-one-table-for-each-sample)
    - [21. Combining contig-level coverage and taxonomy into one table for each sample](#21-combining-contig-level-coverage-and-taxonomy-into-one-table-for-each-sample)
    - [22. Generating normalized, gene- and contig-level coverage summary tables of KO-annotations and taxonomy across samples](#22-generating-normalized-gene--and-contig-level-coverage-summary-tables-of-ko-annotations-and-taxonomy-across-samples)
    - [23. **M**etagenome-**A**ssembled **G**enome (MAG) recovery](#23-metagenome-assembled-genome-mag-recovery)
    - [24. Generating MAG-level functional summary overview](#24-generating-mag-level-functional-summary-overview)
  - [**Read-based Feature Table Decontamination**]()
    - [25. R environment setup](#25-r-environment-setup)
      - [25a. Load libraries](#25a-load-libraries)
      - [25b. Define custom functions](#25b-define-custom-functions)
      - [25c. Set global variables](#25c-set-global-variables)
      - [25d. Import kaiju taxonomy data](#25d-import-kaiju-taxonomy-data)
      - [25e. Import kraken2 taxonomy data](#25e-import-kraken2-taxonomy-data)
      - [25f. Import sample metadata](#25f-import-sample-metadata)
    - [26. Feature table decontamination](#26-read-based-processing-feature-table-decontamination)
      - [26a. Taxonomy Filtering](#26a-taxonomy-filtering)
      - [26b. Decontamination with decontam](#26b-decontamination-with-decontam)
        - [26b.i. Setup variables](#26bi-setup-variables)
        - [26b.ii. Identify prevalence of contaminant sequences](#26bii-identify-prevalence-of-contaminant-sequences)
        - [26b.iii. Decontaminated taxonomy plots](#26biii-decontaminated-taxonomy-plots)
  - [**Assembly-based Feature Table Decontamination**]()
    [27. ](#27-assembly-based-processing-decontamination)


<br>

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
|glue| 1.8.0 | [https://cran.r-project.org/web/packages/glue/index.html](https://cran.r-project.org/web/packages/glue/index.html) |
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
model="hac"
input_directory=/path/to/pod5/or/fast5/data
kit_name=SQK-RPB004

dorado basecaller "hac" ${input_directory} \
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

- basecalled.bam (basecalled nanopore data in bam format, output from [step 1](#1-basecalling))

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

- /path/to/fastq/output/ (directory containing spilt fastq files)[step 2a](#2a-split-fastq))

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

- /path/to/raw_data/sample.fastq.gz (raw reads, output from [Step 2](#2-demultiplexing))

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

- sample_filtered.fastq (raw reads, output from [Step 2](#2-demultiplexing))

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

### 6. Assemble Contaminants

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

### 7. Contaminant Removal

#### 7a. Build Contaminant Index and Map Reads

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

- /path/to/contaminant_assembly/assembly.fasta (Contaminant assembly, output from [Step 6](#6-assemble-contaminants))
- /path/to/trimmed_reads/sample_trimmed.fastq (Filtered and trimmed reads, output from [Step 5a](#5a-trim-filtered-data))

**Output Data**

- sample.sam (Reads aligned to contaminant assembly)

#### 7b. Sort and Index Contaminant Alignments
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

- sample.sam (Reads aligned to contaminant assembly, output from [Step 7a](#7a-build-contaminant-index-and-map-reads))

**Output Data:**

- sample_sorted.bam (sorted mapping to contaminant assembly)
- sample_sorted.bam.bai (index of sorted mapping to contaminant assembly)

#### 7c. Gather Contaminant Mapping Metrics

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

- sample_sorted.bam (sorted mapping to contaminant assembly, output from [Step 7b](#7b-sort-and-index-contaminant-alignments))
- sample_sorted.bam.bai (index of sorted mapping to contaminant assembly, output from [Step 7b](#7b-sort-and-index-contaminant-alignments))

**Output Data:**

- sample_flagstats.txt (SAM FLAG counts)
- sample_stats.txt (comprehensive alignment statistics)
- sample_idxstats.txt (contig alignment summary statistics)

#### 7d. Generate Decontaminated Read Files
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

- sample_sorted.bam (sorted mapping to contaminant assembly, output from [Step 7b](#7b-sort-and-index-contaminant-alignments))

**Output Data:**

- sample_blank_removed.fastq.gz (blank removed reads in fastq format)

#### 7e. Contaminant Removal QC

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

- sample_blank_removed.fastq.gz (blank removed reads, output from [Step 7d](#7d-generate-decontaminated-read-files))

**Output Data:**

- **/path/to/noblank_nanoplot_output/sample_NanoPlot-report.html** (NanoPlot html summary)
- /path/to/noblank_nanoplot_output/sample_NanoPlot_<date>_<time>.log (NanoPlot log file)
- /path/to/noblank_nanoplot_output/sample_NanoStats.txt (text file containing basic statistics)


#### 7f. Compile Contaminant Removal QC

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

- /path/to/noblank_nanoplot_output/*NanoStats.txt (NanoPlot output data, output from [Step 7d](#7d-generate-decontaminated-read-files))

**Output Data:**

- **noblank_multiqc_report/noblank_multiqc.html** (multiqc output html summary)
- **noblank_multiqc_report/noblank_multiqc_data.zip** (zip archive containing multiqc output data)

<br>

---

### 8. Host Removal

#### 8a. Build or download host database

##### 8a.i. Download from URL

```bash

  # Downloading and unpacking database from ${host_url}
  wget -O host.tar.gz --timeout=3600 --tries=0 --continue  host_url

  mkdir kraken2_host_db/ && \
  tar -zxvf -C kraken2_host_db/ && \
  rm -rf  host.tar.gz # Cleaning up
```

**Parameter Definitions:**

- `--timeout` - network timeout in seconds
- `--tries` - number of times to retry the download
- `--continue` - continue getting a partially downloaded file (if it exists)
- `host_url` - positional argument specifying the URl for the host database

**Output Data:**

- kraken2_host_db/ - Kraken2 database directory


##### 8a.ii. Build from custom reference

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


##### 8a.iii. Build from host name

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

- `host_name` - host database name (one of )

**Output Data:**

- kraken2_host_db/ - Kraken2 database directory

#### 8b. Remove host reads
```bash
kraken2 --db kraken2_host_db/ --gzip-compressed --threads NumberOfThreads --use-names \
        --output sample-kraken2-output.txt --report sample-kraken2-report.tsv \
        --unclassified-out sample_host_removed.fastq sample_blank_removed.fastq.gz
gzip sample_host_removed.fastq
```

**Parameter Definitions:**

- `--db` - specifies the directory holding the kraken2 database files created in [Step 8a](#8a-build-or-download-host-database)
- `--gzip-compressed` - specifies the input fastq files are gzip-compressed
- `--threads` - number of parallel processing threads to use
- `--use-names` - specifies adding taxa names in addition to taxon IDs
- `--output` - specifies the name of the kraken2 read-based output file (one line per read)
- `--report` - specifies the name of the kraken2 report output file (one line per taxa, with number of reads assigned to it)
- `--unclassified-out` - name of output file of reads that were not classified i.e non-host reads.
- `sample_blank_removed.fastq.gz` - positional argument specifying the input read file

**Input Data:**

- sample_blank_removed.fastq.gz (gzipped blank removed fastq file, output from [Step 7d](#7d-generate-decontaminated-read-files))

**Output Data:**

- sample-kraken2-output.txt (kraken2 read-based output file (one line per read))
- sample-kraken2-report.tsv (kraken2 report output file (one line per taxa, with number of reads assigned to it))
- **sample_host_removed.fastq.gz** (host-read removed, gzipped fastq file)

<br>

---

## Read-based Processing

### 9. R Environment Setup

> Taxonomy bar plot, heatmaps and feature decontamination with decontam is performed in R.

#### 9a. Load libraries

```R
library(decontam)
library(phyloseq)
library(tidyverse)
library(glue)
library(pheatmap)
library(pavian)
```

#### 9b. Define Custom Functions

##### get_last_assignment()
<details>
  <summary>retrieves the last taxonomy assignment from a taxonomy string</summary>

  ```R
  get_last_assignment <- function(taxonomy_string, split_by=';', remove_prefix=NULL){
    # A function to get the last taxonomy assignment from a taxonomy string 
    split_names <- strsplit(x =  taxonomy_string , split = split_by) %>% 
      unlist()
    
    level_name <- split_names[[length(split_names)]]
    
    if(level_name == "_"){
      return(taxonomy_string)
    }
    
    if(!is.null(remove_prefix)){
      level_name <- gsub(pattern = remove_prefix, replacement = '', x = level_name)
    }
    
    return(level_name)
  }
  ```
  **Function Parameter Definitions:**
  - `taxonomy_string` - a character string containing a list of taxonomy assignments
  - `split_by=` - a character string containing a regular expression used to split the `taxonomy_string`
  - `remove_prefix=` - a character string containing a regular expression to be matched and removed, default=`NULL`

  **Returns:** the last taxonomy assignment listed in the `taxonomy_string`
</details>

##### mutate_taxonomy()
<details>
  <summary>ensure that the taxonomy column is named "taxonomy" and aggregate duplicates to ensure that taxonomy names are unique</summary>

  ```R
  mutate_taxonomy <- function(df, taxonomy_column="taxonomy"){
    
    # make sure that the taxonomy column is always named taxonomy
    col_index <- which(colnames(df) == taxonomy_column)
    colnames(df)[col_index] <- 'taxonomy'
    df <- df %>% dplyr::mutate(across( where(is.numeric), \(x) tidyr::replace_na(x,0)  ) )%>% 
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

  **Returns:** a dataframe with unique taxonomy names stored in a column named "taxonomy"

</details>

##### process_kaiju_table()
<details>
  <summary>reformat kaiju output table</summary>

  ```R
process_kaiju_table <- function(file_path, taxon_col="taxon_name"{
  
 abs_abun_df <-  read_delim(file = file_path,
                            delim = "\t",
                            col_names = TRUE) %>% 
             select(sample, reads, taxonomy=!!sym(taxon_col)) %>%
             pivot_wider(names_from = "sample", values_from = "reads",
                             names_sort = TRUE) %>%
             mutate_taxonomy
  
   
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
process_kraken_table <- function(reports_dir){

reports <- read_reports(reports_dir)

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

species_names <- species_table[,"species"]
rownames(species_table) <- species_names

return(species_table)

}
```
  **Function Parameter Definitions:**
  - `reports_dir` - path to a directory containing kraken2 reports 

  **Returns:** a kraken species count dataframe with samples and species as columns and rowsx, respectively.

</details>



##### filter_rare()
<details>
  <summary>filter out rare and non_microbial taxonomy assignments</summary>

  ```R
  filter_rare <- function(species_table, non_microbial, threshold=1){
    
    clean_tab_count  <-  species_table %>% 
      filter(str_detect(Species, non_microbial, negate = TRUE))  
    
    clean_tab <- clean_tab_count %>% 
      mutate( across( where(is.numeric)  , \(x) (x/sum(x, na.rm = TRUE))*100 ) )
    
    rownames(clean_tab) <- clean_tab$Species
    clean_tab  <- clean_tab[,-1] 
    
    
    # Get species with relative abundance less than 1% in all samples
    rare_species <- map(clean_tab, .f = \(col) rownames(clean_tab)[col < threshold])
    rare <- Reduce(intersect, rare_species)
    
    rownames(clean_tab_count) <- clean_tab_count$Species
    clean_tab_count  <- clean_tab_count[,-1] 
    
    abund_table <- clean_tab_count[!(rownames(clean_tab_count) %in% rare), ]
    
    return(abund_table)
  }
  ```
  **Function Parameter Definitions:**
  - `species_table` - the dataframe to filter
  - `non_microbial` - a character vector denoting the string used to identify a species as non-microbial
  - `threshold=` - abundance threshold used to determine if the relative abundance is rare, value denotes a percentage between 0 and 100.

  **Returns:** a dataframe with rare and non_microbrial assignemnts removed
</details>


##### make_plot()
<details>
  <summary>create bar plot of relative abundance</summary>

  ```R
  # Make bar plot
  make_plot <- function(abund_table, metadata, colors2use, publication_format){
    
    abund_table_wide <- abund_table %>% 
        as.data.frame() %>% 
        rownames_to_column("Sample_ID") %>% 
        inner_join(metadata) %>% 
        select(!!!colnames(metadata), everything()) %>% 
        mutate(Sample_ID = Sample_ID %>% str_remove("barcode"))
        
    abund_table_long <- abund_table_wide  %>%
        pivot_longer(-colnames(metadata), 
                    names_to = "Species",
                    values_to = "relative_abundance")
      
    p <- ggplot(abund_table_long, mapping = aes(x=Sample_ID, y=relative_abundance, fill=Species)) +
         geom_col() +
         scale_fill_manual(values = colors2use) + 
         labs(x=NULL, y="Relative Abundance (%)") + 
         publication_format

    return(p)
  }
  ```
  **Function Parameter Definitions:**
  - `abund_table` - a dataframe containing the data to plot
  - `metadata` - a vector of strings specifying the data to include in the plot
  - `colors2use` - a vector of strings specifying a custom color palette for coloring plots
  - `publication_format` - a ggplot::theme object specifying the custom theme for plotting

  **Returns:** a ggplot bar plot

</details>


##### run_decontam()
<details>
  <summary>Feature table decoxntamination with decontam</summary>

  ```R
run_decontam <- function(feature_table, metadata, contam_threshold=0.1, prev_col=NULL, freq_col=NULL){

            # feature_table - [MATRIX] feature matrix to decontaminate with sample names as column and features as row
            # prev_col - [BOOLEAN] column in metadata to be use for prevalence based analysis
            # freq_col - [NUMERIC] column in metadata to be use for frequency based analysis 
            sub_metadata <- metadata[colnames(feature_table),]
            # Modify NTC concentration
            # Often times the user may set the NTC concentration to zero because they think nothing is 
            # and should be in the negative control but decontam fails if the value is zero.
            # to prevent decontam from failing we use a very small concentration value
            # 0.0000001
            if(!is.null(freq_col)){

            sub_metadata <- sub_metadata %>% 
              mutate(!!freq_col:=map_dbl(!!sym(freq_col), .f= function(conc) { 
                                              if(conc == 0) return(0.0000001) else return(conc) 
                                              } )
                    )
            sub_metadata[, freq_col] <- as.numeric(sub_metadata[,freq_col])

            }

            # Create phyloseq object
            ps <- phyloseq(otu_table(feature_table, taxa_are_rows = TRUE),
                        sample_data(sub_metadata))


           # Run Decontam 
            if(!is.null(freq_col) && !is.null(prev_col)){   

            # Run decontam in both prevalence and frequency modes
            contamdf <- isContaminant(ps, neg=prev_col, conc=freq_col, threshold=contam_threshold) # threshold

            }else if(!is.null(freq_col)){
             
            # Run decontam in frequency mode
            contamdf <- isContaminant(ps, conc=freq_col, threshold=contam_threshold) # threshold

            }else if(!is.null(prev_col)){

            # Run decontam in prevalence mode
            contamdf <- isContaminant(ps, conc=freq_col, threshold=contam_threshold)
            
            }else{

               cat("Both freq_col and prev_col cannot be set tdo NULL\n")
               cat("please supply either one or both column names your metadata")
               cat("for frequency and prevalence based analysis, respectively\n")
               stop()
            }

                                 
         return(contamdf)


}

```
**Function Parameter Definitions:**
  - `metadata` - a vector of strings specifying the data to include in the plot
  - `feature_table` -  feature matrix to decontaminate with sample names as column and features as row
  - `prev_col` - a boolean column in metadata to be use for prevalence based analysis
  - `freq_col` - a numeric column in metadata to be use for frequency based analysis
  - `contam_threshold` -  the probability threshold below which (strictly less than) the null-hypothesis 
                         (not a contaminant) should be rejected in favor of the alternate hypothesis (contaminant).


**Returns:** a dataframe of detailed decontam results

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
# remove white colors
custom_palette <- custom_palette[-c(21:23,
                                    grep(pattern = "white|snow|azure|gray|#FFFAFAFF|aliceblue",
                                         x = custom_palette, 
                                         ignore.case = TRUE)
                                   )
                                ]

# subplots grouping variable
facets <- c("Sample_Type","input_conc_ng", "lambda_spike")
```
**Input Data:** 

*No input data required*

**Output Data:**

- `publication_format` (a ggplot::theme object specifying the custom theme for plotting)
- `custom_palette` (a vector of strings specifying a custom color palette for coloring plots)
- `facets` (a vector of strings listing subplot grouping variables )

<br>

---

### 10. Taxonomic profiling using kaiju

#### 10a. Build kaiju database

```bash
# Make directory that will hold all the download kaiju database
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

- None

**Output Data:**

- kaiju-db/nr_euk/kaiju_db_nr_euk.fmi (fmi file)
- kaiju-db/nodes.dmp (nodes file)
- kaiju-db/names.dmp (names file)


#### 10b. Kaiju Taxonomic Classification

```bash
kaiju -f kaiju-db/nr_euk/kaiju_db_nr_euk.fmi -t kaiju-db/nodes.dmp \
    -z NumberOfThreads \
    -E 1e-05 \
    -i /path/to/decontaminated_reads/sample_host_removed.fastq.gz \
    -o sample_kaiju.out
```

**Parameter Definitions:**

- `-f` - specifies path to the Kaiju database (.fmi) file
- `-t` - specifies path to the Kaiju nodes.dmp file
- `-z` - number of parallel processing threads to use
- `-E` - specifies the minimum E-value in Greedy mode (default: 0.01)
- `-i` - specifies path to the input file
- `-o` - specifies the name of output file

**Input Data:**

- kaiju-db/nr_euk/kaiju_db_nr_euk.fmi (fmi file, output from [Step 9a](#9a-build-kaiju-database))
- kaiju-db/nodes.dmp (nodes file, output from [Step 9a](#9a-build-kaiju-database))
- sample_host_removed.fastq.gz (gzipped decontaminated reads fastq file, output from [Step 8](#8-host-removal))

**Output Data:**

- sample_kaiju.out (kaiju output file)

#### 10c. Compile kaiju taxonomy results

```bash
# Merge kaiju reports to one table at the species level
  kaiju2table -t nodes.dmp -n names.dmp -p -r species \
              -o merged_kaiju_table.tsv *_kaiju.out

# Covert the file names to sample names
sed -i -E 's/.+\/(.+)_kaiju\.out/\1/g' merged_kaiju_table.tsv && \
sed -i -E 's/file/sample/' merged_kaiju_table.tsv
```

**Parameter Definitions:**

- `-n` - specifies path to the Kaiju names.dmp file
- `-t` - specifies path to the Kaiju nodes.dmp file
- `-r` - specifies taxonomic rank, must be one of: phylum, class, order, family, genus, species
- `-o` - specifies the name of krona formatted kaiju output file
- `*_kaiju.out` - positional argument specifying the path to the Kaiju output file (output from [Step 9ai](#9ai-read-taxonomic-classification-using-kaiju))

**Input Data:**

- kaiju-db/nodes.dmp (nodes file, output from [Step 9a](#9a-build-kaiju-database))
- kaiju-db/names.dmp (names file, output from [Step 9a](#9a-build-kaiju-database))
- *kaiju.out (kaiju report files, output from [Step 9b](#9b-kaiju-taxonomic-classification))

**Output Data:**

- **merged_kaiju_table.tsv** (Compiled kaiju table at the species taxon level)

#### 10d. Convert kaiju output to krona format

```bash
kaiju2krona -u -n kaiju-db/names.dmp -t kaiju-db/nodes.dmp \
            -i sample_kaiju.out -o sample.krona
```

**Parameter Definitions:**

- `-u` - include count for unclassified reads in output
- `-n` - specifies path to the Kaiju names.dmp file
- `-t` - specifies path to the Kaiju nodes.dmp file
- `-i` - specifies path to the Kaiju output file (output from [Step 9b](#9b-kaiju-taxonomic-classification))
- `-o` - specifies the name of krona formatted kaiju output file

**Input Data:**
- kaiju-db/nodes.dmp (nodes file, output from [Step 9a](#9a-build-kaiju-database))
- kaiju-db/names.dmp (names file, output from [Step 9a](#9a-build-kaiju-database))
- sample_kaiju.out (kaiju output file, output from [Step 9ai](#9ai-read-taxonomic-classification-using-kaiju))

**Output Data:**

- sample.krona (krona formatted kaiju output)

#### 10e. Compile kaiju krona report

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
- `${KTEXT_FILES[*]}` - a array positional arguement with the follow content: 
                     sample_1.krona,sample_1 sample_2.krona,sample_2 .. sample_n.krona,sample_n.

**Input Data:**
*.krona (all sample .krona formatted files, output from [Step 9e](#9e-convert-kaiju-output-to-krona-format)) 

                      
**Output Data:**

- kaiju-report.html (compiled krona html report output)


#### 10f. Create kaiju species count table

```R
library(tidyverse)
feature_table <- process_kaiju_table (file_path="merged_kaiju_table.tsv")
write_csv(x = feature_table, file = "kaiju_species_table.csv")
```

**Parameter Definitions:**

- `file_path` - path to compiled kaiju table at the species taxon level
- `x`  - feature table dataframe to write to file
- `file` - path to where to write kaiju count table per sample.

**Input Data:**

- merged_kaiju_table.tsv (Compiled kaiju table at the species taxon level, from [Step 9c](#9c-compile-kaiju-taxonomy-results))

**Output Data:**

- kaiju_species_table.csv (kaiju species count table in csv format)


#### 10g. Read-in tables

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
- kaiju_species_table.csv (path to kaiju species taable from [step 10f](#10f-create-kaiju-species-count-table))

**Output Data:**

- `metadata` - a dataframe of sample-wise metadata
- `species_table` - a dataframe of species count per sample
---

#### 10h. Taxonomy barplots

```R

```

<br>

---

### 11. Taxonomic Profiling using Kraken2

#### 11a. Download kraken2 database

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

- `O` - name of file to download the url content to.
- `--timeout=3600` - specifies the network timeout to seconds seconds
- `--tries=0` - retry downdload infinitely.
- `--continue` -  continue getting a partially-downloaded file
- `*_URL` - position arguement specifying the url to download a particular resource from.


**Input Data:**

- `INSPECT_URL=` - url specifying the location of kraken2 inspect file
- `LIRARY_REPORT_URL=` -  url specifying the location of kraken2 library report file
- `MD5_URL=` -  url specifying the location of md5 file of kraken database
- `DB_URL=` - url specifying the location of the main kraken database archive in .tar.gz format

**Output Data:**

- kraken2-db/  (a directory containing kraken 2 database files)

#### 11b. Taxonomic Classification

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

- kraken2-db/ (a direcory containing kraken 2 database files, output from [Step 10a](#10a-download-kraken2-database))
- sample_host_removed.fastq.gz (gzipped reads fastq file, output from [Step 7d](#7d-generate-decontaminated-read-files))

**Output Data:**

- sample-kraken2-output.txt (kraken2 read-based output file (one line per read))
- sample-kraken2-report.tsv (kraken2 report output file (one line per taxa, with number of reads assigned to it))

#### 11c. Convert Kraken2 output to Krona format

```bash
kreport2krona.py --report-file sample-kraken2-report.tsv  --output sample.krona
```

**Parameter Definitions:**

- `--output` - specifies the name of the krona output file
- `--report-file` - specifies the name of the input kraken2 report file

**Input Data:**

- sample-kraken2-report.tsv (kraken report, output from [Step 10b](#10b-taxonomic-classification)

**Output Data:**

- sample.krona (krona formatted kraken2 output)


#### 11d. Compile kraken2 krona report

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
- `${KTEXT_FILES[*]}` - a array positional arguement with the follow content: 
                     sample_1.krona,sample_1 sample_2.krona,sample_2 .. sample_n.krona,sample_n.

**Input Data:**

- *.krona (all sample .krona formatted files, output from [Step 10d](#10d-convert-kraken2-output-to-krona-format)) 

                      
**Output Data:**

- kraken-report.html (compiled krona html report output)

#### 11e. Create kraken species count table

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

#### 11f. Read-in tables

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

- `file` - path to input tables
- `delim` - file delimiter 

**Input Data:**

- metadata_file  (path to sample-wise metadata file)
- kraken_species_table.csv (path to kraken species taable)

**Output Data:**

- `metadata` - a dataframe of sample-wise metadata
- `species_table` - a dataframe


#### 11g. Taxonomy barplots

```R

```

<br>

---

## Assembly-based processing

### 12. Sample assembly

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

- sample_host_removed.fastq.gz (decontaminated raw data in fastq format, output from [Step 8](#8-host-removal))

**Output Data**

- sample_assembly.fasta (sample assembly)
- sample_flye.log (log file)

<br>

---

### 13. Polish assembly

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

- /path/to/decontaminated_raw_data/sample_host_removed.fastq.gz (decontaminated raw data in fastq format, output from [Step 8](#8-host-removal))
- /path/to/assemblies/sample_assembly.fasta (sample assembly, output from [Step 11](#11-sample-assembly))

**Output Data:**

- sample_polished.fasta (polished sample assembly)

---

### 14. Renaming contigs and summarizing assemblies

#### 14a. Renaming contig headers

```bash
bit-rename-fasta-headers -i sample-1_polished.fasta -w c_sample-1 -o sample-1_assembly.fasta
```

**Parameter Definitions:**  

- `-i` – input fasta file
- `-w` – wanted header prefix (a number will be appended for each contig), starts with a “c_” to ensure they won’t start with a number which can be problematic
- `-o` – output fasta file


**Input Data:**

- sample-1_polished.fasta (polished assembly file from [step 12](#12-polish-assembly))

**Output files:**

- **sample-1-assembly.fasta** (contig-renamed assembly file)


#### 14b. Summarizing assemblies

```bash
bit-summarize-assembly -o assembly-summaries_GLmetagenomics.tsv *-assembly.fasta
```

**Parameter Definitions:**  

- `-o` – output summary table

- `*-assembly.fasta` - multiple input assemblies provided as positional arguments

**Input Data:**

- *-assembly.fasta (contig-renamed assembly files from [step 13a](#13a-renaming-contig-headers))

**Output files:**

- **assembly-summaries_GLmetagenomics.tsv** (table of assembly summary statistics)

<br>

---

### 15. Gene prediction
```bash
prodigal -a sample-1-genes.faa -d sample-1-genes.fasta -f gff -p meta -c -q \
         -o sample-1-genes.gff -i sample-1-assembly.fasta
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

- sample-1-assembly.fasta (contig-renamed assembly file from [step 5a](#5a-renaming-contig-headers))

**Output Data:**

- sample-1-genes.faa (gene-calls amino-acid fasta file)
- sample-1-genes.fasta (gene-calls nucleotide fasta file)
- **sample-1-genes.gff** (gene-calls in general feature format)

<br>

#### 15a. Remove line wraps in gene prediction output
```bash
bit-remove-wraps sample-1-genes.faa > sample-1-genes.faa.tmp 2> /dev/null
mv sample-1-genes.faa.tmp sample-1-genes.faa

bit-remove-wraps sample-1-genes.fasta > sample-1-genes.fasta.tmp 2> /dev/null
mv sample-1-genes.fasta.tmp sample-1-genes.fasta
```

**Input Data:**

- sample-1-genes.faa (gene-calls amino-acid fasta file, output from [Step 14](#14-gene-prediction))
- sample-1-genes.fasta (gene-calls nucleotide fasta file, output from [Step 14](#14-gene-prediction))

**Output Data:**

- **sample-1-genes.faa** (gene-calls amino-acid fasta file with line wraps removed)
- **sample-1-genes.fasta** (gene-calls nucleotide fasta file with line wraps removed)

<br>

---

### 16. Functional annotation
> **Notes**  
> The annotation process overwrites the same temporary directory by default. So if running multiple processses at a time, it is necessary to specify a specific temporary directory with the `--tmp-dir` argument as shown below.


#### 16a. Downloading reference database of HMM models (only needs to be done once)

```bash
curl -LO ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz
curl -LO ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz
tar -xzvf profiles.tar.gz
gunzip ko_list.gz 
```

#### 16b. Running KEGG annotation

```bash
exec_annotation -p profiles/ -k ko_list --cpu NumberOfThreads -f detail-tsv -o sample-1-KO-tab.tmp \
                --tmp-dir sample-1-tmp-KO --report-unannotated sample-1-genes.faa 
```

**Parameter Definitions:**

- `-p` – specifies the directory holding the downloaded reference HMMs
- `-k` – specifies the downloaded reference KO  (Kegg Orthology) terms 
- `--cpu` – specifies the number of searches to run in parallel
- `-f` – specifies the output format
- `-o` – specifies the output file name
- `--tmp-dir` – specifies the temporary directory to write to (needed if running more than one process concurrently, see Notes above)
- `--report-unannotated` – specifies to generate an output for each entry
- `sample-1-genes.faa` – the input file is specified as a positional argument 


**Input Data:**

- sample-1-genes.faa (amino-acid fasta file, from [step 6](#6-gene-prediction))
- profiles/ (reference directory holding the KO HMMs)
- ko_list (reference list of KOs to scan for)

**Output Data:**

- sample-1-KO-tab.tmp (table of KO annotations assigned to gene IDs)


#### 16c. Filtering output to retain only those passing the KO-specific score and top hits

```bash
bit-filter-KOFamScan-results -i sample-1-KO-tab.tmp -o sample-1-annotations.tsv

# removing temporary files
rm -rf sample-1-tmp-KO/ sample-1-KO-annots.tmp
```

**Parameter Definitions:**  

- `-i` – specifies the input table
- `-o` – specifies the output table


**Input Data:**

- sample-1-KO-tab.tmp (table of KO annotations assigned to gene IDs from [step 7b](#7b-running-kegg-annotation))

**Output Data:**

- sample-1-annotations.tsv (table of KO annotations assigned to gene IDs)

<br>

---

### 17. Taxonomic classification

#### 17a. Pulling and un-packing pre-built reference db (only needs to be done once)

```bash
wget tbb.bio.uu.nl/bastiaan/CAT_prepare/CAT_prepare_20200618.tar.gz
tar -xvzf CAT_prepare_20200618.tar.gz
```

#### 17b. Running taxonomic classification

```bash
CAT contigs -c sample-1-assembly.fasta -d CAT_prepare_20200618/2020-06-18_database/ \
            -t CAT_prepare_20200618/2020-06-18_taxonomy/ -p sample-1-genes.faa \
            -o sample-1-tax-out.tmp -n NumberOfThreads -r 3 --top 4 --I_know_what_Im_doing --no-stars
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

- sample-1-assembly.fasta (assembly file from [step 5a](#5a-renaming-contig-headers))
- sample-1-genes.faa (gene-calls amino-acid fasta file from [step 6](#6-gene-prediction))

**Output Data:**

- sample-1-tax-out.tmp.ORF2LCA.txt (gene-calls taxonomy file)
- sample-1-tax-out.tmp.contig2classification.txt (contig taxonomy file)

#### 17c. Adding taxonomy info from taxids to genes

```bash
CAT add_names -i sample-1-tax-out.tmp.ORF2LCA.txt -o sample-1-gene-tax-out.tmp \
              -t CAT_prepare_20200618/2020-06-18_taxonomy/ --only_official --exclude-scores
```

**Parameter Definitions:**  

- `-i` – specifies the input taxonomy file
- `-o` – specifies the output file 
- `-t` – specifies the CAT reference taxonomy database
- `--only_official` – specifies to add only standard taxonomic ranks
- `--exclude-scores` - specifies to exclude bit-score support scores in the lineage

**Input Data:**

- sample-1-tax-out.tmp.ORF2LCA.txt (gene-calls taxonomy file from [step 8b](#8b-running-taxonomic-classification))

**Output Data:**

- sample-1-gene-tax-out.tmp (gene-calls taxonomy file with lineage info added)



#### 17d. Adding taxonomy info from taxids to contigs

```bash
CAT add_names -i sample-1-tax-out.tmp.contig2classification.txt -o sample-1-contig-tax-out.tmp \
              -t CAT-ref/2020-06-18_taxonomy/ --only_official --exclude-scores
```

**Parameter Definitions:**  

- `-i` – specifies the input taxonomy file
- `-o` – specifies the output file 
- `-t` – specifies the CAT reference taxonomy database
- `--only_official` – specifies to add only standard taxonomic ranks
- `--exclude-scores` - specifies to exclude bit-score support scores in the lineage


**Input Data:**

- sample-1-tax-out.tmp.contig2classification.txt (contig taxonomy file from [step 8b](#8b-running-taxonomic-classification))

**Output Data:**

- sample-1-contig-tax-out.tmp (contig taxonomy file with lineage info added)


#### 17e. Formatting gene-level output with awk and sed

```bash
awk -F $'\t' ' BEGIN { OFS=FS } { if ( $3 == "lineage" ) { print $1,$3,$5,$6,$7,$8,$9,$10,$11 } \
    else if ( $2 == "ORF has no hit to database" || $2 ~ /^no taxid found/ ) \
    { print $1,"NA","NA","NA","NA","NA","NA","NA","NA" } else { n=split($3,lineage,";"); \
    print $1,lineage[n],$5,$6,$7,$8,$9,$10,$11 } } ' sample-1-gene-tax-out.tmp | \
    sed 's/no support/NA/g' | sed 's/superkingdom/domain/' | sed 's/# ORF/gene_ID/' | \
    sed 's/lineage/taxid/'  > sample-1-gene-tax-out.tsv
```

#### 17f. Formatting contig-level output with awk and sed

```bash
awk -F $'\t' ' BEGIN { OFS=FS } { if ( $2 == "classification" ) { print $1,$4,$6,$7,$8,$9,$10,$11,$12 } \
    else if ( $2 == "no taxid assigned" ) { print $1,"NA","NA","NA","NA","NA","NA","NA","NA" } \
    else { n=split($4,lineage,";"); print $1,lineage[n],$6,$7,$8,$9,$10,$11,$12 } } ' sample-1-contig-tax-out.tmp | \
    sed 's/no support/NA/g' | sed 's/superkingdom/domain/' | sed 's/^# contig/contig_ID/' | \
    sed 's/lineage/taxid/' > sample-1-contig-tax-out.tsv

  # clearing intermediate files
rm sample-1*.tmp*
```

**Input Data:**

- sample-1-gene-tax-out.tmp (gene-calls taxonomy file with lineage info added from [step 8c](#8c-adding-taxonomy-info-from-taxids-to-genes))
- sample-1-contig-tax-out.tmp (contig taxonomy file with lineage info added from [step 8d](#8d-adding-taxonomy-info-from-taxids-to-contigs))


**Output Data:**

- sample-1-gene-tax-out.tsv (gene-calls taxonomy file with lineage info added reformatted)
- sample-1-contig-tax-out.tsv (contig taxonomy file with lineage info added reformatted)

<br>

---

### 18. Read-Mapping

#### 18a. Align Reads to Sample Assembly

```bash
minimap2 -a -x map-ont \
        -t NumberOfThreads \
        sample_assembly.fasta sample_host_removed.fastq.gz \
        > sample.sam  2> sample-mapping-info.txt | 
```

**Parameter Definitions:**

- `-t` - number of parallel processing threads to use
- `-a` – output in SAM format
- `-x map-ont` - specifies preset for mapping Nanopore reads to a reference

**Input Data**

- /path/to/assemblies/sample_assembly.fasta (Sample assembly, output from [Step 13a](#13a-renaming-contig-headers))
- /path/to/trimmed_reads/sample_host_removed.fastq.gz (Filtered and trimmed reads, output from [Step 8](#8-host-removal))

**Output Data**

- sample.sam (Reads aligned to contaminant assembly)

#### 18b. Sort and Index Assembly Alignments

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

- sample.sam (Reads aligned to sample assembly, output from [Step 13c](#13c-read-mapping))

**Output Data:**

- sample_sorted.bam (sorted mapping to sample assembly)
- sample_sorted.bam.bai (index of sorted mapping to sample assembly)

<br>

---

### 19. Getting coverage information and filtering based on detection
> **Notes**  
> “Detection” is a metric of what proportion of a reference sequence recruited reads (see [here](http://merenlab.org/2017/05/08/anvio-views/#detection)). Filtering based on detection is one way of helping to mitigate non-specific read-recruitment.

#### 19a. Filtering coverage levels based on detection

```bash
  # pileup.sh comes from the bbduk.sh package
pileup.sh -in sample-1.bam fastaorf=sample-1-genes.fasta outorf=sample-1-gene-cov-and-det.tmp \
          out=sample-1-contig-cov-and-det.tmp
```

**Parameter Definitions:**  

- `-in` – the input bam file
- `fastaorf=` – input gene-calls nucleotide fasta file
- `outorf=` – the output gene-coverage tsv file
- `out=` – the output contig-coverage tsv file


#### 19b. Filtering gene coverage based on requiring 50% detection and parsing down to just gene ID and coverage

```bash
grep -v "#" sample-1-gene-cov-and-det.tmp | awk -F $'\t' ' BEGIN { OFS=FS } { if ( $10 <= 0.5 ) $4 = 0 } \
     { print $1,$4 } ' > sample-1-gene-cov.tmp

cat <( printf "gene_ID\tcoverage\n" ) sample-1-gene-cov.tmp > sample-1-gene-coverages.tsv
```

Filtering contig coverage based on requiring 50% detection and parsing down to just contig ID and coverage:
```bash
grep -v "#" sample-1-contig-cov-and-det.tmp | awk -F $'\t' ' BEGIN { OFS=FS } { if ( $5 <= 50 ) $2 = 0 } \
     { print $1,$2 } ' > sample-1-contig-cov.tmp

cat <( printf "contig_ID\tcoverage\n" ) sample-1-contig-cov.tmp > sample-1-contig-coverages.tsv

  # removing intermediate files

rm sample-1-*.tmp
```

**Input Data:**

- sample-1.bam (mapping file from [step 9b](#9b-performing-mapping-conversion-to-bam-and-sorting))
- sample-1-genes.fasta (gene-calls nucleotide fasta file from [step 6](#6-gene-prediction))

**Output Data:**

- sample-1-gene-coverages.tsv (table with gene-level coverages)
- sample-1-contig-coverages.tsv (table with contig-level coverages)

<br>

---

### 20. Combining gene-level coverage, taxonomy, and functional annotations into one table for each sample
> **Notes**  
> Just uses `paste`, `sed`, and `awk`, all are standard in any Unix-like environment.  

```bash
paste <( tail -n +2 sample-1-gene-coverages.tsv | sort -V -k 1 ) <( tail -n +2 sample-1-annotations.tsv | sort -V -k 1 | cut -f 2- ) \
      <( tail -n +2 sample-1-gene-tax-out.tsv | sort -V -k 1 | cut -f 2- ) > sample-1-gene-tab.tmp

paste <( head -n 1 sample-1-gene-coverages.tsv ) <( head -n 1 sample-1-annotations.tsv | cut -f 2- ) \
      <( head -n 1 sample-1-gene-tax-out.tsv | cut -f 2- ) > sample-1-header.tmp

cat sample-1-header.tmp sample-1-gene-tab.tmp > sample-1-gene-coverage-annotation-and-tax.tsv

  # removing intermediate files
rm sample-1*tmp sample-1-gene-coverages.tsv sample-1-annotations.tsv sample-1-gene-tax-out.tsv
```

**Input Data:**

- sample-1-gene-coverages.tsv (table with gene-level coverages from [step 10b](#10b-filtering-gene-coverage-based-on-requiring-50-detection-and-parsing-down-to-just-gene-id-and-coverage))
- sample-1-annotations.tsv (table of KO annotations assigned to gene IDs from [step 7c](#7c-filtering-output-to-retain-only-those-passing-the-ko-specific-score-and-top-hits))
- sample-1-gene-tax-out.tsv (gene-level taxonomic classifications from [step 8f](#8f-formatting-contig-level-output-with-awk-and-sed))


**Output Data:**

- **sample-1-gene-coverage-annotation-and-tax.tsv** (table with combined gene coverage, annotation, and taxonomy info)

<br>

---

### 21. Combining contig-level coverage and taxonomy into one table for each sample
> **Notes**  
> Just uses `paste`, `sed`, and `awk`, all are standard in any Unix-like environment.  

```bash
paste <( tail -n +2 sample-1-contig-coverages.tsv | sort -V -k 1 ) \
      <( tail -n +2 sample-1-contig-tax-out.tsv | sort -V -k 1 | cut -f 2- ) > sample-1-contig.tmp

paste <( head -n 1 sample-1-contig-coverages.tsv ) <( head -n 1 sample-1-contig-tax-out.tsv | cut -f 2- ) \
      > sample-1-contig-header.tmp
      
cat sample-1-contig-header.tmp sample-1-contig.tmp > sample-1-contig-coverage-and-tax.tsv

  # removing intermediate files
rm sample-1*tmp sample-1-contig-coverages.tsv sample-1-contig-tax-out.tsv
```

**Input Data:**

- sample-1-contig-coverages.tsv (table with contig-level coverages from [step 10b](#10b-filtering-gene-coverage-based-on-requiring-50-detection-and-parsing-down-to-just-gene-id-and-coverage))
- sample-1-contig-tax-out.tsv (contig-level taxonomic classifications from [step 8f](#8f-formatting-contig-level-output-with-awk-and-sed))


**Output Data:**

- **sample-1-contig-coverage-and-tax.tsv** (table with combined contig coverage and taxonomy info)

<br>

---

### 22. Generating normalized, gene- and contig-level coverage summary tables of KO-annotations and taxonomy across samples

> **Notes**  
> * To combine across samples to generate these summary tables, we need the same "units". This is done for annotations based on the assigned KO terms, and all non-annotated functions are included together as "Not annotated". It is done for taxonomic classifications based on taxids (full lineages included in the table), and any not classified are included together as "Not classified". 
> * The values we are working with are coverage per gene (so they are number of bases recruited to the gene normalized by the length of the gene). These have been normalized by making the total coverage of a sample 1,000,000 and setting each individual gene-level coverage its proportion of that 1,000,000 total. So basically percent, but out of 1,000,000 instead of 100 to make the numbers more friendly. 

#### 22a. Generating gene-level coverage summary tables

```bash
bit-GL-combine-KO-and-tax-tables *-gene-coverage-annotation-and-tax.tsv -o Combined
```

**Parameter Definitions:**  

- `*-gene-coverage-annotation-and-tax.tsv` - positional arguments specifying the input tsv files, can be provided as a space-delimited list of files, or with wildcards like above

- `-o` – specifies the output prefix


**Input Data:**

- *-gene-coverage-annotation-and-tax.tsv (tables with combined gene coverage, annotation, and taxonomy info generated for individual samples from [step 11](#11-combining-gene-level-coverage-taxonomy-and-functional-annotations-into-one-table-for-each-sample))

**Output Data:**

- **Combined-gene-level-KO-function-coverages-CPM_GLmetagenomics.tsv** (table with all samples combined based on KO annotations; normalized to coverage per million genes covered)
- **Combined-gene-level-taxonomy-coverages-CPM_GLmetagenomics.tsv** (table with all samples combined based on gene-level taxonomic classifications; normalized to coverage per million genes covered)
- **Combined-gene-level-KO-function-coverages_GLmetagenomics.tsv** (table with all samples combined based on KO annotations)
- **Combined-gene-level-taxonomy-coverages_GLmetagenomics.tsv** (table with all samples combined based on gene-level taxonomic classifications)


#### 22b. Generating contig-level coverage summary tables

```bash
bit-GL-combine-contig-tax-tables *-contig-coverage-and-tax.tsv -o Combined
```
**Parameter Definitions:**  

- `*-contig-coverage-and-tax.tsv` - positional arguments specifying the input tsv files, can be provided as a space-delimited list of files, or with wildcards like above

- `-o` – specifies the output prefix


**Input Data:**

- *-contig-coverage-annotation-and-tax.tsv (tables with combined contig coverage, annotation, and taxonomy info generated for individual samples from [step 12](#12-combining-contig-level-coverage-and-taxonomy-into-one-table-for-each-sample))

**Output Data:**

- **Combined-contig-level-taxonomy-coverages-CPM_GLmetagenomics.tsv** (table with all samples combined based on contig-level taxonomic classifications; normalized to coverage per million genes covered)
- **Combined-contig-level-taxonomy-coverages_GLmetagenomics.tsv** (table with all samples combined based on contig-level taxonomic classifications)

<br>

---

### 23. **M**etagenome-**A**ssembled **G**enome (MAG) recovery

#### 23a. Binning contigs

```bash
jgi_summarize_bam_contig_depths --outputDepth sample-1-metabat-assembly-depth.tsv --percentIdentity 97 --minContigLength 1000 --minContigDepth 1.0  --referenceFasta sample-1-assembly.fasta sample-1.bam

metabat2  --inFile sample-1-assembly.fasta --outFile sample-1 --abdFile sample-1-metabat-assembly-depth.tsv -t NumberOfThreads

mkdir sample-1-bins
mv sample-1*bin*.fasta sample-1-bins
zip -r sample-1-bins.zip sample-1-bins
```

**Parameter Definitions:**  

-  `--outputDepth` – specifies the output depth file
-  `--percentIdentity` – minimum end-to-end percent identity of a mapped read to be included
-  `--minContigLength` – minimum contig length to include
-  `--minContigDepth` – minimum contig depth to include
-  `--referenceFasta` – the assembly fasta file generated in step 5a
-  `sample-1.bam` – final positional arguments are the bam files generated in step 9
-  `--inFile` - the assembly fasta file generated in step 5a
-  `--outFile` - the prefix of the identified bins output files
-  `--abdFile` - the depth file generated by the previous `jgi_summarize_bam_contig_depths` command
-  `-t` - number of parallel processing threads to use


**Input Data:**

- sample-1-assembly.fasta (assembly fasta file created in [step 5a](#5a-renaming-contig-headers))
- sample-1.bam (bam file created in [step 9b](#9b-performing-mapping-conversion-to-bam-and-sorting))

**Output Data:**

- **sample-1-metabat-assembly-depth.tsv** (tab-delimited summary of coverages)
- sample-1-bins/sample-1-bin\*.fasta (fasta files of recovered bins)
- **sample-1-bins.zip** (zip file containing fasta files of recovered bins)

#### 23b. Bin quality assessment
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

- sample-1-bins/sample-1-bin\*.fasta (bin fasta files generated in [step 14a](#14a-binning-contigs))

**Output Data:**

- **bins-overview_GLmetagenomics.tsv** (tab-delimited file with quality estimates per bin)
- checkm-output-dir (directory holding detailed checkm outputs)

#### 23c. Filtering MAGs

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

- bins-overview_GLmetagenomics.tsv (tab-delimited file with quality estimates per bin from [step 14b](#14b-bin-quality-assessment))

**Output Data:**

- checkm-MAGs-overview.tsv (tab-delimited file with quality estimates per MAG)
- MAGs/\*.fasta (directory holding high-quality MAGs)
- **\*-MAGs.zip** (zip files containing directories of high-quality MAGs)


#### 23d. MAG taxonomic classification
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

- MAGs/\*.fasta (directory holding high-quality MAGs from [step 14c](#14c-filtering-mags))

**Output Data:**

- gtdbtk-output-dir/gtdbtk.\*.summary.tsv (files with assigned taxonomy and info)

#### 23e. Generating overview table of all MAGs

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

- assembly-summaries_GLmetagenomics.tsv (table of assembly summary statistics from [step 5b](#5b-summarizing-assemblies))
- MAGs/\*.fasta (directory holding high-quality MAGs from [step 14c](#14c-filtering-mags))
- checkm-MAGs-overview.tsv (tab-delimited file with quality estimates per MAG from [step 14c](#14c-filtering-mags))
- gtdbtk-output-dir/gtdbtk.\*.summary.tsv (directory of files with assigned taxonomy and info from [step 14d](#14d-mag-taxonomic-classification))

**Output Data:**

- **MAGs-overview_GLmetagenomics.tsv** (a tab-delimited overview of all recovered MAGs)


<br>

---

### 24. Generating MAG-level functional summary overview

#### 24a. Getting KO annotations per MAG
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

- \*-gene-coverage-annotation-and-tax.tsv (sample gene-coverage-annotation-and-tax.tsv file generated in [step 11](#11-combining-gene-level-coverage-taxonomy-and-functional-annotations-into-one-table-for-each-sample))
- MAGs/\*.fasta (directory holding high-quality MAGs from [step 14c](#14c-filtering-mags))

**Output Data:**

- **MAG-level-KO-annotations_GLmetagenomics.tsv** (tab-delimited table holding MAGs and their KO annotations)


#### 24b. Summarizing KO annotations with KEGG-Decoder

```bash
KEGG-decoder -v interactive -i MAG-level-KO-annotations_GLmetagenomics.tsv -o MAG-KEGG-Decoder-out_GLmetagenomics.tsv
```

**Parameter Definitions:**  

- `-v interactive` – specifies to create an interactive html output
- `-i` – specifies the input MAG-level-KO-annotations_GLmetagenomics.tsv file generated in [step 15a](#15a-getting-ko-annotations-per-mag)
- `-o` – specifies the output table

**Input Data:**

- MAG-level-KO-annotations_GLmetagenomics.tsv (tab-delimited table holding MAGs and their KO annotations, generated in [step 15a](#15a-getting-ko-annotations-per-mag))

**Output Data:**

- **MAG-KEGG-Decoder-out_GLmetagenomics.tsv** (tab-delimited table holding MAGs and their proportions of genes held known to be required for specific pathways/metabolisms)

- **MAG-KEGG-Decoder-out_GLmetagenomics.html** (interactive heatmap html file of the above output table)

<br>

---

## Read-based Feature Table Decontamination
> Feature table decontamination is performed in R.  

### 25. R Environment Setup

#### 25a. Load libraries

```R
library(decontam)
library(phyloseq)
library(tidyverse)
library(DT)
library(plotly)
library(glue)
library(pheatmap)
library(pavian)
```

#### 25b. Define Custom Functions

##### get_last_assignment()
<details>
  <summary>retrieves the last taxonomy assignment from a taxonomy string</summary>

  ```R
  get_last_assignment <- function(taxonomy_string, split_by=';', remove_prefix=NULL){
    # A function to get the last taxonomy assignment from a taxonomy string 
    split_names <- strsplit(x =  taxonomy_string , split = split_by) %>% 
      unlist()
    
    level_name <- split_names[[length(split_names)]]
    
    if(level_name == "_"){
      return(taxonomy_string)
    }
    
    if(!is.null(remove_prefix)){
      level_name <- gsub(pattern = remove_prefix, replacement = '', x = level_name)
    }
    
    return(level_name)
  }
  ```
  **Function Parameter Definitions:**
  - `taxonomy_string` - a character string containing a list of taxonomy assignments
  - `split_by=` - a character string containing a regular expression used to split the `taxonomy_string`
  - `remove_prefix=` - a character string containing a regular expression to be matched and removed, default=`NULL`

  **Returns:** the last taxonomy assignment listed in the `taxonomy_string`
</details>

##### mutate_taxonomy()
<details>
  <summary>ensure that the taxonomy column is named "taxonomy" and aggregate duplicates to ensure that taxonomy names are unique</summary>

  ```R
  mutate_taxonomy <- function(df, taxonomy_column="taxonomy"){
    
    # make sure that the taxonomy column is always named taxonomy
    col_index <- which(colnames(df) == taxonomy_column)
    colnames(df)[col_index] <- 'taxonomy'
    df <- df %>% dplyr::mutate(across( where(is.numeric), \(x) tidyr::replace_na(x,0)  ) )%>% 
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

  **Returns:** a dataframe with unique taxonomy names stored in a column named "taxonomy"

</details>

##### process_kaiju_table()
<details>
  <summary>reformat kaiju output table</summary>

  ```R
  process_kaiju_table <- function(file_path, taxon_col="taxon_path",
                                  kingdom=NULL, remove_non_microbial = TRUE){
    
    kaija_table <- read_delim(file = file_path,
                              delim = "\t",
                              col_names = TRUE)
    
    if(remove_non_microbial){
      # Remove non-microbial and unclassified assignments in this case Metazoa for animal assignments
      non_microbial_indices <- grep(pattern = "unclassified|assigned|Metazoa|Chordata|Nematoda|Arthropoda|Annelida|Brachiopoda|Mollusca|Cnidaria|Streptophyta",
                                    x = kaija_table[[taxon_col]])
      
      if(!is_empty(non_microbial_indices)){
        kaija_table <- kaija_table[-non_microbial_indices,]
      }
      
    }
    
    if(!is.null(kingdom)){
      kingdom_indices <- grep(pattern = kingdom ,
                              x = kaija_table[[taxon_col]])
      if(!is_empty(kingdom_indices)){
        kaija_table <- kaija_table[kingdom_indices,]
      }
    }
    
    
    abs_abun_df <- pivot_wider(data = kaija_table %>% dplyr::select(sample,reads,taxonomy=!!sym(taxon_col)), 
                              names_from = "sample", values_from = "reads",
                              names_sort = TRUE) %>% mutate_taxonomy
    
    rel_abun_df <- pivot_wider(data = kaija_table %>% dplyr::select(sample,percent,taxonomy=!!sym(taxon_col)), 
                              names_from = "sample", values_from = "percent",
                              names_sort = TRUE) %>% mutate_taxonomy
    
    # Set the taxon names as row names, drop the taxonomy column and convert to a matrix
    rownames(abs_abun_df) <- abs_abun_df[,"taxonomy"]
    rownames(rel_abun_df) <- rel_abun_df[,"taxonomy"]
    
    abs_abun_df <- abs_abun_df[,-(which(colnames(abs_abun_df) == "taxonomy"))]
    rel_abun_df <- rel_abun_df[,-(which(colnames(rel_abun_df) == "taxonomy"))]
    
    abs_abun_matrix <- as.matrix(abs_abun_df)
    rel_abun_matrix <- as.matrix(rel_abun_df)
    
    final_tables <- list("relative_table"=rel_abun_matrix,
                        "abundance_table"=abs_abun_matrix)
    return(final_tables)
    
  }
  ```
  **Function Parameter Definitions:**
  - `file_path` - file path to the tab-delimited kaiju output table file
  - `taxon_col=`- name of the taxon column in the input data file, default="taxon_path"
  - `kingdom=` - a character string containing a regular expression used to filter for specific kingdoms, default=`NULL`
  - `remove_non_microbial=` - a boolean specifying whether or not to remove non-microbial and unclassified assuments, default=`TRUE`

  **Returns:** a dataframe with reformated kaiju output

</details>

##### create_dt()
<details>
  <summary>create an HTML widget to display rectangular data (`matrix` or `dataframe`) using the DataTables Javascript library</summary>

```R
create_dt <- function(table2show, caption=NULL) {
  DT::datatable(table2show,
                rownames = FALSE, # remove row numbers
                filter = "top", # add filter on top of columns
                extensions = "Buttons", # add download buttons
                caption=caption,
                options = list(
                  autoWidth = TRUE,
                  dom = "Blfrtip", # location of the download buttons
                  buttons = c("copy", "csv", "excel", "pdf", "print"), # download buttons
                  pageLength = 5, # show first 5 entries, default is 10
                  order = list(0, "asc") # order the title column by ascending order
                ),
                escape = FALSE # make URLs clickable) 
  )
}
```
**Function Parameter Definitions:**
- `table2show` - a `matrix` or `dataframe` containing tabular data to display
- `caption=` - a character vector to use as the caption for the table

</details>

##### filter_rare()
<details>
  <summary>filter out rare and non_microbial taxonomy assignments</summary>

  ```R
  filter_rare <- function(species_table, non_microbial, threshold=1){
    
    clean_tab_count  <-  species_table %>% 
      filter(str_detect(Species, non_microbial, negate = TRUE))  
    
    clean_tab <- clean_tab_count %>% 
      mutate( across( where(is.numeric)  , \(x) (x/sum(x, na.rm = TRUE))*100 ) )
    
    rownames(clean_tab) <- clean_tab$Species
    clean_tab  <- clean_tab[,-1] 
    
    
    # Get species with relative abundance less than 1% in all samples
    rare_species <- map(clean_tab, .f = \(col) rownames(clean_tab)[col < threshold])
    rare <- Reduce(intersect, rare_species)
    
    rownames(clean_tab_count) <- clean_tab_count$Species
    clean_tab_count  <- clean_tab_count[,-1] 
    
    abund_table <- clean_tab_count[!(rownames(clean_tab_count) %in% rare), ]
    
    return(abund_table)
  }
  ```
  **Function Parameter Definitions:**
  - `species_table` - the dataframe to filter
  - `non_microbial` - a character vector denoting the string used to identify a species as non-microbial
  - `threshold=` - abundance threshold used to determine if the relative abundance is rare, value denotes a percentage between 0 and 100.

  **Returns:** a dataframe with rare and non_microbrial assignemnts removed
</details>


##### make_plot()
<details>
  <summary>create bar plot of relative abundance</summary>

  ```R
  # Make bar plot
  make_plot <- function(abund_table, metadata, colors2use, publication_format){
    
    abund_table_wide <- abund_table %>% 
        as.data.frame() %>% 
        rownames_to_column("Sample_ID") %>% 
        inner_join(metadata) %>% 
        select(!!!colnames(metadata), everything()) %>% 
        mutate(Sample_ID = Sample_ID %>% str_remove("barcode"))
        
    abund_table_long <- abund_table_wide  %>%
        pivot_longer(-colnames(metadata), 
                    names_to = "Species",
                    values_to = "relative_abundance")
      
    p <- ggplot(abund_table_long, mapping = aes(x=Sample_ID, y=relative_abundance, fill=Species)) +
         geom_col() +
         scale_fill_manual(values = colors2use) + 
         labs(x=NULL, y="Relative Abundance (%)") + 
         publication_format

    return(p)
  }
  ```
  **Function Parameter Definitions:**
  - `abund_table` - a dataframe containing the data to plot
  - `metadata` - a vector of strings specifying the data to include in the plot
  - `colors2use` - a vector of strings specifying a custom color palette for coloring plots
  - `publication_format` - a ggplot::theme object specifying the custom theme for plotting

  **Returns:** a ggplot bar plot

</details>

##### get_colors2use()
<details>
  <summary>get colors to use in plots</summary>

  ```R
  get_colors2use <- function(species, expected_microbes, microbe_colors, custom_palette){
    
    unexpected_microbes <- setdiff(species, expected_microbes)
    
    start <- length(species)+1
    end <-  length(species) + length(unexpected_microbes)
    unexpected_microbes_colors <-  custom_palette[start:end]
    names(unexpected_microbes_colors) <- unexpected_microbes
    colors2use <- append(microbe_colors,unexpected_microbes_colors)
    return(colors2use)
    
  }
  ```
  **Function Parameter Definitions:**
  - `species` - a vector specifying the list of species that will use this color palette, used to set the number of colors in the palette
  - `expected_microbes` - the list of microbe species that were expected in the data
  - `microbe_colors` - colors assigned to the expected microbes
  - `custom_palette` - a vector of strings specifying a custom color palette

  **Returns:** a vector of strings specifying the color palette to use for the input species list

</details>

##### process_taxonomy()
<details>
  <summary>process taxonomy assignment table</summary>

  ```R
  process_taxonomy <- function(taxonomy, prefix='\\w__') {
    # Function to process a taxonopmy assignment table
    #1. ~ taxonomy is a string specifying the taxonomic assignment file name
    #2 prefix ~ is a regular expression specifying the characters to remove
    # from the taxon names  '\\w__'  for greengenes and 'D_\\d__' for SILVA
    
    
    taxonomy <- apply(X = taxonomy, MARGIN = 2, FUN = as.character) 
    
    #taxonomy[,'species'] <- paste(taxonomy[,'genus'],taxonomy[,'species'])
    # replace NAa with Other and delete the D_num__ prefix from the taxonomy names
    for (rank in colnames(taxonomy)) {
      #delete the taxonomy prefix
      taxonomy[,rank] <- gsub(pattern = prefix, x = taxonomy[, rank],
                              replacement = '')
      indices <- which(is.na(taxonomy[,rank]))
      taxonomy[indices, rank] <- rep(x = "Other", times=length(indices)) 
      #replace empty cell
      indices <- which(taxonomy[,rank] == "")
      taxonomy[indices,rank] <- rep(x = "Other", times=length(indices))
    }
    taxonomy <- apply(X = taxonomy,MARGIN = 2,
                      FUN =  gsub,pattern = "_",replacement = " ") %>% 
      as.data.frame(stringAsfactor=F)
    return(taxonomy)
  }
  ```
  **Function Parameter Definitions:**
  - `taxonomy` - string specifying the taxonomy assignment file path
  - `prefix=` - a regular expression specifying the characters to remove from the taxon names. `\\w__` for greengenes and `D_\\d__` for SILVA

  **Returns:** a dataframe containing the modified taxonomy values

</details>

##### format_taxonomy_table()
<details>
  <summary>format a taxonomy assignment table by appending a suffix to a known name</summary>

  ```R
  format_taxonomy_table <- function(taxonomy=taxonomy.m,stringToReplace="Other",
                                    suffix=";Other") {
    
    for (taxa_index in seq_along(taxonomy)) {
      #indices <- which(taxonomy[,taxa_index] == stringToReplace)
      
      indices <- grep(x = taxonomy[,taxa_index], pattern = stringToReplace)
      
      taxonomy[indices,taxa_index] <- 
        paste0(taxonomy[indices,taxa_index-1],
              rep(x = suffix, times=length(indices)))
      
    }
    return(taxonomy)
  }
  ```
  **Function Parameter Defintions:**
  - `taxonomy=` - a dataframe of the unformatted taxonomy table
  - `stringToReplace=` - a string specifying value to replace in the taxonomy known name
  - `suffix=` - a string containing the suffix to append to the known name

  **Returns:** a dataframe with the formatted taxonomy table
</details>

##### fix_names()
<details>
  <summary>regularize taxonomy names</summary>

  ```R
  fix_names<- function(taxonomy,stringToReplace,suffix){
    #1~ taxonomy is a taxonomy dataframe with taxonomy ranks as column names
    #2~ stringToReplace is a vector of regex strings specifying what to replace
    #3~ suffix is a string specifying the replacement value
    
    
    for(index in seq_along(stringToReplace)){
      taxonomy <- format_taxonomy_table(taxonomy = taxonomy,
                                        stringToReplace=stringToReplace[index], 
                                        suffix=suffix[index])
    }
    return(taxonomy)
  }
  ```
  **Function Parameter Defintions:**
  - `taxonomy=` - a dataframe containing the taxonomy table
  - `stringToReplace=` - a vector of regular expression strings specifying the values to replace
  - `suffix=` - a string specifying the replacement value

  **Returns:** a dataframe containing a taxonomy table with regularized names

</details>

##### read_input_table()
<details>
  <summary>read a tab-delimited file into a dataframe</summary>

  ```R
  # Function to read an input table into a dataframe
  read_input_table <- function(file_name){
    
    # Get depth from file name
    df <- read_delim(file = file_name, delim = "\t", comment = "#")
    return(df)
    
  }
  ```
  **Function Parameter Defintions:**
  - `filename=` - a tab-delimited file containing the data to import. Lines starting with `#` will be ignored.

  **Returns:** a dataframe containing a the data from the provided file

</details>

##### read_contig_table()
<details>
  <summary>read the assembly-based contig annotation table into a dataframe</summary>

  ```R
   # Read Assembly-based contig annotation table
  read_contig_table <- function(file_name, sample_names){
    
    df <- read_input_table(file_name)

    taxonomy_table <- df %>%
      select(domain:species) %>%
      mutate(domain=replace_na(domain, "Unclassified"))
    
    counts_table <- df %>% select(!!sample_names)

    taxonomy_table  <- process_taxonomy(taxonomy_table)
    taxonomy_table <- fix_names(taxonomy_table, "Other", ";_")

    df <- bind_cols(taxonomy_table, counts_table)
    
    return(df)
  }
  ```
  **Function Parameter Defintions:**
  - `file_name` - a tab-delimited file containing the assembly-based contig annotation table
  - `sample_names` - a list of sample names to use in the dataframe. Filters out samples in the file that are not present in the sample_names list.
  **Returns:** a dataframe containing a the data from the provided file


</details>



#### 25c. Set global variables

```R

kraken_taxonomy_outdir <- "kraken2_taxonomy/"
kaiju_taxonomy_outdir <- "kaiju_taxonomy/"

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
# remove white colors
custom_palette <- custom_palette[-c(21:23,
                                    grep(pattern = "white|snow|azure|gray|#FFFAFAFF|aliceblue",
                                         x = custom_palette, 
                                         ignore.case = TRUE)
                                   )
                                ]
# Define expected microbes to use for filtering
expected_microbes <- c("Pseudomonas aeruginosa", "Salmonella enterica",
                       "Limosilactobacillus fermentum", "Lactobacillus fermentum", "Staphylococcus aureus",
                       "Enterococcus faecalis", "Escherichia coli",
                       "Listeria monocytogenes", "Bacillus subtilis", "Bacillus spizizenii",
                       "Saccharomyces cerevisiae", "Cryptococcus neoformans")
orig_expected_microbes <- c("Pseudomonas aeruginosa", "Salmonella enterica",
                       "Limosilactobacillus fermentum", "Staphylococcus aureus",
                       "Enterococcus faecalis", "Escherichia coli",
                       "Listeria monocytogenes", "Bacillus spizizenii",
                       "Saccharomyces cerevisiae", "Cryptococcus neoformans")
orig_expected_microbes <- c(sort(orig_expected_microbes), "Escherichia phage Lambda")

# Define expected microbe color palette
microbe_colors <- custom_palette[1:length(orig_expected_microbes)]
names(microbe_colors) <- orig_expected_microbes

# Define human associated microbes
human_associated_microbes <- c("Staphylococcus epidermedis", "Staphylococcus hominis", "Cutibacterium acnes",
                               "Staphylococcus haemolyticus", "Malassezia", "Corynebacterium", "Micrococcus",
                               "Hoylesella shahii", "Streptococcus mitis",
                               "Eubacterium saphenum", "Lawsonella clevelandensis")

# subplots grouping variable
facets_kaiju <- c("Sample_Type","input_conc_ng", "lambda_spike")
facets_kraken2 <- c("Sample_Type","input_conc_ng")

assembly_directory <- here('data/assembly_annotation/kraken2/medaka')
assembly_based_summary <- glue("{assembly_directory}/Assembly-based-processing-overview_GLmetagenomics.tsv")
bin_overview <- glue("{assembly_directory}/bins-overview_GLmetagenomics.tsv")
mag_overview <- glue("{assembly_directory}/MAGs-overview_GLmetagenomics.tsv")
assembly_summary <- glue("{assembly_directory}/assembly-summaries_GLmetagenomics.tsv")
contig_taxonomy <- glue("{assembly_directory}/Combined-contig-level-taxonomy-coverages-CPM_GLmetagenomics.tsv")
gene_function <- glue("{assembly_directory}/Combined-gene-level-KO-function-coverages-CPM_GLmetagenomics.tsv")
gene_taxonomy <- glue("{assembly_directory}/Combined-gene-level-taxonomy-coverages-CPM_GLmetagenomics.tsv")

metdata_file <- here("mapping/metadata.txt")
plot_dir <- here('results/')

analysis_name <- "kraken-medaka"
sample_type <- "assembly"


```
**Input Data:** 

*No input data required*

**Output Data:**

- `publication_format` (a ggplot::theme object specifying the custom theme for plotting)
- `custom_palette` (a vector of strings specifying a custom color palette for coloring plots)
- `expected_microbes` (a vector of strings listing microbes that may be found in the samples)
- `orig_expected_microbes` (a vector of strings listing microbes that may be found in the samples plus "Escherichia phage Lambda")
- `microbe_colors` (a vector of strings specifying the custom color palette to use for coloring the `orig_expected_microbes`)
- `human_associated_microbes` (a vector of strings listing microbes that are known to be found in humans)
- `facets_kaiju` (a vector of strings listing subplot grouping variables for kaiju data)
- `facets_kraken2` (a vector of strings listing subplot grouping variables for kraken2 data)
- `kraken_taxonomy_outdir` (a path to the output folder for read taxonomy output based on kraken2 processing)
- `kaiju_taxonomy_outdir` (a path to the output folder for read taxonomy output based on kaiju processing)

#### 25d. Import Kaiju Taxonomy Data

```R
kaiju_table <- "/path/to/kaiju_read_taxonomy/merged_kaiju_summary_species.tsv"
feature_table <- process_kaiju_table(kaiju_table, taxon_col="taxon_name", remove_non_microbial = FALSE)$abundance_table

# Create Species table (species raw read count by barcode)
feature_table %>% as.data.frame %>% 
  rownames_to_column("Species") %>%
  pivot_longer(-Species, names_to = "Barcode", values_to = "Reads") %>% 
  write_delim(species_csv, delim=',')

## The number of reads classified at the species level
colSums(feature_table) %>%
  enframe(name = "Barcode", value = "Number of reads") %>%
  write_delim("{kaiju_taxonomy_outdir}species_counts{assay_suffix}.csv", delim=",")

species_table <- feature_table

```
**Input Data:**
- `kaiju_table` (the merged kaiju summary data at the species taxon level, output from [Step 9f](#9f-compile-kaiju-taxonomy-results))
- `kaiju_taxonomy_outdir` (a path to the output folder for read taxonomy output based on kaiju processing)
- `assay_suffix` (standard GeneLab assay suffix to use in output files)

**Output Data:**
- `species_table` (dataframe of relative abundance data)
- **kaiju_taxonomy/species_counts_GLMetagenomics.csv** (Number of reads classified for each species)


##### 25e. Import Kraken2 Taxonomy Data

```R
kraken_reports_dir <- "/path/to/read_taxonomy/kraken2_output/"

# import kraken2 reports
reports <- pavian::read_reports(kraken_reports_dir)

# create taxonomy overview
summary_table  <- pavian::summarize_reports(reports)
rownames(summary_table) <- rownames(summary_table) %>% str_split("-") %>% map_chr(\(x) pluck(x, 1))
summary_table %>% rownames_to_column("Sample_ID") %>% write_delim('{kraken_taxonomy_outdir}kraken_taxonomy_overview.csv', delim=',')

samples <- names(reports) %>% str_split("-") %>% map_chr(\(x) pluck(x, 1))
merged_reports  <- pavian::merge_reports2(reports, col_names = samples)
taxonReads <- merged_reports$taxonReads
cladeReads <- merged_reports$cladeReads
tax_data <- merged_reports[["tax_data"]]

#Create species table
species_table <- tax_data %>% 
  bind_cols(cladeReads) %>%
  filter(taxRank %in% c("U","S")) %>% 
  select(-contains("tax")) %>%
  zero_if_na() %>% 
  filter(name != 0) %>%  # drop unknown taxonomies
  group_by(name) %>% 
  summarise(across(everything(), sum)) %>% 
  ungroup() %>% 
  as.data.frame()

species_names <- species_table[,"name"]
rownames(species_table) <- species_names

taxonomy_col <- match("name", colnames(species_table))
species_table <- species_table[,-taxonomy_col]

species_table <- apply(X = species_table, MARGIN = 2, FUN = as.numeric)
rownames(species_table) <- species_names

# calculate total number of reads for each sample
colSums(species_table) %>%
  enframe(name = "Sample", value = "Number of reads") %>%
  write_delim("{kraken_taxonomy_outdir}species_counts{assay_suffix}.csv", delim=",")
```

**Input Data:**
- `kraken_reports` (the per-sample kraken reports, output from [Step 10a](#10a-taxonomic-classification))
- `kraken_taxonomy_outdir` (a path to the output folder for read taxonomy output based on kraken2 processing)
- `assay_suffix` (standard GeneLab assay suffix to use in output files)


**Output Data:**
- `species_table` (a dataframe of species raw read counts by barcode)
- **kraken_taxonomy/species_counts_GLMetagenomics.csv** (a dataframe of per-sample read counts)
- **kraken_taxonomy/kraken_taxonomy_overview.csv** (Comma-separated table containing a summary of Kraken2 taxonomy classification)


#### 25f. Import Sample Metadata

```R
# define input files
metadata_file <- "/path/to/metadata.txt"

# Import metadata
metadata <- read_delim(metdata_file , delim = "\t") %>% as.data.frame()
row.names(metadata) <- metadata$Sample_ID
```

**Input Data:** 

- `metadata_file` (a file containing sample metadata for the study, columns are: Sample_ID (string), Sample_Type (string), input_conc_ng (float), lambda_spike ('no'/'yes'), Sample_or_Control (string))

**Output Data:**

- `metadata` (a dataframe containing sample metadata for the study with the sampleIDs as the row names)

--- 

### 26. Feature table decontamination

The read-based feature table decontamination and taxonomy QC are performed using the same functions for both kraken2 and kaiju generated taxonomies.

#### 26a. Taxonomy filtering

```R
# with unclassified data
output_dir <- "{taxonomy_type}_taxonomy/"
abundance_threshold <- 0.5

species <- species_table %>% as.data.frame %>% 
  rownames_to_column("Species") %>% pull(Species) %>% unique()
colors2use <- get_colors2use(species, orig_expected_microbes, microbe_colors, custom_palette)

abund_table <- species_table %>% 
               as.data.frame %>% 
               mutate( across(everything(), \(x) (x/sum(x, na.rm = TRUE))*100 ) ) %>% 
               rownames_to_column("Species") 
  
rownames(abund_table) <- abund_table$Species
  
abund_table <- abund_table[, -match(x = "Species", colnames(abund_table))] %>% t

# excluding unclassified and host reads
non_microbial <- "Unclassified|unclassified|Homo sapien"

# Get species with relative abundance greater than 0.5 in all the samples
clean_tab <- species_table %>% 
  as.data.frame %>% 
  rownames_to_column("Species") 

abund_table <- filter_rare(clean_tab, non_microbial, threshold=abundance_threshold)
species <- rownames(abund_table)
colors2use <- get_colors2use(species, orig_expected_microbes, microbe_colors, custom_palette)

species_abund_table <- abund_table %>% 
                    as.data.frame %>% 
                   mutate( across( where(is.numeric)  , \(x) (x/sum(x, na.rm = TRUE))*100 ) )

abund_table <- species_abund_table %>% t

# Without human-associated microbes
unwanted <- str_c(c(non_microbial, human_associated_microbes), collapse = "|")
clean_tab2 <- filter_rare(clean_tab, unwanted, threshold=abundance_threshold)
clean_tab2 <- clean_tab2   %>% 
  mutate( across( where(is.numeric)  , \(x) (x/sum(x, na.rm = TRUE))*100 ) )
abund_table <- clean_tab2 %>% t
species <- rownames(clean_tab2)
colors2use <- get_colors2use(species, orig_expected_microbes, microbe_colors, custom_palette)

p <- make_plot(abund_table, metadata, colors2use, publication_format) + 
  facet_wrap(facets, scales = "free_x", nrow=1)

p$data %>% mutate(run="Ultra Low", host_read_removal="kraken", taxonomy="{taxonomy_type}") %>% write_delim(file="{output_dir}/{taxonomy_type}_no_unwanted{assay_suffix}.tsv", delim = "\t")

# Expected microbes alone 
non_microbial <- "Unclassifed|unclassified|Homo sapien"

clean_tab2 <- clean_tab %>% 
  filter(str_detect(Species, non_microbial, negate = TRUE))  %>% 
    filter(str_detect(Species, str_c(expected_microbes, collapse = "|"))) %>%  #select only the expected microbes
  mutate( across( where(is.numeric)  , \(x) (x/sum(x, na.rm = TRUE))*100 ) )

rownames(clean_tab2) <- clean_tab2$Species
clean_tab2  <- clean_tab2[,-1] 
abund_table <- clean_tab2 %>% t
species <- rownames(clean_tab2)
colors2use <- get_colors2use(species, orig_expected_microbes, microbe_colors, custom_palette)

p <- make_plot(abund_table, metadata, colors2use, publication_format) + 
  facet_wrap(facets, scales = "free_x", nrow=1)

p$data %>% mutate(run="Ultra Low", host_read_removal="kraken", taxonomy="{taxonomy_type}") %>% write_delim(file="{output_dir}{taxonomy_type}_expected{assay_suffix}.tsv", delim = "\t")

# Without Unclassified and host reads alone

# Get species with relative abundance greater than 1 in all the samples
clean_tab2 <- clean_tab %>% 
  as.data.frame %>% 
  filter(str_detect(Species, non_microbial, negate = TRUE))  %>% 
  mutate( across( where(is.numeric)  , \(x) (x/sum(x, na.rm = TRUE))*100 ) )

rownames(clean_tab2) <- clean_tab2$Species
clean_tab2  <- clean_tab2[,-1] 
abund_table <- clean_tab2 %>% t
species <- rownames(clean_tab2)
colors2use <- get_colors2use(species, orig_expected_microbes, microbe_colors, custom_palette)

p <- make_plot(abund_table, metadata, colors2use, publication_format) + 
  facet_wrap(facets, scales = "free_x", nrow=1)

#Without removing taxonomies with relative abundance less than 0.5%
p$data %>% mutate(run="Ultra Low", host_read_removal="kraken", taxonomy="{taxonomy_type}") %>% write_delim(file="{output_dir}{taxonomy_type}_no_filt{assay_suffix}.tsv", delim = "\t")

# Filter out unclassified, human reads and rare species

# Rare species here are classified as species with a relative abundance less than 0.5% across
# all samples.

# Get species with relative abundance greater than 0.5 in all the samples
abund_table <- filter_rare(clean_tab, non_microbial, threshold=abundance_threshold)
species <- rownames(abund_table)
colors2use <- get_colors2use(species, orig_expected_microbes, microbe_colors, custom_palette)

species_abund_table <- abund_table %>% 
                    as.data.frame %>% 
                   mutate( across( where(is.numeric)  , \(x) (x/sum(x, na.rm = TRUE))*100 ) )

abund_table <- species_abund_table %>% t

p <- make_plot(abund_table, metadata, colors2use, publication_format) + 
  facet_wrap(facets, scales = "free_x", nrow=1)

p$data %>% mutate(run="Ultra Low", host_read_removal="kraken", taxonomy="{taxonomy_type}") %>% write_delim(file="{output_dir}{taxonomy_type}_filtered{assay_suffix}.tsv", delim = "\t")
```

**Parameter Definitions:**
- `abundance_threshold` - threshold for defining rare species, default=0.5
- `taxonomy_type` - string specify which tool was used to create the input taxonomy, either `kaiju` or `kraken`
- `assay_suffix` - string specifying an assay suffix to use for output file nameing in this dataset (default: GLMetagenomics)


**Input Data:**
- `species_table` (dataframe of relative abundance data, from [Step 24d](#24d-import-kaiju-taxonomy-data) if using kaiju taxonomies or [Step 24e](#24e-import-kraken2-taxonomy-data) is using kraken taxonomies)
- `facets` (a vector of strings listing subplot grouping variables for either kaiju or kraken data, from [Step 24c](#24c-set-global-variables))


**Output Data**
- `species_abund_table` (a dataframe containing filtered realtive abundance values)
- **<kraken|kaiju>_taxonomy/<kraken|kaiju>_expected_GLMetagenomics.tsv** ()
- **<kraken|kaiju>_taxonomy/<kraken|kaiju>_no_filt_GLMetagenomics.tsv** ()
- **<kraken|kaiju>_taxonomy/<kraken|kaiju>_filtered_GLMetagenomics.tsv** ()

---

#### 26b. Decontamination with Decontam

##### 26b.i. Setup variables
```R
feature_table <- species_abund_table #species_table
sub_metadata <- metadata[colnames(feature_table),]
# Modify NTC concentration
sub_metadata <- sub_metadata %>% 
  mutate(input_conc_ng=map2_dbl(Sample_Type, input_conc_ng,
                                .f= function(type, conc) { 
                                  if(conc == 0) return(0.0000001) else return(conc) 
                                  } )
         )
sub_metadata$input_conc_ng <- as.numeric(sub_metadata$input_conc_ng)
ps <- phyloseq(otu_table(feature_table, taxa_are_rows = TRUE),
            sample_data(sub_metadata))
```

**Input Data:**
- `species_abund_table` (a dataframe containing filtered relative abundance values, from [Step ](#25a-taxonomy-filtering))

**Output Data:**
- `ps` (phyloseq object of the relative abundance values with NTC metadata added)

##### 26b.ii. Identify prevalence of contaminant sequences
The prevalence (presence/absence across samples) of each sequence feature in 
true positive samples is compared to the prevalence in negative controls to 
identify contaminants.

```R
contam_threshold <- 0.1
output_dir <- "{taxonomy_type}_taxonomy_decontam/"
# In our phyloseq object, "Sample_or_Control" is the sample variable that holds 
# the negative control information. We’ll summarize that data as a logical 
# variable, with TRUE for control samples, as that is the form required by isContaminant
sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control_Sample"
contamdf <- isContaminant(ps, neg="is.neg", conc="input_conc_ng", threshold=contam_threshold) # threshold

#### Create contaminant table
contamdf %>%
  mutate( across( where(is.numeric), \(x) round(x, digits = 2) ) ) %>%
  rownames_to_column("Species") %>% 
  write_delim(file="{output_dir}{taxonomy_type}_contaminant_table{assay_suffix}.tsv", delim = "\t")

table(contamdf$contaminant)

contamdf %>% filter(contaminant == TRUE) %>% 
  write_delim(file="{output_dir}{taxonomy_type}_filtered_contaminant_table{assay_suffix}.tsv", delim = "\t")


isExpected <- str_detect(rownames(contamdf), pattern = str_c(expected_microbes, collapse = "|"))
contamdf[isExpected,] %>%
  select(-p.freq) %>%
  mutate( across( where(is.numeric), \(x) round(x, digits = 3) ) ) %>% 
  write_delim(file="{output_dir}{taxonomy_type}_contaminant_table_expected_microbes{assay_suffix}.tsv", delim = "\t")
```

**Parameter Defintitions:**
- `contam_threshold` - probability threshold below which the null hypothesis (not a contaminant) should be rejected in favor of the alternate hypothesis (contaminant) (default: 0.1)
- `taxonomy_type` - string specify which tool was used to create the input taxonomy, either `kaiju` or `kraken`
- `assay_suffix` - string specifying an assay suffix to use for output file nameing in this dataset (default: GLMetagenomics)


**Input Data:**
- `ps` (phyloseq object of the relative abundance values with NTC metadata added, from [Step ](#25bi-setup-variables))

**Output Data:**
- `contam_df` (dataframe of contaminant table)
- **<kaiju|kraken>_taxonomy_decontam/<kaiju|kraken>_contaminant_table_GLMetagenomics.tsv** (tab-delimited table of classification information for all input sequences)
- **<kaiju|kraken>_taxonomy_decontam/<kaiju|kraken>_filtered_contaminant_table_GLMetagenomics.tsv** (tab-delimited table of classification information for all sequences identified as contaminants)
- **<kaiju|kraken>_taxonomy_decontam/<kaiju|kraken>_contaminant_table_expected_microbes_GLMetagenomics.tsv** (tab-delimited table of classification information for expected microbes)

##### 26b.iii. Decontaminated taxonomy plots

```R
output_dir <- "{taxonomy_type}_taxonomy_decontam/"
contaminants <- contamdf %>%
  as.data.frame %>%
  rownames_to_column("Species") %>%
  filter(contaminant == TRUE) %>% pull(Species)
species <- species_abund_table  %>% 
  as.data.frame %>% 
  rownames_to_column("Species") %>%
  filter(str_detect(Species, pattern = str_c(contaminants, collapse = "|"), negate = TRUE)) %>%
  pull(Species) %>%
  unique()
colors2use <- get_colors2use(species, orig_expected_microbes, microbe_colors, custom_palette)

abund_table <- species_abund_table %>% 
                    as.data.frame  %>% 
                    rownames_to_column("Species") %>% 
                    filter(str_detect(Species, 
                                      pattern = str_c(contaminants,
                                                      collapse = "|"),
                                      negate = TRUE)) %>%
                    mutate( across( where(is.numeric)   , \(x) (x/sum(x, na.rm = TRUE))*100 ) )
  
rownames(abund_table) <- abund_table$Species
  
abund_table <- abund_table[, -match(x = "Species", colnames(abund_table))] %>% t
  
abund_table_wide <- abund_table %>% 
    as.data.frame() %>% 
    rownames_to_column("Sample_ID") %>% 
    inner_join(metadata) %>% 
    select(!!!colnames(metadata), everything()) %>% 
    mutate(Sample_ID = Sample_ID %>% str_remove("barcode"))
    
  
abund_table_long <- abund_table_wide  %>%
    pivot_longer(-colnames(metadata), 
                 names_to = "Species",
                 values_to = "relative_abundance")
  
p <- ggplot(abund_table_long, mapping = aes(x=Sample_ID, 
                                              y=relative_abundance, fill=Species)) +
    geom_col() +
    scale_fill_manual(values = colors2use) + 
    labs(x=NULL, y="Relative Abundance (%)") + 
    publication_format + 
  facet_wrap(facets, scales = "free_x", nrow=1)

#### Taxonomy plot without contaminants

# Taxonomy plot after contaminant removal at a set threshold of 0.1
# ggsave(filename = "results/species_plot.png", plot = p,
#          device = "png", width = 10, height = 6, units = "in", dpi = 300)
ggplotly(p) %>% saveWidget(file = "{output_dir}{taxonomy_type}_taxonomy_plots_no_contam{assay_suffix}.html")
```
**Parameter Definitions:**
- `taxonomy_type` - string specify which tool was used to create the input taxonomy, either `kaiju` or `kraken`
- `assay_suffix` - string specifying an assay suffix to use for output file nameing in this dataset (default: GLMetagenomics)

**Input Data:**
- `species_abund_table` ()
- `contam_df` ()

**Output Data:**
- **<kaiju|kraken>_taxonomy_decontam/<kaiju|kraken>_taxonomy_plots_no_contam_GLMetagenomics.html** (Plot of taxonomies for decontaminated data)

---

### 27. Assembly-based processing decontamination
Medaka assembly annotation of kraken decontaminated low biomass samples
Quality filtered and trimmed reads were decontaminated (host (human) reads filtered out) using kraken2. Assembly of the clean reads was performed using metaflye followed by polishing with medaka. The polished assembly was annotated using our standard assembly annotation pipeline with prodigal used to predict genes, CAT used for taxonomy assignment of genes and contigs and KOFamScan for genes functional annotation.  

#### 27a. 

```R

# assembly based summary

overview_table <-  read_input_table(assembly_summary)

col_names <- names(overview_table) %>% str_remove_all("-assembly")

overview_table <- overview_table %>%
  set_names(col_names) %>% 
  select(!!sort(col_names)) %?%

write_delim(overview_table, "")

# Contig annotation
sample_order <- col_names[-1] %>% sort()
contig_table <-  read_contig_table(contig_taxonomy, sample_order)

# expected microbes
species_contig_table %>% 
  filter( str_detect( species, str_c(expected_microbes, collapse = "|") ) ) %>% 
  knitr::kable()

contig.m <- species_contig_table %>%
  group_by(species) %>%
  summarise(across(everything(), sum)) %>%
  as.data.frame()

rownames(contig.m) <- contig.m[['species']]
contig.m <- contig.m[,-match("species", colnames(contig.m))] %>% as.matrix()

colours <- colorRampPalette(c('white','red'))(255)
colnames(contig.m) <-  colnames(contig.m) %>% str_remove_all("barcode")

main_groups <- c("Sample_Type", "input_conc_ng", "lambda_spike" )
annotation_colors <- list(Sample_Type= c(Even="red", NTC="blue"),
                          input_conc_ng = c("0.005"="gold4", "0.010"="darkred",
                                            "0.000"="pink", "0.100"="green"),
                          lambda_spike = c(no="black", yes="purple"))
col_annotation <- as.data.frame(metadata)[,rev(main_groups), drop=FALSE]

# Heatmap of normalized counts
pheatmap(mat = contig.m,
         cluster_cols = FALSE, 
         cluster_rows = FALSE, 
         col = colours, 
         angle_col = 0, 
         display_numbers = TRUE,
         fontsize = 12,
         number_format = "%.0f",
         annotation_col = col_annotation,
        annotation_colors = annotation_colors)

# Gene Annotation
gene_taxonomy_table <-  read_contig_table(gene_taxonomy, sample_order)

species_gene_table <- gene_taxonomy_table %>%
  select(species, !!sample_order) %>% 
  group_by(species) %>% 
  summarise(across(everything(), sum)) 
create_dt( species_gene_table %>%  mutate(across(where(is.numeric), \(x) round(x,digits = 2))) )

# expected microbes
species_gene_table %>% 
  filter( str_detect( species, str_c(expected_microbes, collapse = "|") ) ) %>% 
  knitr::kable()

# heatmap of normalized counts
gene.m <- species_gene_table %>% as.data.frame()

rownames(gene.m) <- gene.m[['species']]
gene.m <- gene.m[,-match("species", colnames(gene.m))] %>% as.matrix()
colnames(gene.m) <-  colnames(gene.m) %>% str_remove_all("barcode")

# with unclassified
pheatmap(mat = gene.m,
         cluster_cols = FALSE, 
         cluster_rows = FALSE, 
         col = colours, 
         angle_col = 0, 
         display_numbers = TRUE,
         fontsize = 12, 
         number_format = "%.0f",
        annotation_col = col_annotation,
        annotation_colors = annotation_colors)

# without unclassified
heatmap(mat = gene.m[-match("Unclassified;_;_;_;_;_;_", rownames(gene.m)),],
         cluster_cols = FALSE, 
         cluster_rows = FALSE, 
         col = colours, 
         angle_col = 0, 
         display_numbers = TRUE,
         fontsize = 12, 
         number_format = "%.0f", 
         annotation_col = col_annotation,
        annotation_colors = annotation_colors)

# abundant taxa without unclassified
taxa <- rowSums(gene.m) %>% sort()
abund_taxa <- taxa[taxa>2000] %>% names

abund_gene.m <- gene.m[abund_taxa,]
pheatmap(mat = abund_gene.m[-match("Unclassified;_;_;_;_;_;_", rownames(abund_gene.m)),],
         cluster_cols = FALSE, 
         cluster_rows = FALSE, 
         col = colours, 
         angle_col = 0, 
         display_numbers = TRUE,
         fontsize = 12, 
         number_format = "%.0f",
         annotation_col = col_annotation,
        annotation_colors = annotation_colors)

