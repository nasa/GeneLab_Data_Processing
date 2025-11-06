# Bioinformatics pipeline for amplicon Illumina sequencing data  

> **This page holds an overview and instructions for how GeneLab processes Illumina amplicon sequencing datasets. Exact processing commands for specific datasets that have been released are available in the [GLDS_Processing_Scripts](../GLDS_Processing_Scripts) sub-directory and/or are provided with their processed data in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/).**  

---

**Date:** May 5, 2025  
**Revision:** C  
**Document Number:** GL-DPPD-7104  

**Submitted by:**  
Olabiyi Obayomi, Alexis Torres, and Michael D. Lee (GeneLab Data Processing Team)

**Approved by:**  
Samrawit Gebre (OSDR Project Manager)  
Danielle Lopez (OSDR Deputy Project Manager)  
Jonathan Galazka (OSDR Project Scientist)  
Amanda Saravia-Butler (GeneLab Science Lead)  
Barbara Novak (GeneLab Data Processing Lead)  

---

## Updates from previous version

Software Updates and Changes:

| Program      | Previous Version | New Version   |
|:-------------|:-----------------|:--------------|
| FastQC       | 0.11.9           | 0.12.1        |
| MultiQC      | 1.9              | 1.27.1        |
| Cutadapt     | 2.3              | 5.0           |
| R-base       | 4.1.1            | 4.4.2         |
| DADA2        | 1.20.0           | 1.34.0        |
| DECIPHER     | 2.20.0           | 3.2.0         |
| biomformat   | 1.20.0           | 1.34.0        |
| ANCOMBC      | N/A              | 2.8.0         |
| broom        | N/A              | 1.0.7         |
| DescTools    | N/A              | 0.99.59       |
| DESeq2       | N/A              | 1.46.0        |
| dp_tools     | N/A              | 1.3.8         |
| FSA          | N/A              | 0.9.6         |
| ggdendro     | N/A              | 0.2.0         |
| ggrepel      | N/A              | 0.9.6         |
| ggplot2      | N/A              | 3.5.1         |
| glue         | N/A              | 1.8.0         |
| hexbin       | N/A              | 1.28.3        |
| mia          | N/A              | 1.14.0        |
| phyloseq     | N/A              | 1.50.0        |
| rcolorbrewer | N/A              | 1.1.3         |
| taxize       | N/A              | 0.10.0        |
| tidyverse    | N/A              | 2.0.0         |
| vegan        | N/A              | 2.6-10        |
| vsn          | N/A              | 3.74.0        |
| patchwork    | N/A              | 1.3.0         |
| rstatix      | N/A              | 0.7.2         |
| multcompView | N/A              | 0.1-10        |
| scales       | N/A              | 1.4.0         |
| dendextend   | N/A              | 1.19.0        |

- Added new processing steps in R to generate processed data outputs for alpha and beta diversity, taxonomic summary plots, and differential abundance:
  - Alpha Diversity Analysis ([Step 7](#7-alpha-diversity-analysis))
  - Beta Diversity Analysis ([Step 8](#8-beta-diversity-analysis))
  - Group-wise and Sample-wise Taxonomic Summary Plots ([Step 9](#9-taxonomy-plots))
  - Differential Abundance Testing ([Step 10](#9-differential-abundance-analysis)) with 
    ANCOMBC 1 ([Step 10a](#10a-ancombc-1)), ANCOMBC 2 ([Step 10b](#10b-ancombc-2)), and Deseq2 ([Step 10c](#10c-deseq2))
- Assay-specific suffixes were added where needed for OSDR ("_GLAmpSeq")
- Updated [DECIPHER](https://www2.decipher.codes/data/Downloads/TrainingSets/) reference files to the following:
  - ITS UNITE: "UNITE\_v2024\_April2024.RData" 
  - SILVA SSU r138: "SILVA\_SSU\_r138\_2\_2024.RData"
  - PR2 v4.13: "PR2\_v4\_13\_March2021.RData"
- Added persistent reference links to DECIPHER databases on Figshare and replaced reference links to 
  DECIPHER [website](https://www2.decipher.codes/data/Downloads/TrainingSets/)  
  - [SILVA SSU r138](https://figshare.com/ndownloader/files/52846199)
  - [UNITE v2024](https://figshare.com/ndownloader/files/52846346)
  - [PR2 v4.13](https://figshare.com/ndownloader/files/46241917)

---

# Table of contents  

- [**Software used**](#software-used)
- [**Reference databases used**](#reference-databases-used)
- [**General processing overview with example commands**](#general-processing-overview-with-example-commands)
  - [**1. Raw Data QC**](#1-raw-data-qc)
    - [1a. Raw Data QC](#1a-raw-data-qc)
    - [1b. Compile Raw Data QC](#1b-compile-raw-data-qc)
  - [**2. Trim Primers**](#2-trim-primers)
  - [**3. Quality Filtering**](#3-quality-filtering)
  - [**4. Filtered Data QC**](#4-filtered-data-qc)
    - [4a. Filtered Data QC](#4a-filtered-data-qc)
    - [4b. Compile Filtered Data QC](#4b-compile-filtered-data-qc)
  - [**5. Calculate Error model, Apply DADA2 Algorithm, Assign Taxonomy, and Create Output Tables**](#5-calculate-error-model-apply-dada2-algorithm-assign-taxonomy-and-create-output-tables)
    - [5a. Learning the Error Rates](#5a-learning-the-error-rates)
    - [5b. Inferring Sequences](#5b-inferring-sequences)
    - [5c. Merging Forward and Reverse Reads; Not Needed if Data are Single-End](#5c-merging-forward-and-reverse-reads-not-needed-if-data-are-single-end)
    - [5d. Generating Sequence Table with Counts per Sample](#5d-generating-sequence-table-with-counts-per-sample)
    - [5e. Removing Putative Chimeras](#5e-removing-putative-chimeras)
    - [5f. Assigning Taxonomy](#5f-assigning-taxonomy)
    - [5g. Generating and Writing Standard Outputs](#5g-generating-and-writing-standard-outputs)
  - [**6. Amplicon Seq Data Analysis Set Up**](#6-amplicon-seq-data-analysis-set-up)
    - [6a. Create Sample Runsheet](#6a-create-sample-runsheet)
    - [6b. R Environment Set Up](#6b-r-environment-set-up)
      - [6b.i. Load Libraries](#6bi-load-libraries)
      - [6b.ii. Define Custom Functions](#6bii-define-custom-functions)
      - [6b.iii. Set Variables](#6biii-set-variables)
      - [6b.iv. Read-in Input Tables](#6biv-read-in-input-tables)
      - [6b.v. Preprocessing](#6bv-preprocessing)
  - [**7. Alpha Diversity Analysis**](#7-alpha-diversity-analysis)
    - [7a. Rarefaction Curves](#7a-rarefaction-curves)
    - [7b. Richness and Diversity Estimates](#7b-richness-and-diversity-estimates)
    - [7c. Plot Richness and Diversity Estimates](#7c-plot-richness-and-diversity-estimates)
    - [7d. Package Alpha Diversity Plots](#7d-package-alpha-diversity-plots)
  - [**8. Beta Diversity Analysis**](#8-beta-diversity-analysis)
  - [**9. Taxonomy Plots**](#9-taxonomy-plots)
  - [**10. Differential Abundance Testing**](#10-differential-abundance-testing)
    - [10a. ANCOMBC 1](#10a-ancombc-1)
    - [10b. ANCOMBC 2](#10b-ancombc-2)
    - [10c. DESeq2 ](#10c-deseq2)

---

# Software used  

|Program|Version|Relevant Links|
|:------|:-----:|:-------------|
|FastQC|0.12.1|[https://www.bioinformatics.babraham.ac.uk/projects/fastqc/](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)|
|MultiQC|1.27.1|[https://multiqc.info/](https://multiqc.info/)|
|Cutadapt|5.0|[https://cutadapt.readthedocs.io/en/stable/](https://cutadapt.readthedocs.io/en/stable/)|
|R-base|4.4.2|[https://www.r-project.org/](https://www.r-project.org/)|
|DADA2|1.34.0|[https://www.bioconductor.org/packages/release/bioc/html/dada2.html](https://www.bioconductor.org/packages/release/bioc/html/dada2.html)|
|DECIPHER|3.2.0|[https://bioconductor.org/packages/release/bioc/html/DECIPHER.html](https://bioconductor.org/packages/release/bioc/html/DECIPHER.html)|
|biomformat|1.34.0|[https://github.com/joey711/biomformat](https://github.com/joey711/biomformat)|
|ANCOMBC|2.8.0|[https://github.com/FrederickHuangLin/ANCOMBC](https://github.com/FrederickHuangLin/ANCOMBC)|
|broom|1.0.7|[https://CRAN.R-project.org/package=broom](https://CRAN.R-project.org/package=broom)|
|DescTools|0.99.59|[https://andrisignorell.github.io/DescTools/](https://andrisignorell.github.io/DescTools/)|
|DESeq2|1.46.0|[https://bioconductor.org/packages/release/bioc/html/DESeq2.html](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)|
|FSA|0.9.6|[https://CRAN.R-project.org/package=FSA](https://CRAN.R-project.org/package=FSA)|
|ggdendro|0.2.0|[https://CRAN.R-project.org/package=ggdendro](https://CRAN.R-project.org/package=ggdendro)|
|ggrepel|0.9.6|[https://CRAN.R-project.org/package=ggrepel](https://CRAN.R-project.org/package=ggrepel)|
|ggplot2|3.5.1|[https://CRAN.R-project.org/package=ggplot2](https://CRAN.R-project.org/package=ggplot2)|
|glue|1.8.0|[https://glue.tidyverse.org/](https://glue.tidyverse.org/)|
|hexbin|1.28.3|[https://CRAN.R-project.org/package=hexbin](https://CRAN.R-project.org/package=hexbin)|
|mia|1.14.0|[https://github.com/microbiome/mia](https://github.com/microbiome/mia)|
|phyloseq|1.50.0|[https://bioconductor.org/packages/release/bioc/html/phyloseq.html](https://bioconductor.org/packages/release/bioc/html/phyloseq.html)|
|rcolorbrewer|1.1.3|[https://CRAN.R-project.org/package=RColorBrewer](https://CRAN.R-project.org/package=RColorBrewer)|
|taxize|0.10.0|[https://docs.ropensci.org/taxize/](https://docs.ropensci.org/taxize/)|
|tidyverse|2.0.0|[https://CRAN.R-project.org/package=tidyverse](https://CRAN.R-project.org/package=tidyverse)|
|vegan|2.6-10|[https://cran.r-project.org/package=vegan](https://cran.r-project.org/package=vegan)|
|vsn|3.74.0|[https://bioconductor.org/packages/release/bioc/html/vsn.html](https://bioconductor.org/packages/release/bioc/html/vsn.html)|
|patchwork|1.3.0|[https://CRAN.R-project.org/package=patchwork](https://CRAN.R-project.org/package=patchwork)|
|rstatix|0.7.2|[https://CRAN.R-project.org/package=rstatix](https://CRAN.R-project.org/package=rstatix)|
|multcompView|0.1-10|[https://CRAN.R-project.org/package=multcompView](https://CRAN.R-project.org/package=multcompView)|
|scales|1.4.0|[https://CRAN.R-project.org/package=scales](https://CRAN.R-project.org/package=scales)|
|dendextend|1.19.0|[https://CRAN.R-project.org/package=dendextend](https://CRAN.R-project.org/package=dendextend)|

# Reference databases used
<update figshare links once the updated DBs are downloaded>
  
|Program used|Database|DECIPHER Link|GeneLab Figshare Link|GeneLab Download Date|
|:-----------|:------:|:------------|--------------------:|--------------------:|
|DECIPHER| SILVA SSU r138_2 | [https://www2.decipher.codes/data/Downloads/TrainingSets/SILVA_SSU_r138_2_2024.RData](https://www2.decipher.codes/data/Downloads/TrainingSets/SILVA_SSU_r138_2_2024.RData) |[SILVA_SSU_r138_2_2024.RData](https://figshare.com/ndownloader/files/52846199)| 03/06/2025 |
|DECIPHER| UNITE v2024 | [https://www2.decipher.codes/data/Downloads/TrainingSets/UNITE_v2024_April2024.RData](https://www2.decipher.codes/data/Downloads/TrainingSets/UNITE_v2024_April2024.RData) | [UNITE_v2024_April2024.RData](https://figshare.com/ndownloader/files/52846346)| 03/06/2025 |
|DECIPHER| PR2 v4.13 | [https://www2.decipher.codes/data/Downloads/TrainingSets/PR2_v4_13_March2021.RData](https://www2.decipher.codes/data/Downloads/TrainingSets/PR2_v4_13_March2021.RData) | [PR2_v4_13_March2021.RData](https://figshare.com/ndownloader/files/46241917)| 05/10/2024 |
---

# General processing overview with example commands  

> Exact processing commands for specific datasets are available in the [GLDS_Processing_Scripts](../GLDS_Processing_Scripts) sub-directory of this repository, and/or are provided with their processed data in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/).
>
> Output files listed in **bold** below are included with each Amplicon Seq processed dataset in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/).

---

## 1. Raw Data QC  

<br>

### 1a. Raw Data QC  

```bash
fastqc -o raw_fastqc_output *.fastq.gz
```

**Parameter Definitions:**

* `-o` – the output directory to store results
* `*.fastq.gz` – the input reads are specified as a positional argument, and can be given all at once with wildcards like this, or as individual arguments with spaces in between them

**Input Data:**

* \*fastq.gz (raw reads)

**Output Data:**

* \*fastqc.html (FastQC output html summary)
* \*fastqc.zip (FastQC output data)


<br>  

### 1b. Compile Raw Data QC  

```bash
multiqc --interactive -n raw_multiqc_<tech_type>_GLAmpSeq \
 -o /path/to/raw_multiqc/output/raw_multiqc_<tech_type>_GLAmpSeq_report \
 /path/to/directory/containing/raw_fastqc/files

zip -q -r raw_multiqc_<tech_type>_GLAmpSeq_report.zip raw_multiqc_<tech_type>_GLAmpSeq_report
```

**Parameter Definitions:**
**multiqc**
- `--interactive` – force reports to use interactive plots
- `-n` – prefix name for output files
- `-o` – the output directory to store results
- `/path/to/directory/containing/raw_fastqc/files` – the directory holding the output data from the FastQC run, provided as a positional argument

**zip**
- `-q` – quiet operation
- `-r` - recurse into directories
- `raw_multiqc_<tech_type>_GLAmpSeq_report.zip` – positional argument naming the zip output file
- `raw_multiqc_<tech_type>_GLAmpSeq_report` – positional argument naming the input folder to package

**Input Data:**

* \*fastqc.zip (FastQC output data, output from [Step 1a](#1a-raw-data-qc))

**Output Data:**

* **raw_multiqc_<tech_type>_GLAmpSeq_report.zip** (zip containing the following)
  * **raw_multiqc_<tech_type>_GLAmpSeq.html** (multiqc output html summary)
  * **raw_multiqc_<tech_type>_GLAmpSeq_data** (directory containing multiqc output data)

<br>  

---

## 2. Trim Primers  

The location and orientation of primers in the data is important to understand in deciding how to do this step. `cutadapt` has many options for primer identification and removal, which are described in detail in the [cutadapt adapter type documentation](https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types).  

The following example commands show how it was done for some samples of [GLDS-200](https://osdr.nasa.gov/bio/repo/data/studies/OSD-200), which was 2x250 sequencing of the 16S gene using these primers:  
* forward: 5'-GTGCCAGCMGCCGCGGTAA-3'  
* reverse: 5'-GGACTACVSGGGTATCTAAT-3'  

Due to the size of the target amplicon and the type of sequencing done here, both forward and reverse primers are expected to be on each of the forward and reverse reads. It therefore takes "linked" primers as input for forward and reverse reads, specified in the example command below by the `...` between them. It also expects that the primers start at the first position of the reads ("anchored"), specified with the leading `^` characters in the example command below.  

The following website is useful for reverse complementing primers and dealing with degenerate bases appropriately: [http://arep.med.harvard.edu/labgc/adnan/projects/Utilities/revcomp.html](http://arep.med.harvard.edu/labgc/adnan/projects/Utilities/revcomp.html)  

```bash
cutadapt -a ^GTGCCAGCMGCCGCGGTAA...ATTAGATACCCSBGTAGTCC -A ^GGACTACVSGGGTATCTAAT...TTACCGCGGCKGCTGGCAC \
         -o sample1_<tech_type>_GLAmpSeq_R1_trimmed.fastq.gz \
         -p sample1_<tech_type>_GLAmpSeq_R2_trimmed.fastq.gz \
         sample1_<tech_type>_GLAmpSeq_R1_raw.fastq.gz sample1_<tech_type>_GLAmpSeq_R2_raw.fastq.gz \
         --discard-untrimmed
```

**Parameter Definitions:**

*	`-a` – specifies the primers and orientations expected on the forward reads (when primers are linked as noted above)
*	`-A` – specifies the primers and orientations expected on the reverse reads (when primers are linked as noted above)
*	`-o` – specifies file path/name of forward, primer-trimmed reads
*	`-p` – specifies file path/name of reverse, primer-trimmed reads
*	`sample1_<tech_type>_GLAmpSeq_R1_raw.fastq.gz` – this and following "R2" file are positional arguments specifying the forward and reverse reads, respectively, for input
*	`--discard-untrimmed` – this filters out those reads where the primers were not found as expected

**Input Data:**

* \*fastq.gz (raw reads)

**Output Data:**

* **\*trimmed.fastq.gz** (trimmed reads)
* **trimmed-read-counts_<tech_type>_GLAmpSeq.tsv** (per sample read counts before and after trimming)
* **cutadapt_<tech_type>_GLAmpSeq.log** (log file of standard output and error from cutadapt)

<br>

---

## 3. Quality Filtering
> The following is run in an R environment.  

Specific settings required will depend on the dataset being processing. These include parameters such as `truncLen`, which might depend on the target amplicon and its size, and `maxEE` which might depend on the quality of the sequencing run. For instance, when working with ITS data, it may be preferable to omit using the `truncLen` parameter if the target amplified region is expected to vary to lengths greater than the read size. More information on these parameters can be found at these sites:  
* [https://benjjneb.github.io/dada2/tutorial.html](https://benjjneb.github.io/dada2/tutorial.html)  
* [https://astrobiomike.github.io/amplicon/dada2_workflow_ex](https://astrobiomike.github.io/amplicon/dada2_workflow_ex)  


The following is an example from a [GLDS-200](https://osdr.nasa.gov/bio/repo/data/studies/OSD-200) sample that used paired-end 2x250 sequencing with the following 16S primers:  
* forward: 5'-GTGCCAGCMGCCGCGGTAA-3'
* reverse: 5'- GGACTACVSGGGTATCTAAT-3'

```R
filtered_out <- filterAndTrim(fwd="sample1_<tech_type>_GLAmpSeq_R1_trimmed.fastq.gz",
                              filt="sample1_<tech_type>_GLAmpSeq_R1_filtered.fastq.gz",
                              rev="sample1_<tech_type>_GLAmpSeq_R2_trimmed.fastq.gz", 
                              filt.rev="sample1_<tech_type>_GLAmpSeq_R1_filtered.fastq.gz",
                              truncLen=c(220, 160), maxN=0, maxEE=c(2,2),
                              truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)
```

**Parameter Definitions:**

*	`filtered_out <-` – specifies the variable that will store the summary results within in our R environment
*	`filterAndTrim()` – the DADA2 function we are calling, with the following parameters set within it
*	`fwd=` – specifying the path to the forward reads, here "sample1_<tech_type>_GLAmpSeq_R1_trimmed.fastq.gz"
*	`filt=` – specifying the path to where the output forward reads will be written
*	`rev=` – specifying the path to the reverse reads, here "sample1_<tech_type>_GLAmpSeq_R2_trimmed.fastq.gz"; only applicable if paired-end
*	`filt.rev=` – specifying the path to where the output reverse reads will be written; only applicable if paired-end
*	`truncLen=c(220, 160)` – specifying the forward reads to be truncated at 220 bp, and the reverse to be truncated at 160 bps (note that this parameter also functions as a minimum-length filter); would only have 1 value if not paired-end
*	`maxN=0` – setting the maximum allowed Ns to 0, any reads with an N will be filtered out
*	`maxEE=c(2,2)` – setting maximum expected error allowed to 2 for each forward and reverse read; would only have 1 value if not paired-end
*	`truncQ=2` – looking from the lower-quality end of each read, truncate at the first base with a quality score lower than 2
*	`rm.phix=TRUE` – filter out reads with exact kmers matching the PhiX genome
*	`compress=TRUE` – gzip-compress the output filtered reads
*	`multithread=TRUE` – determine number of cores available and run in parallel when possible (can also take an integer specifying the number to run)

**Input Data:**

* \*.trimmed.fastq.gz (primer-trimmed reads, output from [Step 2](#2-trim-primers))

**Output Data:**

* **\*filtered.fastq.gz** (filtered reads)
* **filtered-read-counts_<tech_type>_GLAmpSeq.tsv** (a tab-separated file containing per sample read counts before and after filtering)

<br>

---

## 4. Filtered Data QC

<br>

### 4a. Filtered Data QC
```bash
fastqc -o filtered_fastqc_output/ *filtered.fastq.gz
```

**Parameter Definitions:**

*	`-o` – the output directory to store results  
*	`*filtered.fastq.gz` – the input reads are specified as a positional argument, and can be given all at once with wildcards like this, or as individual arguments with spaces in between them  

**Input Data:**

* \*fastq.gz (filtered reads)

**Output Data:**

* \*fastqc.html (FastQC output html summary)
* \*fastqc.zip (FastQC output data)

<br>

### 4b. Compile Filtered Data QC
```bash
multiqc --interactive -n filtered_multiqc_<tech_type>_GLAmpSeq \
 -o /path/to/filtered_multiqc/output/filtered_multiqc_<tech_type>_GLAmpSeq_report \
 /path/to/directory/containing/filtered_fastqc/files

zip -q -r filtered_multiqc_<tech_type>_GLAmpSeq_report.zip filtered_multiqc_<tech_type>_GLAmpSeq_report
```

**Parameter Definitions:**
**multiqc**
- `--interactive` – force reports to use interactive plots
- `-n` – prefix name for output files
- `-o` – the output directory to store results
- `/path/to/directory/containing/filtered_fastqc/files` – the directory holding the output data from the FastQC run, provided as a positional argument

**zip**
- `-q` – quiet operation
- `-r` - recurse into directories
- `filtered_multiqc_<tech_type>_GLAmpSeq_report.zip` – positional argument naming the zip output file
- `filtered_multiqc_<tech_type>_GLAmpSeq_report` – positional argument naming the input folder to package

**Input Data:**

* \*fastqc.zip (FastQC output data, output from [Step 4a](#4a-filtered-data-qc))

**Output Data:**

* **filtered_multiqc_<tech_type>_GLAmpSeq_report.zip** (zip containing the following)
  * **filtered_multiqc_<tech_type>_GLAmpSeq_report.html** (multiqc output html summary)
  * **filtered_multiqc_<tech_type>_GLAmpSeq_data** (directory containing multiqc output data)

<br>

---

## 5. Calculate Error Mdel, Apply DADA2 Algorithm, Assign Taxonomy, and Create Output Tables
> The following is run in an R environment.  

These example commands as written assume paired-end data, with notes included on what would be different if working with single-end data. The taxonomy reference database used below is an example only, suitable for the example 16S dataset ([GLDS-200](https://osdr.nasa.gov/bio/repo/data/studies/OSD-200)) used here. Other taxonomy references databases designed for DECIPHER can be found here: [https://www2.decipher.codes/data/Downloads/TrainingSets/](https://www2.decipher.codes/data/Downloads/TrainingSets/)  

<br>

### 5a. Learning the Error Rates
```R
## Forward error rates ##
forward_errors <- learnErrors(fls="sample1_<tech_type>_GLAmpSeq_R1_filtered.fastq.gz", multithread=TRUE)

## Reverse error rates (skip if single-end data) ##
reverse_errors <- learnErrors(fls="sample1_<tech_type>_GLAmpSeq_R2_filtered.fastq.gz", multithread=TRUE)
```

**Parameter Definitions:**  

*	`learnErrors()` – the DADA2 function we are calling, with the following parameters set within it
*	`fls=` – specifies the path to the filtered reads (either forward or reverse)
*	`multithread=TRUE` – determine number of cores available and run in parallel when possible (can also take an integer specifying the number of cores to use)

**Input Data:**

* \*filtered.fastq.gz (filtered reads, output from [Step 3](#3-quality-filtering))

**Output Data:**

* `forward_errors` (a named list containing a numeric matrix with the forward error rates)
* `reverse_errors` (a named list containing a numeric matrix with the reverse error rates (only for paired-end data))

<br>

### 5b. Inferring Sequences
```R
## Inferring forward sequences ##
forward_seqs <- dada(derep="sample1_<tech_type>_GLAmpSeq_R1_filtered.fastq.gz", err=forward_errors, pool="pseudo", multithread=TRUE)

## Inferring reverse sequences (skip if single-end)##
reverse_seqs <- dada(derep="sample1_<tech_type>_GLAmpSeq_R2_filtered.fastq.gz", err=reverse_errors, pool="pseudo", multithread=TRUE)
```

**Parameter Definitions:**  

* `dada()` – the DADA2 function we are calling, with the following parameters set within it
* `derep=` – the path to the filtered reads (either forward or reverse)
* `err=` – the object holding the error profile for the inferred reads (either forward or reverse)
* `pool="pseudo"` – setting the method of incorporating information from multiple samples, "pseudo" instructs the algorithm to perform pseudo-pooling between individually processed samples
* `multithread=TRUE` – determine number of cores available and run in parallel when possible (can also take an integer specifying the number of cores to use)

**Input Data:**

* \*filtered.fastq.gz (filtered reads, output from [Step 3](#3-quality-filtering))
* `forward_errors` (a named list containing a numeric matrix with the forward error rates, output from [Step 5a](#5a-learning-the-error-rates))
* `reverse_errors` (a named list containing a numeric matrix with the reverse error rates, output from [Step 5a](#5a-learning-the-error-rates) (only for paired-end))

**Output Data:**

* `forward_seqs` (a dada-class object containing the forward-read inferred sequences)
* `reverse_seqs` (a dada-class object containing the reverse-read inferred sequences (only for paired-end))

<br>

### 5c. Merging Forward and Reverse Reads; Skip if Data are Single-End
```R
merged_contigs <- mergePairs(dadaF=forward_seqs, derepF="sample1_<tech_type>_GLAmpSeq_R1_filtered.fastq.gz", dadaR=reverse_seqs, derepR="sample1_<tech_type>_GLAmpSeq_R2_filtered.fastq.gz")
```

**Parameter Definitions:** 

* `merged_contigs <-` – specifies the variable that will store the results within in our R environment
* `mergePairs()` – the DADA2 function we are calling, with the following parameters set within it
* `dadaF=` – specifies the object holding the forward-read inferred sequences
* `derepF=` – specifies the path to the filtered forward reads
* `dadaR=` – specifies the object holding the reverse-read inferred sequences
* `derepR=` – specifies the path to the filtered reverse reads

**Input Data:**

* \*filtered.fastq.gz (filtered reads, output from [Step 3](#3-quality-filtering))
* `forward_seqs` (a dada-class object containing forward-read inferred sequences, output from [Step 5b](#5b-inferring-sequences))
* `reverse_seqs` (a dada-class object containing reverse-read inferred sequences, output from [Step 5b](#5b-inferring-sequences))

**Output Data:**

* `merged_contigs` (a dataframe containing the merged contigs)

<br>

### 5d. Generating Sequence Table with Counts per Sample
```R
seqtab <- makeSequenceTable(merged_contigs)
```

**Parameter Definitions:**  

* `seqtab <-` - specifies the variable that will store the results within our R environment
* `makeSequenceTable()` - the DADA2 function we call calling, with either `merged_contigs` for paired-end data (as in this example) or `forward_seqs` for single-end data as input.

**Input Data:**

* `merged_contigs` or `forward_seqs` (a variable containing the merged contigs for paired-end data, output from [Step 5c](#5c-merging-forward-and-reverse-reads-skip-if-data-are-single-end), or the forward-read inferred sequences for single-end data, output from [Step 5b](#5b-inferring-sequences))

**Output Data:**

* `seqtab` (a named integer matrix containing the sequence table)

<br>

### 5e. Removing putative chimeras
```R
seqtab.nochim <- removeBimeraDenovo(unqs=seqtab, method="consensus", multithread=TRUE)
```

**Parameter Definitions:**  

* `seqtab.nochim <-` – specifies the variable that will store the results within in our R environment
* `removeBimeraDenovo()` – the DADA2 function we are calling, with the following parameters set within it
* `unqs=` – specifying the `seqtab` object created above
* `method=` – specifying the method for putative-chimera identification and removal, "consensus" instructs the function to check the samples in the sequence table independently for bimeras and make a consensus decision on each sequence variant 
* `multithread=TRUE` – determine number of cores available and run in parallel when possible (can also take an integer specifying the number to run)

**Input Data:**

* `seqtab` (a named integer matrix containing the sequence table, output from [Step 5d](#5d-generating-sequence-table-with-counts-per-sample))

**Output Data:**

* `seqtab.nochim` (a named integer matrix containing the sequence table filtered to exclude putative chimeras)

<br>

### 5f. Assigning Taxonomy

```R
## Creating a DNAStringSet object from the ASVs: ##
dna <- DNAStringSet(getSequences(seqtab.nochim))

## Downloading the reference R taxonomy object: ##
download.file(url = "https://figshare.com/ndownloader/files/52846199", 
             destfile = "SILVA_SSU_r138_2_2024.RData", 
             method = "libcurl", 
             headers = c("User-Agent" = "Mozilla/5.0"))

## Loading taxonomy object: ##
load("SILVA_SSU_r138_2_2024.RData")

## Classifying sequences:
tax_info <- IdTaxa(test=dna, trainingSet=trainingSet, strand="both", processors=NULL)
```

**Parameter Definitions:** 

- `download.file()` - the R utils function used to download the taxonomy database file
  - `url=` - reference database URL address to download
  - `destfile=` - local path/name for the downloaded file
  - `method=` - specifies the download method to use
  - `headers=` - HTTP headers to pass with the download request
- `IdTaxa()` - the DECIPHER function used to classify the sequences 
  - `test=dna` - DNAStringSet object holding sequences to classify
  - `trainingSet=trainingSet` - specifies the reference database to use
  - `strand="both"` - specifies to check taxonomy assignment in both orientations
  - `processors=NULL` - specifies the number of processors to use, `NULL` indicates to use all available cores or an integer may be provided to manually specify the number to use

**Input Data:**

* `seqtab.nochim` (a named integer matrix containing the filtered sequence table, output from [Step 5e](#5e-removing-putative-chimeras))
* `trainingSet` (a variable provided in the RData object containing the reference database, SILVA_SSU_r138_2_2024.RData)

**Output Data:**

* `tax_info` (the DECIPHER Taxa object containing assigned taxons)

<br>

### 5g. Generating and Writing Standard Outputs

```R
## Giving sequences more manageable names (e.g. ASV_1, ASV_2, …,): ##
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

## Making then writing a fasta of final ASV seqs: ##
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs_<tech_type>_GLAmpSeq.fasta")

## Making then writing a count table: ##
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)

write.table(asv_tab, "counts_<tech_type>_GLAmpSeq.tsv", sep="\t", quote=F, col.names=NA)

## Creating table of taxonomy and setting any that are unclassified as "NA": ##
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
tax_tab <- t(sapply(tax_info, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(tax_tab) <- ranks
rownames(tax_tab) <- gsub(pattern=">", replacement="", x=asv_headers)

write.table(tax_tab, "taxonomy_<tech_type>_GLAmpSeq.tsv", sep = "\t", quote=F, col.names=NA)

## Generating then writing biom file format: ##
biom_object <- make_biom(data=asv_tab, observation_metadata=tax_tab)
write_biom(biom_object, "taxonomy-and-counts_<tech_type>_GLAmpSeq.biom")

## Making a combined taxonomy and count table ##
tax_and_count_tab <- merge(tax_tab, asv_tab)
write.table(tax_and_count_tab, "taxonomy-and-counts_<tech_type>_GLAmpSeq.tsv", sep="\t", quote=FALSE, row.names=FALSE)
```
```bash
zip -j -q taxonomy-and-counts_<tech_type>_GLAmpSeq.biom.zip taxonomy-and-counts_<tech_type>_GLAmpSeq.biom
```

**Parameter Definitions:**
**zip**
- `-j` - junk (don't record) directory names
- `-q` – quiet operation
- `taxonomy-and-counts_<tech_type>_GLAmpSeq.biom.zip` – positional argument naming the zip output file
- `taxonomy-and-counts_<tech_type>_GLAmpSeq.biom` – positional argument naming the input folder to package

**Input Data:**

* `seqtab.nochim` (a named integer matrix containing the filtered sequence table, output from [Step 5e](#5e-removing-putative-chimeras))
* `tax_info` (the DECIPHER Taxa object containing assigned taxons, output from [Step 5f](#5f-assigning-taxonomy))

**Output Data:**

* **ASVs_<tech_type>_GLAmpSeq.fasta** (a fasta file containing the inferred sequences)
* **counts_<tech_type>_GLAmpSeq.tsv** (a tab-separated file containing the sample feature count table)
* **taxonomy_<tech_type>_GLAmpSeq.tsv** (a tab-separated file containing the taxonomy table)
* **taxonomy-and-counts_<tech_type>_GLAmpSeq.tsv** (a tab-separated file containing the combined taxonomy and count table)
* **taxonomy-and-counts_<tech_type>_GLAmpSeq.biom.zip** (a zip package containing the biom-formatted file)
  * **taxonomy-and-counts_<tech_type>_GLAmpSeq.biom** (a biom-formatted file containing the count and taxonomy table)
* **read-count-tracking_<tech_type>_GLAmpSeq.tsv** (a tab-separated file containing the read counts at each processing step)

<br>

---

## 6. Amplicon Seq Data Analysis Set Up  
  
<br>

### 6a. Create Sample Runsheet

> Note: Rather than running the command below to create the runsheet needed for processing, the runsheet may also be created manually by following the examples for [Paired-end](https://github.com/nasa/GeneLab_AmpliconSeq_Workflow/blob/main/examples/runsheet/PE_file.csv) and [Single-end](https://github.com/nasa/GeneLab_AmpliconSeq_Workflow/blob/main/examples/runsheet/SE_file.csv) samples. When creating this table manually, the most important columns for the analyses below are:

* `sample_id` - column with unique sample names.
* `groups`    - column with the groups/treatments that each sample belong to. This column is used for comparison.

```bash
### Download the *ISA.zip file from the OSDR ###
dpt-get-isa-archive \
 --accession OSD-###

### Parse the metadata from the *ISA.zip file to create a sample runsheet ###
dpt-isa-to-runsheet --accession OSD-### \
 --config-type amplicon \
 --config-version Latest \
 --isa-archive *ISA.zip
```

**Parameter Definitions:**

* `--accession OSD-###` - OSD accession ID (replace ### with the OSD number being processed), used to retrieve the urls for the ISA archive and raw reads hosted on the OSDR
* `--config-type` - instructs the script to extract the metadata required for Amplicon Sequencing data processing from the ISA archive
* `--config-version` - specifies the `dp-tools` configuration version to use, a value of `Latest` will specify the most recent version
* `--isa-archive` - specifies the *ISA.zip file for the respective OSD dataset, downloaded in the `dpt-get-isa-archive` command

**Input Data:**

* *No input data required, but the OSD accession ID needs to be indicated, which is used to download the respective ISA archive*

**Output Data:**

* *ISA.zip (compressed ISA directory containing Investigation, Study, and Assay (ISA) metadata files for the respective OSD dataset, used to define sample groups - the *ISA.zip file is located in the [OSDR](https://osdr.nasa.gov/bio/repo/) under 'Files' -> 'Study Metadata Files')

* **{OSD-Accession-ID}_amplicon_v{version}_runsheet.csv** (a comma-separated sample metadata file containing sample group information, version denotes the dp_tools schema used to specify the metadata to extract from the ISA archive)
    > NOTE: if there are multiple valid Amplicon Sequencing assays in the dataset, then multiple runsheets will be generated (1 for each assay). The runsheet filenames will also include the value from the "Parameter Value[Library Selection]" column before the runsheet version ({OSD-Accession-ID}\_amplicon_{library_selection}_v{version}_runsheet.csv). For example, for OSD-268, which has both "16S" and "ITS" assays, two files are generated: OSD-268_amplicon_16S_v2_runsheet.csv and OSD-268_amplicon_ITS_v2_runsheet.csv.

<br>

> The remainder of this document is performed in R.  


### 6b. R Environment Set Up


#### 6b.i. Load Libraries

```R
library(vegan)
library(phyloseq)
library(glue)
library(FSA)
library(multcompView)
library(rstatix)
library(patchwork)
library(RColorBrewer)
library(DESeq2)
library(ggdendro)
library(broom)
library(ggrepel)
library(tools)
library(ANCOMBC)
library(DescTools)
library(taxize)
library(mia)
library(utils)
library(scales)
library(tidyverse)
library(vsn)
library(hexbin)
```

#### 6b.ii. Define Custom Functions

#### calculate_text_size()
<details>
  <summary>calculates text size for plotting based on number of samples and minimum text size</summary>

  ```R
  calculate_text_size <- function(num_samples, start_samples = 25, min_size = 3) {
    max_size <- 11  # Maximum size for up to start_samples
    slope <- -0.15
    
    if (num_samples <= start_samples) {
      return(max_size)
    } else {
      # Calculate the current size with the hard coded slope
      current_size = max_size + slope * (num_samples - start_samples)
      
      # Ensure the size doesn't go below the minimum
      return(max(current_size, min_size))
    }
  }
  ```

  **Function Parameter Definitions:**
  - `num_samples=` - the number of samples to plot
  - `start_samples=25` - the number of samples to start with, used to specify text size
  - `min_size=3` - the minimum text size for plotting  

  **Returns:** an integer representing the calculated text size
</details>

#### expandy()
<details>
  <summary>wrapper around `ggplot2:expand_limits` that expands a plot's y-limit based on a specified set of y-values</summary>

  ```R
  expandy <- function(vec, ymin=NULL) {
    # vec [NUMERIC] - a numeric vector of y values.

    max.val <- max(vec, na.rm=TRUE) + 0.1
    
    expand_limits(y=c(ymin, max.val))
  }
  ```

  **Function Parameter Definitions:**
  - `vec=` - a numeric vector containing y-values
  - `ymin=` - the minimum y-limit

</details>

#### transform_phyloseq()
<details>
  <summary>create a phyloseq object with the appropriate sample count transformation depending on the supplied transformation method ('rarefy' or 'vst')</summary>

  ```R
  transform_phyloseq <- function( feature_table, metadata, method, rarefaction_depth=500){

    # Rarefaction
    if(method == 'rarefy'){
      # Create phyloseq object
      ASV_physeq <- phyloseq(otu_table(feature_table, taxa_are_rows = TRUE),
                            sample_data(metadata))

      # Get the count for every sample sorted in ascending order
      seq_per_sample <- colSums(feature_table) %>% sort()
      # Minimum sequences/count value
      depth <- min(seq_per_sample)

      # Error if the number of sequences per sample left after filtering is 
    # insufficient for diversity analysis
    if(max(seq_per_sample) < 100){
      
      warning_file <- glue("{beta_diversity_out_dir}{output_prefix}beta_diversity_failure{assay_suffix}.txt")
      writeLines(
        text = glue("The maximum sequence count per sample ({max(seq_per_sample)}) is less than 100.
Therefore, beta diversity analysis with rarefaction cannot be performed. Check VST method normalization instead."),
        con = warning_file
      )
      return(NULL)   # stop rarefaction branch, but don't kill script
    }

      # Loop through the sequences per sample and return the count
      # nearest to the minimum required rarefaction depth
      for (count in seq_per_sample) {
        # Get the count equal to rarefaction_depth or nearest to it
        if(count >= rarefaction_depth) {
          depth <- count
          break
        }
      }

      # Error if the depth that ends up being used is also less than 100
    if(depth < 100){
      
      warning_file <- glue("{beta_diversity_out_dir}{output_prefix}beta_diversity_failure{assay_suffix}.txt")
      writeLines(
        text = glue("The rarefaction depth being used in the analysis ({depth}) is less than 100.
Therefore, beta diversity analysis with rarefaction cannot be performed. Check VST method normalization instead."),
        con = warning_file
      )
      return(NULL)   # stop rarefaction branch, but don't kill script
    } 
    
    #Warning if rarefaction depth is between 100 and 500
    if (depth > 100 && depth < 500) {
      warning(glue("Rarefaction depth ({depth}) is between 100 and 500.
Beta diversity results may be unreliable."))
    }

      #----- Rarefy sample counts to even depth per sample
      ps <- rarefy_even_depth(physeq = ASV_physeq,
                              sample.size = depth,
                              rngseed = 1,
                              replace = FALSE,
                              verbose = FALSE)

    # ---- Group check ----
    survived_samples <- sample_names(ps)
    remaining_groups <- unique(metadata[rownames(metadata) %in% survived_samples, groups_colname])
    
    if(length(remaining_groups) < 2){
      warning_file <- glue("{beta_diversity_out_dir}{output_prefix}beta_diversity_failure{assay_suffix}.txt")
      writeLines(
        text = glue("Not enough groups remain after rarefaction at {depth} (only {length(remaining_groups)}). Skipping beta diversity with rarefaction."),
        con = warning_file
      )
      return(NULL)  # stop analysis, like depth failure
    }

    # Write rarefaction depth used into file
    depth_file <- glue("{beta_diversity_out_dir}{output_prefix}rarefaction_depth{assay_suffix}.txt")
    writeLines(
      text = as.character(depth),
      con = depth_file
    )
    
    # Variance Stabilizing Transformation
    }else if(method == "vst"){

      # Using deseq
      # Keep only ASVs with at least 1 count
      feature_table <- feature_table[rowSums(feature_table) > 0, ]
      # Add +1 pseudocount for VST for vst transformation
      feature_table <- feature_table + 1

      # Make the order of samples in metadata match the order in feature table
      metadata <- metadata[colnames(feature_table),]

      # Create VST normalized counts matrix
      # ~1 means no design
      deseq_counts <- DESeqDataSetFromMatrix(countData = feature_table,
                                            colData = metadata,
                                            design = ~1)
      deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
      vst_trans_count_tab <- assay(deseq_counts_vst)

      # Making a phyloseq object with our transformed table
      vst_count_phy <- otu_table(object = vst_trans_count_tab, taxa_are_rows = TRUE)
      sample_info_tab_phy <- sample_data(metadata)
      ps <- phyloseq(vst_count_phy, sample_info_tab_phy)
    }else{
      stop("Please supply a valid normalization method, either 'rarefy' or 'vst' ")
    }
    return(ps)
  }
  ```
  **Function Parameter Definitions:**
  - `feature_table=` - a dataframe containing feature/ASV counts with samples as columns and features as rows
  - `metadata=` - a dataframe containing sample metadata with samples as row names and sample info as columns
  - `method=` - a string specifying the transformation to use: either "rarefy" (rarefaction) or "vst" (variance stabilizing transformation)
  - `rarefaction_depth=500` - the minimum number of reads to simulate during rarefaction
  
  **Returns:** a phyloseq object created from the feature_table, metadata, and specified transformation method
</details>

#### make_dendrogram()
<details>
  <summary>compute hierarchical clustering and create a dendrogram</summary>

  ```R
  make_dendrogram <- function(dist_obj, metadata, groups_colname,
                            group_colors, legend_title){

    # Hierarchical Clustering
    sample_clust <- hclust(d = dist_obj, method = "ward.D2")
    
    # Extract clustering data for plotting
    hcdata <- dendro_data(sample_clust, type = "rectangle")
    segment_data <- segment(hcdata) # sepcifications for tree structure
    label_data <- label(hcdata) %>%
      left_join(metadata %>% 
                  rownames_to_column("label")) # Labels are sample names

    # Plot dendrogram
    dendrogram <- ggplot() +
      # Plot tree
      geom_segment(data = segment_data, 
                  aes(x = x, y = y, xend = xend, yend = yend) 
      ) +
      # Add sample text labels to tree
      geom_text(data = label_data , 
                aes(x = x, y = y, label = label, 
                    color = !!sym(groups_colname) , hjust = 0), 
                size = 4.5, key_glyph = "rect") +
      scale_color_manual(values = group_colors) +
      coord_flip() +
      scale_y_reverse(expand = c(0.2, 0)) +
      labs(color = legend_title) +
      theme_dendro() +
      guides(colour = guide_legend(override.aes = list(size = 5)))+
      theme(legend.key = element_rect(fill=NA),
            text = element_text(face = 'bold'),
            legend.title = element_text(size = 12, face='bold'),
            legend.text = element_text(face = 'bold', size = 11))
    
    return(dendrogram)
  }
  ```
  **Function Parameter Definitions:**
  - `dist_obj=` - a distance object of class 'dist' holding the calculated distance (Euclidean, Bray-Curtis, etc.) between samples
  - `metadata=` - a dataframe containing sample metadata with samples as row names and sample info as columns
  - `groups_colname=` - name of the column in the metadata dataframe to use for specifying sample groups
  - `legend_title=` - legend title to use for plotting

  **Returns:** a dendrogram plot
</details>

#### run_stats()
<details>
  <summary>run variance and adonis (analysis of variance using distance matrices) tests</summary>

  ```R
  run_stats <- function(dist_obj, metadata, groups_colname){

    # Retrieve sample names from the dist object
    samples <- attr(dist_obj, "Label")
    # subset metadata to contain ony samples in the dist_obj
    metadata <- metadata[samples,]

    # Run variance test and present the result in a nicely formatted table / dataframe
    variance_test <- betadisper(d = dist_obj,
                                group = metadata[[groups_colname]]) %>%
      anova() %>% # MAke results anova-like
      broom::tidy() %>%  # make the table 'tidy'
      mutate(across(where(is.numeric), ~round(.x, digits = 2))) # round-up numeric columns

    # Run Adonis test
    adonis_res <- adonis2(formula = dist_obj ~ metadata[[groups_colname]])

    adonis_test <- adonis_res %>%
      broom::tidy() %>% # Make tidy table
      mutate(across(where(is.numeric), ~round(.x, digits = 2))) # round-up numeric columns

    # Return a named list with the variance and adonis test results
    return(list(variance = variance_test, adonis = adonis_test))
  }
  ```
  **Function Parameter Definitions:**
  - `dist_obj=` - distance object of class 'dist' holding the calculated distance (Euclidean, Bray-Curtis, etc.) between samples
  - `metadata=` - a dataframe containing sample metadata with samples as row names and sample info as columns
  - `groups_colname=` - string specifying the name of the column in the metadata dataframe to use for specifying sample groups

  **Returns:** a named list containing the variance and adonis test results as dataframes
</details>

#### plot_pcoa()
<details>
  <summary>generate a Principle Coordinate Analysis (PCoA) plot using phyloseq::ordinate function</summary>

  ```R
  plot_pcoa <- function(ps, stats_res, distance_method,
                        groups_colname, group_colors, legend_title,
                        addtext=FALSE) {

    # Generating a PCoA with phyloseq
    pcoa <- ordinate(physeq = ps, method = "PCoA", distance = distance_method)
    eigen_vals <- pcoa$values$Eigenvalues

    # Calculate the percentage of variance
    percent_variance <- eigen_vals / sum(eigen_vals) * 100

    # Retrieving plot labels
    r2_value <- stats_res$adonis[["R2"]][1]
    prf_value <- stats_res$adonis[["p.value"]][1]
    label_PC1 <- sprintf("PC1 [%.1f%%]", percent_variance[1])
    label_PC2 <- sprintf("PC2 [%.1f%%]", percent_variance[2])

    # Retrieving pcoa vectors
    vectors_df <- pcoa$vectors %>%
                    as.data.frame() %>%
                    rownames_to_column("samples")
    # Creating a dataframe for plotting
    plot_df <- sample_data(ps) %>%
                as.matrix() %>%
                as.data.frame() %>%
                rownames_to_column("samples") %>%
                select(samples, !!groups_colname) %>%
                right_join(vectors_df, join_by("samples"))
    # Plot pcoa
    p <- ggplot(plot_df, aes(x=Axis.1, y=Axis.2,
                            color=!!sym(groups_colname),
                            label=samples)) +
      geom_point(size=1)

    # Add text
    if(addtext){
      p <- p + geom_text(show.legend = FALSE,
                        hjust = 0.3, vjust = -0.4, size = 4)
    }

  # Add annotations to pcoa plot
  p <-  p +  labs(x = label_PC1, y = label_PC2, color = legend_title) +
      coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) +
      scale_color_manual(values = group_colors) +
      theme_bw() + theme(text = element_text(size = 15, face="bold"),
                        legend.direction = "vertical",
                        legend.justification = "center",
                        legend.title = element_text(hjust=0.1)) +
      annotate("text", x = Inf, y = -Inf,
              label = paste("R2:", toString(round(r2_value, 3))),
              hjust = 1.1, vjust = -2, size = 4)+
      annotate("text", x = Inf, y = -Inf,
              label = paste("Pr(>F)", toString(round(prf_value,4))),
              hjust = 1.1, vjust = -0.5, size = 4) + ggtitle("PCoA")

    return(p)
  }
  ```

  **Function Parameter Definitions:**
  - `ps=` - phyloseq object constructed from feature, taxonomy, and metadata tables
  - `stats_res=` - named list containing the variance and adonis test results as dataframes, generated using [run_stats()](#run_stats)
  - `distance_method=` - string specifying the method used to calculate the distance between samples; values can be "euclidean" (Euclidean distance) or "bray" (Bray-Curtis dissimilarity)
  - `groups_colname=` - string specifying the name of the column in the metadata dataframe to use for specifying sample groups
  - `group_colors=` - named character vector of colors for each group in `groups_colname`
  - `legend_title=` - string specifying the legend title to use for plotting
  - `addtext=FALSE` - boolean value specifying if the sample labels should be added to the pcoa plot

  **Returns:** a PCoA plot
</details>

#### remove_rare_features()
<details>
  <summary>filter out rare features from a feature table by occurrence in a fraction of samples depending on the supplied cut-off</summary>

  ```R
  remove_rare_features <- function(feature_table, cut_off_percent=0.75){
    
    # Filter by occurrence in a fraction of samples
    # Define a cut-off for determining what's rare
    cut_off <- cut_off_percent * ncol(feature_table)
    # Get the occurrence for each feature
    feature_occurrence <- rowSums(feature_table > 0)
    # Get names of the abundant features
    abund_features <- names(feature_occurrence[feature_occurrence >= cut_off])
    # Remove rare features
    abun_features.m <- feature_table[abund_features,]
    return(abun_features.m)
  }
  ```
**Function Parameter Definitions:**
  - `feature_table=` - dataframe containing feature/ASV counts with samples as columns and features as rows
  - `cut_off_percent=0.75` - decimal value between 0.001 and 1 specifying the fraction of the total number of samples to determine the most abundant features; by default it removes features that are not present in 3/4 of the total number of samples

  **Returns:** a dataframe of feature/ASV counts filtered to include only features passing the specified threshold
</details>

#### process_taxonomy()
<details>
  <summary>reformat taxonomy table to regularize naming/formatting. Removes extraneous prefixes/suffixes from taxonomy names; replaces empty values with "Other".</summary>

  ```R
  process_taxonomy <- function(taxonomy, prefix='\\w__') {
    
    # Ensure that all columns are of character data type
    taxonomy <- apply(X = taxonomy, MARGIN = 2, FUN = as.character) 
    
    # Loop over every column (rank i.e. domain to species) amd make the necessary edits
    for (rank in colnames(taxonomy)) {
      # Delete the taxonomy prefix
      taxonomy[,rank] <- gsub(pattern = prefix, x = taxonomy[, rank],
                              replacement = '')
      # Delete _number at the end of taxonomy names inserted by the new version of DECIPHER
      taxonomy[,rank] <- gsub(pattern ="_[0-9]+$", x = taxonomy[, rank], replacement = '')

      indices <- which(is.na(taxonomy[,rank]))
      taxonomy[indices, rank] <- rep(x = "Other", times=length(indices)) 
      # replace empty cell with the string 'Other'
      indices <- which(taxonomy[,rank] == "")
      taxonomy[indices,rank] <- rep(x = "Other", times=length(indices))
    }
    # Replace _ with space
    taxonomy <- apply(X = taxonomy,MARGIN = 2,
                      FUN = gsub,pattern = "_",replacement = " ") %>% 
      as.data.frame(stringAsfactor=F)
    return(taxonomy)
  }
  ```
  **Function Parameter Definitions:**
  - `taxonomy=` - dataframe of ASV taxonomy assignments to be processed
  - `prefix='\\w__'` - a regular expression specifying the characters to remove from taxon names; use '\\w__' for greengenes and 'D_\\d__' for SILVA

  **Returns:** a dataframe containing the reformatted ASV taxonomy assignments
</details>

#### format_taxonomy_table()
<details>
  <summary>reformat taxonomy table by appending a suffix to a known name in the previous cell</summary>

  ```R
  format_taxonomy_table <- function(taxonomy, stringToReplace="Other", suffix=";Other") {

    for (taxa_index in seq_along(taxonomy)) {
      
      indices <- grep(x = taxonomy[,taxa_index], pattern = stringToReplace)
      
      taxonomy[indices,taxa_index] <- 
        paste0(taxonomy[indices,taxa_index-1],
              rep(x = suffix, times=length(indices)))
      
    }
    return(taxonomy)
  }
  ```
  **Function Parameter Definitions:**
  - `taxonomy=` - dataframe of ASV taxonomy assignments to be processed
  - `stringToReplace="Other"` - specifies the string to replace, "Other" is used by default
  - `suffix=";Other"` - specifies the replacement string, ";Other" is used by default

  **Returns:** a dataframe containing the reformatted ASV taxonomy assignments
</details>

#### fix_names()
<details>
  <summary>reformat taxonomy table by appending a set of suffixes to a set of known names</summary>

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
  **Custom Functions Used:**
  - [format_taxonomy_table()](#format_taxonomy_table)

  **Function Parameter Definitions:**
  - `taxonomy=` - dataframe of ASV taxonomy assignments to be processed
  - `stringToReplace=` - a vector of regex strings specifying what to replace in the taxonomy dataframe
  - `suffix=` - a vector of regex strings specifying the replacement strings to use

  **Returns:** a dataframe containing the reformatted ASV taxonomy assignments
</details>

#### make_feature_table()
<details>
  <summary>generate taxon level count matrix based on a taxonomy table and an existing feature table, only retains taxa found in at least one sample</summary>

  ```R    
  make_feature_table <- function(count_matrix,taxonomy,
                                taxon_level, samples2keep=NULL){

    feature_counts_df <- data.frame(taxon_level=taxonomy[,taxon_level],
                                    count_matrix, check.names = FALSE,
                                    stringsAsFactors = FALSE)

    feature_counts_df <- aggregate(.~taxon_level,data = feature_counts_df,
                                  FUN = sum)
    rownames(feature_counts_df) <- feature_counts_df[,"taxon_level"]
    feature_table <- feature_counts_df[,-1]
    # Retain only taxa found in at least one sample
    taxa2keep <- rowSums(feature_table) > 0
    feature_table <- feature_table[taxa2keep,]

    if(!is.null(samples2keep)){
      feature_table <- feature_table[,samples2keep]
      # Retain only taxa found in at least one sample
      taxa2keep <- rowSums(feature_table) > 0
      feature_table <- feature_table[taxa2keep,]
    }
    return(feature_table)
  }
  ```
  **Function Parameter Definitions:**
  - `count_matrix=` - a dataframe containing ASVs or OTUs and their respective counts
  - `taxonomy=` - dataframe of ASV taxonomy assignments to be processed
  - `taxon_level=` - string defining the taxon levels, i.e. domain to species
  - `samples2keep=NULL` - a list of sample strings specifying the samples to keep; default value of NULL keeps all samples
  
  **Returns:** a dataframe containing a taxon level count matrix filtered for taxa found in at least one sample
</details>

#### group_low_abund_taxa()
<details>
  <summary>group rare taxa or return a table with the rare taxa</summary>

  ```R
  group_low_abund_taxa <- function(abund_table, threshold=0.05,
                                  rare_taxa=FALSE) {

    # Initialize an empty vector that will contain the indices for the
    # low abundance columns/ taxa to group
    taxa_to_group <- c()
    #initialize the index variable of species with low abundance (taxa/columns)
    index <- 1

    # Loop over every column or taxa then check to see if the max abundance is less than the set threshold
    # if true, save the index in the taxa_to_group vector variable

    for (column in seq_along(abund_table)){
      if(max(abund_table[, column], na.rm = TRUE) < threshold ){
        taxa_to_group[index] <- column
        index = index + 1
      }
    }
    if(is.null(taxa_to_group)){
      message(glue::glue("Rare taxa were not grouped. please provide a higher threshold than {threshold} for grouping rare taxa, only numbers are allowed."))
    return(abund_table)
    }

    if(rare_taxa){
      abund_table <- abund_table[, taxa_to_group, drop=FALSE]
    }else{
      #remove the low abundant taxa or columns
      abundant_taxa <-abund_table[, -(taxa_to_group), drop=FALSE]
      #get the rare taxa
      # rare_taxa <-abund_table[, taxa_to_group]
      rare_taxa <- subset(x = abund_table, select = taxa_to_group)
      #get the proportion of each sample that makes up the rare taxa
      rare <- rowSums(rare_taxa)
      #bind the abundant taxa to the rae taxa
      abund_table <- cbind(abundant_taxa, rare)
      #rename the columns i.e the taxa
      colnames(abund_table) <- c(colnames(abundant_taxa), "Rare")
    }

    return(abund_table)
  }
  ```
  **Function Parameter Definitions:**
  - `abund_table=` - relative abundance matrix with taxa as columns and samples as rows
  - `threshold=0.05` - a number between 0.001 and 1 specifying the threshold for grouping rare taxa
  - `rare_taxa=FALSE` - boolean specifying if only rare taxa should be returned, if set to TRUE then a table with only the rare taxa will be returned
  
  **Returns:** a dataframe containing a relative abundance matrix with taxa as columns and samples as rows
</details>

#### collapse_samples()
<details>
  <summary>collapse samples in a feature table with a user-defined function based on group in metadata</summary>

  ```R
  collapse_samples <- function(taxon_table, metadata, group, fun=sum, convertToRelativeAbundance=FALSE){

    common.ids <- intersect(rownames(taxon_table), rownames(metadata))
    metadata <- droplevels(metadata[common.ids, , drop=FALSE])
    taxon_table <- taxon_table[common.ids, , drop=FALSE]
    taxon_table <- cbind(subset(x = metadata, select=group), taxon_table)

    taxon_table <- aggregate(reformulate(termlabels = group, response = '.'),
                            data = taxon_table, FUN = fun)
    rownames(taxon_table) <- taxon_table[, 1]
    taxon_table <- taxon_table[,-1]
    if(convertToRelativeAbundance){
      taxon_table <- t(apply(X = taxon_table, MARGIN = 1, FUN = function(x) x/sum(x)))
    }

    final <- list(taxon_table,metadata)
    names(final) <- c("taxon_table","metadata")
    return(final)
  }
  ```
  **Function Parameter Definitions:**
  - `taxon_table=` - a dataframe of feature counts with samples as rows and features (e.g. ASVs or OTUs) as columns
  - `metadata=` - a dataframe containing sample metadata with samples as row names and sample info as columns
  - `group=` - a string specifying the column in the metadata dataframe containing the sample groups that will be used to collapse the samples
  - `fun=sum` - a string specifying the R function to apply when collapsing the samples; 'sum' is used by default
  - `convertToRelativeAbundance=FALSE` - a boolean specifying whether or not the value in the taxon table should be converted to per sample relative abundance values
  
  **Returns:** a named list containing two dataframes: `taxon_table` and `metadata`
    - a dataframe of aggregated feature counts by group
    - a dataframe containing group specific metadata for the aggregated feature count
</details>

#### get_ncbi_ids()
<details>
  <summary>retrieve NCBI taxonomy id for a given taxonomy name using taxize</summary>

  ```R
  taxize_options(ncbi_sleep = 0.8)
  get_ncbi_ids <- function(taxonomy, target_region){
    
    if(target_region == "ITS"){
      search_string <- "fungi"
    }else if(target_region == "18S"){
      search_string <- "eukaryote"
    }else{
      search_string <- "bacteria"
    }
    
    uid <- get_uid(taxonomy, division_filter = search_string)
    tax_ids <- uid[1:length(uid)]
    return(tax_ids)
    
  }
  ```
  **Function Parameter Definitions:**
  - `taxonomy=` - string specifying the taxonomy name that will be used to search for the respective NCBI ID
  - `target_region=` - amplicon target region to analyze; options are "16S", "18S", or "ITS"

  **Returns:** an integer of NCBI taxonomic identifiers; if none is found, returns NA
</details>

#### find_bad_taxa()
<details>
  <summary>error handling function for ANCOMBC::ancombc2</summary>

  ```R
  find_bad_taxa <- function(cnd){
    split_res <- strsplit(conditionMessage(cnd), "\n")

    if(split_res == "replacement has 0 rows, data has 1" || 
      split_res == "All taxa contain structural zeros") { 
      
      return(
        list(res=data.frame(taxon=split_res, lfc=NA, se=NA,
                            W=NA, p=NA, q=NA, diff=NA, pass_ss=NA))
      )
    }
    
    bad_taxa <- split_res[[c(1L, 2L)]]
    bad_taxa <- .subset2(strsplit(bad_taxa, ", "), 1L)
    return(bad_taxa)
  }
  ```
  **Function Parameter Definitions:**
  - `cnd=` - specifies the error condition to catch when running the ANCOMBC::ancombc2 function

  **Returns:** a list containing one empty dataframe named 'res' with the same format as the ANCOMBC::ancombc2 primary result
</details>

#### ancombc2()
<details>
  <summary>wrapper around ANCOMBC::ancombc2() function that adds error handling</summary>

  ```R
  ancombc2 <- function(data, ...) {

    tryCatch(
      ANCOMBC::ancombc2(data = data, ...),
      error = function(cnd) {
              
        res <- find_bad_taxa(cnd)
        if( is.data.frame(res[[1]]) ){
          # Returns a manually created empty data.frame
          return(res)
        }else{
          # Returns the names of the bad taxa to exclude from further analysis
          bad_taxa <- res # renaming for readability
        }
        
        # Second error catcher in case it fails the first one 
        tryCatch(
          ANCOMBC::ancombc2(data = data[!rownames(data) %in% bad_taxa, ], ...),
          
          error = function(cnd) {
            # Returns a manually created empty data.frame
            find_bad_taxa(cnd)
          })   
      }
    )
  }
  ```
  **Custom Functions Used:**
  - [find_bad_taxa()](#find_bad_taxa)

  **Function Parameter Definitions:**
  - `data=` - specifies the treeSummarizedExperiment containing the feature, taxonomy and metadata to be analyzed using ancombc2
  - `...` - Other arguments passed on to ancombc2

  **Returns:** an ancombc2 result, or an empty result as returned by [find_bad_taxa()](#find_bad_taxa)
</details>

#### gm_mean()
<details>
  <summary>calculates the geometric mean</summary>

  ```R
    gm_mean <- function(x, na.rm=TRUE) {
      exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
    }
  ```
  **Function Parameter Definitions:**
  - `x=` - a numeric vector specifying the values to calculate the geometirc mean on
  - `na.rm=TRUE` - boolean specifying if NAs should be removed prior to calculating the geometirc mean; default is TRUE
  
  **Returns:** a numeric value representing the geometric mean
</details>

#### plotSparsity()
<details>
  <summary>Plots ASV sparsity. A modification of DESeq2::plotSparsity to generate a ggplot.</summary>

  ```R
    plotSparsity <- function (x, normalized = TRUE, feature="ASV", ...) {
  
  
      if (is(x, "DESeqDataSet")) {

            x <- counts(x, normalized = normalized)
      }
  
      rs <- MatrixGenerics::rowSums(x)
      rmx <- apply(x, 1, max)

      # Prepare plot dataframe
      df <- data.frame(rs=rs, rmx=rmx) %>% 
               mutate(x=rs, y=rmx/rs) %>% 
               filter(x>0)
    
      # Plot
      ggplot(data = df, aes(x=x, y=y), ...) +
        geom_point(size=3) +
        scale_x_log10() + 
        scale_y_continuous(limits = c(0,1)) +
        theme_bw() +
        labs(title = glue("Concentration of {feature} counts over total sum of {feature} counts"),
             x=glue("Sum of counts per {feature}"),
             y=glue("Max {feature} count / Sum of {feature} counts")) + 
        theme(axis.text = element_text(face = "bold", size = 12),
              axis.title = element_text(face = "bold", size = 14),
              title = element_text(face = "bold", size = 14))

  }
    
  ```
  **Function Parameter Definitions:**
  - `x=` -  a matrix or DESeqDataSet to plot
  - `normalized=TRUE` - boolean specifying whether to normalize the counts from a DESeqDataSEt; default is TRUE
  - `feature=` - a string specifying which feature type ("ASV", "OTU", "gene" etc.) is being plotted; Default is "ASV"
  - `...=` - any named argument(s) that can be passed to the ggplot2::ggplot function.

  **Returns:** a sparsity plot of type ggplot
</details>


<br>

#### 6b.iii. Set Variables

```R
# Define a custom palette for plotting
custom_palette <- c("#1F78B4","#33A02C","#FB9A99","#E31A1C","#6A3D9A",
                    "#FDBF6F", "#FF7F00","#CAB2D6","#FF00FFFF", "#B15928",
                    "#000000","#FFC0CBFF", "#A6CEE3", "#8B864EFF","#F0027F",
                    "#666666","#1B9E77", "#E6AB02","#A6761D","#FFFF00FF",
                    "#00FFFFFF", "#FFFF99", "#B2182B","#FDDBC7","#D1E5F0",
                    "#B2DF8A","#CC0033","#FF00CC","#330033", "#999933",
                    "#FF9933", "#FFFAFAFF",colors()) 

# Remove white colors
pattern_to_filter <- "white|snow|azure|gray|#FFFAFAFF|aliceblue"
custom_palette <- custom_palette[-c(21:23, grep(pattern = pattern_to_filter,
                                                x = custom_palette, 
                                                ignore.case = TRUE))]
# Custom theme for plotting
publication_format <- theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.ticks.length=unit(-0.15, "cm"),
        axis.text.x=element_text(margin=ggplot2::margin(t=0.5,r=0,b=0,l=0,unit ="cm")),
        axis.text.y=element_text(margin=ggplot2::margin(t=0,r=0.5,b=0,l=0,unit ="cm")), 
        axis.title = element_text(size = 18,face ='bold.italic', color = 'black'), 
        axis.text = element_text(size = 16,face ='bold', color = 'black'),
        legend.position = 'right', 
        legend.title = element_text(size = 15,face ='bold', color = 'black'),
        legend.text = element_text(size = 14,face ='bold', color = 'black'),
        strip.text = element_text(size = 14,face ='bold', color = 'black'))
```

**Input Data:** 

*No input data required*

**Output Data:**

* `custom_palette` (a vector of strings specifying a custom color palette for coloring plots)
* `publication_format` (a ggplot::theme object specifying the custom theme for plotting)

<br>

#### 6b.iv. Read-in Input Tables

```R
diff_abund_out_dir <- "differential_abundance/"
if(!dir.exists(diff_abund_out_dir)) dir.create(diff_abund_out_dir)
assay_suffix <- "_GLAmpSeq"
output_prefix <- ""
custom_palette <- {COLOR_VECTOR}
groups_colname <- "groups"
sample_colname <- "Sample Name"
metadata_file <- file.path("{OSD-Accession-ID}_amplicon_v{version}_runsheet.csv")
features_file <- file.path("counts_<tech_type>_GLAmpSeq.tsv")
taxonomy_file <- file.path("taxonomy_<tech_type>_GLAmpSeq.tsv")

# Read-in metadata and convert from tibble to dataframe
metadata <- read_csv(file = metadata_file) %>% as.data.frame()
# Set row names
row.names(metadata) <- metadata[[sample_colname]]
# Write out Sample Table
write_csv(x = metadata %>% select(!!sym(sample_colname), !!sym(groups_colname)),
          file = glue("{diff_abund_out_dir}{output_prefix}SampleTable_<tech_type>{assay_suffix}.csv"))

# Delete sample column since the rownames now contain sample names
metadata[,sample_colname] <- NULL
# Get unique group names
group_column_values <- metadata %>% pull(!!sym(groups_colname))
group_levels <- unique(group_column_values)

# Write out table listing contrasts used for all differential abundance methods
# Get pairwise combinations
pairwise_comp.m <- utils::combn(group_levels, 2)
# Create comparison names  
comparisons <- paste0("(", pairwise_comp.m[2,], ")v(", pairwise_comp.m[1,], ")")
names(comparisons) <- comparisons
# Create contrasts table
contrasts_df <- data.frame(
  " " = c("1", "2"),
  rbind(pairwise_comp.m[2,], pairwise_comp.m[1,]) %>% as.data.frame() %>% setNames(comparisons),
  check.names = FALSE
)
write_csv(x = contrasts_df,
          file = glue("{diff_abund_out_dir}{output_prefix}contrasts_<tech_type>{assay_suffix}.csv"))

# Add colors to metadata that equals the number of groups
num_colors <- length(group_levels)
palette <- 'Set1'
number_of_colors_in_palette <- 9
if(num_colors <= number_of_colors_in_palette){
   colors <- RColorBrewer::brewer.pal(n = num_colors, name = palette)
}else{
  colors <- custom_palette[1:num_colors]
}

# ------ Metadata ----- #
# Assign color names to each group
group_colors <- setNames(colors, group_levels)
metadata <- metadata %>%
  mutate(color = map_chr(!!sym(groups_colname),
                         function(group) { group_colors[group] }
                         ) 
        ) # assign group specific colors to each row in metadata

# Retrieve sample names
sample_names <- rownames(metadata)
deseq2_sample_names <- make.names(sample_names, unique = TRUE)

# Subset metadata to contain only the groups and color columns
sample_info_tab <- metadata %>%
  select(!!groups_colname, color) %>% # select groups and color columns
  arrange(!!sym(groups_colname)) # metadata by groups column

# Retrieves unique colors
values <- sample_info_tab %>% pull(color) %>% unique()

# ---- Import Feature or ASV table ---- #
feature_table <- read.table(file = features_file, header = TRUE,
                              row.names = 1, sep = "\t", 
                              check.names = FALSE)

# ---- Import Taxonomy table ---- #
taxonomy_table <- read.table(file = taxonomy_file, header = TRUE,
                              row.names = 1, sep = "\t"),
                              check.names = FALSE)
```

**Input Data:**

* `diff_abund_out_dir` (a string specifying the path to the output folder for the differential abundance results, default is "differential_abundance/")
* `assay_suffix` (a string specifying the suffix to be added to output files; default is the Genelab assay suffix, "_GLAmpSeq")
* `output_prefix` (a string specifying an additional prefix to be added to the output files; default is no additional prefix, "")
* `groups_colname` (a string specifying the name of the column in the metadata table containing the group names)
* `sample_colname` (a string specifying the name of the column in the metadata table containing the sample names)
* `custom_palette` (a vector of strings specifying a custom color palette for coloring plots, output from [6b.iii. Set Variables](#6biii-set-variables))
* {OSD-Accession-ID}_amplicon_v{version}_runsheet.csv (a comma-separated sample metadata file containing sample group information, output from [Step 6a](#6a-create-sample-runsheet))
*	counts_<tech_type>_GLAmpSeq.tsv (a tab-separated file containing sample feature counts table (i.e. ASV or OTU table), output from [Step 5g](#5g-generating-and-writing-standard-outputs))
* taxonomy_<tech_type>_GLAmpSeq.tsv (a tab-separated file containing feature taxonomy table containing ASV taxonomy assignments, output from [Step 5g](#5g-generating-and-writing-standard-outputs))

**Output Data:**

* `metadata` (a dataframe containing the sample metadata, with samples as row names and sample info as columns)
* `feature_table` (a dataframe containing the sample feature counts table (i.e. ASV or OTU table) from the input counts file)
* `taxonomy_table` (a dataframe containing ASV taxonomy assignments from the input taxonomy file)
* `sample_info_tab` (a dataframe containing a subset of the metadata dataframe with only the groups and color columns)
* `values` (a character vector of unique color values for each group)
* `sample_names` (a character vector of sample names)
* `deseq2_sample_names` (a character vector of unique sample names)
* `group_colors` (a named character vector of colors for each group)
* `group_levels` (a character vector of unique group names)
* **differential_abundance/SampleTable_<tech_type>_GLAmpSeq.csv** (a comma-separated file containing a table with two columns: "Sample Name" and "groups"; the output_prefix denotes the method used to compute the differential abundance)
* **differential_abundance/contrasts_<tech_type>_GLAmpSeq.csv** (a comma-separated file listing all pairwise group comparisons)

<br>

#### 6b.v. Preprocessing
Filters the feature and taxonomy tables to include only features that (a) pass the specified prevalence and library count thresholds and (b) are not from Chloroplast or Mitochondrial Organelle contamination.

```R
feature_table <- {DATAFRAME} # from step [Read-in Input Tables]
taxonomy_table <- {DATAFRAME} # from step [Read-in Input Tables]
target_region <- "16S" # 16S, 18S, or ITS
remove_rare <- FALSE # TRUE OR FALSE
prevalence_cutoff <- 0
library_cutoff <- 0


if(remove_rare){
  
  # Remove samples with less than library-cutoff
  message(glue("Dropping samples with less than {library_cutoff} read counts"))
  feature_table <- feature_table[,colSums(feature_table) >= library_cutoff]
  # Remove rare ASVs
  message(glue("Dropping features with prevalence less than {prevalence_cutoff * 100}%"))
  feature_table <- remove_rare_features(feature_table,
                                        cut_off_percent = prevalence_cutoff)
}

# Preprocess ASV and taxonomy tables

message(glue("There are {sum(is.na(taxonomy_table$domain))} features without 
           taxonomy assignments. Dropping them..."))

# Dropping features that couldn't be assigned taxonomy
# For beta and alpha diversity only, unassigned ASVs are not dropped in DA analyses
taxonomy_table <- taxonomy_table[-which(is.na(taxonomy_table$domain)),]

# Handle case where no domain was assigned but a phylum was.
if(all(is.na(taxonomy_table$domain))){
  
  if(target_region == "ITS"){
    taxonomy_table$domain <- "Fungi"
  }else if(target_region == "18S"){
    taxonomy_table$domain <- "Eukaryotes"
  }else{
    taxonomy_table$domain <- "Bacteria"
  }
  
}

# Removing Chloroplast and Mitochondria Organelle DNA contamination
asvs2drop <- taxonomy_table %>%
               unite(col="taxonomy",domain:species) %>%
               filter(str_detect(taxonomy, "[Cc]hloroplast|[Mm]itochondria")) %>%
               row.names()
taxonomy_table <- taxonomy_table[!(rownames(taxonomy_table) %in% asvs2drop),]

# Clean taxonomy names
feature_names <- rownames(taxonomy_table)
taxonomy_table <- process_taxonomy(taxonomy_table)
rownames(taxonomy_table) <- feature_names
taxonomy_table <- fix_names(taxonomy_table, "Other", ";_")

# Subset tables 

# Get features common to the taxonomy and feature table 
common_ids <- intersect(rownames(feature_table), rownames(taxonomy_table))

# Subset the feature and taxonomy tables to contain 
# only features found in both tables
feature_table <- feature_table[common_ids,]
taxonomy_table <- taxonomy_table[common_ids,]

# drop samples with zero sequence counts
samples2keep <-  colnames(feature_table)[colSums(feature_table) > 0]

feature_table <- feature_table[, samples2keep]
metadata <- metadata[samples2keep,]
```
**Custom Functions Used:**

* [remove_rare_features()](#remove_rare_features)
* [process_taxonomy()](#process_taxonomy)
* [fix_names()](#fix_names)

**Parameter Definitions:**

* `remove_rare` - boolean specifying if rare features and samples should be filtered out based on the `prevalence_cutoff` and `library_cutoff` cutoff thresholds, respectively, prior to analysis; default is FALSE 
* `prevalence_cutoff` - a decimal between 0 and 1 specifying the proportion of samples required to contain a taxon in order to keep the taxon when `remove_rare` is set to TRUE; default is 0, i.e. do not exclude any taxon / feature
* `library_cutoff` - a numerical value specifying the number of total counts a sample must have across all features to be retained when `remove_rare` is set to TRUE; default is 0, i.e. no samples will be dropped
* `target_region` - a string specifying the amplicon target region; options are either "16S", "18S", or "ITS"

**Input Data:**

* `feature_table` (a dataframe containing sample feature counts (i.e. ASV or OTU table), output from [6b.iv. Read-in Input Tables](#6biv-read-in-input-tables))
* `taxonomy_table` (a dataframe containing ASV taxonomy assignments, output from [6b.iv. Read-in Input Tables](#6biv-read-in-input-tables))

**Output Data:**

* `feature_table` (a dataframe containing a filtered subset of the samples feature counts (i.e. ASV or OTU table) after removing features that do not meet the filtering thresholds or that belong to Chroroplast or Mitochondrial organelles, and samples with zero sequence counts)
* `taxonomy_table` (a dataframe containing a filtered subset of the feature taxonomy table after removing ASV taxonomy assignments for features that do not meet the filtering thresholds or that belong to Chroroplast or Mitochondrial organelles)

<br>

---

## 7. Alpha Diversity Analysis

Alpha diversity examines the variety and abundance of taxa within individual samples. Rarefaction curves are utilized to 
visually represent this diversity, plotting the number of unique sequences (ASVs) identified against the total number of 
sequences sampled, offering a perspective on the saturation and completeness of sampling. Metrics like Observed features 
estimates and Shannon diversity indices are employed to quantify the richness (total number of unique sequences) and 
diversity (combination of richness and evenness) within these samples, respectively.

### 7a. Rarefaction Curves
```R
# Create output directory if it doesn't already exist
alpha_diversity_out_dir <- "alpha_diversity/"
if(!dir.exists(alpha_diversity_out_dir)) dir.create(alpha_diversity_out_dir)
sample_info_tab <- {DATAFRAME} 
feature_table <- {DATAFRAME}
taxonomy_table <- {DATAFRAME}
group_colors <- {NAMED_VECTOR} 
groups_colname <- "groups"
rarefaction_depth <- 500
legend_title <- "Groups"
assay_suffix <- "_GLAmpSeq"
output_prefix <- ""

# Create phyloseq object
ASV_physeq <- phyloseq(otu_table(feature_table, taxa_are_rows = TRUE),
                       tax_table(as.matrix(taxonomy_table)), 
                       sample_data(sample_info_tab))

seq_per_sample <- colSums(feature_table) %>% sort()

# Get rarefaction depth
# minimum value
depth <- min(seq_per_sample)

# Error if the number of sequences per sample left after filtering is
# insufficient for diversity analysis
if(max(seq_per_sample) < 100){
 
  warning_file <- glue("{alpha_diversity_out_dir}{output_prefix}alpha_diversity_failure_<tech_type>{assay_suffix}.txt")
  writeLines(
    text = glue("The maximum sequence count per sample ({max(seq_per_sample)}) is less than 100.
Therefore, alpha diversity analysis cannot be performed."),
    con = warning_file
  )
  quit(status = 0)
}

for (count in seq_per_sample) {

  if(count >= rarefaction_depth) {
    depth <- count
    break
    }

}

# Error if the depth that ends up being used is also less than 100
if(depth < 100){
 
  warning_file <- glue("{alpha_diversity_out_dir}{output_prefix}alpha_diversity_failure_<tech_type>{assay_suffix}.txt")
  writeLines(
    text = glue("The rarefaction depth being used in the analysis ({depth}) is less than 100.
Therefore, alpha diversity analysis cannot be performed."),
    con = warning_file
  )
  quit(status = 0)
} 

#Warning if rarefaction depth is between 100 and 500
if (depth > 100 && depth < 500) {

  warning(glue("Rarefaction depth ({depth}) is between 100 and 500.
Alpha diversity results may be unreliable."))

}

# -------------------- Rarefy sample counts to even depth per sample
ps.rarefied <- rarefy_even_depth(physeq = ASV_physeq, 
                                 sample.size = depth,
                                 rngseed = 1, 
                                 replace = FALSE, 
                                 verbose = FALSE)

# Write rarefaction depth used into file to be used in protocol
depth_file <- glue("{alpha_diversity_out_dir}{output_prefix}rarefaction_depth_<tech_type>{assay_suffix}.txt")
writeLines(
  text = as.character(depth),
  con = depth_file
)


# ------------------- Rarefaction curve
# Calculate a rough estimate of the step sample step size for plotting.
# This is meant to keep plotting time constant regardless of sample depth
step <- (50*depth)/1000

p <- rarecurve(t(otu_table(ps.rarefied)) %>% as.data.frame(),
               step = step,
               col = sample_info_tab[["color"]], 
               lwd = 2, ylab = "ASVs", cex=0.5,
               label = FALSE, tidy = TRUE)


sample_info_tab_names <- sample_info_tab %>% rownames_to_column("Site")

p <- p %>% left_join(sample_info_tab_names, by = "Site")

# Sample rarefaction curves

rareplot <- ggplot(p, aes(x = Sample, y = Species,
                          group = Site, color = !!sym(groups_colname))) + 
  geom_line() + 
  scale_color_manual(values = group_colors) +
  labs(x = "Number of Sequences", y = "Number of ASVs", color = legend_title) + 
  theme_bw() +
  theme(legend.position = "right",
        text = element_text(face = 'bold', size = 15),
        legend.text = element_text(face = 'bold', size = 14),
        legend.direction = "vertical",
        legend.justification = "center",
        legend.box.just = "center",
        legend.title = element_text(size = 15, face='bold'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(t = 10, r = 20, b = 10, l = 10, unit = "pt"))

ggsave(filename = glue("{alpha_diversity_out_dir}/{output_prefix}rarefaction_curves_<tech_type>{assay_suffix}.png"),
       plot=rareplot, width = 14, height = 8.33, dpi = 300, limitsize = FALSE)
```

**Input Data:**

* `alpha_diversity_out_dir` (a string specifying the path to the output folder for the alpha diversity results, default is "alpha_diversity/")
* `rarefaction_depth` (an integer specifying the minimum number of reads to simulate during rarefaction for alpha diversity estimation)
* `groups_colname` (a string specifying the name of the column in the `sample_info_tab` table containing the group names)
* `legend_title` (a string specifying the legend title for plotting)
* `assay_suffix` (a string specifying the suffix to be added to output files; default is the Genelab assay suffix, "_GLAmpSeq")
* `output_prefix` (a string specifying an additional prefix to be added to the output files; default is no additional prefix, "")
* `sample_info_tab` (a dataframe containing a subset of the metadata dataframe with only the groups and color columns, output from [6b.iv. Read-in Input Tables](#6biv-read-in-input-tables))
* `feature_table` (a dataframe containing a filtered subset of the samples feature dataframe (i.e. ASV), output from [6b.v. Preprocessing](#6bv-preprocessing))
* `taxonomy_table` (a dataframe containing a filtered subset of the feature taxonomy dataframe with ASV taxonomy assignments, output from [6b.v. Preprocessing](#6bv-preprocessing))
* `group_colors` (a named character vector of colors for each group, output from [6b.iv. Read-in Input Tables](#6biv-read-in-input-tables))

**Output Data:**
> NOTE: If alpha diversity analysis couldn't be perfomed due to insufficient sequence counts per samples, a failure file (alpha_diversity_failure_<tech_type>_GLAmpSeq.txt) will be generated instead of the below output, and the subsequent alpha diversity steps (7b, 7c, and 7d) will be skipped.

* `ps.rarefied` (a phyloseq object of the sample features (i.e. ASV) with feature counts derived from the `feature_table`, resampled such that all samples have the same library size)
* rarefaction_curves_<tech_type>_GLAmpSeq.png (plot containing the rarefaction curves for each sample)
* alpha_diversity/rarefaction_depth_<tech_type>_GLAmpSeq.txt (rarefaction depth value used in alpha analysis)

<br>

### 7b. Richness and Diversity Estimates
```R
metadata <- {DATAFRAME}
groups_colname <- "groups"
assay_suffix <- "_GLAmpSeq"
output_prefix <- ""

# ------------------  Richness and diversity estimates  ------------------#

# Statistics table
diversity_metrics <- c("Observed", "Chao1", "Shannon", "Simpson")
names(diversity_metrics) <- diversity_metrics
diversity.df <- estimate_richness(ps.rarefied,
                                  measures = diversity_metrics) %>%
               select(-se.chao1) %>%
               rownames_to_column("samples")

merged_table <- metadata %>%
  rownames_to_column("samples") %>%
  inner_join(diversity.df)

diversity_stats <- map_dfr(.x = diversity_metrics, function(metric){
  

  number_of_groups <- merged_table[,groups_colname] %>% unique() %>% length()
  
  if (number_of_groups < 2){
    warning_file <- glue("{alpha_diversity_out_dir}{output_prefix}alpha_diversity_warning.txt")
    original_groups <- length(unique(metadata[[groups_colname]]))
    writeLines(
      text = glue("Group count information:
Original number of groups: {original_groups}
Number of groups after filtering: {number_of_groups}

There are less than two groups to compare, hence, pairwise comparisons cannot be performed.
Please ensure that your metadata contains two or more groups to compare..."),
      con = warning_file
    )
    quit(status = 0)
  }else if(number_of_groups == 2){
    
  df <- data.frame(y=merged_table[,metric], x=merged_table[,groups_colname]) %>%
      wilcox_test(y~x) %>% 
      adjust_pvalue(method = "bonferroni") %>%
      select(group1, group2, W=statistic, p, p.adj) %>% 
      mutate(Metric=metric) %>% 
    add_significance(p.col='p.adj', output.col = 'p.signif') %>% 
    select(Metric,group1, group2, W, p, p.adj, p.signif)
      
  }else{
    
  res <- dunnTest(merged_table[,metric],merged_table[,groups_colname])
  
  df <- res$res %>%
    separate(col = Comparison, into = c("group1", "group2"), sep = " - ") %>% 
    mutate(Metric=metric)  %>% 
    rename(p=P.unadj, p.adj=P.adj) %>% 
    add_significance(p.col='p.adj', output.col = 'p.signif') %>% 
    select(Metric,group1, group2, Z, p, p.adj, p.signif)
  
  }

  return(df)
})

# Write diversity statistics table to file
write_csv(x = diversity_stats, 
            file = glue("{alpha_diversity_out_dir}/{output_prefix}statistics_table_<tech_type>{assay_suffix}.csv"))

# Get different letters indicating statistically significant group comparisons for every diversity metric
comp_letters <- data.frame(group = group_levels)
colnames(comp_letters) <- groups_colname

walk(.x = diversity_metrics, function(metric = .x) {

  sub_comp <- diversity_stats %>% filter(Metric == metric)

  sanitize <- function(x) gsub("-", "_", x)
  g1 <- sanitize(sub_comp$group1)
  g2 <- sanitize(sub_comp$group2)

  safe_names <- paste(g1, g2, sep = "-")
  orig_names <- paste(sub_comp$group1, sub_comp$group2, sep = "-")
  safe_to_orig <- setNames(orig_names, safe_names)

  p_values <- setNames(sub_comp$p.adj, safe_names)

  letters <- multcompView::multcompLetters(p_values)$Letters
  names(letters) <- safe_to_orig[names(letters)]

  letters_df <- enframe(letters,
                        name = groups_colname,
                        value = glue("{metric}_letter"))

  comp_letters <<- comp_letters %>% left_join(letters_df)
})

# Summary table
diversity_table <- metadata %>%
  rownames_to_column("samples") %>%
  inner_join(diversity.df) %>%
  group_by(!!sym(groups_colname)) %>%
  summarise(N = n(), across(Observed:Simpson,
                   .fns = list(mean = mean, se = se),
                   .names = "{.col}_{.fn}")) %>%
  mutate(across(where(is.numeric), ~round(.x, digits = 2))) %>%
  left_join(comp_letters) %>%
  mutate(Observed = glue("{Observed_mean} ± {Observed_se}{Observed_letter}"),
         Chao1 = glue("{Chao1_mean} ± {Chao1_se}{Chao1_letter}"),
         Shannon = glue("{Shannon_mean} ± {Shannon_se}{Shannon_letter}"),
         Simpson = glue("{Simpson_mean} ± {Simpson_se}{Simpson_letter}")
         ) %>%
  select(-contains("_"))

# Write diversity summary table to file 
write_csv(x = diversity_table,
            file = glue("{alpha_diversity_out_dir}/{output_prefix}summary_table_<tech_type>{assay_suffix}.csv"))
```

**Input Data:**

* `groups_colname` (a string specifying the name of the column in the metadata table containing the group names)
* `assay_suffix` (a string specifying the suffix to be added to output files; default is the Genelab assay suffix, "_GLAmpSeq")
* `output_prefix` (a string specifying an additional prefix to be added to the output files; default is no additional prefix, "")
* `metadata` (a dataframe containing the sample metadata, with samples as row names and sample info as columns, output from [6b.iv. Read-in Input Tables](#6biv-read-in-input-tables))
* `ps.rarefied` (a phyloseq object of the sample features (i.e. ASV) with feature counts, resampled such that all samples have the same library size, output from [7a. Rarefaction Curves](#7a-rarefaction-curves))
* `group_levels` (a character vector of unique group names, output from [6b.iv. Read-in Input Tables](#6biv-read-in-input-tables))

**Output Data:**

* **alpha_diversity/statistics_table_<tech_type>_GLAmpSeq.csv** (a comma-separated table containing the z-score, p-value, and adjusted p-value statistics for each pairwise comparison for all metrics evaluated, Observed, Chao1, Shannon, and Simpson)
* **alpha_diversity/summary_table_<tech_type>_GLAmpSeq.csv** (a comma-separated table containing the sample number and mean +/- standard error of each metric (Observed, Chao1, Shannon, and Simpson) for each group)

<br>

### 7c. Plot Richness and Diversity Estimates

```R
sample_info_tab <- {DATAFRAME}
metadata <- {DATAFRAME}
groups_colname <- "groups"
legend_title <- "Groups"
assay_suffix <- "_GLAmpSeq"
output_prefix <- ""

# ------------------ Make richness by sample dot plots ---------------------- #

number_of_samples <- length(rownames(sample_info_tab))
richness_sample_label_size <- calculate_text_size(number_of_samples)
metrics2plot <- c("Observed", "Shannon")
names(metrics2plot) <- metrics2plot

samples_order <- metadata %>% arrange(!!sym(groups_colname)) %>% rownames()

richness_by_sample <- plot_richness(ps.rarefied, color = groups_colname,
                                    measures = metrics2plot)

richness_by_sample <- ggplot(richness_by_sample$data %>% 
                                mutate(samples = factor(samples, 
                                                        levels=samples_order)),
                              aes(x=samples, y=value, colour = !!sym(groups_colname))) +
  geom_point() +
  geom_errorbar(aes(ymin=value-se, ymax = value+se),
                width=0.2, position=position_dodge(0.9)) +
  facet_wrap(~variable, scales = "free_y") +
  scale_color_manual(values = group_colors) + 
  theme_bw() +labs(x = NULL, color = legend_title, y="Alpha Diversity Measure") +
  theme(
    text = element_text(face = 'bold', size = 15),
    legend.text = element_text(face = 'bold', size = 14),
    legend.position = "bottom",
    legend.direction = "vertical",
    legend.justification = "center",
    legend.box.just = "center",
    legend.title = element_text(face = 'bold', size = 15, hjust = 0.09),
    axis.text.x = element_text(angle = 90,
                               size = richness_sample_label_size,
                               vjust = 0.5,  # Vertically center the text
                               hjust = 1),
    axis.ticks.length=unit(-0.15, "cm"),
    strip.text = element_text(size = 14,face ='bold')
  )

# Save sample plot
ggsave(filename = glue("{alpha_diversity_out_dir}/{output_prefix}richness_and_diversity_estimates_by_sample_<tech_type>{assay_suffix}.png"),
       plot=richness_by_sample, width = 14, height = 8.33, 
       dpi = 300, units = "in", limitsize = FALSE)

# ------------------- Make richness by group box plots ----------------------- #
richness_by_group <- plot_richness(ps.rarefied, x = groups_colname, 
                                   color = groups_colname,
                                   measures = metrics2plot)

p <- map(.x = metrics2plot, .f = function(metric){
    
  p <- ggplot(richness_by_group$data %>% filter(variable == metric), 
              aes(x=!!sym(groups_colname), y=value, fill=!!sym(groups_colname)) 
              ) +
    geom_point() +
    geom_boxplot() +
    scale_fill_manual(values = group_colors) + 
    theme_bw() + labs(fill = legend_title, x = NULL, y= metric) +
    theme(
      text = element_text(size = 15, face = 'bold'),
      legend.text = element_text(face = 'bold', size = 14),
      legend.position = "right",
      legend.direction = "vertical",
      legend.justification = "center",
      legend.box.just = "center",
      legend.title = element_text(face = 'bold', size = 15),
      axis.text.x = element_blank(),
      axis.ticks.length=unit(-0.15, "cm"),
      strip.text = element_text(size = 14,face ='bold')
    ) 
  
  
  summary_table <- p$data %>%
    select(!!sym(groups_colname), value) %>% 
    group_by(!!sym(groups_colname)) %>%
    summarise(max=max(value), range=max(value)-min(value)) %>%
    left_join(comp_letters %>%
                select(!!sym(groups_colname), label= !!sym( glue("{metric}_letter") ) 
                       ) 
              )
  text_size <- 6
  
  # Calculate a constant to add to the max value of each group
  # to determine where each group text will be added 
  toAdd <- if_else(condition = max(summary_table$range) <= 5,
                   true = min(summary_table$range),
                   false = (median(summary_table$range) - min(summary_table$range)) / 20
                    )
  
  # Add text to plot
  p + geom_text(data=summary_table,
                mapping = aes(y=max+toAdd, label=label, fontface = "bold"),
                size = text_size)
})

richness_by_group <- wrap_plots(p, ncol = 2, guides = 'collect') + 
                     plot_annotation(caption = "If letters are shared between two groups, then they are not significantly different (q-value > 0.05)",
                                     theme = theme(plot.caption = element_text(face = 'bold.italic')) 
                                     )

# Save group plot
width <- 3.6 * length(group_levels)
ggsave(filename = glue("{output_prefix}richness_and_diversity_estimates_by_group_<tech_type>{assay_suffix}.png"),
       plot=richness_by_group, width = width, 
       height = 8.33, dpi = 300, units = "in",
       path = alpha_diversity_out_dir)

```
**Custom Functions Used:**

* [calculate_text_size()](#calculate_text_size)

**Input Data:**

* `groups_colname` (a string specifying the name of the column in the metadata table containing the group names)
* `legend_title` (a string specifying the legend title for plotting)
* `assay_suffix` (a string specifying the suffix to be added to output files; default is the Genelab assay suffix, "_GLAmpSeq")
* `output_prefix` (a string specifying an additional prefix to be added to the output files; default is no additional prefix, "")
* `metadata` (a dataframe containing the sample metadata, with samples as row names and sample info as columns, output from [6b.iv. Read-in Input Tables](#6biv-read-in-input-tables))
* `sample_info_tab` (a dataframe containing a subset of the metadata dataframe with only the groups and color columns, output from [6b.iv. Read-in Input Tables](#6biv-read-in-input-tables))
* `ps.rarefied` (a phyloseq object of the sample features (i.e. ASV) with feature counts, resampled such that all samples have the same library size, output from [7a. Rarefaction Curves](#7a-rarefaction-curves))
* `group_levels` (a character vector of unique group names, output from [6b.iv. Read-in Input Tables](#6biv-read-in-input-tables))
* `group_colors` (a named character vector of colors for each group, output from [6b.iv. Read-in Input Tables](#6biv-read-in-input-tables))


**Output Data:**

* richness_and_diversity_estimates_by_sample_<tech_type>_GLAmpSeq.png (dot plots containing richness and diversity estimates for each sample)
* richness_and_diversity_estimates_by_group_<tech_type>_GLAmpSeq.png (box plots containing richness and diversity estimates for each group)

<br>

### 7d. Package alpha diversity plots

```bash
zip -q alpha_diversity_plots_<tech_type>_GLAmpSeq.zip *.png
```

**Parameter Definitions:**

- `-q` – quiet operation
- `alpha_diversity_plots_<tech_type>_GLAmpSeq.zip` – positional argument naming the zip output file
- `*.png` – positional argument naming the input file(s) to package

**Input Data:**

* rarefaction_curves_<tech_type>_GLAmpSeq.png (plot containing the rarefaction curves for each sample, output from [Step 7a](#7a-rarefaction-curves))
* richness_and_diversity_estimates_by_sample_<tech_type>_GLAmpSeq.png (dot plots containing richness and diversity estimates for each sample, output from [Step 7c](#7c-plot-richness-and-diversity-estimates))
* richness_and_diversity_estimates_by_group_<tech_type>_GLAmpSeq.png (box plots containing richness and diversity estimates for each group, output from [Step 7c](#7c-plot-richness-and-diversity-estimates))

**Output Data:**

* **alpha_diversity/alpha_diversity_plots_<tech_type>_GLAmpSeq.zip**

<br>  

---

## 8. Beta Diversity Analysis

Beta diversity measures the variation in species composition between different samples or environments. A common practice in working with a new dataset is to generate some exploratory visualizations like ordinations and hierarchical clusterings. These give us a quick overview of how our samples relate to each other and can be a way to check for problems like batch effects.

Two normalization methods are supported before performing hierarchical clustering: variance stabilizing transformation (VST) and rarefaction. After rarefaction, the default Bray-Curtis dissimilarity can be used to generate dissimilarity matrices for hierarchical clustering. VST, however, generates negative values which are incompatible with calculating Bray-Curtis dissimilarity. For VST transformed data, Euclidean distance is used instead.

```R
beta_diversity_out_dir <- "beta_diversity/"
if(!dir.exists(beta_diversity_out_dir)) dir.create(beta_diversity_out_dir)
metadata <- {DATAFRAME}
feature_table <- {DATAFRAME} 
group_colors <- {NAMED_VECTOR} 
groups_colname <- "groups"
rarefaction_depth <- 500
legend_title <- "Groups"
assay_suffix <- "_GLAmpSeq"
output_prefix <- ""
distance_methods <- c("euclidean", "bray")
normalization_methods <- c("vst", "rarefy")

options(warn=-1) # ignore warnings

# Run the analysis
walk2(.x = normalization_methods, .y = distance_methods,
      .f = function(normalization_method, distance_method){
  
# Create transformed phyloseq object
ps <- transform_phyloseq(feature_table, metadata, 
                        method = normalization_method,
                        rarefaction_depth = rarefaction_depth)

# Skip downstream analysis when normalization by rarefaction fails
if (is.null(ps)) {
  message(glue("{normalization_method} failed. Skipping downstream analysis."))
  return(NULL)
}
  
  # ---------Clustering and dendrogram plotting

  # Extract normalized count table
  count_tab <- otu_table(ps)

  # VSD validation check point
  if(normalization_method == "vst"){
  
  # Visualize the sd vs the rank of the mean plot.
  mead_sd_plot <- meanSdPlot(t(count_tab))$gg + 
    theme_bw() +
    labs(title = "MEAN-SD VST validation Plot",
         x="Rank Of Mean",
         y="Standard Deviation") + 
    theme(axis.text = element_text(face = "bold", size = 12),
          axis.title = element_text(face = "bold", size = 14),
          title = element_text(face = "bold", size = 14))
  
  # Save VSD validation plot
  ggsave(filename = glue("{beta_diversity_out_dir}/{output_prefix}vsd_validation_plot_<tech_type>{assay.suffix}.png"),
         plot = mead_sd_plot, width = 14, height = 10, 
         dpi = 300, units = "in", limitsize = FALSE)
  }


  # Calculate distance between samples
  dist_obj <- vegdist(t(count_tab), method = distance_method)

  # Make dendrogram
  dendrogram <- make_dendrogram(dist_obj, metadata, groups_colname,
                              group_colors, legend_title)

  # Save dendrogram
  ggsave(filename = glue("{beta_diversity_out_dir}/{output_prefix}{distance_method}_dendrogram_<tech_type>{assay_suffix}.png"),
       plot = dendrogram, width = 14, height = 10, 
       dpi = 300, units = "in", limitsize = FALSE)

  #---------------------------- Run stats
  # Checking homogeneity of variance and comparing groups using adonis test

  stats_res <- run_stats(dist_obj, metadata, groups_colname)
  write_csv(x = stats_res$variance, 
            file = glue("{beta_diversity_out_dir}/{output_prefix}{distance_method}_variance_table_<tech_type>{assay_suffix}.csv"))

  write_csv(x = stats_res$adonis, 
            file = glue("{beta_diversity_out_dir}/{output_prefix}{distance_method}_adonis_table_<tech_type>{assay_suffix}.csv"))

  #---------------------------- Make PCoA
  # Unlabeled PCoA plot
  ordination_plot_u <- plot_pcoa(ps, stats_res, distance_method, 
                                 groups_colname, group_colors, legend_title) 
  ggsave(filename=glue("{beta_diversity_out_dir}/{output_prefix}{distance_method}_PCoA_without_labels_<tech_type>{assay_suffix}.png"),
       plot=ordination_plot_u, width = 14, height = 8.33, 
       dpi = 300, units = "in", limitsize = FALSE)

  # Labeled PCoA plot
  ordination_plot <- plot_pcoa(ps, stats_res, distance_method,
                               groups_colname, group_colors, legend_title,
                               addtext=TRUE) 
  ggsave(filename=glue("{beta_diversity_out_dir}/{output_prefix}{distance_method}_PCoA_w_labels_<tech_type>{assay_suffix}.png"),
       plot=ordination_plot, width = 14, height = 8.33, 
       dpi = 300, units = "in", limitsize = FALSE)

})
```
```bash
zip -q bray_curtis_plots_<tech_type>_GLAmpSeq.zip bray*.png

zip -q euclidean_distance_plots_<tech_type>_GLAmpSeq.zip euclidean*.png
```

**Custom Functions Used in R:**

* [transform_phyloseq()](#transform_phyloseq)
* [make_dendrogram()](#make_dendrogram)
* [run_stats()](#run_stats)
* [plot_pcoa()](#plot_pcoa)

**Parameter Definitions:**

**zip**
- `-q` – quiet operation
- `*_plots_<tech_type>_GLAmpSeq.zip` – positional argument naming the zip output file
- `*.png` – positional argument naming the input file(s) to package

**Input Data:**

* `rarefaction_depth` (an integer specifying the minimum number of reads to simulate during rarefaction)
* `groups_colname` (a string specifying the name of the column in the metadata table containing the group names)
* `legend_title` (a string specifying the legend title for plotting)
* `assay_suffix` (a string specifying the suffix to be added to output files; default is the Genelab assay suffix, "_GLAmpSeq")
* `output_prefix` (a string specifying an additional prefix to be added to the output files; default is no additional prefix, "")
* `normalization_methods` (a string vector specifying the method(s) to use for normalizing sample counts; "vst" (variance stabilizing transform) and "rarefy" (rarefaction) are supported)
* `distance_methods` (a string vector specifying the method(s) to use to calculate the distance between samples; "vst" transformed data uses "euclidean" (Euclidean distance) and "rarefy" transformed data uses "bray" (Bray-Curtis distance))
* `metadata` (a dataframe containing the sample metadata, with samples as row names and sample info as columns, output from [6b.iv. Read-in Input Tables](#6biv-read-in-input-tables))
* `feature_table` (a dataframe containing a filtered subset of the samples feature dataframe (i.e. ASV), output from [6b.v. Preprocessing](#6bv-preprocessing))
* `group_colors` (a named character vector of colors for each group, output from [6b.iv. Read-in Input Tables](#6biv-read-in-input-tables))

**Output Data:**
> NOTE: If beta diversity analysis using rarefaction couldn't be perfomed due to insufficient sequence counts per samples or not enough groups surviving rarefaction, a failure file (beta_diversity_failure_<tech_type>\_GLAmpSeq.txt) will be generated instead of the expected bray* output and the rarefaction_depth_<tech_type>_GLAmpSeq.txt. All euclidean* output files and vsd_validation_plot\_<tech_type>_GLAmpSeq.png will still be generated.

* **beta_diversity/<distance_method>\_plots_<tech_type>_GLAmpSeq.zip** (zip containing the following)
  * **<distance_method>\_dendrogram_<tech_type>_GLAmpSeq.png** (dendrogram(s) of the specified distance, Euclidean or Bray-Curtis, - based hierarchical clustering of the samples, colored by experimental groups)
  * **<distance_method>_PCoA_without_labels\_<tech_type>_GLAmpSeq.png** (Principle Coordinates Analysis plots of VST transformed and rarefy transformed ASV counts for Euclidean and Bray-Curtis distance methods, respectively, without sample labels)
  * **<distance_method>_PCoA_w_labels\_<tech_type>_GLAmpSeq.png** (Principle Coordinates Analysis plots of VST transformed and rarefy transformed ASV counts for Euclidean and Bray-Curtis distance methods, respectively, with sample labels)
* **beta_diversity/<distance_method>_adonis_table\_<tech_type>_GLAmpSeq.csv** (comma-separated table(s) containing the degrees of freedom (df), sum of squares (SumOfSqs), coefficient of determination (R^2), F-statistic (statistic), and p-value for the model (variation explained by experimental groups) and residual (unexplained variation) sources of variation (terms) for the specified distance analysis, Euclidean or Bray-Curtis)
* **beta_diversity/<distance_method>_variance_table\_<tech_type>_GLAmpSeq.csv** (comma-separated table(s) containing the degrees of freedom (df), sum of squares (sumsq), mean square (meansq), F-statistic (statistic), and p-value for the groups (variation explained by experimental groups) and residual (unexplained variation) sources of variation (terms) for the specified distance analysis, Euclidean or Bray-Curtis)
* **beta_diversity/vsd_validation_plot\_<tech_type>_GLAmpSeq.png** (VST transformation validation diagnostic plot)
* beta_diversity/rarefaction_depth\_<tech_type>_GLAmpSeq.txt (rarefaction depth value used in beta analysis)

<br>

---

## 9. Taxonomy Plots

Taxonomy summaries provide insights into the composition of microbial communities at various taxonomy levels.

```R
taxonomy_plots_out_dir <- "taxonomy_plots/"
if(!dir.exists(taxonomy_plots_out_dir)) dir.create(taxonomy_plots_out_dir)
metadata <- {DATAFRAME} 
feature_table <- {DATAFRAME}
taxonomy_table <- {DATAFRAME}
custom_palette <- {COLOR_VECTOR}
publication_format <- {GGPLOT_THEME}
groups_colname <- "groups"
assay_suffix <- "_GLAmpSeq"
output_prefix <- ""

# -------------------------Prepare feature tables -------------------------- #
# For ITS and 18S datasets the taxonomy columns may also contain kingdom and division taxonomy levels
# which will break the code. To avoid this, we only plot the phylum to species levels.
taxon_levels <- c("phylum", "class", "order", "family", "genus", "species") # Plot only phylum to species
names(taxon_levels) <- taxon_levels
taxon_tables <- map(.x = taxon_levels,
                    .f = make_feature_table,
                    count_matrix = feature_table,
                    taxonomy = taxonomy_table)

# ----------------------- Sample abundance plots -------------------------- #
group_rare <- TRUE
samples_order <- metadata %>% arrange(!!sym(groups_colname)) %>% rownames()
dont_group <- c("phylum")
# In percentage
# phylum 1%, class 3%, order 3%, family 8%, genus 8% and species 9%
thresholds <- c(phylum=1,class=3, order=3, family=8, genus=8, species=9)
# Convert from wide to long format
relAbundance_tbs_rare_grouped <- map2(.x = taxon_levels,
                                     .y = taxon_tables, 
                                     .f = function(taxon_level=.x,
                                                   taxon_table=.y){
                                      
                                       print(taxon_level)
                                       taxon_table <- apply(X = taxon_table, MARGIN = 2,
                                                            FUN = function(x) x/sum(x)) * 100
                                       
                                       
                                       taxon_table <- as.data.frame(taxon_table %>% t())
                                       if(group_rare && !(taxon_level %in% dont_group)){
                                         
                                         taxon_table <- group_low_abund_taxa(taxon_table %>%
                                                                               as.data.frame(check.names=FALSE,
                                                                                             stringAsFactor=FALSE),
                                                                             threshold = thresholds[taxon_level])
                                         
                                       }
                                       taxon_table$samples <- rownames(taxon_table)
                                       
                                       
                                       # Change data frame from wide to long format
                                       taxon_table <- taxon_table %>% 
                                         pivot_longer(cols = -samples, names_to = taxon_level, values_to = "relativeAbundance")
                                       taxon_table$samples <- factor(x = taxon_table$samples, 
                                                                     levels = samples_order)
                                       return(taxon_table)
                                     })

x_lab <- "Samples"
y_lab <- "Relative abundance (%)"
x <- 'samples'
y <- "relativeAbundance"
facet_by <- reformulate(groups_colname)
number_of_samples <- length(samples_order)


if(number_of_samples >=  30 ){

    plot_width <- 0.6 * number_of_samples

}else{
 
    plot_width <- 14
}

# Make sample plots
walk2(.x = relAbundance_tbs_rare_grouped, .y = taxon_levels, 
                           .f = function(relAbundance_tb, taxon_level){
                             
                             df <- relAbundance_tb %>%
                               left_join(metadata %>% rownames_to_column("samples"))
                             
                          p <- ggplot(data = df, mapping = aes(x= !!sym(x), y=!!sym(y) )) +
                               geom_col(aes(fill = !!sym(taxon_level) )) + 
                               facet_wrap(facet_by, scales = "free", 
                               nrow = 1, labeller = label_wrap_gen(width=10)) +
                               publication_format +
                               labs(x = x_lab , y = y_lab, fill= tools::toTitleCase(taxon_level)) + 
                               scale_fill_manual(values = custom_palette) +
                               theme(axis.text.x=element_text(
                                 margin=margin(t=0.5,r=0,b=0,l=0,unit ="cm"),
                                 angle = 90, 
                                 hjust = 0.5, vjust = 0.5)) + 
                               labs(x=NULL)
                          
                          ggsave(filename = glue("{taxonomy_plots_out_dir}/{output_prefix}samples_{taxon_level}_<tech_type>{assay_suffix}.png"),
                                 plot=p, width = plot_width, height = 8.5, dpi = 300, limitsize = FALSE)
                          
                           })

# ------------------------ Group abundance plots ----------------------------- #
# In percentage
# phylum 1% and 2% for class to species.
thresholds <- c(phylum=1,class=2, order=2, family=2, genus=2, species=2)

# Convert from wide to long format for every treatment group of interest
group_rare <- TRUE # should rare taxa be grouped ?
maximum_number_of_taxa <- 500 # If the number of taxa is more than this then rare taxa will be grouped anyway. 

group_relAbundance_tbs <- map2(.x = taxon_levels, .y = taxon_tables, 
                                     .f = function(taxon_level=.x, taxon_table=.y){
                                       
                                       taxon_table <- as.data.frame(taxon_table %>% t()) 
                                       taxon_table <- (collapse_samples(taxon_table = taxon_table,
                                                                        metadata = metadata, group = groups_colname,
                                                                        convertToRelativeAbundance = TRUE)$taxon_table * 100 ) %>%
                                         as.data.frame(check.names=FALSE)
                                       
                                       if(ncol(taxon_table) > maximum_number_of_taxa){
                                         group_rare <- TRUE
                                       }
                                       
                                       if(group_rare && !(taxon_level %in% dont_group)){
                                         taxon_table <- group_low_abund_taxa(taxon_table %>%
                                                                               as.data.frame(check.names=FALSE,
                                                                                             stringAsFactor=FALSE),
                                                                             threshold = thresholds[taxon_level])
                                         group_rare <- FALSE
                                       }
                                       
                                       taxon_table[,groups_colname] <- rownames(taxon_table)
                                       
                                       
                                       # Change from wide to long format
                                       taxon_table <- taxon_table %>% 
                                         pivot_longer(cols = -!!sym(groups_colname),
                                                      names_to = taxon_level,
                                                      values_to = "relativeAbundance")
                              
                                       return(taxon_table)
                                       
                                     })

# Make bar plots
y_lab <- "Relative abundance (%)"
y <- "relativeAbundance"
number_of_groups <- length(group_levels)
plot_width <- 2.5 * number_of_groups

# Cap the maximum plot width to 50 regardless of the number of groups
if(plot_width >  50 ){
  
  plot_width <- 50
}

walk2(.x = group_relAbundance_tbs, .y = taxon_levels, 
                           .f = function(relAbundance_tb=.x, taxon_level=.y){
                             
                             p <- ggplot(data =  relAbundance_tb %>%
                                           mutate(X = str_wrap(!!sym(groups_colname),
                                                             width = 10) # wrap long group names
                                                  ), 
                                        mapping = aes(x = X , y = !!sym(y))) +
                               geom_col(aes(fill = !!sym(taxon_level))) + 
                               publication_format +
                               theme(axis.text.x=element_text(
                                 margin=margin(t=0.5,r=0,b=0,l=0,unit ="cm"),
                                 angle = 0, 
                                 hjust = 0.5, vjust = 0.5)) + 
                               labs(x = NULL , y = y_lab, fill = tools::toTitleCase(taxon_level)) + 
                               scale_fill_manual(values = custom_palette)
                             ggsave(filename = glue("{taxonomy_plots_out_dir}/{output_prefix}groups_{taxon_level}_<tech_type>{assay_suffix}.png"),
                                    plot=p, width = plot_width, height = 10, dpi = 300, limitsize = FALSE)
                           })
```
```bash
zip -q sample_taxonomy_plots_<tech_type>_GLAmpSeq.zip samples_<tech_type>_GLAmpSeq.png

zip -q group_taxonomy_plots_<tech_type>_GLAmpSeq.zip groups_<tech_type>_GLAmpSeq.png
```

**Custom Functions Used in R**

* [make_feature_table()](#make_feature_table)
* [group_low_abund_taxa()](#group_low_abund_taxa)
* [collapse_samples()](#collapse_samples)

**Parameter Definitions:**

**zip**
- `-q` – quiet operation
- `*_taxonomy_plots_<tech_type>_GLAmpSeq.zip` – positional argument naming the zip output file
- `*.png` – positional argument naming the input file(s) to package

**Input Data:**

* `groups_colname` (a string specifying the name of the column in the metadata table containing the group names)
* `assay_suffix` (a string specifying the suffix to be added to output files; default is the Genelab assay suffix, "_GLAmpSeq")
* `output_prefix` (a string specifying an additional prefix to be added to the output files; default is no additional prefix, "")
* `metadata` (a dataframe containing the sample metadata, with samples as row names and sample info as columns, output from [6b.iv. Read-in Input Tables](#6biv-read-in-input-tables))
* `feature_table` (a dataframe containing a filtered subset of the samples feature dataframe (i.e. ASV), output from [6b.v. Preprocessing](#6bv-preprocessing))
* `taxonomy_table` (a dataframe containing a filtered subset of the feature taxonomy dataframe with ASV taxonomy assignments, output from [6b.v. Preprocessing](#6bv-preprocessing))
* `custom_palette` (a vector of strings specifying a custom color palette for coloring plots, output from [6b.iii. Set Variables](#6biii-set-variables))
* `publication_format` (a ggplot::theme object specifying the custom theme for plotting, output from [6b.iii. Set Variables](#6biii-set-variables))
* `group_levels` (a character vector of unique group names, output from [6b.iv. Read-in Input Tables](#6biv-read-in-input-tables))


**Output Data:**

* **taxonomy_plots/sample_taxonomy_plots_<tech_type>_GLAmpSeq.zip** (zip containing the following)
  * **samples_<taxon_level>_<tech_type>_GLAmpSeq.png** (barplots of the relative abundance of the specified taxon level for each sample)
* **taxonomy_plots/group_taxonomy_plots_<tech_type>_GLAmpSeq.zip** (zip containing the following)
  * **groups_<taxon_level>_<tech_type>_GLAmpSeq.png** (barplots of the relative abundance of the specified taxon level for each group)

Where `taxon_level` is all of phylum, class, order, family, genus, and species.

> Please note that the species plot can be misleading as short amplicon sequences can't be used to accurately predict species.

<br>

---


## 10. Differential Abundance Testing

Differential abundance testing aims to uncover specific taxa that exhibit notable variations across different conditions, complemented by visualizations like volcano plots to illustrate these disparities and their implications on ASV abundance and overall microbial community dynamics. ANCOMBC 1, ANCOMBC 2, and DESeq2 provide 3 different methods for calculating differential abundance. ANCOMBC 1 and 2 were specifically designed to handle the compositional nature of microbiome data. ANCOMBC 2 is an improved version of ANCOMBC 1 particularly for datasets with high sparsity, small sample sizes, or longitudinal and correlated experimental designs. This pipeline also implements DESeq2 because it is a popular choice for differential abundance. DESeq2 assumes a negative binomial model and can have issues with sparse data, which is frequently true in microbiome datasets. Two diagnostic plots (VST validation and ASV sparsity plots) help assess whether DESeq2 is appropriate for a given dataset. The VST validation plot assesses whether VST is successfully stabilizing variance. The ASV sparsity plot helps users determine if their data is too sparse to use DESeq2 to assess differential abundance.

> In general, we recommend using the ANCOMBC 2 differential abundance data, although all differential abundance data outputs should be evaluated in the context of the question(s) the user seeks to answer with these data.


### 10a. ANCOMBC 1

```R
# Create output directory if it doesn't already exist
diff_abund_out_dir <- "differential_abundance/ancombc1/"
if(!dir.exists(diff_abund_out_dir)) dir.create(diff_abund_out_dir, recursive = TRUE)
metadata <- {DATAFRAME}
feature_table <- {DATAFRAME}
taxonomy_table <- {DATAFRAME}
feature <- "ASV"
groups_colname <- "groups"
samples_column <- "Sample Name"
assay_suffix <- "_GLAmpSeq"
target_region <- "16S" # "16S", "18S" or "ITS"
output_prefix <- ""
prevalence_cutoff <- 0
library_cutoff <- 0
remove_struc_zero <- FALSE
threads <- 5

# # Get long asv taxonomy names and clean
species <- taxonomy_table %>%
  unite(species,domain:species,sep = ";") %>% 
pull %>% str_replace_all("Other", "_")

taxonomy_table <- fix_names(taxonomy_table, "Other", ";_")

taxonomy_table[,"species"] <- species

# Create phyloseq object from feature, taxonomy and sample metadata tables
ps <- phyloseq(otu_table(feature_table, taxa_are_rows = TRUE),
               sample_data(metadata),
               tax_table(as.matrix(taxonomy_table)))

# Convert phyloseq to tree summarized experiment object
tse <-  mia::makeTreeSummarizedExperimentFromPhyloseq(ps)

# Get unique group comparison as a matrix
pairwise_comp.m <- utils::combn((metadata[,group] %>% unique %>% sort), 2)
pairwise_comp_df <- pairwise_comp.m %>% as.data.frame 
# Name the columns in the pairwise matrix as group1vgroup2
colnames(pairwise_comp_df) <- map_chr(pairwise_comp_df,
                                      \(col) str_c(col, collapse = "v"))
comparisons <- colnames(pairwise_comp_df)
names(comparisons) <- comparisons


# ---------------------- Run ANCOMBC 1 ---------------------------------- #
set.seed(123)
final_results_bc1 <- map(pairwise_comp_df, function(col){
  
  group1 <- col[1]
  group2 <- col[2]
  
  # Subset the treeSummarizedExperiment object to contain only samples
  # in group1 and group2 
  tse_sub <- tse[, tse[[groups_colname]] %in% c(group1, group2)]
  
  # Note that by default, levels of a categorical variable in R are sorted 
  # alphabetically. 
  # Changing the reference group by reordering the factor levels
  tse_sub[[groups_colname]] <- factor(tse_sub[[groups_colname]] , levels = c(group1, group2))

  # Run ancombc (uses default parameters unless specified in pipeline Parameter Definitions)
  tryCatch({
    out <- ancombc(data = tse_sub, 
                    formula = groups_colname, 
                    p_adj_method = "fdr", prv_cut = prevalence_cutoff,
                    lib_cut = library_cutoff, 
                    group = groups_colname , struc_zero = remove_struc_zero,
                    neg_lb = TRUE, 
                    conserve = TRUE,
                    n_cl = threads, verbose = TRUE)
    
    # ------ Set data frame names ---------# 
    # lnFC 
    lfc <- out$res$lfc %>%
      as.data.frame() %>% 
      select(-contains("Intercept")) %>% 
      set_names(
        c("taxon",
          glue("Lnfc_({group2})v({group1})"))
      )
    
    # SE
    se <- out$res$se %>%
      as.data.frame() %>% 
      select(-contains("Intercept")) %>%
      set_names(
        c("taxon",
          glue("Lnfc.SE_({group2})v({group1})"))
      )
    
    # W    
    W <- out$res$W %>%
      as.data.frame() %>% 
      select(-contains("Intercept")) %>%
      set_names(
        c("taxon",
          glue("Stat_({group2})v({group1})"))
      )
    
    # p_val
    p_val <- out$res$p_val %>%
      as.data.frame() %>% 
      select(-contains("Intercept")) %>%
      set_names(
        c("taxon",
          glue("P.value_({group2})v({group1})"))
      )
    
    # q_val
    q_val <- out$res$q_val %>%
      as.data.frame() %>% 
      select(-contains("Intercept")) %>% 
      set_names(
        c("taxon",
          glue("Q.value_({group2})v({group1})"))
      )
    
    # Diff_abn
    diff_abn <- out$res$diff_abn %>%
      as.data.frame() %>% 
      select(-contains("Intercept")) %>%
      set_names(
        c("taxon",
          glue("Diff_({group2})v({group1})"))
      )
    
    # Merge the dataframes to one results dataframe
    res <- lfc %>%
      left_join(se) %>%
      left_join(W) %>% 
      left_join(p_val) %>% 
      left_join(q_val) %>% 
      left_join(diff_abn)
    
    return(res)
  }, error = function(e) {
    # Create log message
    log_msg <- c(
      "\nANCOMBC1 analysis failed for comparison: ", group1, " vs ", group2,
      "\nError: ", e$message,
      "\n\nDiagnostics:",
      paste("- Number of taxa after filtering:", nrow(taxonomy_table)),
      paste("- Number of samples in group", group1, ":", sum(tse_sub[[group]] == group1)),
      paste("- Number of samples in group", group2, ":", sum(tse_sub[[group]] == group2)),
      "\nPossibly insufficient data for ANCOMBC1 analysis. Consider adjusting filtering parameters or group assignments."
    )
    
    # Write to log file
    writeLines(log_msg, 
              file.path(diff_abund_out_dir, 
                       glue("{output_prefix}ancombc1_failure.txt")))
    
    # Print to console and quit
    message(log_msg)
    quit(status = 0)
  })
})

# ------------ Create merged stats pairwise dataframe ----------------- #
# Initialize the merged stats dataframe to contain the taxon column for joining
merged_stats_df <- final_results_bc1[[names(final_results_bc1)[1]]] %>%
  as.data.frame() %>% select(taxon)

# Loop over the results of every comparison and join it the pre-existing 
# stats table
walk(comparisons[names(final_results_bc1)], .f = function(comparison){
  
  # Get comparison specific statistics
  df <- final_results_bc1[[comparison]] %>% as.data.frame()
  
  # Merge it to the pre-existing statistics table
  merged_stats_df <<- merged_stats_df %>%
    dplyr::full_join(df, by = join_by("taxon"))
  
})

# Sort ASVs in ascending order
merged_stats_df <- merged_stats_df %>% 
  rename(!!feature := taxon) %>%
  mutate(!!feature := SortMixed(!!sym(feature)))

# ------ Get comparison names
# Since all column groups i.e. lnFC, pval, W, etc. have the same
# suffixes as comparison names, we only need to extract the comparison names
# from one of them. Here we extract them from the "lnFC" prefixed columns
comp_names <- merged_stats_df %>% 
  select(starts_with("Lnfc_", ignore.case = FALSE)) %>%
  colnames() %>% str_remove_all("Lnfc_")
names(comp_names) <- comp_names

# -------------- Make volcano plots ------------------ #
volcano_plots <- map(comp_names, function(comparison){
  
  # Construct column names for columns to be selected
  comp_col  <- c(
    glue("Lnfc_{comparison}"),
    glue("Lnfc.SE_{comparison}"),
    glue("Stat_{comparison}"),
    glue("P.value_{comparison}"),
    glue("Q.value_{comparison}"),
    glue("Diff_{comparison}")
  )

  sub_res_df <- merged_stats_df %>% 
    select(!!feature, all_of(comp_col)) %>% drop_na()
  colnames(sub_res_df) <- str_replace_all(colnames(sub_res_df),
                                          pattern = "(.+)_.+", 
                                          replacement = "\\1")

  # Set default pvalue and plot dimensions
  p_val <- 0.1
  plot_width_inches <- 11.1
  plot_height_inches <- 8.33
  
  # Retrieve a vector of the 2 groups being compared
  groups_vec <- comparison %>%
    str_replace_all("\\)v\\(", ").vs.(") %>%  # replace ')v(' with ').vs.(' to enhance accurate groups splitting 
    str_remove_all("\\(|\\)") %>%  # remove brackets
    str_split(".vs.") %>% unlist  # split groups to list then convert to a vector
  
  group1 <- groups_vec[1]
  group2 <- groups_vec[2]
  
  ###### Long x-axis label adjustments ##########
  x_label <- glue("ln Fold Change\n\n( {group1} vs {group2} )")
  label_length <- nchar(x_label)
  max_allowed_label_length <- plot_width_inches * 10
  
  # Construct x-axis label with new line breaks if was too long
  if (label_length > max_allowed_label_length){
    x_label <- glue("ln Fold Change\n\n( {group1} \n vs \n {group2} )")
  }

  # Make plot
  p <- ggplot(sub_res_df %>% mutate(diff = Q.value <= p_val), 
              aes(x=Lnfc, y=-log10(Q.value), 
                  color=diff, label=!!sym(feature))) +
    geom_point(alpha=0.7, size=2) +
    scale_color_manual(values=c("TRUE"="red", "FALSE"="black"),
                       labels=c(paste0("qval > ", p_val), 
                                paste0("qval \u2264 ", p_val))) +
    geom_hline(yintercept = -log10(p_val), linetype = "dashed") +
    ggrepel::geom_text_repel(show.legend = FALSE) +
    expandy(-log10(sub_res_df$Q.value)) + # Expand plot y-limit
    coord_cartesian(clip = 'off') +
    scale_y_continuous(oob = scales::oob_squish_infinite) + # prevent clipping of infinite values
    labs(x= x_label, y="-log10(Q-value)", 
         title = "Volcano Plot", color=NULL,
         caption = glue("dotted line: q-value = {p_val}")) + 
    theme_bw() +
    theme(legend.position="top", legend.key = element_rect(colour=NA),
          plot.caption = element_text(face = 'bold.italic'))
  # Save plot
  file_name <- glue("{output_prefix}{comparison %>% str_replace_all('[:space:]+','_')}_volcano_<tech_type>{assay.suffix}.png")
  ggsave(filename = file_name,
         plot = p, device = "png", width = plot_width_inches,
         height = plot_height_inches, units = "in",
         dpi = 300, path = diff_abund_out_dir)

  return(p)
})

# ------------------- Add NCBI id to feature, i.e. ASV -------------- #
# Get the best/least possible taxonomy name for the ASVs
tax_names <- map_chr(str_replace_all(taxonomy_table$species, ";_","") %>%
                       str_split(";"),
                     function(row) row[length(row)])

df <- data.frame(ASV=rownames(taxonomy_table), best_taxonomy=tax_names)
colnames(df) <- c(feature, "best_taxonomy")

# Pull NCBI IDS for unique taxonomy names
# Filter out unannotated entries before querying NCBI
valid_taxonomy <- df$best_taxonomy %>% unique() %>% setdiff("_")
df2_valid <- data.frame(best_taxonomy = valid_taxonomy) %>%
  mutate(NCBI_id=get_ncbi_ids(best_taxonomy, target_region),
         .after = best_taxonomy)

# Add unannotated entries with NA NCBI_id
df2_invalid <- data.frame(best_taxonomy = "_", NCBI_id = NA)
df2 <- rbind(df2_valid, df2_invalid)

df <- df %>%
  left_join(df2, join_by("best_taxonomy")) %>% 
  right_join(merged_stats_df)

# Manually creating a normalized table because normalized 
# tables differ by comparison
normalized_table <- as.data.frame(feature_table + 1) %>%
  rownames_to_column(feature) %>%
  mutate(across( where(is.numeric), log ) )

# Create a missing values / NAs dataframe of samples that were dropped
# due to prefiltering steps (prevalence and library cut offs filtering)
# proir to running ANCOMBC
samples <- metadata[[samples_column]]
samplesdropped <- setdiff(x = samples, y = colnames(normalized_table)[-1])
missing_df <- data.frame(ASV=normalized_table[[feature]],
                         matrix(data = NA, 
                                nrow = nrow(normalized_table),
                                ncol = length(samplesdropped)
                         )
)
colnames(missing_df) <- c(feature,samplesdropped)

# Create mean and standard deviation table
group_levels <- metadata[, groups_colname] %>% unique() %>% sort()
group_means_df <- normalized_table[feature]
walk(group_levels, function(group_level){

  mean_col <- glue("Group.Mean_({group_level})")
  std_col <- glue("Group.Stdev_({group_level})")
  
  # Samples that belong to the current group
  Samples <- metadata %>%
    filter(!!sym(groups_colname) == group_level) %>%
    pull(!!sym(samples_column))
  # Samples that belong to the current group that are in the normalized table
  Samples <- intersect(colnames(normalized_table), Samples)
  
  # Calculate means and standard deviations for the current comparison
  temp_df <- normalized_table %>% select(!!feature, all_of(Samples)) %>% 
    rowwise() %>%
    mutate(!!mean_col := mean(c_across(where(is.numeric)), na.rm = TRUE),
           !!std_col := sd(c_across(where(is.numeric)), na.rm = TRUE) ) %>% 
    select(!!feature,!!sym(mean_col), !!sym(std_col))
  
  # Merge the current comparison's means and stdandard deviations
  # to previous ones
  group_means_df <<- group_means_df %>% left_join(temp_df)
  
})

# Append missing sample columns to the normalized table
normalized_table <- normalized_table %>%
  left_join(missing_df, by = feature) %>%
  select(!!feature, all_of(samples))

# Compute globally/ASV normalized means and standard deviations 
All_mean_sd <- normalized_table %>%
  rowwise() %>%
  mutate(All.mean=mean(c_across(where(is.numeric)), na.rm = TRUE),
         All.stdev=sd(c_across(where(is.numeric)), na.rm = TRUE) ) %>% 
  select(!!feature, All.mean, All.stdev)

# Merge the taxonomy table to the stats table
merged_df <- df %>%
  left_join(taxonomy_table %>%
              as.data.frame() %>%
              rownames_to_column(feature)) %>% 
  select(!!feature, domain:species,everything())

# Merge all the pre-combined dataframes in the desired order
merged_df <- merged_df %>%
  select(!!sym(feature):NCBI_id) %>%
  left_join(normalized_table, by = feature) %>%
  left_join(merged_df) %>%
  left_join(All_mean_sd) %>%
  left_join(group_means_df, by = feature) %>%
  mutate(across(where(is.matrix), as.numeric))

# Write out results of differential abundance using ANCOMBC 1
output_file <- glue("{diff_abund_out_dir}/{output_prefix}ancombc1_differential_abundance_<tech_type>{assay_suffix}.csv")
# Write combined table to file but before that drop
# all columns of inferred differential abundance by ANCOMBC 
write_csv(merged_df %>%
            select(-starts_with("Diff_")),
          output_file)

```
```bash
zip -q ancombc1_volcano_plots_<tech_type>_GLAmpSeq.zip *_volcano_<tech_type>_GLAmpSeq.png
```

**Custom Functions Used:**

* [expandy()](#expandy)
* [get_ncbi_ids()](#get_ncbi_ids)
* [fix_names()](#fix_names)
  
**Parameter Definitions:**

**R script**
* `ancombc()` - ANCOMBC::ancombc function (*using the following non-default values:*)
  * `data` - TreeSummarizedExperiment object created from `feature_table` input data
  * `formula` - a string specifying the variable in the metadata to use for the fixed effects formula (e.g. group names), set by `groups_colname` input data
  * `prv_cut` - fraction between 0 and 1 specifying the taxon prevalence cut-off, set by `prevalence_cutoff` input data
  * `lib_cut` - a numerical threshold for filtering samples based on library sizes, set by `library_cutoff` input data
  * `group` - the name of the group variable in the metadata, set by `groups_colname` input data
  * `struc_zero` - logical value indicating whether or not group-wise structural zeros should be detected, set by `remove_struc_zero` input data
  * `n_cl` - the number of processes to run in parallel, set by `threads` input data
  * `p_adj_method` - a string specifying the p-value adjustment method for multiple comparisons testing, set to "fdr" to standardize the multiple comparisons method across all three differential abundance methods
  * `neg_lb` - logical value specifying whether to classify a taxon as a structural zero using its asymptotic lower bound, set to "TRUE"
  * `conserve` - logical value indicating where or not a conservative variance estimator should be used for the test statistic, set to "TRUE"
  * `verbose` - logical value specifying whether or not to generate verbose output

**zip**
- `-q` – quiet operation
- `ancombc1_volcano_plots_<tech_type>_GLAmpSeq.zip` – positional argument naming the zip output file
- `*_volcano_<tech_type>_GLAmpSeq.png` – positional argument naming the input file(s) to package


**Input Data:**

* `feature` (a string specifying the feature type, i.e. "ASV" or "OTU")
* `groups_colname` (a string specifying the name of the column in the metadata table containing the group names)
* `samples_column` (a string specifying the name of the column in the metadata table containing the sample names)
* `assay_suffix` (a string specifying the suffix to be added to output files; default is the Genelab assay suffix, "_GLAmpSeq")
* `output_prefix` (a string specifying an additional prefix to be added to the output files; default is no additional prefix, "")
* `threads` (a number specifying the number of cpus to use for parallel processing)
* `prevalence_cutoff` (a decimal between 0 and 1 specifying the proportion of samples required to contain a taxon in order to keep the taxon when `remove_rare` (set in [Step 6b.iv. Preprocessing](#6bv-preprocessing)) is set to TRUE; default is 0, i.e. do not exclude any taxon / feature)
* `library_cutoff` (a numerical variable specifying the number of total counts a sample must have across all features to be retained when `remove_rare` (set in [Step 6b.iv. Preprocessing](#6bv-preprocessing)) is set to TRUE; default is 0, i.e. no samples will be dropped)
* `target_region` (a string specifying the amplicon target region; options are either "16S", "18S", or "ITS")
* `remove_struc_zero` (a boolean variable specifying whether or not structural zeros (a.k.a ASVs with zero count in at least one group) should be removed; default is FALSE i.e. structural zeros won't be removed)
* `metadata` (a dataframe containing the sample metadata, with samples as row names and sample info as columns, output from [6b.iv. Read-in Input Tables](#6biv-read-in-input-tables))
* `feature_table` (a dataframe containing a filtered subset of the samples feature dataframe (i.e. ASV), output from [6b.v. Preprocessing](#6bv-preprocessing))
* `taxonomy_table` (a dataframe containing a filtered subset of the feature taxonomy dataframe with ASV taxonomy assignments, output from [6b.v. Preprocessing](#6bv-preprocessing))

<br>

**Output Data:**

* **differential_abundance/ancombc1/ancombc1_volcano_plots_<tech_type>_GLAmpSeq.zip** (zip containing the following)
  * **(\<group1\>)v(\<group2\>)\_volcano_<tech_type>_GLAmpSeq.png** (volcano plots for each pariwise comparison)
* **differential_abundance/ancombc1/ancombc1_differential_abundance_<tech_type>_GLAmpSeq.csv** (a comma-separated ANCOM-BC1 differential abundance results table containing the following columns:
  - ASV (identified ASVs)
  - taxonomic assignment columns
  - NCBI identifier for the best taxonomic assignment for each ASV 
  - Normalized abundance values for each ASV for each sample
  - For each pairwise group comparison:
    - natural log of the fold change (Lnfc)
    - standard error for the lnFC (Lnfc.SE)
    - test statistic from the primary result (Stat)
    - P-value (P.value)
    - Adjusted p-value (Q.value)
  - All.mean (mean across all samples)
  - All.stdev (standard deviation across all samples) 
  - For each group:
    - Group.Mean_(group) (mean within group)
    - Group.Stdev_(group) (standard deviation within group))

<br>

---

### 10b. ANCOMBC 2

```R
diff_abund_out_dir <- "differential_abundance/ancombc2/"
if(!dir.exists(diff_abund_out_dir)) dir.create(diff_abund_out_dir, recursive = TRUE)
metadata <- {DATAFRAME} 
feature_table <- {DATAFRAME}
taxonomy_table <- {DATAFRAME}
feature <- "ASV"
target_region <- "16S" # "16S" , "18S" or "ITS"
groups_colname <- "groups"
samples_column <- "Sample Name"
assay_suffix <- "_GLAmpSeq"
output_prefix <- ""
prevalence_cutoff <- 0  # from [Step 6b.v. Preprocessing]
library_cutoff <- 0  # from [Step 6b.v. Preprocessing]
remove_struc_zero <- FALSE
threads <- 5

# # Get long asv taxonomy names and clean
species <- taxonomy_table %>%
  unite(species,domain:species,sep = ";") %>% 
pull %>% str_replace_all("Other", "_")

taxonomy_table <- fix_names(taxonomy_table, "Other", ";_")

taxonomy_table[,"species"] <- species

# Create phyloseq object
ps <- phyloseq(otu_table(feature_table, taxa_are_rows = TRUE),
               sample_data(metadata),
               tax_table(as.matrix(taxonomy_table)))

# Convert phyloseq to tree summarized experiment object
tse <- mia::makeTreeSummarizedExperimentFromPhyloseq(ps)

# Getting the reference group and making sure that it is the reference 
# used in the analysis
group_levels <- metadata[, groups_colname] %>% unique() %>% sort()
ref_group <- group_levels[1] # the first group is used as the reference group by default
tse[[groups_colname]] <- factor(tse[[groups_colname]] , levels = group_levels)


# ---------------------- Run ANCOMBC2 ---------------------------------- #
# Run ancombc2 (uses default parameters unless specified in pipeline Parameter Definitions)
output <- ancombc2(data = tse,
                   fix_formula = groups_colname,
                   p_adj_method = "fdr",
                   prv_cut = prevalence_cutoff, 
                   lib_cut = library_cutoff, s0_perc = 0.05,
                   group = groups_colname, struc_zero = remove_struc_zero,
                   n_cl = threads, verbose = TRUE,
                   pairwise = TRUE, 
                   iter_control = list(tol = 1e-5, max_iter = 20,
                                       verbose = FALSE),
                   mdfdr_control = list(fwer_ctrl_method = "fdr", B = 100), 
                   lme_control = NULL, trend_control = NULL)

# For 2-group comparisons, use res instead of mapping across pairwise results in res_pair
is_two_group <- length(unique(tse[[group]])) == 2

# Create new column names - the original column names given by ANCOMBC are
# difficult to understand
tryCatch({
  # Check if this is a 2-group comparison (using res instead of res_pair)
  if(is_two_group) {
    # For 2-group comparisons, use the group-specific columns
    group_cols <- colnames(output$res)[grepl(paste0("^[a-zA-Z_]+_", group), colnames(output$res))]
    if(length(group_cols) > 0) {
      # Extract group name from the first group-specific column
      group_name <- str_replace(group_cols[1], paste0("^[a-zA-Z_]+_", group), "")
      # Create comparison name
      comparison_name <- glue("({group_name})v({ref_group})")
      
      new_colnames <- c(
        feature,  # Keep the feature column name
        glue("Lnfc_{comparison_name}"),
        glue("Lnfc.SE_{comparison_name}"),
        glue("Stat_{comparison_name}"),
        glue("P.value_{comparison_name}"),
        glue("Q.value_{comparison_name}"),
        glue("Diff_{comparison_name}"),
        glue("Passed_ss_{comparison_name}")
      )
    } else {
      stop("Could not identify group-specific column for 2-group comparison")
    }
  } else {
    # Multi-group comparisons
    new_colnames <- map_chr(output$res_pair  %>% colnames, 
                            function(colname) {
                              # Columns comparing a group to the reference group
                              if(str_count(colname,group) == 1){
                                str_replace_all(string=colname, 
                                                pattern=glue("(.+)_{group}(.+)"),
                                                replacement=glue("\\1_(\\2)v({ref_group})")) %>% 
                                str_replace(pattern = "^lfc_", replacement = "Lnfc_") %>% 
                                str_replace(pattern = "^se_", replacement = "Lnfc.SE_") %>% 
                                str_replace(pattern = "^W_", replacement = "Stat_") %>%
                                str_replace(pattern = "^p_", replacement = "P.value_") %>%
                                str_replace(pattern = "^q_", replacement = "Q.value_") %>%
                                str_replace(pattern = "^diff_", replacement = "Diff_") %>%
                                str_replace(pattern = "^passed_ss_", replacement = "Passed_ss_")
                                
                              # Columns with normal two groups comparison
                              } else if(str_count(colname,group) == 2){
                                
                                str_replace_all(string=colname, 
                                                pattern=glue("(.+)_{group}(.+)_{group}(.+)"),
                                                replacement=glue("\\1_(\\2)v(\\3)")) %>% 
                                str_replace(pattern = "^lfc_", replacement = "Lnfc_") %>% 
                                str_replace(pattern = "^se_", replacement = "Lnfc.SE_") %>% 
                                str_replace(pattern = "^W_", replacement = "Stat_") %>%
                                str_replace(pattern = "^p_", replacement = "P.value_") %>%
                                str_replace(pattern = "^q_", replacement = "Q.value_") %>%
                                str_replace(pattern = "^diff_", replacement = "Diff_") %>%
                                str_replace(pattern = "^passed_ss_", replacement = "Passed_ss_")
                                
                                # Feature/ ASV column 
                              } else{
                                
                                return(colname)
                              }
                            } )
  }
}, error = function(e) {
  writeLines(c("ANCOMBC2 script failed at res_pair processing:", e$message,
              "\n\nDiagnostics:",
              paste("- Number of taxa after filtering:", nrow(taxonomy_table)),
              paste("- Number of groups:", length(unique(tse[[group]]))),
              paste("- Sample sizes per group:"),
              paste("  ", paste(names(table(tse[[group]])), "=", table(tse[[group]]), collapse="\n  ")),
              "\nPossibly insufficient data for ANCOMBC2 analysis. Consider adjusting filtering parameters or group assignments."), 
            file.path(diff_abund_out_dir, glue("{output_prefix}ancombc2_failure.txt")))
  quit(status = 0)
})

# Change the column named taxon to the feature name e.g. ASV
new_colnames[match("taxon", new_colnames)] <- feature


# Rename columns
if(is_two_group) {
  # For 2-group comparisons, we need to select the group-specific columns and rename them
  # The columns are named like "lfc_groupsGround Control", "se_groupsGround Control", etc.
  
  group_specific_cols <- colnames(output$res)[grepl(paste0("^[a-zA-Z_]+_", group), colnames(output$res))]
  
  # Create a new data frame with the selected columns
  paired_stats_df <- output$res %>%
    select(taxon, all_of(group_specific_cols)) %>%
    set_names(new_colnames)
} else {
  # Multi-group comparisons
  paired_stats_df <- output$res_pair  %>%  set_names(new_colnames)
}

# Get the unique comparison names 
uniq_comps <- str_replace_all(new_colnames, ".+_(\\(.+\\))", "\\1") %>% unique()
uniq_comps <- uniq_comps[-match(feature, uniq_comps)]

# ------ Sort columns by group comparisons --------#
# Create a data frame containing only the feature/ASV column
res_df <- paired_stats_df[1] 
walk(uniq_comps, function(comp){
  
  # Get the results for a comparison
  temp_df <- paired_stats_df %>% select(!!sym(feature), contains(comp))
  
  # Merge the current comparison to previous comparisons by feature/ASV id
  res_df <<- res_df %>% left_join(temp_df)
})

# --------- Add NCBI id to feature ---------------#

# Get the best taxonomy assigned to each ASV
tax_names <- map_chr(str_replace_all(taxonomy_table$species, ";_","") %>%
                       str_split(";"),
                     function(row) row[length(row)])

df <- data.frame(ASV=rownames(taxonomy_table), best_taxonomy=tax_names)
colnames(df) <- c(feature, "best_taxonomy")

# Querying NCBI...
# Pull NCBI IDS for unique taxonomy names
# Filter out unannotated entries before querying NCBI
valid_taxonomy <- df$best_taxonomy %>% unique() %>% setdiff("_")
df2_valid <- data.frame(best_taxonomy = valid_taxonomy) %>%
  mutate(NCBI_id=get_ncbi_ids(best_taxonomy, target_region),
         .after = best_taxonomy)

# Add unannotated entries with NA NCBI_id
df2_invalid <- data.frame(best_taxonomy = "_", NCBI_id = NA)
df2 <- rbind(df2_valid, df2_invalid)

df <- df %>%
  left_join(df2, join_by("best_taxonomy")) %>% 
  right_join(res_df)

# Retrieve the normalized table
normalized_table <- output$bias_correct_log_table %>%
  rownames_to_column(feature) %>%
  mutate(across(where(is.numeric), ~replace_na(.x, replace=0)))

# Create a missing values / NAs dataframe of samples that were dropped
# due to prefiltering steps (prevalence and library cut offs filtering)
# proir to running ANCOMBC2
samples <- metadata[[samples_column]]
samplesdropped <- setdiff(x = samples, y = colnames(normalized_table)[-1])
missing_df <- data.frame(ASV=normalized_table[[feature]],
           matrix(data = NA, 
                  nrow = nrow(normalized_table),
                  ncol = length(samplesdropped)
                  )
           )
colnames(missing_df) <- c(feature, samplesdropped)

group_means_df <- normalized_table[feature]
walk(group_levels, function(group_level){

  mean_col <- glue("Group.Mean_({group_level})")
  std_col <- glue("Group.Stdev_({group_level})")
  
  # Samples that belong to the current group
  Samples <- metadata %>%
    filter(!!sym(groups_colname) == group_level) %>%
    pull(!!sym(samples_column))
  # Samples that belong to the current group that are in the normalized table
  Samples <- intersect(colnames(normalized_table), Samples)
  
  temp_df <- normalized_table %>% select(!!feature, all_of(Samples)) %>% 
    rowwise() %>%
    mutate(!!mean_col := mean(c_across(where(is.numeric)), na.rm = TRUE),
           !!std_col := sd(c_across(where(is.numeric)), na.rm = TRUE) ) %>% 
    select(!!feature,!!sym(mean_col), !!sym(std_col))
  
  group_means_df <<- group_means_df %>% left_join(temp_df)
  
})

# Append missing samples columns to normalized table
normalized_table <- normalized_table %>%
  left_join(missing_df, by = feature) %>% 
  select(!!feature, all_of(samples))

# Calculate global mean and standard deviation
All_mean_sd <- normalized_table %>%
  rowwise() %>%
  mutate(All.mean=mean(c_across(where(is.numeric)), na.rm = TRUE),
         All.stdev=sd(c_across(where(is.numeric)), na.rm = TRUE) ) %>% 
  select(!!feature, All.mean, All.stdev)

# Append the taxonomy table to the ncbi and stats table
merged_df <- df %>%
  left_join(taxonomy_table %>%
              as.data.frame() %>%
              rownames_to_column(feature)) %>% 
  select(!!feature,domain:species,everything())

# Combine tables in the desired order
merged_df <- merged_df %>%
  select(!!sym(feature):NCBI_id) %>%
  left_join(normalized_table, by = feature) %>%
  left_join(merged_df) %>% 
  left_join(All_mean_sd) %>% 
  left_join(group_means_df, by = feature)

# Writing out results of differential abundance using ANCOMBC2...
output_file <- glue("{diff_abund_out_dir}{output_prefix}ancombc2_differential_abundance_<tech_type>{assay_suffix}.csv")
# Write out merged stats table but before that 
# drop ANCOMBC inferred differential abundance columns
write_csv(merged_df %>%
            select(-starts_with("diff_")),
          output_file)

# ---------------------- Visualization --------------------------------------- #
# ------------ Make volcano ---------------- #
volcano_plots <- map(uniq_comps, function(comparison){
  
  comp_col <- c(
    glue("Lnfc_{comparison}"),
    glue("Lnfc.SE_{comparison}"),
    glue("Stat_{comparison}"),
    glue("P.value_{comparison}"),
    glue("Q.value_{comparison}"),
    glue("Diff_{comparison}"),
    glue("Passed_ss_{comparison}")
  )
  
  sub_res_df <- res_df %>% 
    select(!!feature, all_of(comp_col))
  colnames(sub_res_df) <- str_replace_all(colnames(sub_res_df),
                                          pattern = "(.+)_.+", 
                                          replacement = "\\1")
  # Set default qvalue and plot dimensions.
  p_val <- 0.1
  plot_width_inches <- 11.1
  plot_height_inches <- 8.33
  
  groups_vec <- comparison %>%
    str_replace_all("\\)v\\(", ").vs.(") %>% 
    str_remove_all("\\(|\\)") %>%
    str_split(".vs.") %>% unlist
  
  group1 <- groups_vec[1]
  group2 <- groups_vec[2]
  
  ######Long x-axis label adjustments##########
  x_label <- glue("ln Fold Change\n\n( {group1} vs {group2} )")
  label_length <- nchar(x_label)
  max_allowed_label_length <- plot_width_inches * 10
  
  # Construct x-axis label with new line breaks if was too long
  if (label_length > max_allowed_label_length){
    x_label <- glue("ln Fold Change\n\n( {group1} \n vs \n {group2} )")
  }
  #######################################
  
  p <- ggplot(sub_res_df %>% mutate(diff = Q.value <= p_val),
              aes(x=Lnfc, y=-log10(Q.value), color=diff, label=!!sym(feature))) +
    geom_point(alpha=0.7, size=2) +
    scale_color_manual(values=c("TRUE"="red", "FALSE"="black"),
                       labels=c(paste0("qval > ", p_val), 
                                paste0("qval \u2264 ", p_val))) +
    geom_hline(yintercept = -log10(p_val), linetype = "dashed") +
    ggrepel::geom_text_repel(show.legend = FALSE) + 
    expandy(-log10(sub_res_df$Q.value)) + # Expand plot y-limit
    coord_cartesian(clip = 'off') +
    scale_y_continuous(oob = scales::oob_squish_infinite) + # prevent clipping of infinite values
    labs(x= x_label, y="-log10(Q-value)", 
         title = "Volcano Plot", color=NULL,
         caption = glue("dotted line: q-value = {p_val}")) + 
    theme_bw() +
    theme(legend.position="top", legend.key = element_rect(colour=NA),
          plot.caption = element_text(face = 'bold.italic'))

  # Save plot
  file_name <-  glue("{output_prefix}{comparison %>% str_replace_all('[:space:]+','_')}_volcano<_tech_type>{assay.suffix}.png")
  ggsave(filename = file_name,
         plot = p, device = "png",
         width = plot_width_inches,
         height = plot_height_inches,
         units = "in", dpi = 300, path = diff_abund_out_dir)
  
  return(p)
})
```
```bash
zip -q ancombc2_volcano_plots_<tech_type>_GLAmpSeq.zip *_volcano_<tech_type>_GLAmpSeq.png
```

**Custom Functions Used in R**

* [expandy()](#expandy)
* [get_ncbi_ids()](#get_ncbi_ids)
* [fix_names()](#fix_names)
* [ancombc2()](#ancombc2) (the wrapper function that calls ANCOMBC::ancombc2)

**Parameter Definitions:**

**R script**
* `ancombc2()` - ANCOMBC::ancombc2 function (*pipeline uses default values unless defined below*)
  * `data` - a TreeSummarizedExperiment object created from `feature_table` input data
  * `fix_formula` - a string specifying the variable in the metadata to use for the fixed effects formula (e.g. group names), set by `groups_colname` input data
  * `prv_cut` - fraction between 0 and 1 specifying the taxon prevalence cut-off, set by `prevalence_cutoff` input data
  * `lib_cut` - a numerical threshold for filtering samples based on library sizes, set by `library_cutoff` input data
  * `group` - the name of the group variable in the metadata, set by `groups_colname` input data
  * `struc_zero` - logical value indicating whether or not group-wise structural zeros should be detected, set by `remove_struc_zero` input data
  * `n_cl` - specifies the number of processes to run in parallel, set to `threads` input data
  * `p_adj_method` - a string specifying the p-value adjustment method for multiple comparisons testing, set to "fdr" to standardize the multiple comparisons method across all three differential abundance methods.
  * `pairwise` - logical value indicating whether or not to perform the pairwise directional test, set to "TRUE" to compute all pairwise comparisons
  * `iter_control` - a named list of control parameters for the iterative MLE or RMEL algorithm
    * `tol` - iteration convergence tolerance, set to "1e-5" to match ancombc
  * `mdfdr_control` - a named list of control parameters for mixed directional false discovery rate (mdFDR)
    * `fwer_ctrl_method` - family-wise error controlling procedure, set to 'fdr' to match p_adj_method
  * `lme_control` - a named list of control parameters for mixed model fitting, set to 'NULL' to disable
  * `verbose` - logical value specifying whether or not to generate verbose output

**zip**
- `-q` – quiet operation
- `ancombc2_volcano_plots_<tech_type>_GLAmpSeq.zip` – positional argument naming the zip output file
- `*_volcano_<tech_type>_GLAmpSeq.png` – positional argument naming the input file(s) to package

**Input Data:**

* `feature` (a string specifying the feature type, i.e. "ASV" or "OTU")
* `groups_colname` (a string specifying the name of the column in the metadata table containing the group names)
* `samples_column` (a string specifying the name of the column in the metadata table containing the sample names)
* `assay_suffix` (a string specifying the suffix to be added to output files; default is the Genelab assay suffix, "_GLAmpSeq")
* `output_prefix` (a string specifying an additional prefix to be added to the output files; default is no additional prefix, "")
* `threads` (a number specifying the number of cpus to use for parallel processing)
* `prevalence_cutoff` (a decimal between 0 and 1 specifying the proportion of samples required to contain a taxon in order to keep the taxon when `remove_rare` (set in [Step 6b.iv. Preprocessing](#6bv-preprocessing)) is set to TRUE; default is 0, i.e. do not exclude any taxon/feature)
* `library_cutoff` (a numerical value specifying the number of total counts a sample must have across all features to be retained when `remove_rare` (set in [Step 6b.iv. Preprocessing](#6bv-preprocessing)) is set to TRUE; default is 0, i.e. no samples will be dropped)
* `target_region` (a string specifying the amplicon target region; options are either "16S", "18S", or "ITS")
* `remove_struc_zero` (a boolean value specifying whether or not structural zeros (a.k.a ASVs with zero count in at least one group) should be removed; default is FALSE i.e. structural zeros won't be removed)
* `metadata` (a dataframe containing the sample metadata, with samples as row names and sample info as columns, output from [6b.iv. Read-in Input Tables](#6biv-read-in-input-tables))
* `feature_table` (a dataframe containing a filtered subset of the samples feature dataframe (i.e. ASV), output from [6b.v. Preprocessing](#6bv-preprocessing))
* `taxonomy_table` (a dataframe containing a filtered subset of the feature taxonomy dataframe with ASV taxonomy assignments, output from [6b.v. Preprocessing](#6bv-preprocessing))

**Output Data:**

* **differential_abundance/ancombc2/ancombc2_volcano_plots_<tech_type>_GLAmpSeq.zip** (zip containing the following)
  * **(\<group1\>)v(\<group2\>)\_volcano_<tech_type>_GLAmpSeq.png** (volcano plots for each pariwise comparison)
* **differential_abundance/ancombc2/ancombc2_differential_abundance_<tech_type>_GLAmpSeq.csv** (a comma-separated ANCOM-BC2 differential abundance results table containing the following columns:
  - ASV (identified ASVs)
  - taxonomic assignment columns
  - NCBI identifier for the best taxonomic assignment for each ASV 
  - Normalized abundance values for each ASV for each sample
  - For each pairwise group comparison:
    - natural log of the fold change (Lnfc)
    - standard error for the lnFC (Lnfc.SE)
    - test statistic from the primary result (Stat)
    - P-value (P.value)
    - Adjusted p-value (Q.value)
  - All.mean (mean across all samples)
  - All.stdev (standard deviation across all samples) 
  - For each group:
    - Group.Mean_(group) (mean within group)
    - Group.Stdev_(group) (standard deviation within group))

<br>

---

### 10c. DESeq2

```R
# Create output directory if it doesn't already exist
diff_abund_out_dir <- "differential_abundance/deseq2/"
if(!dir.exists(diff_abund_out_dir)) dir.create(diff_abund_out_dir, recursive = TRUE)
metadata <- {DATAFRAME}
feature_table <- {DATAFRAME}
taxonomy_table <- {DATAFRAME}
feature <- "ASV"
samples_column <- "Sample Name"
groups_colname <- "groups" 
assay_suffix <- "_GLAmpSeq"
target_region <- "16S" # "16S", "18S" or "ITS"
output_prefix <- ""

# Get long asv taxonomy names and clean
species <- taxonomy_table %>%
  unite(species,domain:species,sep = ";") %>% 
pull %>% str_replace_all("Other", "_")

taxonomy_table <- fix_names(taxonomy_table, "Other", ";_")

taxonomy_table[,"species"] <- species

# Create phyloseq object from the feature, taxonomy and metadata tables 
ASV_physeq <- phyloseq(otu_table(feature_table, taxa_are_rows = TRUE),
                       tax_table(as.matrix(taxonomy_table)),
                       sample_data(metadata))
# Convert the phyloseq object to a deseq object 
deseq_obj <- phyloseq_to_deseq2(physeq = ASV_physeq,
                                design = reformulate(groups_colname))

# Add pseudocount if any 0 count samples are present
if (sum(colSums(counts(deseq_obj)) == 0) > 0) {
  # Add pseudo count of 1
  count_data <- counts(deseq_obj) + 1 
  # Make a columns of integer type
  count_data <- as.matrix(apply(count_data, 2, as.integer))
  rownames(count_data) <- rownames(counts(deseq_obj))
  colnames(count_data) <- colnames(counts(deseq_obj))
  counts(deseq_obj) <- count_data
}

# ---------------------- Run DESeq ---------------------------------- #
# https://rdrr.io/bioc/phyloseq/src/inst/doc/phyloseq-mixture-models.R 
deseq_modeled <- tryCatch({
  # Attempt to run DESeq, if error occurs then attempt an alternative 
  # size factor estimation method
  DESeq(deseq_obj)
}, error = function(e) {
  message("Error encountered in DESeq, applying alternative \
          method for size factor estimation...")
  
  geoMeans <- apply(counts(deseq_obj), 1, gm_mean)
  
  # Apply the alternative size factor estimation method
  deseq_obj <- estimateSizeFactors(deseq_obj, geoMeans=geoMeans)
  
  # Call DESeq again with alternative geom mean size est
  tryCatch({
    DESeq(deseq_obj)
  }, error = function(e2) {

    writeLines(c("Error:", e2$message,
                "\nUsing gene-wise estimates as final estimates instead of standard curve fitting."), 
              file.path(diff_abund_out_dir, glue("{output_prefix}deseq2_warning.txt")))
    
    # Use gene-wise estimates as final estimates
    deseq_obj <- estimateDispersionsGeneEst(deseq_obj)
    dispersions(deseq_obj) <- mcols(deseq_obj)$dispGeneEst
    # Continue with testing using nbinomWaldTest
    nbinomWaldTest(deseq_obj)
  })
})


# Make ASV Sparsity plot
sparsity_plot <- plotSparsity(deseq_modeled) 
ggsave(filename = glue("{diff_abund_out_dir}/{output_prefix}asv_sparsity_plot_<tech_type>{assay.suffix}.png"),
       plot = sparsity_plot, width = 14, height = 10, dpi = 300, units = "in")

# Get unique group comparison as a matrix
pairwise_comp.m <- utils::combn((metadata[,group] %>% unique %>% sort), 2)
pairwise_comp_df <- pairwise_comp.m %>% as.data.frame 
# Set the colnames as group1vgroup2
colnames(pairwise_comp_df) <- map_chr(pairwise_comp_df,
                                      \(col) str_c(col, collapse = "v"))
comparisons <- colnames(pairwise_comp_df)
names(comparisons) <- comparisons

# Retrieve statistics table
merged_stats_df <- data.frame(ASV=rownames(feature_table))
colnames(merged_stats_df) <- feature

walk(pairwise_comp_df, function(col){
  
  group1 <- col[1]
  group2 <- col[2]

# Retrieve the statistics table for the current pair and rename the columns
df <- results(deseq_modeled, contrast = c(group, group2, group1)) %>%
  data.frame() %>%
  rownames_to_column(feature) %>% 
  set_names(c(feature ,
              glue("baseMean_({group2})v({group1})"),
              glue("Log2fc_({group2})v({group1})"),
              glue("lfcSE_({group2})v({group1})"), 
              glue("Stat_({group2})v({group1})"), 
              glue("P.value_({group2})v({group1})"),
              glue("Adj.p.value_({group2})v({group1})") 
            ))
            
  merged_stats_df <<- merged_stats_df %>% 
                          dplyr::left_join(df, join_by(!!feature))
})

# ---------------------- Add NCBI id to feature, i.e. ASV
# Get the best / lowest possible taxonomy assignment for the features, i.e. ASVs
tax_names <- map_chr(str_replace_all(taxonomy_table$species, ";_","") %>%
                       str_split(";"),
                     function(row) row[length(row)])

df <- data.frame(ASV=rownames(taxonomy_table), best_taxonomy=tax_names)
colnames(df) <- c(feature, "best_taxonomy")

# Pull NCBI IDS for unique taxonomy names
# Filter out unannotated entries before querying NCBI
valid_taxonomy <- df$best_taxonomy %>% unique() %>% setdiff("_")
df2_valid <- data.frame(best_taxonomy = valid_taxonomy) %>%
  mutate(NCBI_id=get_ncbi_ids(best_taxonomy, target_region),
         .after = best_taxonomy)

# Add unannotated entries with NA NCBI_id
df2_invalid <- data.frame(best_taxonomy = "_", NCBI_id = NA)
df2 <- rbind(df2_valid, df2_invalid)

# -------- Retrieve deseq normalized table from the deseq model
normalized_table <- counts(deseq_modeled, normalized=TRUE) %>% 
                        as.data.frame() %>%
                        rownames_to_column(feature)

# Creating a dataframe of samples that were dropped because they didn't
# meet our cut-off criteria
samples <- metadata[[samples_column]]
samplesdropped <- setdiff(x = samples, y = colnames(normalized_table)[-1])
missing_df <- data.frame(ASV=normalized_table[[feature]],
                         matrix(data = NA, 
                                nrow = nrow(normalized_table),
                                ncol = length(samplesdropped)
                         )
)
colnames(missing_df) <- c(feature,samplesdropped)

# Calculate mean and standard deviation of all ASVs for each group in 
# a dataframe called group_means_df
group_levels <- metadata[, groups_colname] %>% unique() %>% sort()
group_means_df <- normalized_table[feature]
walk(group_levels, function(group_level){
  
  # Initializing mean and std column names
  mean_col <- glue("Group.Mean_({group_level})")
  std_col <- glue("Group.Stdev_({group_level})")
  
  # Get a vector of samples that belong to the current group
  Samples <- metadata %>%
    filter(!!sym(groups_colname) == group_level) %>%
    pull(!!sym(samples_column))
  # Retain only samples that belong to the current group that are in the normalized table
  Samples <- intersect(colnames(normalized_table), Samples)
  
  # Calculate the means and standard deviations for the current group 
  temp_df <- normalized_table %>% select(!!feature, all_of(Samples)) %>% 
    rowwise() %>%
    mutate(!!mean_col := mean(c_across(where(is.numeric)), na.rm = TRUE),
           !!std_col := sd(c_across(where(is.numeric)), na.rm = TRUE) ) %>% 
    select(!!feature,!!sym(mean_col), !!sym(std_col))
  
  group_means_df <<- group_means_df %>% left_join(temp_df)
  
})

# Append mean, standard deviation and missing samples to the normalized table
normalized_table <- normalized_table %>%
  left_join(missing_df, by = feature) %>% 
  select(!!feature, all_of(samples))

# Calculate mean global means and standard deviations
All_mean_sd <- normalized_table %>%
  rowwise() %>%
  mutate(All.mean=mean(c_across(where(is.numeric)), na.rm = TRUE),
         All.stdev=sd(c_across(where(is.numeric)), na.rm = TRUE) ) %>% 
  select(!!feature, All.mean, All.stdev)

# Add taxonomy
merged_df <- df  %>% # statistics table
  left_join(taxonomy_table %>%
              as.data.frame() %>%
              rownames_to_column(feature)) %>%  # append taxonomy table
  select(!!feature, domain:species,everything())  # select columns of interest

# Merge all prepared tables in the desired order
merged_df <- merged_df %>%
  select(!!sym(feature):NCBI_id) %>%  # select only the features and NCBI ids
  left_join(normalized_table, by = feature) %>%  # append the normalized table
  left_join(merged_df) %>%  # append the stats table
  left_join(All_mean_sd) %>%  # append the global/ASV means and stds
  left_join(group_means_df, by = feature) %>% # append the group means and stds
  mutate(across(where(is.matrix), as.numeric)) # convert meatrix columns to numeric columns

# Defining the output file
output_file <- glue("{diff_abund_out_dir}/{output_prefix}deseq2_differential_abundance_<tech_type>{assay_suffix}.csv")
# Writing out results of differential abundance using DESeq2
# after dropping baseMean columns
write_csv(merged_df %>%
            select(-starts_with("baseMean_")), 
          output_file)

# ------------------------- Make volcano plots ------------------------ #
# Loop over group pairs and make a volcano comparing the pair 
walk(pairwise_comp_df, function(col){
  
  group1 <- col[1]
  group2 <- col[2]
  
  # Setting plot dimensions 
  plot_width_inches <- 11.1
  plot_height_inches <- 8.33
  p_val <- 0.1 # logfc cutoff
  
  # Retrieve data for plotting
  deseq_res <- results(deseq_modeled, contrast = c(group, group2, group1))
  volcano_data <- as.data.frame(deseq_res)
  volcano_data <- volcano_data[!is.na(volcano_data$padj), ]
  volcano_data$significant <- volcano_data$padj <= p_val
  
   ######Long x-axis label adjustments##########
  x_label <- glue("Log2 Fold Change\n\n( {group2} vs {group1} )")
  label_length <- nchar(x_label)
  max_allowed_label_length <- plot_width_inches * 10
  
  # Construct x-axis label with new line breaks if was too long
  if (label_length > max_allowed_label_length){
    x_label <- glue("Log2 Fold Change\n\n( {group2} \n vs \n {group1} )")
  }
  #######################################
  
  # ASVs promoted in space on right, reduced on left
  p <- ggplot(volcano_data %>% 
                as.data.frame() %>% 
                rownames_to_column(feature),
              aes(x = log2FoldChange, y = -log10(padj), 
                  color = significant, label = !!sym(feature)) 
              ) +
    geom_point(alpha=0.7, size=2) +
    geom_hline(yintercept = -log10(p_val), linetype = "dashed") +
    scale_color_manual(values=c("black", "red"), 
                       labels=c(paste0("padj > ", p_val), 
                                paste0("padj \u2264 ", p_val))) +
    ggrepel::geom_text_repel(show.legend = FALSE) + 
    expandy(-log10(volcano_data$padj)) + # Expand plot y-limit
    coord_cartesian(clip = 'off') +
    scale_y_continuous(oob = scales::oob_squish_infinite) +
    theme_bw() +
    labs(title = "Volcano Plot",
         x = x_label,
         y = "-Log10 P-value",
         color = NULL,
         caption = glue("dotted line: padj = {p_val}")) +
    theme(legend.position="top", legend.key = element_rect(colour=NA),
          plot.caption = element_text(face = 'bold.italic'))
  
  # --- Save Plot
  # Replace space in group name with underscore 
  group1 <- str_replace_all(group1, "[:space:]+", "_")
  group2 <- str_replace_all(group2, "[:space:]+", "_")
  ggsave(filename = glue("{output_prefix}({group2})v({group1})_volcano_<tech_type>{assay_suffix}.png"),
         plot = p,
         width = plot_width_inches, 
         height = plot_height_inches, 
         dpi = 300, 
         path = diff_abund_out_dir)
})
```
```bash
zip -q deseq2_volcano_plots_<tech_type>_GLAmpSeq.zip *_volcano_<tech_type>_GLAmpSeq.png
```

**Custom Functions Used in R:**

* [expandy()](#expandy)
* [get_ncbi_ids()](#get_ncbi_ids)
* [fix_names()](#fix_names)
* [gm_mean()](#gm_mean)
* [plotSparsity()](#plotSparsity)

**Parameter Definitions:**

**R script**
* *pipeline uses default values for `DESeq()` analysis* 

**zip**
- `-q` – quiet operation
- `deseq2_volcano_plots_<tech_type>_GLAmpSeq.zip` – positional argument naming the zip output file
- `*_volcano_<tech_type>_GLAmpSeq.png` – positional argument naming the input file(s) to package

**Input Data:**

* `feature` (a string specifying the feature type, i.e. "ASV" or "OTU")
* `groups_colname` (a string specifying the name of the column in the metadata table containing the group names)
* `samples_column` (a string specifying the name of the column in the metadata table containing the sample names)
* `assay_suffix` (a string specifying the suffix to be added to output files; default is the Genelab assay suffix, "_GLAmpSeq")
* `output_prefix` (a string specifying an additional prefix to be added to the output files; default is no additional prefix, "")
* `target_region` (a string specifying the amplicon target region; options are either "16S", "18S", or "ITS")
* `metadata` (a dataframe containing the sample metadata, with samples as row names and sample info as columns, output from [6b.iv. Read-in Input Tables](#6biv-read-in-input-tables))
* `feature_table` (a dataframe containing a filtered subset of the samples feature dataframe (i.e. ASV), output from [6b.v. Preprocessing](#6bv-preprocessing))
* `taxonomy_table` (a dataframe containing a filtered subset of the feature taxonomy dataframe with ASV taxonomy assignments, output from [6b.v. Preprocessing](#6bv-preprocessing))

**Output Data:**

* **differential_abundance/deseq2/deseq2_volcano_plots_<tech_type>_GLAmpSeq.zip** (zip containing the following)
  * **(\<group1\>)v(\<group2\>)_volcano\_<tech_type>_GLAmpSeq.png** (volcano plots for each pariwise comparison)
* **differential_abundance/deseq2/deseq2_differential_abundance_<tech_type>_GLAmpSeq.csv** (a comma-separated DESeq2 differential abundance results table containing the following columns:
  - ASV (identified ASVs)
  - taxonomic assignment columns
  - NCBI identifier for the best taxonomic assignment for each ASV 
  - Normalized abundance values for each ASV for each sample
  - For each pairwise group comparison:
    - log2 of the fold change (Log2fc)
    - standard error for the log2FC (lfcSE)
    - test statistic from the primary result (Stat)
    - P-value (P.value)
    - Adjusted p-value (Adj.p.value)
  - All.mean (mean across all samples)
  - All.stdev (standard deviation across all samples) 
  - For each group:
    - Group.Mean_(group) (mean within group)
    - Group.Stdev_(group) (standard deviation within group))
* **differential_abundance/deseq2/asv_sparsity_plot_<tech_type>_GLAmpSeq.png** (a diagnostic plot of ASV sparsity to be used to assess if running DESeq2 is appropriate)
<br>

---
