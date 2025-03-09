# Bioinformatics pipeline for amplicon Illumina sequencing data  

> **This page holds an overview and instructions for how GeneLab processes Illumina amplicon sequencing datasets. Exact processing commands for specific datasets that have been released are available in the [GLDS_Processing_Scripts](../GLDS_Processing_Scripts) sub-directory and/or are provided with their processed data in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/).**  

---

**Date:** March XX, 2024  
**Revision:** B  
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
| MultiQC      | 1.9              | 1.19          |
| Cutadapt     | 2.3              | 4.6           |
| R-base       | 4.1.1            | 4.4.1         |
| DADA2        | 1.20.0           | 1.30.0        |
| DECIPHER     | 2.20.0           | 2.30.0        |
| biomformat   | 1.20.0           | 1.30.0        |
| ANCOMBC      | N/A              | 2.6.0         |
| broom        | N/A              | 1.0.7         |
| DescTools    | N/A              | 0.99.57       |
| DESeq2       | N/A              | 1.42.0        |
| FSA          | N/A              | 0.9.5         |
| ggdendro     | N/A              | 0.2.0         |
| ggrepel      | N/A              | 0.9.6         |
| glue         | N/A              | 1.8.0         |
| mia          | N/A              | 1.12.0        |
| phyloseq     | N/A              | 1.46.0        |
| rcolorbrewer | N/A              | 1.1_3         |
| taxize       | N/A              | 0.9.100.1     |
| tidyverse    | N/A              | 2.0.0         |
| vegan        | N/A              | 2.6.4         |

- Added new processing steps in R to generate processed data outputs for alpha and beta diversity, 
  taxonomic summary plots, and differential abundance:
  - Alpha Diversity Analysis ([Step 7](#7-alpha-diversity-analysis))
  - Beta Diversity Analysis ([Step 8](#8-beta-diversity-analysis))
  - Group-wise and Sample-wise Taxonomic Summary Plots ([Step 9](#9-taxonomy-plots))
  - Differential Abundance Testing ([Step 10](#9-differential-abundance-analysis)) with 
    ANCOMBC 1 ([Step 10a](#10a-ancombc-1)), ANCOMBC 2 ([Step 10b](#10b-ancombc-2)), and Deseq2 ([Step 10c](#10c-deseq2))
- Assay-specific suffixes were added where needed for OSDR ("_GLAmpSeq")
- Updated reference files:
  - ITS UNITE: "UNITE\_v2023\_July2023.RData" from [DECIPHER](https://www2.decipher.codes/data/Downloads/TrainingSets/)
- Added persistent reference links to DECIPHER databases on Figshare and replaced reference links to 
  DECIPHER [website]([http://www2.decipher.codes/Classification/TrainingSets/](https://www2.decipher.codes/data/Downloads/TrainingSets/)) 
  - [SILVA SSU r138](https://figshare.com/ndownloader/files/46245217)
  - [UNITE v2023](https://figshare.com/ndownloader/files/49181545)
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
      - [Load Libraries](#load-libraries)
      - [Load Functions](#load-functions)
      - [Set Variables](#set-variables)
      - [Read-in Input Tables](#read-in-input-tables)
      - [Preprocessing](#preprocessing)
  - [**7. Alpha Diversity Analysis**](#7-alpha-diversity-analysis)
  - [**8. Beta Diversity Analysis**](#8-beta-diversity-analysis)
  - [**9. Taxonomy Plots**](#9-taxonomy-plots)
  - [**10. Differential Abundance Testing**](#10-differential-abundance-testing)
    - [10a. ANCOMBC 1](#10a-ancombc-1)
    - [10b. ANCOMBC 2](#10b-ancombc-2)
    - [10c. DESeq2 ](#10c-deseq2)

---

# Software used  

|Program|Version|Relevant Links|
|:------|:------:|:-------------|
|FastQC|0.12.1|[https://www.bioinformatics.babraham.ac.uk/projects/fastqc/](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)|
|MultiQC|1.19|[https://multiqc.info/](https://multiqc.info/)|
|Cutadapt|4.6|[https://cutadapt.readthedocs.io/en/stable/](https://cutadapt.readthedocs.io/en/stable/)|
|R-base|4.4.1|[https://www.r-project.org/](https://www.r-project.org/)|
|DADA2|1.30.0|[https://www.bioconductor.org/packages/release/bioc/html/dada2.html](https://www.bioconductor.org/packages/release/bioc/html/dada2.html)|
|DECIPHER|2.30.0|[https://bioconductor.org/packages/release/bioc/html/DECIPHER.html](https://bioconductor.org/packages/release/bioc/html/DECIPHER.html)|
|biomformat|1.30.0|[https://github.com/joey711/biomformat](https://github.com/joey711/biomformat)|
|ANCOMBC|2.6.0|[https://github.com/FrederickHuangLin/ANCOMBC](https://github.com/FrederickHuangLin/ANCOMBC)|
|broom|1.0.7|[https://CRAN.R-project.org/package=broom](https://CRAN.R-project.org/package=broom)|
|DescTools|0.99.57|[https://andrisignorell.github.io/DescTools/](https://andrisignorell.github.io/DescTools/)|
|DESeq2|1.42.0|[https://bioconductor.org/packages/release/bioc/html/DESeq2.html](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)|
|FSA|0.9.5|[https://CRAN.R-project.org/package=FSA](https://CRAN.R-project.org/package=FSA)|
|ggdendro|0.2.0|[https://CRAN.R-project.org/package=ggdendro](https://CRAN.R-project.org/package=ggdendro)|
|ggrepel|0.9.6|[https://CRAN.R-project.org/package=ggrepel](https://CRAN.R-project.org/package=ggrepel)|
|glue|1.8.0|[https://glue.tidyverse.org/](https://glue.tidyverse.org/)|
|mia|1.12.0|[https://github.com/microbiome/mia](https://github.com/microbiome/mia)|
|phyloseq|1.46.0|[https://bioconductor.org/packages/release/bioc/html/phyloseq.html](https://bioconductor.org/packages/release/bioc/html/phyloseq.html)|
|rcolorbrewer|1.1_3|[https://CRAN.R-project.org/package=RColorBrewer](https://CRAN.R-project.org/package=RColorBrewer)|
|taxize|0.9.100.1|[https://docs.ropensci.org/taxize/](https://docs.ropensci.org/taxize/)|
|tidyverse|2.0.0|[https://CRAN.R-project.org/package=tidyverse](https://CRAN.R-project.org/package=tidyverse)|
|vegan|2.6.4|[https://cran.r-project.org/package=vegan](https://cran.r-project.org/package=vegan)|

# Reference databases used
<update figshare links once the updated DBs are downloaded>
|Program used| Database| DECIPHER Link | GeneLab Figshare Link | GeneLab Download Date|
|:-----|:-----:|:-----|--------:|
|DECIPHER| SILVA SSU r138_2 | [https://www2.decipher.codes/data/Downloads/TrainingSets/SILVA_SSU_r138_2_2024.RData](https://www2.decipher.codes/data/Downloads/TrainingSets/SILVA_SSU_r138_2_2024.RData) |[SILVA_SSU_r138_2019.RData](https://figshare.com/ndownloader/files/46245217)| <insert download date >|
|DECIPHER| UNITE v2024 | [https://www2.decipher.codes/data/Downloads/TrainingSets/UNITE_v2024_April2024.RData](https://www2.decipher.codes/data/Downloads/TrainingSets/UNITE_v2024_April2024.RData) | [UNITE_v2023_July2023.RData](https://figshare.com/ndownloader/files/49181545)| <insert download date >|
|DECIPHER| PR2 v4.13 | [https://www2.decipher.codes/data/Downloads/TrainingSets/PR2_v4_13_March2021.RData](https://www2.decipher.codes/data/Downloads/TrainingSets/PR2_v4_13_March2021.RData) | [PR2_v4_13_March2021.RData](https://figshare.com/ndownloader/files/46241917)| <insert download date >|
---

# General processing overview with example commands  

> Exact processing commands for specific datasets are available in the [GLDS_Processing_Scripts](../GLDS_Processing_Scripts) sub-directory of this repository, and/or are provided with their processed data in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/).
>
> Output files listed in **bold** below are included with each Amplicon Seq processed dataset in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/).

---

## 1. Raw Data QC  

<br>

### 1a. Raw Data QC  

```
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

```
multiqc --interactive -n raw_multiqc_GLAmpSeq -o /path/to/raw_multiqc/output/raw_multiqc_GLAmpSeq_report /path/to/directory/containing/raw_fastqc/files

zip -r raw_multiqc_GLAmpSeq_report.zip raw_multiqc_GLAmpSeq_report
```

**Parameter Definitions:**

- `--interactive` – force reports to use interactive plots
- `-n` – prefix name for output files
- `-o` – the output directory to store results
- `/path/to/directory/containing/raw_fastqc/files` – the directory holding the output data from the FastQC run, provided as a positional argument

**Input Data:**

* \*fastqc.zip (FastQC output data, output from [Step 1a](#1a-raw-data-qc))

**Output Data:**

* **raw_multiqc_GLAmpSeq_report.zip** (zip containing the following)
  * **raw_multiqc_GLAmpSeq.html** (multiqc output html summary)
  * **raw_multiqc_GLAmpSeq_data** (directory containing multiqc output data)

<br>  

---

## 2. Trim Primers  

The location and orientation of primers in the data is important to understand in deciding how to do this step. `cutadapt` has many options for primer identification and removal, which are described in detail in the [cutadapt adapter type documentation](https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types).  

The following example commands show how it was done for some samples of [GLDS-200](https://osdr.nasa.gov/bio/repo/data/studies/OSD-200), which was 2x250 sequencing of the 16S gene using these primers:  
* forward: 5’-GTGCCAGCMGCCGCGGTAA-3’  
* reverse: 5’-GGACTACVSGGGTATCTAAT-3’  

Due to the size of the target amplicon and the type of sequencing done here, both forward and reverse primers are expected to be on each of the forward and reverse reads. It therefore takes “linked” primers as input for forward and reverse reads, specified in the example command below by the `...` between them. It also expects that the primers start at the first position of the reads (“anchored”), specified with the leading `^` characters in the example command below.  

The following website is useful for reverse complementing primers and dealing with degenerate bases appropriately: [http://arep.med.harvard.edu/labgc/adnan/projects/Utilities/revcomp.html](http://arep.med.harvard.edu/labgc/adnan/projects/Utilities/revcomp.html)  

```
cutadapt -a ^GTGCCAGCMGCCGCGGTAA...ATTAGATACCCSBGTAGTCC -A ^GGACTACVSGGGTATCTAAT...TTACCGCGGCKGCTGGCAC \
         -o sample1_R1_trimmed.fastq.gz -p sample1_R2_trimmed.fastq.gz sample1_R1_raw.fastq.gz sample1_R2_raw.fastq.gz \
         --discard-untrimmed
```

**Parameter Definitions:**

*	`-a` – specifies the primers and orientations expected on the forward reads (when primers are linked as noted above)
*	`-A` – specifies the primers and orientations expected on the reverse reads (when primers are linked as noted above)
*	`-o` – specifies file path/name of forward, primer-trimmed reads
*	`-p` – specifies file path/name of reverse, primer-trimmed reads
*	`sample1_R1_raw.fastq.gz` – this and following “R2” file are positional arguments specifying the forward and reverse reads, respectively, for input
*	`--discard-untrimmed` – this filters out those reads where the primers were not found as expected

**Input Data:**

* \*fastq.gz (raw reads)

**Output Data:**

* **\*trimmed.fastq.gz** (trimmed reads)
* **trimmed-read-counts_GLAmpSeq.tsv** (per sample read counts before and after trimming)
* **cutadapt_GLAmpSeq.log** (log file of standard output and error from cutadapt)

<br>

---

## 3. Quality Filtering
> The following is run in an R environment.  

Specific settings required will depend on the dataset being processing. These include parameters such as `truncLen`, which might depend on the target amplicon and its size, and `maxEE` which might depend on the quality of the sequencing run. For instance, when working with ITS data, it may be preferable to omit using the `truncLen` parameter if the target amplified region is expected to vary to lengths greater than the read size. More information on these parameters can be found at these sites:  
* [https://benjjneb.github.io/dada2/tutorial.html](https://benjjneb.github.io/dada2/tutorial.html)  
* [https://astrobiomike.github.io/amplicon/dada2_workflow_ex](https://astrobiomike.github.io/amplicon/dada2_workflow_ex)  


The following is an example from a [GLDS-200](https://osdr.nasa.gov/bio/repo/data/studies/OSD-200) sample that used paired-end 2x250 sequencing with the following 16S primers:  
* forward: 5’-GTGCCAGCMGCCGCGGTAA-3’
* reverse: 5’- GGACTACVSGGGTATCTAAT-3’

```
filtered_out <- filterAndTrim(fwd=“sample1_R1_trimmed.fastq.gz”, filt=“sample1_R1_filtered.fastq.gz”,
                              rev=“sample1_R2_trimmed.fastq.gz”, filt.rev=“sample1_R1_filtered.fastq.gz”,
                              truncLen=c(220, 160), maxN=0, maxEE=c(2,2),
                              truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)
```

**Parameter Definitions:**

*	`filtered_out <-` – specifies the variable that will store the summary results within in our R environment
*	`filterAndTrim()` – the DADA2 function we are calling, with the following parameters set within it
*	`fwd=` – specifying the path to the forward reads, here “sample1_R1_trimmed.fastq.gz”
*	`filt=` – specifying the path to where the output forward reads will be written
*	`rev=` – specifying the path to the reverse reads, here “sample1_R2_trimmed.fastq.gz”; only applicable if paired-end
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
* **filtered-read-counts_GLAmpSeq.tsv** (per sample read counts before and after filtering)

<br>

---

## 4. Filtered Data QC

<br>

### 4a. Filtered Data QC
```
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
```
multiqc --interactive -n filtered_multiqc_GLAmpSeq -o /path/to/filtered_multiqc/output/filtered_multiqc_GLAmpSeq_report /path/to/directory/containing/filtered_fastqc/files

zip -r filtered_multiqc_GLAmpSeq_report.zip filtered_multiqc_GLAmpSeq_report
```

**Parameter Definitions:**

- `--interactive` – force reports to use interactive plots
- `-n` – prefix name for output files
- `-o` – the output directory to store results
- `/path/to/directory/containing/filtered_fastqc/files` – the directory holding the output data from the FastQC run, provided as a positional argument

**Input Data:**

* \*fastqc.zip (FastQC output data, output from [Step 4a](#4a-filtered-data-qc))

**Output Data:**

* **filtered_multiqc_GLAmpSeq_report.zip** (zip containing the following)
  * **filtered_multiqc_GLAmpSeq_report.html** (multiqc output html summary)
  * **filtered_multiqc_GLAmpSeq_data** (directory containing multiqc output data)

<br>

---

## 5. Calculate Error Mdel, Apply DADA2 Algorithm, Assign Taxonomy, and Create Output Tables
> The following is run in an R environment.  

These example commands as written assume paired-end data, with notes included on what would be different if working with single-end data. The taxonomy reference database used below is an example only, suitable for the example 16S dataset ([GLDS-200](https://osdr.nasa.gov/bio/repo/data/studies/OSD-200)) used here. Other taxonomy references databases designed for DECIPHER can be found here: [https://www2.decipher.codes/data/Downloads/TrainingSets/](https://www2.decipher.codes/data/Downloads/TrainingSets/)  

<br>

### 5a. Learning the Error Rates
```R
## Forward error rates ##
forward_errors <- learnErrors(fls=“sample1_R1_filtered.fastq.gz”, multithread=TRUE)

## Reverse error rates (skip if single-end data) ##
reverse_errors <- learnErrors(fls=“sample1_R2_filtered.fastq.gz”, multithread=TRUE)
```

**Parameter Definitions:**  

*	`learnErrors()` – the DADA2 function we are calling, with the following parameters set within it
*	`fls=` – specifies the path to the filtered reads (either forward or reverse)
*	`multithread=TRUE` – determine number of cores available and run in parallel when possible (can also take an integer specifying the number of cores to use)

**Input Data:**

* \*filtered.fastq.gz (filtered reads, output from [Step 3](#3-quality-filtering))

**Output Data:**

* `forward_errors` (variable containing a named list storing a numeric matrix with the forward error rates)
* `reverse_errors` (variable containing a named list storing a numeric matrix with the reverse error rates (only for paired-end data))

<br>

### 5b. Inferring Sequences
```R
## Inferring forward sequences ##
forward_seqs <- dada(derep=“sample1_R1_filtered.fastq.gz”, err=forward_errors, pool=“pseudo”, multithread=TRUE)

## Inferring reverse sequences (skip if single-end)##
reverse_seqs <- dada(derep=“sample1_R2_filtered.fastq.gz”, err=reverse_errors, pool=“pseudo”, multithread=TRUE)
```

**Parameter Definitions:**  

* `dada()` – the DADA2 function we are calling, with the following parameters set within it
* `derep=` – the path to the filtered reads (either forward or reverse)
* `err=` – the object holding the error profile for the inferred reads (either forward or reverse)
* `pool=“pseudo”` – setting the method of incorporating information from multiple samples, "pseudo" instructs the algorithm to perform pseudo-pooling between individually processed samples
* `multithread=TRUE` – determine number of cores available and run in parallel when possible (can also take an integer specifying the number of cores to use)

**Input Data:**

* \*filtered.fastq.gz (filtered reads, output from [Step 3](#3-quality-filtering))
* `forward_errors` (forward error rates, output from [Step 5a](#5a-learning-the-error-rates))
* `reverse_errors` (reverse error rates, output from [Step 5a](#5a-learning-the-error-rates) (only for paired-end))

**Output Data:**

* `forward_seqs` (variable containing a dada-class object storing the forward-read inferred sequences)
* `reverse_seqs` (variable containing a dada-class object storing the reverse-read inferred sequences (only for paired-end))

<br>

### 5c. Merging Forward and Reverse Reads; Skip if Data are Single-End
```R
merged_contigs <- mergePairs(dadaF=forward_seqs, derepF=“sample1_R1_filtered.fastq.gz”, dadaR=reverse_seqs, derepR=“sample1_R2_filtered.fastq.gz”)
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
* `forward_seqs` (forward-read inferred sequences, output from [Step 5b](#5b-inferring-sequences))
* `reverse_seqs` (reverse-read inferred sequences, output from [Step 5b](#5b-inferring-sequences))

**Output Data:**

* `merged_contigs` (variable containing a data.frame storing the merged contigs)

<br>

### 5d. Generating Sequence Table with Counts per Sample
```R
seqtab <- makeSequenceTable(merged_contigs)
```

**Parameter Definitions:**  

* If single-end data, instead of “merged_contigs”, the forward_seqs object would be provided to the `makeSequenceTable()` function here

**Input Data:**

* `merged_contigs` or `forward_seqs` (for paired-end data, the merged contigs, output from [Step 5d](#5d-generating-sequence-table-with-counts-per-sample), for single-end data, the `forward_seqs` object, output from [Step 5b](#5b-inferring-sequences))

**Output Data:**

* `seqtab` (a variable containing a named integer matrix containing the sequence table)

<br>

### 5e. Removing putative chimeras
```R
seqtab.nochim <- removeBimeraDenovo(unqs=seqtab, method=“consensus”, multithread=TRUE)
```

**Parameter Definitions:**  

* `seqtab.nochim <-` – specifies the variable that will store the results within in our R environment
* `removeBimeraDenovo()` – the DADA2 function we are calling, with the following parameters set within it
* `unqs=` – specifying the “seqtab” object created above
* `method=` – specifying the method for putative-chimera identification and removal, "consensus" instructs the function to check the samples in the sequence table independently for bimeras and make a consensus decision on each sequence variant 
* `multithread=TRUE` – determine number of cores available and run in parallel when possible (can also take an integer specifying the number to run)

**Input Data:**

* `seqtab` (sequence table, output from [Step 5d](#5d-generating-sequence-table-with-counts-per-sample))

**Ouptut Data:**

* `seqtab.nochim` (variable containing a named integer matrix containing the sequence table without putative chimeras)

<br>

### 5f. Assigning Taxonomy
<Change download file URL in the command below> 

```R
## Creating a DNAStringSet object from the ASVs: ##
dna <- DNAStringSet(getSequences(seqtab.nochim))

## Downloading the reference R taxonomy object: ##
download.file( url="https://figshare.com/ndownloader/files/46245217", destfile=“SILVA_SSU_r138_2_2024.RData”)

## Loading taxonomy object: ##
load(“SILVA_SSU_r138_2_2024.RData”)

## Classifying sequences:
tax_info <- IdTaxa(test=dna, trainingSet=trainingSet, strand=“both”, processors=NULL)
```

**Parameter Definitions:** 

- `download.file()`
  - `url=` - reference database URL address to download
  - `destfile=` - local path/name for the downloaded file
- `IdTaxa()`
  - `test=dna` - DNAStringSet object holding sequences to classify
  - `trainingSet=trainingSet` - specifies the reference database to use
  - `strand="both"` - specifies to check taxonomy assignment in both orientations
  - `processors=NULL` - specifies the number of processors to use, `NULL` indicates to use all available cores or an integer may be provided to manually specify the number to use

**Input Data:**

* `seqtab.nochim` (variable containing the sequence table without putative chimeras, output from [Step 5e](#5e-removing-putative-chimeras))
* `trainingSet` (variable provided in the RData object holding the reference database, SILVA_SSU_r138_2_2024.RData)

**Output Data:**

* `tax_info` (variable containing the DECIPHER Taxa object containing assigned taxons)

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
write(asv_fasta, "ASVs_GLAmpSeq.fasta")

## Making then writing a count table: ##
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)

write.table(asv_tab, "counts_GLAmpSeq.tsv", sep="\t", quote=F, col.names=NA)

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

write.table(tax_tab, "taxonomy_GLAmpSeq.tsv", sep = "\t", quote=F, col.names=NA)

## Generating then writing biom file format: ##
biom_object <- make_biom(data=asv_tab, observation_metadata=tax_tab)
write_biom(biom_object, "taxonomy-and-counts_GLAmpSeq.biom")

## Making a combined taxonomy and count table ##
tax_and_count_tab <- merge(tax_tab, asv_tab)
write.table(tax_and_count_tab, "taxonomy-and-counts_GLAmpSeq.tsv", sep="\t", quote=FALSE, row.names=FALSE)
```

**Input Data:**

* `seqtab.nochim` (variable containing the sequence table without chimeras, output from [Step 5e](#5e-removing-putative-chimeras))
* `tax_info` (variable containing the DECIPHER Taxa object, output from [Step 5f](#5f-assigning-taxonomy))

**Output Data:**

* **ASVs_GLAmpSeq.fasta** (inferred sequences)
* **counts_GLAmpSeq.tsv** (count table)
* **taxonomy_GLAmpSeq.tsv** (taxonomy table)
* **taxonomy-and-counts_GLAmpSeq.tsv** (combined taxonomy and count table)
* **taxonomy-and-counts_GLAmpSeq.biom** (count and taxonomy table in biom format)
* **read-count-tracking_GLAmpSeq.tsv** (read counts at each processing step)

<br>

---

## 6. Amplicon Seq Data Analysis Set Up  
  
<br>

### 6a. Create Sample Runsheet

> Note: Rather than running the command below to create the runsheet needed for processing, the runsheet may also be created manually by following the examples for [Paired-end](../Workflow_Documentation/NF_AmpIllumina-B/workflow_code/PE_file.csv) and [Single-end](../Workflow_Documentation/NF_AmpIllumina-B/workflow_code/SE_file.csv) samples. When creating this table manually, the most important columns for the analyses below are:

* `sample_id` - column with unique sample names.
* `groups`    - column with the groups/treatments that each sample belong to. This column is used for comparison.

```bash
### Download the *ISA.zip file from the OSDR ###

dpt-get-isa-archive \
 --accession OSD-###

### Parse the metadata from the *ISA.zip file to create a sample runsheet ###

dpt-isa-to-runsheet --accession OSD-### \
 --config-type AmpSeq \
 --config-version Latest \
 --isa-archive *ISA.zip
```

**Parameter Definitions:**

* `--accession OSD-###` - OSD accession ID (replace ### with the OSD number being processed), used to retrieve the urls for the ISA archive and raw reads hosted on the OSDR
* `--config-type` - instructs the script to extract the metadata required for Amplicon Sequencing data processing from the ISA archive
* `--config-version` - specifies the `dp-tools` configuration version to use, a value of `Latest` will specify the most recent version
* `--isa-archive` - specifies the *ISA.zip file for the respective OSD dataset, downloaded in the `dpt-get-isa-archive` command


**Input Data:**

* No input data required but the OSD accession ID needs to be indicated, which is used to download the respective ISA archive 

**Output Data:**

* *ISA.zip (compressed ISA directory containing Investigation, Study, and Assay (ISA) metadata files for the respective OSD dataset, used to define sample groups - the *ISA.zip file is located in the [OSDR](https://osdr.nasa.gov/bio/repo/) under 'Files' -> 'Study Metadata Files')

* **{OSD-Accession-ID}_AmpSeq_v{version}_runsheet.csv** (table containing metadata required for processing, version denotes the dp_tools schema used to specify the metadata to extract from the ISA archive)

<br>

> The remainder of this document is performed in R.  


### 6b. R Environment Set Up


#### Load Libraries

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
```

#### Load Functions

```R
# Function to calculate text size for plotting
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

# A function to expand a plot's y-limit
expandy <- function(vec, ymin=NULL) {
  # vec [NUMERIC] - a numeric vector of y valuesx.

  max.val <- max(vec, na.rm=TRUE) + 0.1
  min.log <- floor(log10(max.val))
  
  # expand_limits(y=c(ymin, ceiling(max.val/10^min.log)*10^min.log))
  expand_limits(y=c(ymin, max.val))
}


# A function to create a phyloseq object with the appropriate
# sample count transformation depending on the supplied transformation method,
# either 'rarefy' or  'vst'
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
    
    # Loop through the sequences per sample and return the count
    # nearest to the minimum required rarefaction depth
    for (count in seq_per_sample) {
      # Get the count equal to rarefaction_depth or nearest to it
      if(count >= rarefaction_depth) {
        depth <- count
        break
      }
      
    }
    
    #----- Rarefy sample counts to even depth per sample
    ps <- rarefy_even_depth(physeq = ASV_physeq, 
                            sample.size = depth,
                            rngseed = 1, 
                            replace = FALSE, 
                            verbose = FALSE)

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

# -----------  A function  Hierarchical Clustering and dendogram plotting
make_dendogram <- function(dist_obj, metadata, groups_colname,
                           group_colors, legend_title){

  # Hierarchical Clustering
  sample_clust <- hclust(d = dist_obj, method = "ward.D2")
  
  # Extract clustering data for plotting
  hcdata <- dendro_data(sample_clust, type = "rectangle")
  segment_data <- segment(hcdata) # sepcifications for tree structure
  label_data <- label(hcdata) %>%
    left_join(metadata %>% 
                rownames_to_column("label")) # Labels are sample names

  # Plot dendogram
  dendogram <- ggplot() +
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
  
  
  return(dendogram)
  
}

# A function to arun variance test and adonis test
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

# Make PCoA
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

# A function to filter out rare features from a feature table depending
# on the supplied cut off.
remove_rare_features <- function(feature_table, cut_off_percent=3/4){
  
  # Filter by occurrence in a fraction of samples
  # Define a cut-off for determining what's rare
  cut_off <- cut_off_percent * ncol(feature_table)
  # Get the occurrence for each feature
  feature_occurence <- rowSums(feature_table > 0)
  # Get names of the abundant features
  abund_features <- names(feature_occurence[feature_occurence >= cut_off])
  # Remove rare features
  abun_features.m <- feature_table[abund_features,]
  return(abun_features.m)
}

 # Function to process a taxonopmy assignment table
process_taxonomy <- function(taxonomy, prefix='\\w__') {
  
  # Ensure that all columns are of character data type
  taxonomy <- apply(X = taxonomy, MARGIN = 2, FUN = as.character) 
  
  # Loop over every column (rank i.e. domain to species) amd make the necessary edits
  for (rank in colnames(taxonomy)) {
    # Delete the taxonomy prefix
    taxonomy[,rank] <- gsub(pattern = prefix, x = taxonomy[, rank],
                            replacement = '')
    indices <- which(is.na(taxonomy[,rank]))
    taxonomy[indices, rank] <- rep(x = "Other", times=length(indices)) 
    # replace empty cell with the string 'Other'
    indices <- which(taxonomy[,rank] == "")
    taxonomy[indices,rank] <- rep(x = "Other", times=length(indices))
  }
  # Replace _ with space
  taxonomy <- apply(X = taxonomy,MARGIN = 2,
                    FUN =  gsub,pattern = "_",replacement = " ") %>% 
    as.data.frame(stringAsfactor=F)
  return(taxonomy)
}

# Function to format a taxonomy assignment table by appending a suffix
# to a known name
format_taxonomy_table <- function(taxonomy, stringToReplace="Other", suffix=";Other") {

  for (taxa_index in seq_along(taxonomy)) {
    
    indices <- grep(x = taxonomy[,taxa_index], pattern = stringToReplace)
    
    taxonomy[indices,taxa_index] <- 
      paste0(taxonomy[indices,taxa_index-1],
             rep(x = suffix, times=length(indices)))
    
  }

  return(taxonomy)
}

fix_names<- function(taxonomy,stringToReplace,suffix){
  
  for(index in seq_along(stringToReplace)){
    taxonomy <- format_taxonomy_table(taxonomy = taxonomy,
                                      stringToReplace=stringToReplace[index], 
                                      suffix=suffix[index])
  }
  return(taxonomy)
}

# A function to generate taxon level count matrix based on a taxonomy table and
# an existing feature table
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


# Function to group rare taxa or return a table with the rare taxa
group_low_abund_taxa <- function(abund_table, threshold=0.05,
                                 rare_taxa=FALSE) {
  
  # Intialize an empty vector that will contain the indices for the
  # low abundance columns/ taxa to group
  taxa_to_group <- c()
  #intialize the index variable of species with low abundance (taxa/columns)
  index <- 1
  
  # Loop over every column or taxa then check to see if the max abundance is less than the set threshold
  # if true, save the index in the taxa_to_group vector variable
  while(TRUE){
    
    for (column in seq_along(abund_table)){
      if(max(abund_table[,column]) < threshold ){
        taxa_to_group[index] <- column
        index = index + 1
      }
    }
    if(is.null(taxa_to_group)){
      threshold <- readline("please provide a higher threshold for grouping rare taxa, only numbers are allowed   ")
      threshold <- as.numeric(threshold)
    }else{
      break
    }
    
  }
  
  
  if(rare_taxa){
    abund_table <- abund_table[,taxa_to_group,drop=FALSE]
  }else{
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


# Function to collapse samples in a feature table with a defined function(fun)
# based on a group in metadata 
collapse_samples <- function(taxon_table, metadata, group, fun=sum, convertToRelativeAbundance=FALSE){

  common.ids <- intersect(rownames(taxon_table),rownames(metadata))
  metadata <- droplevels(metadata[common.ids,,drop=FALSE])
  taxon_table <- taxon_table[common.ids,,drop=FALSE]
  taxon_table <- cbind(subset(x = metadata, select=group),taxon_table)
  
  taxon_table <- aggregate(reformulate(termlabels = group, response = '.'),
                           data = taxon_table, FUN = fun)
  rownames(taxon_table) <- taxon_table[,1]
  taxon_table <- taxon_table[,-1]
  if(convertToRelativeAbundance){
    taxon_table <- t(apply(X = taxon_table, MARGIN = 1, FUN = function(x) x/sum(x)))
  }
  
  final <- list(taxon_table,metadata)
  names(final) <- c("taxon_table","metadata")
  return(final)
}

taxize_options(ncbi_sleep = 0.8)
# A function to retrieve NCBI taxonomy id for a given taxonomy name
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

# Error handling function when running ANCOMBC2
find_bad_taxa <- function(cnd){

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

# A function to run ANCOMBC2 while handling common errors 
ancombc2 <- function(data, ...) {

  tryCatch(
    ANCOMBC::ancombc2(data = data, ...),
    error = function(cnd) {
            
      res  <- find_bad_taxa(cnd)
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

# Geometric mean function used when running DESeq2
gm_mean <- function(x, na.rm=TRUE) {
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

```

#### Set Variables

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
        strip.text =  element_text(size = 14,face ='bold', color = 'black'))
```

**Function Parameter Definitions:**

- `calculate_text_size()`
  - `num_samples=` - the number of samples to plot
  - `start_samples = 25` - default number of samples, used to specify text size
  - `min_size = 3` - minimum text size for plotting
- `expandy()`
  - `vec=` - a numeric vector containing y-values
- `transform_phyloseq()`
  - `feature_table=` - dataframe containing feature / ASV counts with samples as columns and features as rows
  - `metadata=` - dataframe containing sample metadata with samples as row names and sample info as columns
  - `method=` - either "rarefy" or "vst" to specify rarefaction or variance stabilizing transformation, respectively
  - `rarefaction_depth=500` - sample rarefaction to even depth when the Bray Curtis distance method is used
- `make_dendogram()`
  - `dis_obj=` - a distance object holding the calculated distance (euclidean, bry curtid etc.) between samples
  - `metadata=` - dataframe containing sample metadata with samples as row names and sample info as columns
  - `groups_colname=` - name of the column in the metadata dataframe to use for specifying sample groups
  - `legend_title=` - legend title to use for plotting
- `run_stats()`
  - `dis_obj=` - a distance object holding the calculated distance (euclidean, bry curtid etc.) between samples
  - `metadata=` - dataframe containing sample metadata with samples as row names and sample info as columns
  - `groups_colname=` - name of the column in the metadata dataframe to use for specifying sample groups
- `plot_pcoa()`
  - `ps=` - phyloseq object contructed from feature, taxonmy, and metadata tables
  - `stats_res=` - named list generated after running the `run_stats()` function; the list should contain variance and adonis tests dataframes
  - `distance_method=` - method used to calculate the distance between samples; values can be "euclidean" or "bray" for euclidean or Bray Curtis distance, respectively
  - `groups_colname=` - name of the column in the metadata dataframe to use for specifying sample groups
  - `group_colors=` - named character vector of colors for each group in `groups_colname`
  - `legend_title=` - legend title to use for plotting
  - `addtext=FALSE` - boolean to specify if sample labels should be added to the pcoa plot; default=FALSE
- `remove_rare_features()`
  - `feature_table=` - a feature table matrix with samples as columns and features as rows
  - `cut_off_percent=3/4` - the cut-off fraction or decimal between 0.001 to 1 of the total number of samples to determine the most abundant features; by default it removes features that are not present in 3/4 of the total number of samples
- `process_taxonomy()`
  - `taxonomy=` - the taxonomy dataframe to be processed
  - `prefix='\\w__'` - a regular expression specifying the characters to remove from taxon names; use '\\w__'  for greengenes and 'D_\\d__' for SILVA
- `format_taxonomy_table()`
  - `taxonomy=` - the taxonomy dataframe to be formatted
  - `stringToReplace="Other"` - specifies the string to replace, "Other" is used by default
  - `suffix=";Other"` - specifies the replacement string, ";Other" is used by dfault
- `fix_names()`
  - `taxonomy=` - a taxonomy dataframe with taxonomy ranks as column names
  - `stringToReplace=` - a vector of regex strings specifying what to replace in the taxonomy dataframe
  - `suffix=` - a vector of regex strings specifying the replacement strings to use
- `make_feature_table()`
  - `count_matrix=` - matrix containing ASVs or OTUs and their respecitve counts
  - `taxonomy=` - a taxonomy matrix with taxonomy ranks as column names
  - `taxon_level=` - a string defining the taxon levels, i.e. domain to species
  - `samples2keep=NULL` - specifies the samples to evaluate; the default value of NULL instructs the function to only retain taxa found in at least one sample
- `group_low_abund_taxa()`
  - `abund_table=` - relative abundance matrix with taxa as columns and samples as rows
  - `threshold=0.05` - specifies the threshold for filtering out rare taxa; default value is 0.05
  - `rare_taxa=FALSE` - boolean specifying if only rare taxa should be returned, if set to TRUE then a table with only the rare taxa will be returned; default is FALSE
- `collapse_samples()`
  - `taxon_table=` - a matrix count table with samples as rows and features (e.g. ASVs or OTUs) as columns
  - `metadata=` - dataframe containing sample metadata with samples as row names and sample info as columns
  - `group=` - specifies the column in the `metadata` dataframe containing the sample groups that will be used to collapse the samples
  - `fun=sum` - specifies the function to apply when collapsing the samples; sum is used by default
  - `convertToRelativeAbundance=FALSE` - boolean specifying if the value in the taxon table should be converted to per sample relative abundance values; default is FALSE
- `get_ncbi_ids()`
  - `taxonomy=` - string specifying the taxonomy name that will be used to search for the respective NCBI ID
  - `target_region=` - the amplicon target region to analyze; options are "16S", "18S", or "ITS"
- `find_bad_taxa()`
  - `cnd=` - specifies the error condition to catch when running the `ancombc2()` function
- `ancombc2()`
  - `data=` - specifies the treeSummarizedExperiment containing the feature, taxonomy and metdata to be analyzed using ancombc2
- `gm_mean()`
  - `x=` - a numeric vector specifying the values to calculate the geometirc mean on
  - `na.rm=TRUE` - boolean specifying if NAs should be removed prior to calculating the geometirc mean; default is TRUE  

**Input Data:** 

*No input data required*

**Output Data:**

* `custom_palette` (variable containing a character vector defining a custom color palette for coloring plots)
* `publication_format` (variable specifying the custom ggplot theme for plotting)
* *several functions, indicated in the Function Parameter Definitions above, that will be called in the following pipeline steps*  

<br>

#### Read-in Input Tables

```R
custom_palette  <- {COLOR_VECTOR}
groups_colname <- "groups"
sample_colname <- "Sample Name"
metadata_file <- file.path("{OSD-Accession-ID}_AmpSeq_v{version}_runsheet.csv")
features_file <- file.path("counts_GLAmpSeq.tsv")
taxonomy_file <- file.path("taxonomy_GLAmpSeq.tsv")

# Read-in metadata and convert from tibble to dataframe
metadata <- read_csv(file = metadata_file) %>% as.data.frame()
# Set row names
row.names(metadata) <- metadata[[sample_colname]]
# Write out Sample Table
write_csv(x = metadata %>% select(!!sym(samples_column), !!sym(groups_colname)),
          file = glue("{diff_abund_out_dir}{output_prefix}SampleTable{assay_suffix}.csv"))

# Delete sample column since the rownames now contain sample names
metadata[,sample_colname] <- NULL
# Get unique group names
group_column_values <- metadata %>% pull(!!sym(groups_colname))
group_levels <- unique(group_column_values)

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

# Subset metadata to contain on the groups and color columns
sample_info_tab <- metadata %>%
  select(!!groups_colname, color) %>% # select groups and color columns
  arrange(!!sym(groups_colname)) # metadata by groups column

# Retrieves unique colors
values <- sample_info_tab %>% pull(color) %>% unique()

# ---- Import Feature or ASV table ---- #
feature_table <- read.table(file = features_file, header = TRUE,
                            row.names = 1, sep = "\t")

# ---- Import Taxonomy table ---- #
taxonomy_table <-  read.table(file = taxonomy_file, header = TRUE,
                              row.names = 1, sep = "\t")
```

**Parameter Definitions:**

*	`groups_colname` - variable containing the name of the column in the metadata table containing the group names
* `sample_colname` - variable contianing the name of the column in the metadata table containing the sample names

**Input Data:**

* {OSD-Accession-ID}_AmpSeq_v{version}_runsheet.csv (a comma-separated sample metadata file containing sample group information, output from [Step 6a](#6a-create-sample-runsheet))
*	counts_GLAmpSeq.tsv (a tab separated samples feature table (i.e. ASV or OTU table) containing feature counts, output from [Step 5g](#5g-generating-and-writing-standard-outputs))
* taxonomy_GLAmpSeq.tsv (a tab separated feature taxonomy table containing ASV taxonomy assignments, output from [Step 5g](#5g-generating-and-writing-standard-outputs))
* `custom_palette` (a color palette, output from [Set Variables](#set-variables))

**Output Data:**

* `metadata` (variable containing a metadata dataframe with samples as row names and sample info, including groups and group colors, as columns)
* `feature_table` (variable containing a samples feature dataframe (i.e. ASV) with feature counts)
* `taxonomy_table` (variable containing a feature taxonomy dataframe containing ASV taxonomy assignments)
* `sample_info_tab` (variable containing a subset of the metadata dataframe with samples as row names and group names and group colors as columns)
* `values` (variable containing a character vector of unique color values for each group)
* `sample_names` (variable containing a character vector of sample names)
* `deseq2_sample_names` (variable containing a character vector of unique sample names)
* `group_colors` (variable containing a named character vector of colors for each group)
* `group_levels` (variable containing a character vector of unique group names)

<br>

#### Preprocessing

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
taxonomy_table <- taxonomy_table[-which(is.na(taxonomy_table$domain)),]

# Handle case where no domain was assigned but a phylum was.
if(all(is.na(taxonomy$domain))){
  
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
taxonomy_table  <- process_taxonomy(taxonomy_table)
rownames(taxonomy_table) <- feature_names
taxonomy_table <- fix_names(taxonomy_table, "Other", ";_")


# Get long asv taxonomy names and clean
species <- taxonomy_table %>%
  unite(species,domain:species,sep = ";") %>% 
pull %>% str_replace_all("Other", "_")

taxonomy_table <- fix_names(taxonomy_table, "Other", ";_")

taxonomy_table[,"species"] <- species


# Subset tables 

# Get features common to the taxonomy and feature table 
common_ids <- intersect(rownames(feature_table), rownames(taxonomy_table))

# Subset the feature and taxonomy tables to contain 
# only features found in both tables
feature_table <- feature_table[common_ids,]
taxonomy_table <- taxonomy_table[common_ids,]
```
**Parameter Definitions:**

* `remove_rare`       - should rare features and samples be filtered out prior to analysis? If true, rare features and 
                        samples will be removed according to the cutoffs set below.
* `prevalence_cutoff` - If `remove_rare` is true, a numerical fraction between 0 and 1. 
                        Taxa with prevalences(the proportion of samples in which the taxon is present) less than this 
                        will be excluded from the analysis. Default is 0, i.e. do not exclude any taxon / feature.
* `library_cutoff`    - If `remove_rare` is true, a numerical threshold for filtering samples based on library sizes. 
                        Samples with library sizes less than lib_cutoff will be excluded in the analysis. Default is 0, 
                        i.e. no sample will be dropped. To discard samples with read counts less than or equal to 100, 
                        set to 100.
* `target_region`     - amplicon target region. Options are either "16S", "18S", or "ITS"
* `feature_table`     - ASV count table variable from [Read-in Input Tables](#read-in-input-tables)
* `taxonomy_table`    - ASV taxonomy table variable from [Read-in Input Tables](#read-in-input-tables)

**Output Data:**
* `feature_table` (data.frame of sample feature table containing only features in common with taxonomy_table, i.e. ASV counts data.frame)
* `taxonomy_table` (data.frame taxonomy table subset containing features in common with feature_table, i.e. ASV taxonomy data.frame)
* `sample_info_tab` (data.frame of sample information, i.e. a subset of sample metadata)

<br>

---

## 7. Alpha Diversity Analysis

Alpha diversity examines the variety and abundance of taxa within individual samples. Rarefaction curves are utilized to 
visually represent this diversity, plotting the number of unique sequences (ASVs) identified against the total number of 
sequences sampled, offering a perspective on the saturation and completeness of sampling. Metrics like Observed features 
estimates and Shannon diversity indices are employed to quantify the richness (total number of unique sequences) and 
diversity (combination of richness and evenness) within these samples.

> Please note that if you'd like to run the code in this section, make sure that you [load the libraries](#load-libraries) 
and [functions](#load-functions), [read-in input tables](#read-in-input-tables) and [preprocess](#preprocessing) them in R 
by running the lines of code in [section 6](#6-amplicon-seq-data-analysis-set-up) sequentially, particularly those 
in [section 6b](#6b-r-environment-set-up).

```R
# Create output directory if it doesn't already exist
alpha_diversity_out_dir <- "alpha_diversity/"
if(!dir.exists(alpha_diversity_out_dir)) dir.create(alpha_diversity_out_dir)
metadata  <-  {DATAFRAME}
sample_info_tab <-  {DATAFRAME} 
feature_table <- {DATAFRAME}
taxonomy_table <- {DATAFRAME}
group_colors  <- {NAMED_VECTOR} 
groups_colname  <- "groups"
rarefaction_depth <- 500
legend_title  <- "Groups"
assay_suffix  <- "_GLAmpSeq"
output_prefix  <- ""

# Create phyloseq object
ASV_physeq <- phyloseq(otu_table(feature_table, taxa_are_rows = TRUE),
                       tax_table(as.matrix(taxonomy_table)), 
                       sample_data(sample_info_tab))

seq_per_sample <- colSums(feature_table) %>% sort()

# Get rarefaction depth
# minimum value
depth <- min(seq_per_sample)

for (count in seq_per_sample) {
  
  if(count >= rarefaction_depth) {
    depth <- count
    break
    }
  
}

# -------------------- Rarefy sample counts to even depth per sample
ps.rarefied <- rarefy_even_depth(physeq = ASV_physeq, 
                                 sample.size = depth,
                                 rngseed = 1, 
                                 replace = FALSE, 
                                 verbose = FALSE)


# ------------------- Rarefaction curve
# Calculate a rough estimate of the step sample step size for plotting.
# This is meant to keep plotting time constant regardless of sample depth
step <-  (50*depth)/1000

p <- rarecurve(t(otu_table(ps.rarefied)) %>% as.data.frame(),
               step = step,
               col = sample_info_tab[["color"]], 
               lwd = 2, ylab = "ASVs", cex=0.5,
               label = FALSE, tidy = TRUE)


sample_info_tab_names <-  sample_info_tab %>% rownames_to_column("Site")

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

ggsave(filename = glue("{alpha_diversity_out_dir}/{output_prefix}rarefaction_curves{assay_suffix}.png"),
       plot=rareplot, width = 14, height = 8.33, dpi = 300)

# ------------------  Richness and diversity estimates  ------------------#

# Statistics table
diversity_metrics <- c("Observed", "Chao1", "Shannon", "Simpson")
names(diversity_metrics) <- diversity_metrics
diversity.df <- estimate_richness(ps.rarefied, 
                                  measures = diversity_metrics) %>% 
               select(-se.chao1) %>%
               rownames_to_column("samples")
               

merged_table  <- metadata %>%
  rownames_to_column("samples") %>%
  inner_join(diversity.df)

diversity_stats <- map_dfr(.x = diversity_metrics, function(metric){
  
  res <- dunnTest(merged_table[,metric],merged_table[,groups_colname])
  
  df <- res$res %>%
    separate(col = Comparison, into = c("group1", "group2"), sep = " - ") %>% 
    mutate(Metric=metric)  %>% 
    rename(p=P.unadj, p.adj=P.adj) %>% 
    mutate(p.format=round(p,digits = 2))
  
  add_significance(df, p.col='p.adj', output.col = 'p.signif') %>% 
    select(Metric,group1, group2, Z, p, p.adj, p.format, p.signif) %>% 
    mutate(across(where(is.numeric), ~round(.x, digits = 2)))
  
})


# Write diversity statistics table to file
write_csv(x = diversity_stats, 
            file = glue("{alpha_diversity_out_dir}/{output_prefix}statistics_table{assay_suffix}.csv"))

# Get diffrent letters compare groups for every diversity metric
comp_letters <- data.frame(group = group_levels)
colnames(comp_letters) <- groups_colname

walk(.x = diversity_metrics, function(metric=.x){
  
  sub_comp <- diversity_stats %>% filter(Metric == metric)
  p_values <- sub_comp$p.adj # holm p adjusted values by default
  names(p_values) <- paste(sub_comp$group1,sub_comp$group2, sep = "-")
  
  letters_df <-  enframe(multcompView::multcompLetters(p_values)$Letters,
                      name = groups_colname,
                      value = glue("{metric}_letter"))
  comp_letters <<- comp_letters %>% left_join(letters_df)
   
})


# Summary table
diversity_table <- metadata %>%
  rownames_to_column("samples") %>%
  inner_join(diversity.df) %>%
  group_by(!!sym(groups_colname))  %>% 
  summarise(N=n(), across(Observed:Simpson,
                   .fns=list(mean=mean, se=se),
                   .names = "{.col}_{.fn}")) %>% 
  mutate(across(where(is.numeric), ~round(.x, digits = 2))) %>%
  left_join(comp_letters) %>% 
  mutate(Observed = glue("{Observed_mean} ± {Observed_se}{Observed_letter}"),
         Chao1 = glue("{Chao1_mean} ± {Chao1_se}{Chao1_letter}"),
         Shannon = glue("{Shannon_mean} ± {Shannon_se}{Shannon_letter}"),
         Simpson = glue("{Simpson_mean} ± {Simpson_se}{Simpson_letter}")
         ) %>%
  select (-contains("_"))

# Write diversity summary table to file 
write_csv(x = diversity_table, 
            file = glue("{alpha_diversity_out_dir}/{output_prefix}summary_table{assay_suffix}.csv"))


# ------------------ Make richness by sample dot plots ---------------------- #

number_of_samples <- length(rownames(sample_info_tab))
richness_sample_label_size <- calculate_text_size(number_of_samples)
metrics2plot <- c("Observed", "Shannon")
names(metrics2plot) <- metrics2plot

samples_order <- metadata %>% arrange(!!sym(groups_colname)) %>% rownames()

richness_by_sample <- plot_richness(ps.rarefied, color = groups_colname,
                                    measures = metrics2plot)

richness_by_sample <-  ggplot(richness_by_sample$data %>% 
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
    strip.text =  element_text(size = 14,face ='bold')
  )

# Save sample plot
ggsave(filename = glue("{alpha_diversity_out_dir}/{output_prefix}richness_and_diversity_estimates_by_sample{assay_suffix}.png"),
       plot=richness_by_sample, width = 14, height = 8.33, dpi = 300, units = "in")


# ------------------- Make richness by group box plots ----------------------- #
richness_by_group <- plot_richness(ps.rarefied, x = groups_colname, 
                                   color = groups_colname,
                                   measures = metrics2plot)

p <- map(.x = metrics2plot, .f = function(metric){
    
  p <-  ggplot(richness_by_group$data %>%  filter(variable == metric), 
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
      strip.text =  element_text(size = 14,face ='bold')
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

richness_by_group <- wrap_plots(p, ncol = 2, guides =  'collect')  + 
                     plot_annotation(caption = "If letters are shared between two groups, then they are not significantly different (q-value > 0.05)",
                                     theme = theme(plot.caption = element_text(face = 'bold.italic')) 
                                     )

# Save group plot
width <- 3.6 * length(group_levels)
ggsave(filename = glue("{output_prefix}richness_and_diversity_estimates_by_group{assay_suffix}.png"),
       plot=richness_by_group, width = width, 
       height = 8.33, dpi = 300, units = "in",
       path = alpha_diversity_out_dir)

```
**Parameter Definitions:**

* `rarefaction_depth` – minimum rarefaction depth for alpha diversity estimation
* `groups_colname`    - name of group column in metadata to be analyzed
* `legend_title`      - legend title for plotting
* `assay_suffix`      - Genelab assay suffix. Default : "_GLAmpSeq"
* `output_prefix`     - additional prefix to be added to output files. Default: ""

**Input Data:**
* `metadata` (sample metadata with the group/treatment to be analyzed, output from [Preprocessing](#preprocessing))
* `sample_info_tab` (sample information, output from [Preprocessing](#preprocessing))
* `feature_table` (ASV counts table, output from [Preprocessing](#preprocessing))
* `taxonomy_table` (taxonomy table, output from [Preprocessing](#preprocessing))
* `publication_format` (custom ggplot theme, output from [Set Variables](#set-variables))
* `group_levels` (group levels to compare, output from [Read-in Input Tables](#read-in-input-tables))
* `group_colors` (colors for each group, output from [Read-in Input Tables](#read-in-input-tables))

**Output Data:**

* **alpha_diversity/<output_prefix>rarefaction_curves_GLAmpSeq.png** (Rarefaction curves)
* **alpha_diversity/<output_prefix>statistics_table_GLAmpSeq.csv** (Statistics Table)
* **alpha_diversity/<output_prefix>summary_table_GLAmpSeq.csv** (Summary Table)
* **alpha_diversity/<output_prefix>richness_and_diversity_estimates_by_sample_GLAmpSeq.png** (Samples Dot Plot)
* **alpha_diversity/<output_prefix>richness_and_diversity_estimates_by_group_GLAmpSeq.png** (Group boxplot)

<br>

---

## 8. Beta Diversity Analysis

Beta diversity measures the variation in species composition between different samples or environments. A common practice in working with a new dataset is to generate some exploratory visualizations like ordinations and hierarchical clusterings. These give us a quick overview of how our samples relate to each other and can be a way to check for problems like batch effects.

> Please note that if you'd like to run the code in this section, make sure that you [load the libraries](#load-libraries) 
and [functions](#load-functions), [read-in input tables](#read-in-input-tables) and [preprocess](#preprocessing) them in R 
by running the lines of code in [section 6](#6-amplicon-seq-data-analysis-set-up) sequentially, particularly those in [section 6b](#6b-r-environment-set-up).

```R
beta_diversity_out_dir <- "beta_diversity/"
if(!dir.exists(beta_diversity_out_dir)) dir.create(beta_diversity_out_dir)
metadata  <-  {DATAFRAME}  
feature_table <- {DATAFRAME} 
group_colors  <- {NAMED_VECTOR} 
groups_colname  <- "groups"
rarefaction_depth <- 500
legend_title  <- "Groups"
assay_suffix  <- "_GLAmpSeq"
output_prefix  <- ""
distance_methods <- c("euclidean", "bray")
normalization_methods <- c("vst", "rarefy")
legend_title <- NULL

options(warn=-1) # ignore warnings
# Run the analysis
walk2(.x = normalization_methods, .y = distance_methods,
      .f = function(normalization_method, distance_method){
  
# Create transformed phyloseq object
ps <- transform_phyloseq(feature_table, metadata, 
                         method = normalization_method,
                         rarefaction_depth = rarefaction_depth)

# ---------Clustering and dendogram plotting

# Extract normalized count table
count_tab <- otu_table(ps)

# Calculate distance between samples
dist_obj <- vegdist(t(count_tab), method = distance_method)

# Make dendogram
dendogram <- make_dendogram(dist_obj, metadata, groups_colname,
                            group_colors, legend_title)

# Save dendogram
ggsave(filename = glue("{output_prefix}{distance_method}_dendrogram{assay_suffix}.png"),
       plot = dendogram, width = 14, 
       height = 10, dpi = 300,
      units = "in", path = beta_diversity_out_dir)

#---------------------------- Run stats
# Checking homogeneity of variance and comparing groups using adonis test

stats_res <- run_stats(dist_obj, metadata, groups_colname)
write_csv(x = stats_res$variance, 
          file = glue("{beta_diversity_out_dir}/{output_prefix}{distance_method}_variance_table{assay_suffix}.csv"))

write_csv(x = stats_res$adonis, 
          file = glue("{beta_diversity_out_dir}/{output_prefix}{distance_method}_adonis_table{assay_suffix}.csv"))

#---------------------------- Make PCoA
# Unlabeled PCoA plot
ordination_plot_u <- plot_pcoa(ps, stats_res, distance_method, 
                               groups_colname,group_colors, legend_title) 
ggsave(filename=glue("{beta_diversity_out_dir}/{output_prefix}{distance_method}_PCoA_without_labels{assay_suffix}.png"),
       plot=ordination_plot_u, width = 14, height = 8.33, dpi = 300, units = "in")

# Labeled PCoA plot
ordination_plot <- plot_pcoa(ps, stats_res, distance_method,
                             groups_colname, group_colors, legend_title,
                             addtext=TRUE) 
ggsave(filename=glue("{beta_diversity_out_dir}/{output_prefix}{distance_method}_PCoA_w_labels{assay_suffix}.png"),
       plot=ordination_plot, width = 14, height = 8.33, dpi = 300, units = "in")

})
```
**Parameter Definitions:**

* `rarefaction_depth`     – minimum rarefaction depth when using Bray Curtis distance
* `groups_colname`        - name of group column in metadata to be analyzed
* `legend_title`          - legend title for plotting
* `assay_suffix`          - Genelab assay suffix. Default : "_GLAmpSeq"
* `output_prefix`         - additional prefix to be added to output files, Default: ""
* `distance_methods`      - method used to calculate the distance between samples, either "euclidean" and "bray" for 
                            euclidean and Bray Curtis distance, respectively
* `normalization_methods` - method for normalizing sample counts, either "vst" and "rarefy" for variance stabilizing 
                            transformation and rarefaction, respectively

**Input Data:**
* `metadata` (sample metadata with the group/treatment to be analyzed, output from [Preprocessing](#preprocessing))
* `feature_table` (ASV counts table, output from [Preprocessing](#preprocessing))
* `group_colors` (colors for each group, output from [Read-in Input Tables](#read-in-input-tables))

**Output Data:**

* **beta_diversity/<output_prefix><distance_method>_dendrogram_GLAmpSeq.png** (Dendogram)
* **beta_diversity/<output_prefix><distance_method>_adonis_table_GLAmpSeq.csv** (Adonis Stats Table)
* **beta_diversity/<output_prefix><distance_method>_PCoA_without_labels_GLAmpSeq.png** (Unlabeled PCoA)
* **beta_diversity/<output_prefix><distance_method>_PCoA_w_labels_GLAmpSeq.png** (Labeled PCoA)

<br>

---


## 9. Taxonomy Plots

Taxonomy summaries provide insights into the composition of microbial communities at various taxonomic levels.


> Please note that if you'd like to run the code in this section, make sure that you [load the libraries](#load-libraries) 
and [functions](#load-functions), [read-in input tables](#read-in-input-tables) and [preprocess](#preprocessing) them in R 
by running the lines of code in [section 6](#6-amplicon-seq-data-analysis-set-up) sequentially, particularly those in [section 6b](#6b-r-environment-set-up).

```R
taxonomy_plots_out_dir <- "taxonomy_plots/"
if(!dir.exists(taxonomy_plots_out_dir)) dir.create(taxonomy_plots_out_dir)
metadata  <-  {DATAFRAME} 
feature_table <- {DATAFRAME}
taxonomy_table <- {DATAFRAME}
custom_palette  <- {COLOR_VECTOR}
publication_format <- {GGPLOT_THEME}
groups_colname <- "groups"
assay_suffix  <- "_GLAmpSeq"
output_prefix  <- ""


# -------------------------Prepare feature tables -------------------------- #
taxon_levels <- colnames(taxonomy_table)
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
# -1 drops the kingdom level since all the microbes are bacteria
relAbundace_tbs_rare_grouped <- map2(.x = taxon_levels[-1],
                                     .y = taxon_tables[-1], 
                                     .f = function(taxon_level=.x,
                                                   taxon_table=.y){
                                      
                                       
                                       taxon_table <- apply(X = taxon_table, MARGIN = 2,
                                                            FUN = function(x) x/sum(x)) * 100
                                       
                                       
                                       taxon_table <- as.data.frame(taxon_table %>% t())
                                       if(group_rare && !(taxon_level %in% dont_group)){
                                         
                                         taxon_table <- group_low_abund_taxa(taxon_table %>%  
                                                                               as.data.frame(check.names=FALSE,
                                                                                             stringAsFactor=FALSE),
                                                                             threshold =  thresholds[taxon_level])
                                         
                                       }
                                       taxon_table$samples <- rownames(taxon_table)
                                       
                                       
                                       # Change data frame from wide to long format
                                       taxon_table <- taxon_table %>% 
                                         pivot_longer(cols = -samples, names_to = taxon_level, values_to = "relativeAbundance")
                                       taxon_table$samples <- factor(x = taxon_table$samples, 
                                                                     levels = samples_order)
                                       return(taxon_table)
                                     })


x_lab <- "Samples"/
y_lab <- "Relative abundance (%)"  
x <-  'samples'  
y <- "relativeAbundance"
facet_by <- reformulate(groups_colname)
number_of_samples <- length(samples_order)
plot_width <- 0.6 * number_of_samples

# Make sample plots
walk2(.x = relAbundace_tbs_rare_grouped, .y = taxon_levels[-1], 
                           .f = function(relAbundace_tb, taxon_level){
                             
                             df <- relAbundace_tb %>%
                               left_join(metadata %>% rownames_to_column("samples"))
                             
                          p <- ggplot(data =  df, mapping = aes(x= !!sym(x), y=!!sym(y) )) +
                               geom_col(aes(fill = !!sym(taxon_level) )) + 
                               facet_wrap(facet_by, scales = "free", 
                               nrow = 1, labeller = label_wrap_gen()) +
                               publication_format +
                               labs(x = x_lab , y = y_lab, fill= tools::toTitleCase(taxon_level)) + 
                               scale_fill_manual(values = custom_palette) +
                               theme(axis.text.x=element_text(
                                 margin=margin(t=0.5,r=0,b=0,l=0,unit ="cm"),
                                 angle = 90, 
                                 hjust = 0.5, vjust = 0.5)) + 
                               labs(x=NULL)
                          
                          ggsave(filename = glue("{taxonomy_plots_out_dir}/{output_prefix}samples_{taxon_level}{assay_suffix}.png"),
                                 plot=p, width = plot_width, height = 8.5, dpi = 300)
                          
                           })



# ------------------------ Group abundance plots ----------------------------- #
# In percentage
# phylum 1% and 2% for class to species.
thresholds <- c(phylum=1,class=2, order=2, family=2, genus=2, species=2)

# Convert from wide to long format for every treatment group of interest
group_rare <- TRUE # should rare taxa be grouped ?
maximum_number_of_taxa <- 500 # If the number of taxa is more than this then rare taxa will be grouped anyway. 

group_relAbundace_tbs <- map2(.x = taxon_levels[-1], .y = taxon_tables[-1], 
                                     .f = function(taxon_level=.x, taxon_table=.y){
                                       
                                       taxon_table <- as.data.frame(taxon_table %>% t()) 
                                       taxon_table <- (collapse_samples(taxon_table = taxon_table,
                                                                        metadata = metadata, group = groups_colname, 
                                                                        convertToRelativeAbundance = TRUE)$taxon_table * 100 )  %>%  
                                         as.data.frame(check.names=FALSE)
                                       
                                       if(ncol(taxon_table) > maximum_number_of_taxa){
                                         group_rare <- TRUE
                                       }
                                       
                                       if(group_rare){
                                         taxon_table <- group_low_abund_taxa(taxon_table %>%  
                                                                               as.data.frame(check.names=FALSE,
                                                                                             stringAsFactor=FALSE),
                                                                             threshold =  thresholds[taxon_level])
                                         group_rare <- FALSE
                                       }
                                       
                                       taxon_table[,groups_colname] <- rownames(taxon_table)
                                       
                                       
                                       # Change from wide to long format
                                       taxon_table <- taxon_table %>% 
                                         pivot_longer(cols =  -!!sym(groups_colname),
                                                      names_to = taxon_level,
                                                      values_to = "relativeAbundance")
                              
                                       return(taxon_table)
                                       
                                     })


# Make bar plots
y_lab <- "Relative abundance (%)"  
y <- "relativeAbundance"
number_of_groups <- length(group_levels)
plot_width <- 2.5 * number_of_groups
walk2(.x = group_relAbundace_tbs, .y = taxon_levels[-1], 
                           .f = function(relAbundace_tb=.x, taxon_level=.y){
                             
                             p <- ggplot(data =  relAbundace_tb %>%
                                           mutate(X = str_wrap(!!sym(groups_colname),
                                                             width = 10) # wrap long group names
                                                  ), 
                                        mapping = aes(x = X , y = !!sym(y)   )) +
                               geom_col(aes(fill = !!sym(taxon_level))) + 
                               publication_format +
                               theme(axis.text.x=element_text(
                                 margin=margin(t=0.5,r=0,b=0,l=0,unit ="cm"),
                                 angle = 0, 
                                 hjust = 0.5, vjust = 0.5)) + 
                               labs(x = NULL , y = y_lab, fill = tools::toTitleCase(taxon_level)) + 
                               scale_fill_manual(values = custom_palette)
                             ggsave(filename = glue("{taxonomy_plots_out_dir}/{output_prefix}groups_{taxon_level}{assay_suffix}.png"),
                                    plot=p, width = plot_width, height = 10, dpi = 300)
                           })
```
**Parameter Definitions:**

* `groups_colname`     - group column in metadata to be analyzed
* `assay_suffix`       - GeneLab assay suffix, default : "_GLAmpSeq"
* `output_prefix`      - additional prefix to be added to output files . Default: ""

**Input Data:**

* `metadata` (sample metadata with the group/treatment to be analyzed, output from [Preprocessing](#preprocessing))
* `feature_table` (ASV counts table, output from [Preprocessing](#preprocessing))
* `taxonomy_table` (taxonomy information, output from [Preprocessing](#preprocessing))
* `publication_format` (custom ggplot theme, output from [Set Variables](#set-variables))
* `custom_palette` (custom color palette, output from [Set Variables](#set-variables))


**Output Data:**

* **taxonomy_plots/<output_prefix>samples_<taxon_level>_GLAmpSeq.png** (samples barplots)
* **taxonomy_plots/<output_prefix>groups_<taxon_level>_GLAmpSeq.png** (groups barplots)

Where `taxon_level` is all of phylum, class, order, family, genus and species.

> Please note that the species plot can be misleading as short amplicon sequences can't be used to accurately predict species.

<br>

---


## 10. Differential Abundance Testing

Differential abundance testing aims to uncover specific taxa that exhibit notable variations across different conditions, complemented by visualizations like volcano plots to illustrate these disparities and their implications on ASV expression and overall microbial community dynamics. ANCOMBC 1, ANCOMBC 2, and DESeq2 provide 3 different methods for calculating differential abundance.

> Please note that if you'd like to run the code in this section, make sure that you [load the libraries](#load-libraries) 
and [functions](#load-functions), [read-in input tables](#read-in-input-tables) and [preprocess](#preprocessing) them in R 
by running the lines of code in [section 6](#6-amplicon-seq-data-analysis-set-up) sequentially, particularly those in [section 6b](#6b-r-environment-set-up).

### 10a. ANCOMBC 1

```R
# Create output directory if it doesn't already exist
diff_abund_out_dir <- "differential_abundance/ancombc1/"
if(!dir.exists(diff_abund_out_dir)) dir.create(diff_abund_out_dir, recursive = TRUE)
metadata  <-  {DATAFRAME}  
feature_table <- {DATAFRAME}
taxonomy_table <- {DATAFRAME}
publication_format <- {GGPLOT_THEME}
feature <- "ASV"
groups_colname  <- "groups"
samples_column  <- "Sample Name"
assay_suffix  <- "_GLAmpSeq"
target_region <- "16S" # "16S", "18S" or "ITS"
output_prefix  <- ""
prevalence_cutoff <- 0
library_cutoff <- 0
remove_struc_zero <- FALSE
threads <- 5


# Create phyloseq object from feature, taxonomy and sample metadata tables
ps <- phyloseq(otu_table(feature_table, taxa_are_rows = TRUE),
               sample_data(metadata),
               tax_table(as.matrix(taxonomy_table)))

# Convert phyloseq to tree summarized experiment object
tse <-  mia::makeTreeSummarizedExperimentFromPhyloseq(ps)


# Get unique group comparison as a matrix
pairwise_comp.m <- utils::combn(metadata[, groups_colname] %>% unique, 2)
pairwise_comp_df <- pairwise_comp.m %>% as.data.frame 
# Name the columns in the pairwise matrix as group1vgroup2
colnames(pairwise_comp_df) <- map_chr(pairwise_comp_df,
                                      \(col) str_c(col, collapse = "v"))
comparisons <- colnames(pairwise_comp_df)
names(comparisons) <- comparisons

# Write out contrasts table
write_csv(x = pairwise_comp_df,
          file =  glue("{diff_abund_out_dir}{output_prefix}contrasts{assay_suffix}.csv"),
          col_names = FALSE)


# Running ANCOMBC 1
set.seed(123)
final_results_bc1  <- map(pairwise_comp_df, function(col){
  
  group1 <- col[1]
  group2 <- col[2]
  
  # Subset the treeSummarizedExperiment object to contain only samples
  # in group1 and group2 
  tse_sub <-  tse[, tse[[groups_colname]] %in% c(group1, group2)]
  
  # Note that by default, levels of a categorical variable in R are sorted 
  # alphabetically. 
  # Changing the reference group by reordering the factor levels
  tse_sub[[groups_colname]] <- factor(tse_sub[[groups_colname]] , levels = c(group1, group2))
  
  # data - input data. TreeSummarizedExperiment or Phyloseq object
  # assay_name - name of count table in the input data object.
  # tax_level - taxonomy level for aggregation and analysis
  # prv_cut - prevalence cut-off. proportion of samples in which taxon is present.
  # lib_cut - a numerical threshold for filtering samples based on library sizes.
  # p_adj_method - p-value adjustment method for multiple comparisons
  # struc_zero - should group-wise rare taxa be detected
  # neg_lb - whether to classify a taxon as a structural zero using its asymptotic lower bound. i.e.the best the algorithm can possibly achieve 
  # group - name of the group variable in metadata. Only important you'd like to perform global test  can be set to NULL.
  # alpha - significance level
  # n_cl - number of processes to run in parallel
  # global - should a global test be performed to detect significant differences between at least 2 groups (ANOVA-like comparison) 
  # tol - iteration convergence tolerance for the E-M algorithm.
  # max_iter - max iteration 
  # formula - fixed effects formula
  # conserve - should a conservative variance estimator be used for the test statistic? 
  # it is recommended to set to TRUE if your sample size is small and the number of expected differentially abundant taxa is large.
  
  out <-  ancombc(data = tse_sub, assay_name = "counts", 
                  tax_level = NULL, phyloseq = NULL, 
                  formula = groups_colname, 
                  p_adj_method = "fdr", prv_cut = prevalence_cutoff,
                  lib_cut = library_cutoff, 
                  group = groups_colname , struc_zero = remove_struc_zero,
                  neg_lb = TRUE, tol = 1e-5, 
                  max_iter = 100, conserve = TRUE, alpha = 0.05, global = FALSE,
                  n_cl = threads, verbose = TRUE)
  
  # ------ Set data frame names ---------# 
  # LFC 
  lfc <- out$res$lfc %>%
    as.data.frame() %>% 
    select(-contains("Intercept")) %>% 
    set_names(
      c("taxon",
        glue("LnFC_({group2})v({group1})"))
    )
  
  # SE
  se <- out$res$se %>%
    as.data.frame() %>% 
    select(-contains("Intercept")) %>%
    set_names(
      c("taxon",
        glue("lfcSE_({group2})v({group1})"))
    )
  
  # W    
  W <- out$res$W %>%
    as.data.frame() %>% 
    select(-contains("Intercept")) %>%
    set_names(
      c("taxon",
        glue("Wstat_({group2})v({group1})"))
    )
  
  # p_val
  p_val <- out$res$p_val %>%
    as.data.frame() %>% 
    select(-contains("Intercept")) %>%
    set_names(
      c("taxon",
        glue("pvalue_({group2})v({group1})"))
    )
  
  # q_val
  q_val <- out$res$q_val %>%
    as.data.frame() %>% 
    select(-contains("Intercept")) %>% 
    set_names(
      c("taxon",
        glue("qvalue_({group2})v({group1})"))
    )
  
  
  # Diff_abn
  diff_abn <- out$res$diff_abn %>%
    as.data.frame() %>% 
    select(-contains("Intercept")) %>%
    set_names(
      c("taxon",
        glue("diff_({group2})v({group1})"))
    )
  
  # Merge the dataframes to one results dataframe
  res <- lfc %>%
    left_join(se) %>%
    left_join(W) %>% 
    left_join(p_val)  %>% 
    left_join(q_val) %>% 
    left_join(diff_abn)
  
  
  return(res)
  
})



# ------------ Create merged stats pairwise dataframe ----------------- #
# Initialize the merged stats dataframe to contain the taxon column for joining
merged_stats_df <- final_results_bc1[[names(final_results_bc1)[1]]] %>%
  as.data.frame() %>% select(taxon)

# Loop over the results of every comparion and join it the pre-existing 
# stats table
walk(comparisons[names(final_results_bc1)], .f = function(comparison){
  
  # Get comparison specific statistics
  df <-  final_results_bc1[[comparison]] %>% as.data.frame()
  
  # Merge it to the pre-existing statistics table
  merged_stats_df <<- merged_stats_df %>%
    dplyr::full_join(df, by = join_by("taxon"))
  
})

# Sort ASVs in ascending order
merged_stats_df <- merged_stats_df %>% 
  rename(!!feature := taxon) %>%
  mutate(!!feature := SortMixed(!!sym(feature)))


# ------ Get comparison names
# Since all column groups i.e. LnFC, pval, W, etc. have the same
# suffixes as comparison names, we only need to extract the comparion names
# from one of them. Here we extract them from the "LnFC" prefixed columns
comp_names <- merged_stats_df %>% 
  select(starts_with("LnFC")) %>%
  colnames() %>% str_remove_all("LnFC_")
names(comp_names) <- comp_names


# -------------- Make volcano plots ------------------ #
volcano_plots <- map(comp_names, function(comparison){
  
  # Construct column names for columns to be selected
  comp_col  <- c(
    glue("LnFC_{comparison}"),
    glue("lfcSE_{comparison}"),
    glue("Wstat_{comparison}"),
    glue("pvalue_{comparison}"),
    glue("qvalue_{comparison}"),
    glue("diff_{comparison}")
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
    str_remove_all("\\(|\\)") %>% # remove brackets
    str_split(".vs.") %>% unlist # split groups to list then convert to a vector
  
  group1 <- groups_vec[1]
  group2 <- groups_vec[2]
  
  ###### Long x-axis label adjustments ##########
  x_label <- glue("Ln Fold Change\n\n( {group1} vs {group2} )")
  label_length <- nchar(x_label)
  max_allowed_label_length <- plot_width_inches * 10
  
  # Construct x-axis label with new line breaks if was too long
  if (label_length > max_allowed_label_length){
    x_label <- glue("Ln Fold Change\n\n( {group1} \n vs \n {group2} )")
  }

  
  # Make plot
  p <- ggplot(sub_res_df %>% mutate(diff = qvalue <= p_val), 
              aes(x=LnFC, y=-log10(qvalue), 
                  color=diff, label=!!sym(feature))) +
    geom_point(alpha=0.7, size=2) +
    scale_color_manual(values=c("TRUE"="red", "FALSE"="black"),
                       labels=c(paste0("qval > ", p_val), 
                                paste0("qval \u2264 ", p_val))) +
    geom_hline(yintercept = -log10(p_val), linetype = "dashed") +
    ggrepel::geom_text_repel(show.legend = FALSE) +
    expandy(-log10(sub_res_df$qvalue)) + # Expand plot y-limit
    coord_cartesian(clip = 'off') +
    scale_y_continuous(oob = scales::oob_squish_infinite) + # prevent clipping of infinite values
    labs(x= x_label, y="-log10(Q-value)", 
         title = "Volcano Plot", color=NULL,
         caption = glue("dotted line: q-value = {p_val}")) + 
    theme_bw() +
    theme(legend.position="top", legend.key = element_rect(colour=NA),
          plot.caption = element_text(face = 'bold.italic'))
  # Save plot
  file_name <-  glue("{output_prefix}{comparison %>% str_replace_all('[:space:]+','_')}_volcano.png")
  ggsave(filename = file_name,
         plot = p, device = "png", width = plot_width_inches,
         height = plot_height_inches, units = "in",
         dpi = 300, path = diff_abund_out_dir)
  
  
  return(p)
})


# ------------------- Add NCBI id to feature, i.e. ASV -------------- #
# Get the best/least possible taxonomy name for the ASVs
tax_names <- map_chr(str_replace_all(taxonomy_table$species, ";_","")  %>%
                       str_split(";"),
                     function(row) row[length(row)])

df <- data.frame(ASV=rownames(taxonomy_table), best_taxonomy=tax_names)
colnames(df) <- c(feature, "best_taxonomy")

# Pull NCBI IDS for unique taxonomy names
df2 <- data.frame(best_taxonomy = df$best_taxonomy %>%
                    unique()) %>%
  mutate(NCBI_id=get_ncbi_ids(best_taxonomy, target_region),
         .after = best_taxonomy)

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
    mutate(!!mean_col := mean(c_across(where(is.numeric))),
           !!std_col := sd(c_across(where(is.numeric))) ) %>% 
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
  mutate(All.Mean=mean(c_across(where(is.numeric))),
         All.Stdev=sd(c_across(where(is.numeric))) ) %>%
  select(!!feature, All.Mean, All.Stdev)

# Merge the taxonomy table to the stats table
merged_df <- df  %>%
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
  mutate(across(where(is.numeric), ~round(.x, digits=3))) %>%
  mutate(across(where(is.matrix), as.numeric))

# Write out results of differential abundance using ANCOMBC 1
output_file <- glue("{diff_abund_out_dir}/{output_prefix}ancombc1_differential_abundance{assay_suffix}.csv")
# Write combined table to file but before that drop
# all columns of inferred differential abundance by ANCOMBC 
write_csv(merged_df %>%
            select(-starts_with("diff_")),
          output_file)

```

**Parameter Definitions:**

* `feature`            - feature type, i.e. ASV or OTU.
* `groups_colname`     - group column in metadata to be analyzed
* `samples_column`     - specifies the column in metadata containing the sample names in the feature table
* `assay_suffix`       - GeneLab assay suffix, default : "_GLAmpSeq"
* `output_prefix`      - additional prefix to be added to output files . Default: ""
* `threads`            - specifies the number of cpus to use for parallel processing
* `prevalence_cutoff`  - If `remove_rare` is true, a numerical fraction between 0 and 1. Taxa with prevalences (the 
                         proportion of samples in which the taxon is present) less than this will be excluded from 
                         the analysis. Default is 0, i.e. do not exclude any taxa / features.
* `library_cutoff`     - If `remove_rare` is true, a numerical threshold for filtering samples based on library sizes. 
                         Samples with library sizes less than lib_cut will be excluded in the analysis. Default is 0, 
                         i.e. no sample will be dropped. To discard samples with read counts less than or equal to 100, 
                         set to 100.
* `target_region`      - amplicon target region. Options are either "16S", "18S", or "ITS".
* `remove_struc_zero`  - Should structural zeros (a.k.a ASVs with zero count in at least one group) be removed ?
                         Default is FALSE i.e. structural zeros won't be removed.

**Input Data:**

* `metadata` (sample metadata with the group/treatment to be analyzed, output from [Preprocessing](#preprocessing))
* `feature_table` (ASV counts table, output from [Preprocessing](#preprocessing))
* `taxonomy_table` (taxonomy information, output from [Preprocessing](#preprocessing))
* `publication_format` (custom ggplot theme, output from [Set Variables](#set-variables))

**Output Data:**

* **differential_abundance/ancombc1/<output_prefix>(\<group1\>)v(\<group2\>)_volcano.png** (Comparion Volcano Plot)
* **differential_abundance/ancombc1/<output_prefix>ancombc1_differential_abundance_GLAmpSeq.csv** (Statistics Table)
<br>

---

### 10b. ANCOMBC 2

```R
diff_abund_out_dir <- "differential_abundance/ancombc2/"
if(!dir.exists(diff_abund_out_dir)) dir.create(diff_abund_out_dir, recursive = TRUE)
metadata  <-  {DATAFRAME} 
feature_table <- {DATAFRAME}
taxonomy_table <- {DATAFRAME}
publication_format <- {GGPLOT_THEME}
feature <- "ASV"
target_region <- "16S" # "16S" , "18S" or "ITS"
groups_colname  <- "groups"
samples_column  <- "Sample Name"
assay_suffix  <- "_GLAmpSeq"
output_prefix  <- ""
prevalence_cutoff <- 0
library_cutoff <- 0
remove_struc_zero <- FALSE
threads <- 5


# Create phyloseq object
ps <- phyloseq(otu_table(feature_table, taxa_are_rows = TRUE),
               sample_data(metadata),
               tax_table(as.matrix(taxonomy_table)))

# Convert phyloseq to tree summarized experiment object
tse <-  mia::makeTreeSummarizedExperimentFromPhyloseq(ps)

# Getting the reference group and making sure that it is the reference 
# used in the analysis
group_levels <- metadata[, groups_colname] %>% unique() %>% sort()
ref_group <- group_levels[1] # the first group is used as the reference group by default
tse[[groups_colname]] <- factor(tse[[groups_colname]] , levels = group_levels)

# Running ANCOMBC2....
# Run acombc2

  # data - input data. TreeSummarizedExperiment or Phyloseq object
  # assay_name - name of count table in the input data object.
  # tax_level - taxonomy level for aggregation and analysis
  # prv_cut - prevalence cut-off. proportion of samples in which taxon is present.
  # lib_cut - a numerical threshold for filtering samples based on library sizes.
  # p_adj_method - p-value adjustment method for multiple comparisons
  # struc_zero - should group-wise rare taxa be detected
  # neg_lb - whether to classify a taxon as a structural zero using its asymptotic lower bound, i.e.the best the algorithm can possibly achieve 
  # group - name of the group variable in metadata. Only important you'd like to perform global test  can be set to NULL.
  # alpha - significance level
  # n_cl - number of processes to run in parallel
  # global - should a global test be performed to detect significant differences between at least 2 groups (ANOVA-like comparison) 
  # tol - iteration convergence tolerance for the E-M algorithm.
  # max_iter - max iteration 
  # formula - fixed effects formula
  # conserve - should a conservative variance estimator be used for the test statistic? 
  # it is recommended to set to TRUE if your sample size is small and the number of expected differentially abundant taxa is large.

output <- ancombc2(data = tse, assay_name = "counts", tax_level = NULL,
                   fix_formula = groups_colname, rand_formula = NULL,
                   p_adj_method = "fdr", pseudo_sens = TRUE,
                   prv_cut = prevalence_cutoff, 
                   lib_cut = library_cutoff, s0_perc = 0.05,
                   group = groups_colname, struc_zero = remove_struc_zero,
                   neg_lb = FALSE, alpha = 0.05, n_cl = threads, verbose = TRUE,
                   global = TRUE, pairwise = TRUE, 
                   dunnet = TRUE, trend = FALSE,
                   iter_control = list(tol = 1e-5, max_iter = 20,
                                       verbose = FALSE),
                   em_control = list(tol = 1e-5, max_iter = 100),
                   mdfdr_control = list(fwer_ctrl_method = "fdr", B = 100), 
                   lme_control = NULL, trend_control = NULL)



# Create new column names - the original column names given by ANCOMBC are
# difficult to understand
new_colnames <- map_chr(output$res_pair  %>% colnames, 
                        function(colname) {
                          # Columns comparing a group to the reference group
                          if(str_count(colname,groups_colname) == 1){
                            str_replace_all(string=colname, 
                                            pattern=glue("(.+)_{groups_colname}(.+)"),
                                            replacement=glue("\\1_(\\2)v({ref_group})")) %>% 
                            str_replace(pattern = "^lfc_", replacement = "LnFC_") %>% 
                            str_replace(pattern = "^se_", replacement = "lfcSE_") %>% 
                            str_replace(pattern = "^W_", replacement = "Wstat_") %>%
                            str_replace(pattern = "^p_", replacement = "pvalue_") %>%
                            str_replace(pattern = "^q_", replacement = "qvalue_")
                            
                          # Columns with normal two groups comparison
                          } else if(str_count(colname,groups_colname) == 2){
                            
                            str_replace_all(string=colname, 
                                            pattern=glue("(.+)_{groups_colname}(.+)_{groups_colname}(.+)"),
                                            replacement=glue("\\1_(\\2)v(\\3)")) %>% 
                            str_replace(pattern = "^lfc_", replacement = "LnFC_") %>% 
                            str_replace(pattern = "^se_", replacement = "lfcSE_") %>% 
                            str_replace(pattern = "^W_", replacement = "Wstat_") %>%
                            str_replace(pattern = "^p_", replacement = "pvalue_") %>%
                            str_replace(pattern = "^q_", replacement = "qvalue_")
                            
                            # Feature/ ASV column 
                          } else{
                            
                            return(colname)
                          }
                        } )



# Change the column named taxon to the feature name e.g. ASV
new_colnames[match("taxon", new_colnames)] <- feature


# Round numeric values and rename columns
paired_stats_df <- output$res_pair  %>% 
  mutate(across(where(is.numeric), ~round(.x, digits=3))) %>%
  set_names(new_colnames)

# Get the unique comparison names 
uniq_comps <- str_replace_all(new_colnames, ".+_(\\(.+\\))", "\\1") %>% unique()
uniq_comps <- uniq_comps[-match(feature, uniq_comps)]

# Write out contrasts table
uniq_comps %>%
  str_replace_all("\\((.+)\\)v\\((.+)\\)", "\\1.vs.\\2") %>% 
  str_split(".vs.") %>% 
  map(.f = function(comparison) data.frame(comparison)) %>% 
  list_cbind() %>% 
  write_csv(file = glue("{diff_abund_out_dir}{output_prefix}contrasts{assay_suffix}.csv"),
            col_names = FALSE)


# ------ Sort columns by group comparisons --------#
# Create a data frame containing only the feature/ASV column
res_df <- paired_stats_df[1] 
walk(uniq_comps, function(comp){
  
  # Get the results for a comparison
  temp_df <- paired_stats_df %>% select(!!sym(feature), contains(comp))
  
  # Merge the current comparison to previous comparisons by feature/ASV id
  res_df <<- res_df %>% left_join(temp_df)
})

# --------- Add NCBI id to feature  ---------------#

# Get the best taxonomy assigned to each ASV
tax_names <- map_chr(str_replace_all(taxonomy_table$species, ";_","")  %>%
                       str_split(";"),
                     function(row) row[length(row)])

df <- data.frame(ASV=rownames(taxonomy_table), best_taxonomy=tax_names)
colnames(df) <- c(feature, "best_taxonomy")

# Querying NCBI...
# Pull NCBI IDS for unique taxonomy names
df2 <- data.frame(best_taxonomy = df$best_taxonomy %>%
                    unique()) %>%
  mutate(NCBI_id=get_ncbi_ids(best_taxonomy, target_region),
         .after = best_taxonomy)

df <- df %>%
  left_join(df2, join_by("best_taxonomy")) %>% 
  right_join(res_df)


# Retrieve the normalized table
normalized_table <- output$bias_correct_log_table  %>%
  rownames_to_column(feature) %>%
  mutate(across(where(is.numeric), ~replace_na(.x, replace=0)))

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
    mutate(!!mean_col := mean(c_across(where(is.numeric))),
           !!std_col := sd(c_across(where(is.numeric))) ) %>% 
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
  mutate(All.Mean=mean(c_across(where(is.numeric))),
         All.Stdev=sd(c_across(where(is.numeric))) ) %>% 
  select(!!feature, All.Mean, All.Stdev)


# Append the taxonomy table to the ncbi and stats table
merged_df <- df  %>%
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
  left_join(group_means_df, by = feature) %>% 
  mutate(across(where(is.numeric), ~round(.x, digits=3)))


# Writing out results of differential abundance using ANCOMBC2...
output_file <- glue("{diff_abund_out_dir}{output_prefix}ancombc2_differential_abundance{assay_suffix}.csv")
# Write out merged stats table but before that 
# drop ANCOMBC inferred diffrential abundance columns
write_csv(merged_df %>%
            select(-starts_with("diff_")),
          output_file)



# ---------------------- Visualization --------------------------------------- #
# Making volcano plots...
# ------------ Make volcano ---------------- #
volcano_plots <- map(uniq_comps, function(comparison){
  
  comp_col  <- c(
    glue("LnFC_{comparison}"),
    glue("lfcSE_{comparison}"),
    glue("Wstat_{comparison}"),
    glue("pvalue_{comparison}"),
    glue("qvalue_{comparison}"),
    glue("diff_{comparison}"),
    glue("passed_ss_{comparison}")
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
  x_label <- glue("Ln Fold Change\n\n( {group1} vs {group2} )")
  label_length <- nchar(x_label)
  max_allowed_label_length <- plot_width_inches * 10
  
  # Construct x-axis label with new line breaks if was too long
  if (label_length > max_allowed_label_length){
    x_label <- glue("Ln Fold Change\n\n( {group1} \n vs \n {group2} )")
  }
  #######################################
  
  
  
  p <- ggplot(sub_res_df %>% mutate(diff = qvalue <= p_val),
              aes(x=LnFC, y=-log10(qvalue), color=diff, label=!!sym(feature))) +
    geom_point(alpha=0.7, size=2) +
    scale_color_manual(values=c("TRUE"="red", "FALSE"="black"),
                       labels=c(paste0("qval > ", p_val), 
                                paste0("qval \u2264 ", p_val))) +
    geom_hline(yintercept = -log10(p_val), linetype = "dashed") +
    ggrepel::geom_text_repel(show.legend = FALSE) + 
    expandy(-log10(sub_res_df$qvalue)) + # Expand plot y-limit
    coord_cartesian(clip = 'off') +
    scale_y_continuous(oob = scales::oob_squish_infinite) + # prevent clipping of infinite values
    labs(x= x_label, y="-log10(Q-value)", 
         title = "Volcano Plot", color=NULL,
         caption = glue("dotted line: q-value = {p_val}")) + 
    theme_bw() +
    theme(legend.position="top", legend.key = element_rect(colour=NA),
          plot.caption = element_text(face = 'bold.italic'))

  # Save plot
  file_name <-  glue("{output_prefix}{comparison %>% str_replace_all('[:space:]+','_')}_volcano.png")
  ggsave(filename = file_name,
         plot = p, device = "png",
         width = plot_width_inches,
         height = plot_height_inches,
         units = "in", dpi = 300, path = diff_abund_out_dir)
 
  
  return(p)
  
})

```

**Parameter Definitions:**

* `feature`            - feature type, i.e. ASV or OTU.
* `samples_column`     – column in metadata containing the sample names in the feature table
* `groups_colname`     - group column in metadata to be analyzed
* `assay_suffix`       - GeneLab assay suffix, default : "_GLAmpSeq"
* `output_prefix`      - additional prefix to be added to output files . Default: ""
* `threads`            - specifies the number of cpus to use for parallel processing.
* `prevalence_cutoff`  - If `remove_rare` is true, a numerical fraction between 0 and 1. Taxa with prevalences (the 
                         proportion of samples in which the taxon is present) less than this will be excluded from 
                         the analysis. Default is 0, i.e. do not exclude any taxa / features.
* `library_cutoff`     - If `remove_rare` is true, a numerical threshold for filtering samples based on library sizes. 
                         Samples with library sizes less than lib_cut will be excluded in the analysis. Default is 0, 
                         i.e. no sample will be dropped. To discard samples with read counts less than or equal to 100, 
                         set to 100.
* `target_region`      - amplicon target region. Options are either "16S", "18S", or "ITS"
* `remove_struc_zero`  - Should structural zeros (a.k.a ASVs with zero count in at least one group) be removed ?
                         Default is FALSE i.e. structural zeros won't be removed.

**Input Data:**

* `metadata` (sample metadata with the group/treatment to be analyzed, output from [Preprocessing](#preprocessing))
* `feature_table` (ASV counts table, output from [Preprocessing](#preprocessing))
* `taxonomy_table` (taxonomy information, output from [Preprocessing](#preprocessing))
* `publication_format` (custom ggplot theme, output from [Set Variables](#set-variables))

**Output Data:**

* **differential_abundance/ancombc2/<output_prefix>(\<group1\>)v(\<group2\>)_volcano.png** (Comparion Volcano Plot)
* **differential_abundance/ancombc2/<output_prefix>ancombc2_differential_abundance_GLAmpSeq.csv** (Statistics Table)


<br>

---

### 10c. DESeq2

```R
# Create output directory if it doesn't already exist
diff_abund_out_dir <- "differential_abundance/deseq2/"
if(!dir.exists(diff_abund_out_dir)) dir.create(diff_abund_out_dir, recursive = TRUE)
metadata  <-  {DATAFRAME}
feature_table <- {DATAFRAME}
taxonomy_table <- {DATAFRAME}
publication_format <- {GGPLOT_THEME}
feature <- "ASV"
samples_column <- "Sample Name"
groups_colname  <- "groups" 
assay_suffix  <- "_GLAmpSeq"
target_region <- "16S" # "16S", "18S" or "ITS"
output_prefix  <- ""
prevalence_cutoff <- 0
library_cutoff <- 0


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

# Run Deseq
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
  DESeq(deseq_obj)
})


# Get unique group comparison as a matrix
pairwise_comp.m <- utils::combn(metadata[,groups_colname] %>% unique, 2)
pairwise_comp_df <- pairwise_comp.m %>% as.data.frame 
# Set the colnames as group1vgroup2
colnames(pairwise_comp_df) <- map_chr(pairwise_comp_df,
                                      \(col) str_c(col, collapse = "v"))
comparisons <- colnames(pairwise_comp_df)
names(comparisons) <- comparisons

# Write out contrasts table
write_csv(x = pairwise_comp_df,
          file =  glue("{diff_abund_out_dir}{output_prefix}contrasts{assay_suffix}.csv"),
          col_names = FALSE)


# Retrieve statistics table
merged_stats_df <-  data.frame(ASV=rownames(feature_table))
colnames(merged_stats_df) <- feature

walk(pairwise_comp_df, function(col){
  
  group1 <- col[1]
  group2 <- col[2]

# Retrieve the statistics table for the cuurrent pair and rename the columns  
df <- results(deseq_modeled, contrast = c(groups_colname, group1, group2)) %>% # Get stats
  data.frame() %>%
  rownames_to_column(feature) %>% 
  set_names(c(feature ,
              glue("baseMean_({group1})v({group2})"),
              glue("log2FC_({group1})v({group2})"),
              glue("lfcSE_({group1})v({group2})"), 
              glue("stat_({group1})v({group2})"), 
              glue("pvalue_({group1})v({group2})"),
              glue("padj_({group1})v({group2})") 
            )) # rename the columns

            
  merged_stats_df <<- merged_stats_df %>% 
                          dplyr::left_join(df, join_by(!!feature))
})

# ---------------------- Add NCBI id to feature, i.e. ASV

# Get the best / lowest possible taxonomy assignment for the features, i.e. ASVs
tax_names <- map_chr(str_replace_all(taxonomy_table$species, ";_","")  %>%
                       str_split(";"),
                     function(row) row[length(row)])

df <- data.frame(ASV=rownames(taxonomy_table), best_taxonomy=tax_names)
colnames(df) <- c(feature, "best_taxonomy")

# Pull NCBI IDS for unique taxonomy names
df2 <- data.frame(best_taxonomy = df$best_taxonomy %>%
                    unique()) %>%
  mutate(NCBI_id=get_ncbi_ids(best_taxonomy, target_region),
         .after = best_taxonomy)

df <- df %>%
  left_join(df2, join_by("best_taxonomy")) %>% 
  right_join(merged_stats_df)

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
    mutate(!!mean_col := mean(c_across(where(is.numeric))),
           !!std_col := sd(c_across(where(is.numeric))) ) %>% 
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
  mutate(All.Mean=mean(c_across(where(is.numeric))),
         All.Stdev=sd(c_across(where(is.numeric))) ) %>% 
  select(!!feature, All.Mean, All.Stdev)



# Add taxonomy
merged_df <- df  %>% # statistics table
  left_join(taxonomy_table %>%
              as.data.frame() %>%
              rownames_to_column(feature)) %>% # append taxonomy table
  select(!!feature, domain:species,everything()) # select columns of interest

# Merge all prepared tables in the desired order
merged_df <- merged_df %>%
  select(!!sym(feature):NCBI_id) %>%  # select only the features and NCBI ids
  left_join(normalized_table, by = feature) %>% # append the normalized table
  left_join(merged_df) %>% # append the stats table
  left_join(All_mean_sd) %>% # append the global/ASV means and stds
  left_join(group_means_df, by = feature) %>% # append the group means and stds
  mutate(across(where(is.numeric), ~round(.x, digits=3))) %>%  # round numeric columns
  mutate(across(where(is.matrix), as.numeric)) # convert meatrix columns to numeric columns

# Defining the output file
output_file <- glue("{diff_abund_out_dir}/{output_prefix}deseq2_differential_abundance{assay_suffix}.csv")
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
  deseq_res <- results(deseq_modeled, contrast = c(groups_colname, group1, group2))
  volcano_data <- as.data.frame(deseq_res)
  volcano_data <- volcano_data[!is.na(volcano_data$padj), ]
  volcano_data$significant <- volcano_data$padj <= p_val
  
   ######Long x-axis label adjustments##########
  x_label <- glue("Log2 Fold Change\n\n( {group1} vs {group2} )")
  label_length <- nchar(x_label)
  max_allowed_label_length <- plot_width_inches * 10
  
  # Construct x-axis label with new line breaks if was too long
  if (label_length > max_allowed_label_length){
    x_label <- glue("Log2 Fold Change\n\n( {group1} \n vs \n {group2} )")
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
  ggsave(filename = glue("{output_prefix}({group1})v({group2})_volcano.png"),
         plot = p,
         width = plot_width_inches, 
         height = plot_height_inches, 
         dpi = 300, 
         path = diff_abund_out_dir)

})
```

**Parameter Definitions:**

* `feature`            - feature type, i.e. ASV or OTU.
* `samples_column`     – column in metadata containing the sample names in the feature table
* `groups_colname`     - group column in metadata to be analyzed
* `assay_suffix`       - GeneLab assay suffix, default : "_GLAmpSeq"
* `output_prefix`      - additional prefix to be added to output files . Default: ""
* `prevalence_cutoff`  - If `remove_rare` is true, a numerical fraction between 0 and 1. Taxa with prevalences (the 
                         proportion of samples in which the taxon is present) less than this will be excluded from 
                         the analysis. Default is 0, i.e. do not exclude any taxa / features.
* `library_cutoff`     - If `remove_rare` is true, a numerical threshold for filtering samples based on library sizes. 
                         Samples with library sizes less than lib_cut will be excluded in the analysis. Default is 0, 
                         i.e. no sample will be dropped. To discard samples with read counts less than or equal to 100, 
                         set to 100.
* `target_region`      - amplicon target region. Options are "16S", "18S", or "ITS".

**Input Data:**

* `metadata` (sample metadata with the group/treatment to be analyzed, output from [Preprocessing](#preprocessing))
* `feature_table` (ASV counts table, output from [Preprocessing](#preprocessing))
* `taxonomy_table` (taxonomy information, output from [Preprocessing](#preprocessing))
* `publication_format` (custom ggplot theme, output from [Set Variables](#set-variables))


**Output Data:**

* **differential_abundance/deseq2/<output_prefix>(\<group1\>)v(\<group2\>)_volcano.png** (Comparion Volcano Plot)
* **differential_abundance/deseq2/<output_prefix>deseq2_differential_abundance_GLAmpSeq.csv** (Statistics Table)
<br>

---
