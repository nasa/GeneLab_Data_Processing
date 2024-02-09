# Bioinformatics pipeline for amplicon Illumina sequencing data  

> **This page holds an overview and instructions for how GeneLab processes Illumina amplicon datasets. Exact processing commands for specific datasets that have been released are available in the [GLDS_Processing_Scripts](../GLDS_Processing_Scripts) sub-directory and/or are provided with their processed data in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/).**  

---

**Date:** December 14, 2023  
**Revision:** B  
**Document Number:** GL-DPPD-7104-B  

**Submitted by:**  
Alexis Torres and Michael D. Lee (GeneLab Data Processing Team)

**Approved by:**  
Sylvain Costes (GeneLab Project Manager)  
Samrawit Gebre (GeneLab Deputy Project Manager and Acting Genelab Configuration Manager)  
Lauren Sanders (OSDR Project Scientist)  
Amanda Saravia-Butler (GeneLab Data Processing Lead)  

---

## Updates from previous version


- Additional software (R packages) used:
  - vegan
  - tidyverse
  - dendextend
  - ggrepel
  - dplyr
  - RColorBrewer
  - phyloseq
- Inclusion of additional steps and outputs starting from ([step 6](#6-beta-diversity)):
  - Beta Diversity with Hierarchical Clustering ([6a](#6a-hierarchical-clustering)) and Ordination ([6b](#6b-ordination)).
  - Alpha Diversity with Rarefaction Curves ([7a](#7a-rarefaction-curves)) and Richness and Diversity Estimates ([7b](#7b-richness-and-diversity-estimates)).
  - Groupwise and Samplewise Taxonomic Summary Plots ([step 8](#8-taxonomic-summaries)).
  - Differential Abundance Analysis ([step 9](#9-differential-abundance-analysis)) Including Betadisper, Permutational ANOVA ([9a](#9a-betadisper-and-permutational-anova)), DESeq2 ([9b](#9b-differential-abundance-analysis-with-deseq2)) and Volcano Plots ([9c](#9c-volcano-plots)).

<!-- Included R packages -->
- Assay-specific suffixes were added where needed for GeneLab repo ("GLAmpSeq")
- The ITS UNITE reference database used was updated to "UNITE_v2023_July2023.RData", from http://www2.decipher.codes/Classification/TrainingSets/
- Several program versions were updated (all versions listed in [Software used](#software-used) below)

# Table of contents  

- [Software used](#software-used)
- [Reference databases used](#reference-databases-used)
- [General processing overview with example commands](#general-processing-overview-with-example-commands)
  - [1. Raw Data QC](#1-raw-data-qc)
    - [1a. Compile Raw Data QC](#1a-compile-raw-data-qc)
  - [2. Trim Primers](#2-trim-primers)
  - [3. Quality Filtering](#3-quality-filtering)
  - [4. Filtered Data QC](#4-filtered-data-qc)
    - [4a. Compile Filtered Data QC](#4a-compile-filtered-data-qc)
  - [5. Calculate Error model, Apply DADA2 Algorithm, Assign Taxonomy, and Create Output Tables](#5-calculate-error-model-apply-dada2-algorithm-assign-taxonomy-and-create-output-tables)
    - [5a. Learning the Error Rates](#5a-learning-the-error-rates)
    - [5b. Inferring Sequences](#5b-inferring-sequences)
    - [5c. Merging Forward and Reverse Reads; Not Needed if Data are Single-End](#5c-merging-forward-and-reverse-reads-not-needed-if-data-are-single-end)
    - [5d. Generating Sequence Table with Counts per Sample](#5d-generating-sequence-table-with-counts-per-sample)
    - [5e. Removing Putative Chimeras](#5e-removing-putative-chimeras)
    - [5f. Assigning Taxonomy](#5f-assigning-taxonomy)
    - [5g. Generating and Writing Standard Outputs](#5g-generating-and-writing-standard-outputs)
  - [6. Amplicon Seq Data Analysis Set Up](#6-amplicon-seq-data-analysis-set-up)
    - [6a. Create Sample Runsheet](#6a-create-sample-runsheet)
    - [6b. Environment Set Up](#6b-environment-set-up)
  - [7. Beta Diversity](#7-beta-diversity)
    - [7a. Hierarchical Clustering](#7a-hierarchical-clustering)
    - [7b. Ordination](#7b-ordination)
  - [8. Alpha Diversity](#8-alpha-diversity)
    - [8a. Rarefaction Curves](#8a-rarefaction-curves)
    - [8b. Richness and Diversity Estimates](#8b-richness-and-diversity-estimates)
  - [9. Taxonomic Summaries](#9-taxonomic-summaries)
  - [10. Differential Abundance Analysis](#10-differential-abundance-analysis)
    - [10a. Betadisper and Permutational ANOVA](#10a-betadisper-and-permutational-anova)
    - [10b. Differential Abundance Analysis with DESeq2](#10b-differential-abundance-analysis-with-deseq2)
    - [10c. Volcano Plots](#10c-volcano-plots)

---

# Software used  

|Program|Version*|Relevant Links|
|:------|:------:|:-------------|
|FastQC|`0.12.1`|[https://www.bioinformatics.babraham.ac.uk/projects/fastqc/](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)|
|MultiQC|`1.19`|[https://multiqc.info/](https://multiqc.info/)|
|Cutadapt|`4.6`|[https://cutadapt.readthedocs.io/en/stable/](https://cutadapt.readthedocs.io/en/stable/)|
|DADA2|`1.30.0`|[https://www.bioconductor.org/packages/release/bioc/html/dada2.html](https://www.bioconductor.org/packages/release/bioc/html/dada2.html)|
|DECIPHER|`2.30.0`|[https://bioconductor.org/packages/release/bioc/html/DECIPHER.html](https://bioconductor.org/packages/release/bioc/html/DECIPHER.html)|
|biomformat|`1.30.0`|[https://github.com/joey711/biomformat](https://github.com/joey711/biomformat)|
|R-base|`4.3.2`|[https://www.r-project.org/](https://www.r-project.org/)|
|vegan|`2.6.4`|[https://cran.r-project.org/package=vegan](https://cran.r-project.org/package=vegan)|
|tidyverse|`2.0.0`|[https://CRAN.R-project.org/package=tidyverse](https://CRAN.R-project.org/package=tidyverse)|
|dendextend|`1.17.1`|[https://CRAN.R-project.org/package=dendextend](https://CRAN.R-project.org/package=dendextend)|
|ggrepel|`0.9.4`|[https://CRAN.R-project.org/package=ggrepel](https://CRAN.R-project.org/package=ggrepel)|
|dplyr|`1.1.3`|[https://CRAN.R-project.org/package=dplyr](https://CRAN.R-project.org/package=dplyr)|
|rcolorbrewer|`1.1.3`|[https://CRAN.R-project.org/package=RColorBrewer](https://CRAN.R-project.org/package=RColorBrewer)|
|DESeq2|`1.40.2`|[https://bioconductor.org/packages/release/bioc/html/DESeq2.html](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)|
|phyloseq|`1.44.0`|[https://bioconductor.org/packages/release/bioc/html/phyloseq.html](https://bioconductor.org/packages/release/bioc/html/phyloseq.html)|

# Reference databases used

|Program used| Database| Relevant Links|
|:-----|:-----:|--------:|
|DECIPHER| SILVA SSU r138 | [http://www2.decipher.codes/Classification/TrainingSets/SILVA_SSU_r138_2019.RData](http://www2.decipher.codes/Classification/TrainingSets/)|
|DECIPHER| UNITE v2020 | [http://www2.decipher.codes/Classification/TrainingSets/UNITE_v2020_February2020.RData](http://www2.decipher.codes/Classification/TrainingSets/)|

---

# General processing overview with example commands  

> Exact processing commands for specific datasets are available in the [GLDS_Processing_Scripts](../GLDS_Processing_Scripts) sub-directory of this repository, and/or are provided with their processed data in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/).
>
> Output files listed in **bold** below are included with each Amplicon Seq processed dataset in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/).

---

## 1. Raw Data QC  

```
fastqc -o raw_fastqc_output *.fastq.gz
```

**Parameter Definitions:**

* `-o` – the output directory to store results
* `*.fastq.gz` – the input reads are specified as a positional argument, and can be given all at once with wildcards like this, or as individual arguments with spaces in between them

**Input Data:**

* fastq, compressed or uncompressed

**Output Data:**

* fastqc.html (FastQC output html summary)
* fastqc.zip (FastQC output data)


<br>  

### 1a. Compile Raw Data QC  

```
multiqc -o raw_multiqc_output raw_fastqc_output
# this is how it's packaged with our workflow outputs
zip -r raw_multiqc_GLAmpSeq_report.zip raw_multiqc_output
```

**Parameter Definitions:**

*	`-o` – the output directory to store results
*	`raw_fastqc_output/` – the directory holding the output data from the fastqc run, provided as a positional argument

**Input Data:**

* fastqc.zip (FastQC output data)

**Output Data:**

* **raw_multiqc_GLAmpSeq_report.zip** (zip containing the following)
  * **raw_multiqc.html** (multiqc output html summary)
  * **raw_multiqc_data** (directory containing multiqc output data)

<br>  

---

## 2. Trim Primers  

The location and orientation of primers in the data is important to understand in deciding how to do this step. `cutadapt` has many options for primer identification and removal. They are described in detail on their documentation page here: [https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types](https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types)  

The following example commands show how it was done for some samples of [GLDS-200](https://osdr.nasa.gov/bio/repo/data/studies/OSD-200), which was 2x250 sequencing of the 16S gene using these primers:  
* forward: 5’-GTGCCAGCMGCCGCGGTAA-3’  
* reverse: 5’-GGACTACVSGGGTATCTAAT-3’  

Due to the size of the target amplicon and the type of sequencing done here, both forward and reverse primers are expected to be on each of the forward and reverse reads. It therefore takes “linked” primers as input for forward and reverse reads, specified above by the `...` between them. It also expects that the primers start at the first position of the reads (“anchored”), specified with the leading `^` characters.  

The following website is useful for reverse complementing primers and dealing with degenerate bases appropriately: [http://arep.med.harvard.edu/labgc/adnan/projects/Utilities/revcomp.html](http://arep.med.harvard.edu/labgc/adnan/projects/Utilities/revcomp.html)  

```
cutadapt -a ^GTGCCAGCMGCCGCGGTAA...ATTAGATACCCSBGTAGTCC -A ^GGACTACVSGGGTATCTAAT...TTACCGCGGCKGCTGGCAC \
         ## Define what B represents; and define what K represents ##
         -o Primer-trimmed-R1.fq.gz -p Primer-trimmed-R2.fq.gz Input_R1_raw.fastq.gz Input_R2_raw.fastq.gz \
         --discard-untrimmed
```

**Parameter Definitions:**

*	`-a` – specifies the primers and orientations expected on the forward reads (when primers are linked as noted above)

*	`-A` – specifies the primers and orientations expected on the reverse reads (when primers are linked as noted above)

*	`-o` – specifies output of forward, primer-trimmed reads

*	`-p` – specifies output of reverse, primer-trimmed reads

-	`Input_R1_raw.fastq.gz` – this and following “R2” version are positional arguments specifying the forward and reverse reads, respectively, for input

-	`--discard-untrimmed` – this filters out those reads? where the primers were not found as expected

**Input Data:**

* fastq, compressed or uncompressed (original reads)

**Output Data:**

* **trimmed.fastq.gz**, compressed or uncompressed (trimmed reads)
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
filtered_out <- filterAndTrim(fwd=“Primer-trimmed-R1.fq.gz”, filt=“Filtered-R1.fq.gz”,
                              rev=“Primer-trimmed-R2.fq.gz”, filt.rev=“Filtered-R2.fq.gz”,
                              truncLen=c(220, 160), maxN=0, maxEE=c(2,2),
                              truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)
```

**Parameter Definitions:**

*	`filtered_out <-` – specifies the variable that will store the summary results within in our R environment

*	`filterAndTrim()` – the DADA2 function we are calling, with the following parameters set within it

*	`fwd=` – specifying the path to the forward reads, here “Primer-trimmed-R1.fq.gz”

*	`filt=` – specifying the path to where the output forward reads will be written

*	`rev=` – specifying the path to the reverse reads, here “Primer-trimmed-R2.fq.gz”; only applicable if paired-end

*	`filt.rev=` – specifying the path to where the output reverse reads will be written; only applicable if paired-end

*	`truncLen=c(220, 160)` – specifying the forward reads to be truncated at 220 bp, and the reverse to be truncated at 160 bps (note that this parameter also functions as a minimum-length filter); would only have 1 value if not paired-end

*	`maxN=0` – setting the maximum allowed Ns to 0, any reads with an N will be filtered out

*	`maxEE=c(2,2)` – setting maximum expected error allowed to 2 for each forward and reverse read; would only have 1 value if not paired-end

*	`truncQ=2` – looking from the lower-quality end of each read, truncate at the first base with a quality score lower than 2

*	`rm.phix=TRUE` – filter out reads with exact kmers matching the PhiX genome

*	`compress=TRUE` – gzip-compress the output filtered reads

*	`multithread=TRUE` – determine number of cores available and run in parallel when possible (can also take an integer specifying the number to run)

**Input Data:**

* fastq, compressed or uncompressed (primer-trimmed reads)

**Output Data:**

* **filtered.fastq.gz**, compressed or uncompressed (filtered reads)
* **filtered-read-counts_GLAmpSeq.tsv** (per sample read counts before and after filtering)

<br>

---

## 4. Filtered Data QC
```
fastqc -o filtered_fastqc_output/ filtered*.fastq.gz
```

**Parameter Definitions:**

*	`-o` – the output directory to store results  
*	`filtered*.fastq.gz` – the input reads are specified as a positional argument, and can be given all at once with wildcards like this, or as individual arguments with spaces in between them  

**Input Data:**

* fastq, compressed or uncompressed (filtered reads)

**Output Data:**

* fastqc.html (FastQC output html summary)
* fastqc.zip (FastQC output data)

<br>

### 4a. Compile Filtered Data QC
```
multiqc -o filtered_multiqc_output  filtered_fastqc_output
# this is how it's packaged with our workflow outputs
zip -r filtered_multiqc_GLAmpSeq_report.zip filtered_multiqc_output
```

**Parameter Definitions:**

*	`-o` – the output directory to store results
*	`filtered_fastqc_output` – the directory holding the output data from the fastqc run, provided as a positional argument

**Input Data:**

* fastqc.zip (FastQC output data)

**Output Data:**

* **filtered_multiqc_GLAmpSeq_report.zip** (zip containing the following)
  * **filtered_multiqc_report.html** (multiqc output html summary)
  * **filtered_multiqc_data** (directory containing multiqc output data)

<br>

---

## 5. Calculate Error Mdel, Apply DADA2 Algorithm, Assign Taxonomy, and Create Output Tables
> The following is run in an R environment.  

These example commands as written assumes paired-end data, with notes included on what would be different if working with single-end data. The taxonomy reference database used below is as an example only, suitable for the example 16S dataset ([GLDS-200](https://osdr.nasa.gov/bio/repo/data/studies/OSD-200)) used here. But others designed for DECIPHER can be found here: [http://www2.decipher.codes/Downloads.html](http://www2.decipher.codes/Downloads.html)  

<br>

### 5a. Learning the Error Rates
```R
forward_errors <- learnErrors(fls=“Filtered-R1.fq.gz”, multithread=TRUE)
```

**Parameter Definitions:**  

*	`forward_errors <-` – specifies the variable that will store the results within in our R environment

*	`learnErrors()` – the DADA2 function we are calling, with the following parameters set within it

*	`fls=` – the path to the forward filtered reads

*	`multithread=TRUE` – determine number of cores available and run in parallel when possible (can also take an integer specifying the number to run)

```R
reverse_errors <- learnErrors(fls=“Filtered-R2.fq.gz”, multithread=TRUE)
```

**Parameter Definitions:**  

*	same as above, but providing reverse filtered reads; not needed if data are single-end

<br>

### 5b. Inferring Sequences
```R
forward_seqs <- dada(derep=“Filtered-R1.fq.gz”, err=forward_errors, pool=“pseudo”, multithread=TRUE)
```

**Parameter Definitions:**  

*	`forward_seqs <-` – specifies the variable that will store the results within in our R environment

*	`dada()` – the DADA2 function we are calling, with the following parameters set within it

*	`derep=` – the path to the forward filtered reads

*	`err=` – the object holding the error profile for the forward reads created in above step, if not paired-end data, this would be the error-profile object created and the following “reverse_seqs” object would not be created

*	`pool=“pseudo”` – setting the method of incorporating information from multiple samples

*	`multithread=TRUE` – determine number of cores available and run in parallel when possible (can also take an integer specifying the number to run)

<br>

```R
reverse_seqs <- dada(derep=“Filtered-R2.fq.gz”, err=reverse_errors, pool=“pseudo”, multithread=TRUE)
```

**Parameter Definitions:**  

*	same as above, but providing reverse filtered reads and reverse error-profile object; not needed if data are single-end

<br>

### 5c. Merging Forward and Reverse Reads; Not Needed if Data are Single-End
```R
merged_contigs <- mergePairs(dadaF=forward_seqs, derepF=“Filtered-R1.fq.gz”, dadaR=reverse_seqs, derepR=“Filtered-R2.fq.gz”)
```

**Parameter Definitions:** 

*	`merged_contigs <-` – specifies the variable that will store the results within in our R environment

*	`mergePairs()` – the DADA2 function we are calling, with the following parameters set within it

*	`dadaF=` – specifying the object holding the forward-read inferred sequences

*	`derepF=` – specifying the path to the filtered forward reads

*	`dadaR=` – specifying the object holding the reverse-read inferred sequences

*	`derepR=` – specifying the path to the filtered reverse reads

<br>

### 5d. Generating Sequence Table with Counts per Sample
```R
seqtab <- makeSequenceTable(merged_contigs)
```

**Parameter Definitions:**  

*	If single-end data, instead of “merged_contigs”, the forward_seqs object would be provided to the `makeSequenceTable()` function here

<br>

### 5e. Removing putative chimeras
```R
seqtab.nochim <- removeBimeraDenovo(unqs=seqtab, method=“consensus”, multithread=TRUE)
```

**Parameter Definitions:**  

*	`seqtab.nochim <-` – specifies the variable that will store the results within in our R environment

*	`removeBimeraDenovo()` – the DADA2 function we are calling, with the following parameters set within it

*	`unqs=` – specifying the “seqtab” object created above

*	`method=` – specifying the method for putative-chimera identification and removal

*	`multithread=TRUE` – determine number of cores available and run in parallel when possible (can also take an integer specifying the number to run)

<br>

### 5f. Assigning Taxonomy

Creating a DNAStringSet object from the ASVs:
```R
dna <- DNAStringSet(getSequences(seqtab.nochim))
```

Downloading the reference R taxonomy object:
```R
download.file( url=“http://www2.decipher.codes/Classification/TrainingSets/SILVA_SSU_r138_2019.RData”, destfile=“SILVA_SSU_r138_2019.RData”)
```

**Parameter Definitions:**  

*	`download.file()` – the function we are calling, with the following parameters set within it

*	`url=` – specifying the url address of the file to download

*	`destfile=` – specifying the path/name of the file after downloading

<br>

Loading taxonomy object:
```R
load(“SILVA_SSU_r138_2019.RData”)
```

Classifying sequences:
```R
tax_info <- IdTaxa(test=dna, trainingSet=trainingSet, strand=“both”, processors=NULL)
```

**Parameter Definitions:**  

*	`tax_info <-` – specifies the variable that will store the results within in our R environment

*	`IdTaxa()` – the DECIPHER function we are calling, with the following parameters set within it

*	`test=` – specifying the “dna” object created above holding the sequences we want to classify

*	`trainingSet=` – specifying the reference database we downloaded and loaded above

*	`strand=“both”` – specifying to check taxonomy assignment in both orientations 

*	`processors=NULL` – determine number of cores available and run in parallel when possible (can also take an integer specifying the number to run)

<br>

### 5g. Generating and Writing Standard Outputs

Giving sequences more manageable names (e.g. ASV_1, ASV_2, …,):
```R
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}
```

Making then writing a fasta of final ASV seqs:
```R
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs_GLAmpSeq.fasta")
```

Making then writing a count table:
```R
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)

write.table(asv_tab, "counts_GLAmpSeq.tsv", sep="\t", quote=F, col.names=NA)
```

Creating table of taxonomy and setting any that are unclassified as "NA":
```R
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
```

Generating then writing biom file format:
```R
biom_object <- make_biom(data=asv_tab, observation_metadata=tax_tab)
write_biom(biom_object, "taxonomy-and-counts_GLAmpSeq.biom")
```

Making a combined taxonomy and count table
```R
tax_and_count_tab <- merge(tax_tab, asv_tab)
write.table(tax_and_count_tab, "taxonomy-and-counts_GLAmpSeq.tsv", sep="\t", quote=FALSE, row.names=FALSE)
```

**Input Data:**

* fastq, compressed or uncompressed (filtered reads)

**Output Data:**

* `tax_tab` (variable containing the taxonomy table)
* **ASVs_GLAmpSeq.fasta** (inferred sequences)
* **counts_GLAmpSeq.tsv** (count table)
* **taxonomy_GLAmpSeq.tsv** (taxonomy table)
* **taxonomy-and-counts_GLAmpSeq.tsv** (combined taxonomy and count table)
* **taxonomy-and-counts_GLAmpSeq.biom** (count and taxonomy table in biom format)
* **read-count-tracking_GLAmpSeq.tsv** (read counts at each processing step)

<br>

---

## 6. Amplicon Seq Data Analysis Set Up

> The remainder of this document is performed in R.  
  
<br>

### 6a. Create Sample Runsheet

> Note: Rather than running the command below to create the runsheet needed for processing, the runsheet may also be created manually by following the [file specification](../Workflow_Documentation/SW_AmpIllumina-B/examples/runsheet/README.md).

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

- `--accession OSD-###` - OSD accession ID (replace ### with the OSD number being processed), used to retrieve the urls for the ISA archive and raw reads hosted on the OSDR
- `--config-type` - instructs the script to extract the metadata required for Amplicon Sequencing data processing from the ISA archive
- `--config-version` - specifies the `dp-tools` configuration version to use, a value of `Latest` will specify the most recent version
- `--isa-archive` - specifies the *ISA.zip file for the respective OSD dataset, downloaded in the `dpt-get-isa-archive` command


**Input Data:**

- No input data required but the OSD accession ID needs to be indicated, which is used to download the respective ISA archive 

**Output Data:**

- *ISA.zip (compressed ISA directory containing Investigation, Study, and Assay (ISA) metadata files for the respective OSD dataset, used to define sample groups - the *ISA.zip file is located in the [OSDR](https://osdr.nasa.gov/bio/repo/) under 'Files' -> 'Study Metadata Files')

- **{OSD-Accession-ID}_AmpSeq_v{version}_runsheet.csv** (table containing metadata required for processing, version denotes the dp_tools schema used to specify the metadata to extract from the ISA archive)

<br>

### 6b. Environment Set Up

```R
### Import libraries used for processing ###

library(tidyverse)
library(phyloseq)
library(vegan)
library(dendextend)
library(DESeq2)


### Read in the runsheet containing the metadata required for processing ###

runsheet <- read.table(file = "*runsheet.csv", 
                        header = TRUE, row.names = 1, sep = ",")


### Create a sample info table from the runsheet containing the group names and a column of unique colors for each group ###

num_colors <- length(unique(runsheet$groups))
if (num_colors > 9) {
    custom_palette <- colorRampPalette(brewer.pal(9, "Set1"))(num_colors)
    colors <- custom_palette
} else {
    colors <- brewer.pal(num_colors, "Set1")
}
group_colors <- setNames(colors, unique(runsheet$groups))
runsheet <- runsheet %>%
  mutate(!!color_colname := group_colors[.data$groups])
sample_info_tab <- runsheet[, c($groups, $color)]


### Read in the ASV count table, containing the counts of each ASV in each sample ###

count_tab <- read.table(file = "counts_GLAmpSeq.tsv", 
                        header = TRUE, row.names = 1, sep = "\t")
```

**Input Data:**
* \*runsheet.csv (runsheet containing sample metadata required for processing, output from [step 6a](#6a-create-sample-runsheet))  
* counts_GLAmpSeq.tsv (ASV counts table, output from [step 5g](#5g-generating-and-writing-standard-outputs))
  
**Output Data:**
* `count_tab` (variable containing the ASV counts table created from counts.tsv)
* `runsheet` (variable containing sample metadata required for processing)
* `sample_info_tab` (variable containing a subtable of the runsheet, including the 'groups' column and an additional 'color' column with a color for each unique group)
  
<br>

___

## 7. Beta Diversity

Beta diversity measures the variation in species composition between different samples or environments. A common practice in working with a new dataset is to generate some exploratory visualizations like ordinations and hierarchical clusterings. These give us a quick overview of how our samples relate to each other and can be a way to check for problems like batch effects.

Create a DESeq2 object from the counts and the runsheet and apply the Variance Stabilization Transformation (VST):

```R
deseq_counts <- DESeqDataSetFromMatrix(countData = count_tab, 
                                       colData = runsheet, 
                                       design = ~1)
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
vst_counts <- assay(deseq_counts_vst)
```

**Input Data:**
* `count_tab` (variable containing the ASV counts table, output from [step 6b](#6b-environment-set-up))
* `runsheet` (variable containing sample metadata required for processing, output from [step 6b](#6b-environment-set-up))
  
**Output Data:**
* `vst_counts` (variable holding the VST-normalized ASV counts)
  
<br>

### 7a. Hierarchical Clustering

Create a euclidean distance matrix and perform hierarchical clustering.

```R
euc_dist <- dist(t(vst_counts))
euc_clust <- hclust(d = euc_dist, method = "ward.D2")
```

**Parameter Definitions:**  

*	`euc_clust <-` – specifying the variable that will hold the euclidean distance object

*	`hclust()` – the hclust function we are using for hierarchical clustering

*	`d=` – specifying the the input dissimilarity or distance object

* `method=` - specifying the method of clustering to use. "Ward.D2" is one that is commonly used

Create a dendrogram object.

```R
euc_dend <- as.dendrogram(euc_clust, hang = 0.1)
```

**Parameter Definitions:**  

*	`euc_dend <-` – specifying the variable that will hold the dendrogram object

*	`as.dendrogram()` – the dendrogram function we are using to create a dendrogram

*	`euc_clust` – an object that can be converted into a dendrogram

* `hang=` - numeric indicating how the leaves hight should be computed from the height of their parents
  
Color the sample branches by group and plot the dendrogram.

```R
dend_cols <- sample_info_tab$color[order.dendrogram(euc_dend)]
labels_colors(euc_dend) <- dend_cols
png("dendrogram_by_group_GLAmpSeq.png")
euc_dend %>% set("labels_cex", max_cex) %>% plot(ylab = "VST Euc. dist.")
dev.off()
```

**Input Data:**

* `vst_counts` (variable holding the VST-normalized ASV counts, output from [step 7](#7-beta-diversity))
  
**Output Data:**

* **dendrogram_by_group_GLAmpSeq.png** (dendrogram of euclidean distance - based hierarchical clustering of the samples, colored by experimental groups) 
* `euc_dist` (variable containing the samplewise euclidean distance matrix based on VST-normalized counts)

<br>

### 7b. Ordination

Ordination techniques like PCoA help in reducing the dimensionality of the data, allowing us to visualize complex relationships between samples.

Create a physeq object with an OTU table using VST-transformed counts and sample info table. Ordinate the counts using PCoA and euclidean distance.

```R
vst_count_phy <- otu_table(object = vst_counts, taxa_are_rows = TRUE)
sample_info_tab_phy <- sample_data(sample_info_tab)
vst_physeq <- phyloseq(vst_count_phy, sample_info_tab_phy)

vst_pcoa <- ordinate(physeq = vst_physeq, method = "PCoA", distance = "euclidean")
```

**Parameter Definitions:**  

*	`vst_pcoa <-` – specifying the variable that will hold the ordination of the VST-normalized counts

*	`physeq=` – specifying the Phyloseq object that contains the variance stabilized counts and sample metadata

*	`method=` – specifying the ordination method to be used

* `distance=` – specifying the distance metric used for the ordination

**Input Data:**

* `vst_counts` (variable holding the VST-normalized ASV counts, output from [step 7](#7-beta-diversity))
* `sample_info_tab` (variable containing a subtable of the runsheet, including the 'groups' and 'color' columns, output from [step 6b](#6b-environment-set-up))

**Output Data:**

* `vst_physeq` (variable holding the Phyloseq object)
* `vst_pcoa` (variable holding the object containing the coordinates and eigenvalues resulting from the PCoA of the VST-normalized counts)
  
<br>

___

## 8. Alpha Diversity

Alpha diversity examines the variety and abundance of taxa within individual samples. Rarefaction curves are utilized to visually represent this diversity, plotting the number of unique sequences (ASVs) identified against the total number of sequences sampled, offering a perspective on the saturation and completeness of sampling. Metrics like Chao1 richness estimates and Shannon diversity indices are employed to quantify the richness (total number of unique sequences) and diversity (combination of richness and evenness) within these samples.

<br>

### 8a. Rarefaction Curves

```R
rare_curve <- rarecurve(x = t(count_tab), step = 100, col = sample_info_tab$color, 
          lwd = 2, ylab = "ASVs", label = FALSE)
png("rarefaction_curves_GLAmpSeq.png")
plot(rare_curve)
dev.off()
```

**Parameter Definitions:**  

*	`rare_curve <-` – specifies the variable that will store the results within in our R environment

*	`rarecurve()` – the rarefy function we are calling, with the following parameters set within it

*	`x=` - specifies the input data for the rarefaction curve, which should be the transposed counts

*   `step=` - specifies the step size for sample sizes in rarefaction curves

*   `col=` - assigns a color to each line on the rarefaction curve for visual differentiation of sample groups

*   `lwd=` - sets the line width for the curves in the plot

*   `ylab=` - defines the label for the y-axis of the plot

*   `label=` - indicates whether the lines in the plot should be kept

**Input Data:**

* `vst_counts` (variable holding the VST-normalized ASV counts, output from [step 7](#7-beta-diversity))
* `sample_info_tab` (variable containing a subtable of the runsheet, including the 'groups' and 'color' columns, output from [step 6b](#6b-environment-set-up))

**Output Data:**

* **rarefaction_curves_GLAmpSeq.png** (Rarefaction curves plot for all samples)

<br>

### 8b. Richness and Diversity Estimates

```R
count_tab_phy <- otu_table(count_tab, taxa_are_rows = TRUE)
tax_tab_phy <- tax_table(as.matrix(tax_tab))
ASV_physeq <- phyloseq(count_tab_phy, tax_tab_phy, sample_info_tab_phy)
richness_and_diversity_estimates_by_sample <- plot_richness(ASV_physeq, color = "groups", measures = c("Chao1", "Shannon"))
richness_and_diversity_estimates_by_group <- plot_richness(ASV_physeq, x = "groups", color = "groups", measures = c("Chao1", "Shannon"))

ggsave(paste0("richness_and_diversity_estimates_by_sample_GLAmpSeq", ".png"), plot = richness_and_diversity_estimates_by_sample)
ggsave(paste0("richness_and_diversity_estimates_by_group_GLAmpSeq", ".png"), plot = richness_and_diversity_estimates_by_group)
```

**Parameter Definitions:**  

*	`plot_richness()` – the phyloseq function we are calling, with the following parameters set within it

*	`x=` - an optional variable to map to the horizontal axis of the plot

* `color=` - specifies a variable for determining the coloring scheme of the plot

* `measures=` - determines which of the available alpha-diversity measures to include in the plot

**Input Data:**

* `count_tab` (variable containing the ASV counts table, output from [step 6b](#6b-environment-set-up))
* `sample_info_tab` (variable containing a subtable of the runsheet, including the 'groups' and 'color' columns, output from [step 6b](#6b-environment-set-up))  
* `tax_tab` (variable containing the taxonomy table, created in [step 5g](#5g-generating-and-writing-standard-outputs))

**Output Data:**

* **richness_and_diversity_estimates_by_sample_GLAmpSeq.png** (Richness and diversity estimates plot for all samples)
* **richness_and_diversity_estimates_by_group_GLAmpSeq.png** (Richness and diversity estimates plot for all groups)  
* `ASV_physeq` (variable contiaining the Phyloseq object created using an OTU table based on the vst_counts)
  
<br>

___

## 9. Taxonomic Summaries

Taxonomic summaries provide insights into the composition of microbial communities at various taxonomic levels.

```R
proportions_physeq <- transform_sample_counts(ASV_physeq, function(ASV) ASV / sum(ASV))

relative_phyla <- plot_bar(proportions_physeq, x = "groups", fill = "phylum")
relative_classes <- plot_bar(proportions_physeq, x = "groups", fill = "class")

samplewise_phyla <- plot_bar(proportions_physeq, fill = "phylum")
samplewise_classes <- plot_bar(proportions_physeq, fill = "class")

ggsave(filename = "relative_phyla_GLAmpSeq", ".png", plot = relative_phyla)
ggsave(filename = "relative_classes_GLAmpSeq", ".png", plot = relative_classes)
ggsave(filename = "samplewise_relative_phyla_GLAmpSeq", ".png", plot = samplewise_phyla)
ggsave(filename = "samplewise_relative_classes_GLAmpSeq", ".png", plot = samplewise_classes)
```

**Input Data:**

* `ASV_physeq` (variable contiaining the Phyloseq object, output from [step 8b](#8b-richness-and-diversity-estimates))

**Output Data:**

* **relative_phyla_GLAmpSeq.png** (taxonomic summaries plot based on phyla, for all samples)
* **relative_classes_GLAmpSeq.png** (taxonomic summaries plot based on class, for all samples)

* **samplewise_phyla_GLAmpSeq.png** (taxonomic summaries plot based on phyla, for all samples)
* **samplewise_classes_GLAmpSeq.png** (taxonomic summaries plot based on class, for all samples)

<br>

___

## 10. Differential Abundance Analysis

Using Betadisper, permutational ANOVA, and DESeq2, we aim to uncover specific taxa that exhibit notable variations across different conditions, complemented by visualizations like volcano plots to illustrate these disparities and their implications on ASV expression and overall microbial community dynamics.

<br>

### 10a. Betadisper and Permutational ANOVA

Use betadisper to check whether variability of data points in each group is similar.

```R
betadisper(d = euc_dist, group = sample_info_tab$groups) %>% anova()
```

**Parameter Definitions:**  

*	`betadisper()` – the vegan function we are calling, with the following parameters set within it

*	`d=` - specifies the input distance object

* `group=` - specifies the sample grouping information

* `%>% anova()` - Sends the output object from betadisper() to the anova() function to perform the permutational ANOVA test

Use adonis2 to test whether the mean of data differs significantly between groups.

```R
adonis_res <- adonis2(formula = euc_dist ~ sample_info_tab$groups)
```

**Parameter Definitions:**  

*	`adonis_res <-` – specifies the variable that will store the results within in our R environment

*	`adonis2()` – the vegan function we are calling, with the following parameters set within it

*	`formula=` - specifies the model formula or data matrix

Statistics from Adonis2 testing can be incorporated into PCoA visualizations using vst_pcoa which was made earlier.

```R
r2_value <- adonis_res$R2[1]
prf_value <- adonis_res$`Pr(>F)`[1]

label_PC1 <- sprintf("PC1 [%.1f%%]", percent_variance[1])
label_PC2 <- sprintf("PC2 [%.1f%%]", percent_variance[2])

ordination_plot <- plot_ordination(vst_physeq, vst_pcoa, color = "groups") + 
  geom_point(size = 1) + 
  labs(
    col = "Groups", 
    x = label_PC1,
    y = label_PC2
  ) + 
  geom_text(aes(label = rownames(sample_info_tab)), show.legend = FALSE, hjust = 0.3, vjust = -0.4, size = 4) + 
  coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + 
  scale_color_manual(values = unique(sample_info_tab$color),
                     labels = unique(sample_info_tab$groups)) +
  theme_bw() + theme(legend.position = "bottom",  text = element_text(size = 15, ),
                     legend.direction = "vertical",
                     legend.justification = "center",
                     legend.box.just = "center",
                     legend.title.align = 0.5) +
  annotate("text", x = Inf, y = -Inf, label = paste("R2:", toString(round(r2_value, 3))), hjust = 1.1, vjust = -2, size = 4)+
  annotate("text", x = Inf, y = -Inf, label = paste("Pr(>F)", toString(round(prf_value,4))), hjust = 1.1, vjust = -0.5, size = 4)+ ggtitle("PCoA")

ordination_plot_u <- plot_ordination(vst_physeq, vst_pcoa, color = "groups") + 
  geom_point(size = 1) + 
  labs( 
    x = label_PC1,
    y = label_PC2,
    col = "Groups"
  ) +
  coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + 
  scale_color_manual(values = unique(sample_info_tab[[color_colname]]),
                     labels = unique(sample_info_tab$short_group_labels)) +
  theme_bw() + theme(legend.position = "bottom",  text = element_text(size = 15, ),
                     legend.direction = "vertical",
                     legend.justification = "center",
                     legend.box.just = "center",
                     legend.title.align = 0.5) +
  annotate("text", x = Inf, y = -Inf, label = paste("R2:", toString(round(r2_value, 3))), hjust = 1.1, vjust = -2, size = 4)+
  annotate("text", x = Inf, y = -Inf, label = paste("Pr(>F)", toString(round(prf_value,4))), hjust = 1.1, vjust = -0.5, size = 4)+ ggtitle("PCoA")

ggsave(filename=paste0(beta_diversity_out_dir, output_prefix, "PCoA_w_labels_GLAmpSeq", ".png"), plot=ordination_plot)
ggsave(filename=paste0(beta_diversity_out_dir, output_prefix, "PCoA_without_labels_GLAmpSeq", ".png"), plot=ordination_plot_u)

```

**Input Data:**

* `euc_dist` (variable containing the samplewise euclidean distance matrix of transposed VST-normalized counts, output from [step 7a](#7a-hierarchical-clustering))
* `sample_info_tab` (variable containing a subtable of the runsheet, including the 'groups' and 'color' columns, output from [step 6b](#6b-environment-set-up))  
* `vst_physeq` (variable holding the Phyloseq object, output from [step 7b](#7b-ordination))

**Output Data:**

*  **PCoA_w_labels_GLAmpSeq.png** (principle Coordinates Analysis plot of VST transformed ASV counts, with sample labels)
*  **PCoA_without_labels_GLAmpSeq.png** (principle Coordinates Analysis plot of VST transformed ASV counts, without sample labels)

<br>

### 10b. Differential abundance analysis with DESeq2

DESeq2 can be used to identify specific ASVs that have significantly different copy-number counts between sample groups.

```R
deseq_obj <- phyloseq_to_deseq2(physeq = ASV_physeq, design = ~groups)
```

**Parameter Definitions:**  

*	`deseq_obj <-` – specifies the variable that will store the results within in our R environment

*	`phyloseq_to_deseq2()` – the phyloseq function we are calling, with the following parameters set within it

*	`physeq=` - specifies the phyloseq-class object

*   `design=` - a formula specifying the design of the experiment

Run the DESeq() function to normalize for sample read-depth and composition, transform the data, and test for differential abundance between the groups. Save the size-factor-normalized counts.


```R
deseq_modeled <- DESeq(deseq_obj)

write.table(counts(deseq_modeled, normalized=TRUE), file = paste0("normalized_counts_GLAmpSeq.tsv"), sep="\t", row.names=TRUE, quote=FALSE)
```

**Input Data:**

* `ASV_physeq` (variable contiaining the Phyloseq object, output from [step 8b](#8b-richness-and-diversity-estimates))
  
**Output Data:**
* **normalized_counts_GLAmpSeq.tsv** (size factor normalized ASV counts table)
* `deseq_modeled` (variable holding the DESeq2 output data)

<br>

### 10c. Volcano Plots

Define the function for creating the volcano plot and saving the normalized counts for a given contrast.

```R
plot_comparison <- function(group1, group2) {
  
  deseq_res <- results(deseq_modeled, contrast = c("groups", group1, group2))
  norm_tab <- counts(deseq_modeled, normalized = TRUE) %>% data.frame()
  
  volcano_data <- as.data.frame(deseq_res)
  
  p_val <- 0.1
  volcano_data <- volcano_data[!is.na(volcano_data$padj), ]
  volcano_data$significant <- volcano_data$padj <= p_val
  
  p <- ggplot(volcano_data, aes(x=log2FoldChange, y=-log10(padj), color=significant)) +
    geom_point(alpha=0.7, size=2) +
    scale_color_manual(values=c("black", "red"), labels=c(paste0("padj > ", p_val), paste0("padj \u2264 ", p_val))) +
    theme_bw() +
    labs(title="Volcano Plot",
         x=paste("Log2 Fold Change\n(",group1," vs ",group2,")"),
         y="-Log10 P-value",
         color=paste0("")) +
    theme(legend.position="top")
  
  top_points <- volcano_data %>%
    arrange(padj) %>%
    filter(significant) %>%
    head(10)
  
  volcano_plot <- p + geom_text_repel(data=top_points, aes(label=row.names(top_points)), size=3)
  ggsave(filename=paste0("volcano_",
                         gsub(" ", "_", group1),
                         "_vs_",
                         gsub(" ", "_", group2), ".png"),
         plot=volcano_plot,
         width = 11.1, height = 8.33, dpi = 300)
  
  write.csv(deseq_res, file = paste0(gsub(" ", "_", group1),
                        "_vs_",
                        gsub(" ", "_", group2), ".csv"))
}
```

Create volcano plots for all pairwise comparisons.

```R
comparisons <- expand.grid(group1 = unique_groups, group2 = unique_groups)
comparisons <- subset(comparisons, group1 != group2)

apply(comparisons, 1, function(pair) plot_comparison(pair['group1'], pair['group2']))
```
**Parameter Definitions:**  


*	`apply()` – the function we are calling, with the following parameters set within it

*	`physeq=` - the data matrix or array on which the function is to be applied

* `1` – indicates that the function should be applied to each row

* `function(pair) plot_comparison(pair['group1'], pair['group2'])` – an anonymous function which takes a pair of values as input and executes the plot_comparison function on these values

**Input Data:**

* `deseq_modeled` (variable holding the DESeq2 output data, output from [step 10b](#910-differential-abundance-analysis-with-deseq2))
  
**Output Data:**

* **group1_vs_group2.csv** (differential abundance tables for all pairwise contrasts of groups)
* **volcano_group1_vs_group2.png** (volcano plots for all pairwise contrasts of groups)
