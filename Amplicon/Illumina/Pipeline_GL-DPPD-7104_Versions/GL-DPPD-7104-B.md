# Bioinformatics pipeline for amplicon Illumina sequencing data  

> **This page holds an overview and instructions for how GeneLab processes Illumina amplicon datasets. Exact processing commands for specific datasets that have been released are available in the [GLDS_Processing_Scripts](../GLDS_Processing_Scripts) sub-directory and/or are provided with their processed data in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/).**  

---

**Date:** December 30, 2024  
**Revision:** B  
**Document Number:** GL-DPPD-7104-B  

**Submitted by:**  
Olabiyi Obayomi, Alexis Torres, and Michael D. Lee (GeneLab Data Processing Team)

**Approved by:**  
Samrawit Gebre  (GeneLab Project Manager)  
Danielle Lopez (GeneLab Deputy Project Manager)  
Lauren Sanders (OSDR Project Scientist)  
Amanda Saravia-Butler (GeneLab Data Processing Lead)  

---

## Updates from previous version


- Additional software (R packages) used:
  - ANCOMBC
  - broom
  - DESeq2
  - DescTools
  - FSA
  - ggdendro
  - glue
  - ggrepel
  - mia
  - multcompView
  - optparse
  - patchwork
  - phyloseq
  - RColorBrewer
  - rstatix
  - taxize
  - tidyverse
  - tools
  - utils  
  - vegan



- Inclusion of additional steps and outputs starting from ([step 6](#6-amplicon-seq-data-analysis-set-up)):
  - Alpha Diversity Analysis ([step 7](#7-alpha-diversity-analysis)).
  - Beta Diversity Analysis ([step 8](#8-beta-diversity-analysis)).
  - Groupwise and Samplewise Taxonomic Summary Plots ([step 9](#9-taxonomy-plots)).
  - Differential Abundance Testing ([step 10](#9-differential-abundance-analysis)) with ANCOMBC 1 ([10a](#10a-ancombc-1)), ANCOMBC 2 ([10b](#10b-ancombc-2)) and Deseq2 ([10c](#10c-deseq2)).

<!-- Included R packages -->
- Assay-specific suffixes were added where needed for GeneLab repo ("_GLAmpSeq")
- The ITS UNITE reference database used was updated to "UNITE_v2023_July2023.RData", from http://www2.decipher.codes/Classification/TrainingSets/
- Persistent Reference links to RDATA databases on Figshare replaced reference links on DECIPHER's [website](http://www2.decipher.codes/Classification/TrainingSets/) for [SILVA SSU r138](https://figshare.com/ndownloader/files/46245217), [UNITE v2023](https://figshare.com/ndownloader/files/49181545) and [PR2 v4.13](https://figshare.com/ndownloader/files/46241917)
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
    - [6b. R Environment Set Up](#6b-r-environment-set-up)
      - [Load Libraries](#load-libraries)
      - [Load Functions](#load-functions)
      - [Set Variables](#set-variables)
      - [Read-in Input Tables](#read-in-input-tables)
      - [Preprocessing](#preprocessing)
  - [7. Alpha Diversity Analysis](#7-alpha-diversity-analysis)
  - [8. Beta Diversity Analysis](#8-beta-diversity-analysis)
  - [9. Taxonomy Plots](#9-taxonomy-plots)
  - [10. Differential Abundance Testing](#10-differential-abundance-testing)
    - [10a. ANCOMBC 1](#10a-ancombc-1)
    - [10b. ANCOMBC 2](#10b-ancombc-2)
    - [10c. DESeq2 ](#10c-deseq2)

---

# Software used  

|Program|Version*|Relevant Links|
|:------|:------:|:-------------|
|FastQC|`0.12.1`|[https://www.bioinformatics.babraham.ac.uk/projects/fastqc/](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)|
|MultiQC|`1.19`|[https://multiqc.info/](https://multiqc.info/)|
|Cutadapt|`4.6`|[https://cutadapt.readthedocs.io/en/stable/](https://cutadapt.readthedocs.io/en/stable/)|
|R-base|`4.4.1`|[https://www.r-project.org/](https://www.r-project.org/)|
|DADA2|`1.30.0`|[https://www.bioconductor.org/packages/release/bioc/html/dada2.html](https://www.bioconductor.org/packages/release/bioc/html/dada2.html)|
|DECIPHER|`2.30.0`|[https://bioconductor.org/packages/release/bioc/html/DECIPHER.html](https://bioconductor.org/packages/release/bioc/html/DECIPHER.html)|
|biomformat|`1.30.0`|[https://github.com/joey711/biomformat](https://github.com/joey711/biomformat)|
|ANCOMBC|`2.6.0`|[https://github.com/FrederickHuangLin/ANCOMBC](https://github.com/FrederickHuangLin/ANCOMBC)|
|broom|`1.0.7`|[https://CRAN.R-project.org/package=broom](https://CRAN.R-project.org/package=broom)|
|DescTools|`0.99.57`|[https://andrisignorell.github.io/DescTools/](https://andrisignorell.github.io/DescTools/)|
|DESeq2|`1.42.0`|[https://bioconductor.org/packages/release/bioc/html/DESeq2.html](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)|
|FSA|`0.9.5`|[https://CRAN.R-project.org/package=FSA](https://CRAN.R-project.org/package=FSA)|
|ggdendro|`0.2.0`|[https://CRAN.R-project.org/package=ggdendro](https://CRAN.R-project.org/package=ggdendro)|
|ggrepel|`0.9.6`|[https://CRAN.R-project.org/package=ggrepel](https://CRAN.R-project.org/package=ggrepel)|
|glue|`1.8.0`|[https://glue.tidyverse.org/](https://glue.tidyverse.org/)|
|mia|`1.12.0`|[https://github.com/microbiome/mia](https://github.com/microbiome/mia)|
|phyloseq|`1.46.0`|[https://bioconductor.org/packages/release/bioc/html/phyloseq.html](https://bioconductor.org/packages/release/bioc/html/phyloseq.html)|
|rcolorbrewer|`1.1_3`|[https://CRAN.R-project.org/package=RColorBrewer](https://CRAN.R-project.org/package=RColorBrewer)|
|taxize|`0.9.100.1`|[https://docs.ropensci.org/taxize/](https://docs.ropensci.org/taxize/)|
|tidyverse|`2.0.0`|[https://CRAN.R-project.org/package=tidyverse](https://CRAN.R-project.org/package=tidyverse)|
|tools|`4.4.1`|[https://www.R-project.org/](https://www.R-project.org/)|
|utils|`4.4.1`|[https://www.R-project.org/](https://www.R-project.org/)|
|vegan|`2.6.4`|[https://cran.r-project.org/package=vegan](https://cran.r-project.org/package=vegan)|

# Reference databases used

|Program used| Database| Relevant Links|
|:-----|:-----:|--------:|
|DECIPHER| SILVA SSU r138 | [SILVA_SSU_r138_2019.RData](https://figshare.com/ndownloader/files/46245217)|
|DECIPHER| UNITE v2023 | [UNITE_v2023_July2023.RData](https://figshare.com/ndownloader/files/49181545)|
|DECIPHER| PR2 v4.13 | [PR2_v4_13_March2021.RData](https://figshare.com/ndownloader/files/46241917)|
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
  
<br>

### 6a. Create Sample Runsheet

> Note: Rather than running the command below to create the runsheet needed for processing, the runsheet may also be created manually by following the examples for [Paired-end](../Workflow_Documentation/NF_AmpIllumina-B/workflow_code/PE_file.csv) and [Single-end](../Workflow_Documentation/NF_AmpIllumina-B/workflow_code/SE_file.csv) samples.

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
library(tidyverse)
```

#### Load Functions

```R
# Function to calculate text size for plotting
calculate_text_size <- function(num_samples, start_samples = 25, min_size = 3) {
  max_size = 11  # Maximum size for up to start_samples
  slope = -0.15
  
  if (num_samples <= start_samples) {
    return(max_size)
  } else {
    # Calculate the current size with the hard coded slope
    current_size = max_size + slope * (num_samples - start_samples)
    
    # Ensure the size doesn't go below the minimum
    return(max(current_size, min_size))
  }
}

# A function to create a phyloseq object with the appropriate
# sample count transformation depending on the supplied transformation method
# i.e. either 'rarefy' or  'vst'
transform_phyloseq <- function( feature_table, metadata, method, rarefaction_depth=500){
  # feature_table  [DATAFRAME] ~ Feature / ASV count table with samples as columns and features as rows 
  # metadata [DATAFRAME] ~  Samples metadata with samples as row names
  # method [STRING] ~ Distance transformation method to use.
  #                   Either 'rarefy' or 'vst' for rarefaction and variance 
  #                   stabilizing transformation, respectively.
  # rarefaction_depth [INT] ~ Sample rarefaction to even depth when method is 'bray'
  
  if(method == 'rarefy'){
    # Create phyloseq object
    ASV_physeq <- phyloseq(otu_table(feature_table, taxa_are_rows = TRUE),
                           sample_data(metadata))
    
    
    seq_per_sample <- colSums(feature_table) %>% sort()
    # Minimum value
    depth <- min(seq_per_sample)
    
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

# -----------    Hierarchical Clustering and dendogram plotting
make_dendogram <- function(dist_obj, metadata, groups_colname,
                           group_colors, legend_title){
  
  
  sample_clust <- hclust(d = dist_obj, method = "ward.D2")
  
  # Extract clustering data
  hcdata <- dendro_data(sample_clust, type = "rectangle")
  segment_data <- segment(hcdata)
  label_data <- label(hcdata) %>%
    left_join(metadata %>% 
                rownames_to_column("label"))

  # Plot dendogram
  dendogram <- ggplot() +
    geom_segment(data = segment_data, 
                 aes(x = x, y = y, xend = xend, yend = yend)
    ) +
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

# Run variance test and adonis test
run_stats <- function(dist_obj, metadata, groups_colname){
  
  samples <- attr(dist_obj, "Label")
  metadata <- metadata[samples,]
  variance_test <- betadisper(d = dist_obj, 
                              group = metadata[[groups_colname]]) %>%
    anova() %>%
    broom::tidy() %>% 
    mutate(across(where(is.numeric), ~round(.x, digits = 2)))
  
  
  adonis_res <- adonis2(formula = dist_obj ~ metadata[[groups_colname]])
  adonis_test <- adonis_res %>%
    broom::tidy() %>% 
    mutate(across(where(is.numeric), ~round(.x, digits = 2)))
  
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
  
  vectors_df <- pcoa$vectors %>%
                   as.data.frame() %>%
                   rownames_to_column("samples")
  
  plot_df <- sample_data(ps) %>%
               as.matrix() %>%
               as.data.frame() %>%
               rownames_to_column("samples") %>% 
               select(samples, !!groups_colname) %>% 
               right_join(vectors_df, join_by("samples"))
  
  p <- ggplot(plot_df, aes(x=Axis.1, y=Axis.2, 
                           color=!!sym(groups_colname), 
                           label=samples)) +
    geom_point(size=1)

  
  if(addtext){
    p <- p + geom_text(show.legend = FALSE,
                       hjust = 0.3, vjust = -0.4, size = 4)
  }
  
  
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
  
  # feature_table [MATRIX]  feature table matrix with samples as columns and 
  #                         features as rows
  # cut_off_percent [NUMERIC] cut-off fraction  or decimal between 0.001 to 1 
  #                          of the total number of samples to determine the 
  #                          most abundant features. By default it removes 
  #                          features that are not present in 3/4 of the total 
  #                          number of samples
  
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

process_taxonomy <- function(taxonomy, prefix='\\w__') {
  # Function to process a taxonopmy assignment table
  #1. ~ taxonomy is a string specifying the taxonomic assignment file name
  #2 prefix ~ is a regular expression specifying the characters to remove
  # from the taxon names  '\\w__'  for greengenes and 'D_\\d__' for SILVA
  
  
  taxonomy <- apply(X = taxonomy, MARGIN = 2, FUN = as.character) 
  
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

# Function to format a taxonomy assignment table by appending a suffix
# to a known name
format_taxonomy_table <- function(taxonomy=taxonomy.m,stringToReplace="Other",
                                  suffix=";Other") {
  
  for (taxa_index in seq_along(taxonomy)) {
    
    indices <- grep(x = taxonomy[,taxa_index], pattern = stringToReplace)
    
    taxonomy[indices,taxa_index] <- 
      paste0(taxonomy[indices,taxa_index-1],
             rep(x = suffix, times=length(indices)))
    
  }
  return(taxonomy)
}

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

# A function to generate taxon level count matrix based on a taxonomy table and
# an existing feature table
make_feature_table <- function(count_matrix,taxonomy,
                               taxon_level, samples2keep=NULL){

  # EAMPLE:
  # make_feature_table(count_matrix = feature_counts_matrix,
  #                    taxonomy = taxonomy_table, taxon_level = "Phylum")

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
  # abund_table is a relative abundance matrix with taxa as columns and  samples as rows
  #rare_taxa is a boolean specifying if only rare taxa should be returned
  #If set to TRU then a table with only the rare taxa will be returned 
  #intialize an empty vector that will contain the indices for the
  #low abundance columns/ taxa to group
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


# Function to collapse the samples in an oTU table with a defined function(fun)
# based on a group in metadata 
collapse_samples <- function(taxon_table,metadata,group,fun=sum, 
                             convertToRelativeAbundance=FALSE){
  # function to collapse the samples in an oTU table with a defined function(fun)  based on a group in metadata 
  # taxon_table - a matrix count table with samples as rows and features/OTUs as columns
  # metadata - a dataframe to containing the group to collapse samples by. Sample names must be the rownames of the metadata
  # group - an independent factor variable within the metadata to collapse the samples by
  # fun - a function without brackets to apply in order to collapse the samples
  # convertToRelativeAbundance - a boolean set to TRUE OR FALSE if the taxon_table shout be converted to relative abundance
  # default is FALSE
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

# A function to run ANCOMBC2 while handlixnxg commxon 
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
      
      # Second error catcher in case it fails in first one 
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

#### Read-in Input Tables

```R
# Read-in metadata
metadata <- read_csv(file = metadata_file) %>% as.data.frame()
row.names(metadata) <- metadata[[sample_colname]]
metadata[,sample_colname] <- NULL
group_column_values <- metadata %>% pull(!!sym(groups_colname))
group_levels <- unique(group_column_values)

# Add colors to metadata equals to the number of levels
# in the factor groups column
num_colors <- length(group_levels)
palette <- 'Set1'
number_of_colors_in_palette <- 9
if(num_colors <= number_of_colors_in_palette){
   colors <- RColorBrewer::brewer.pal(n = num_colors, name = palette)
}else{
  colors <- custom_palette[1:num_colors]
}

# Metadata
group_colors <- setNames(colors, group_levels)
metadata <- metadata %>%
  mutate(color = map_chr(!!sym(groups_colname),
                         function(group) { group_colors[group] }
                         ) 
        )
sample_names <- rownames(metadata)
deseq2_sample_names <- make.names(sample_names, unique = TRUE)

sample_info_tab <- metadata %>%
  select(!!groups_colname, color) %>%
  arrange(!!sym(groups_colname))


values <- sample_info_tab %>% pull(color) %>% unique()



# Feature or ASV table
feature_table <- read.table(file = features_file, header = TRUE,
                            row.names = 1, sep = "\t")

# Taxonomy table
taxonomy_table <-  read.table(file = taxonomy_file, header = TRUE,
                              row.names = 1, sep = "\t")
```


#### Preprocessing

```R
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
               filter(str_detect(taxonomy, "[Cc]hloroplast|[Mn]itochondria")) %>%
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
# only features found in both table
feature_table <- feature_table[common_ids,]
taxonomy_table <- taxonomy_table[common_ids,]
```


## 7. Alpha Diversity Analysis

Alpha diversity examines the variety and abundance of taxa within individual samples. Rarefaction curves are utilized to visually represent this diversity, plotting the number of unique sequences (ASVs) identified against the total number of sequences sampled, offering a perspective on the saturation and completeness of sampling. Metrics like Chao1 richness estimates and Shannon diversity indices are employed to quantify the richness (total number of unique sequences) and diversity (combination of richness and evenness) within these samples.

```bash
Rscript alpha_diversity.R \
                  --metadata-table amplicon_runsheet.csv \
                  --feature-table counts_GLAmpSeq.tsv \
                  --taxonomy-table taxonomy_GLAmpSeq.tsv \
                  --group groups \
                  --samples-column 'Sample Name' \
                  --rarefaction-depth 500
```
**Parameter Definitions:**

*	`--metadata-table` – specifies the path to a comma separated samples metadata file with the group/treatment to be analyzed
*	`--feature-table` – specifies the path to a tab separated samples feature table i.e. ASV or OTU table
*	`--taxonomy-table` – specifies the path to a feature taxonomy table i.e. ASV taxonomy table
*	`--group` – specifies the group column in metadata to be analyzed
* `--samples-column` – specifies the column in metadata containing the sample names in the feature table
* `--rarefaction-depth`  – specifies the minimum rarefaction depth for alpha diversity estimation

Type `Rscript alpha_diversity.R --help` at the commandline for a full list of available parameters

Content of `alpha_diversity.R`  
```R
# Create output directory if it doesn't already exist
alpha_diversity_out_dir <- "alpha_diversity/"
if(!dir.exists(alpha_diversity_out_dir)) dir.create(alpha_diversity_out_dir)

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
  
  add_significance(df, p.col='p', output.col = 'p.signif') %>% 
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
  p_values <- sub_comp$p
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
                select(groups, label= !!sym( glue("{metric}_letter") ) 
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
                                    mapping = aes(y=max+toAdd, 
                                                  label=label,
                                                  fontface = "bold"),
                                  size = text_size)
})

richness_by_group <- wrap_plots(p, ncol = 2, guides =  'collect')

# Save group plot
width <- 3.6 * length(group_levels)
ggsave(filename = glue("{alpha_diversity_out_dir}/{output_prefix}richness_and_diversity_estimates_by_group{assay_suffix}.png"),
       plot=richness_by_group, width = width, height = 8.33, dpi = 300, units = "in")
```

**Input Data:**

* **amplicon_runsheet.csv** (metadata table - e.g {OSD-Accession-ID}_AmpSeq_v{version}_runsheet.csv )
* **counts_GLAmpSeq.tsv** (count table)
* **taxonomy_GLAmpSeq.tsv** (taxonomy table)

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

```bash
Rscript beta_diversity.R \
                  --metadata-table amplicon_runsheet.csv \
                  --feature-table counts_GLAmpSeq.tsv \
                  --taxonomy-table taxonomy_GLAmpSeq.tsv \
                  --group groups \
                  --samples-column 'Sample Name' \
                  --rarefaction-depth 500
```
**Parameter Definitions:**

*	`--metadata-table` – specifies the path to a comma separated samples metadata file with the group/treatment to be analyzed
*	`--feature-table` – specifies the path to a tab separated samples feature table i.e. ASV or OTU table
*	`--taxonomy-table` – specifies the path to a feature taxonomy table i.e. ASV taxonomy table
*	`--group` – specifies the group column in metadata to be analyzed
* `--samples-column` – specifies the column in metadata containing the sample names in the feature table
* `--rarefaction-depth`  – specifies the minimum rarefaction depth for diversity estimation. Relavant only for Bray Curtis distance calculation between samples 

type `Rscript beta_diversity.R --help` at the commandline for a full list of available parameters

Content of `beta_diversity.R`  

```R

beta_diversity_out_dir <- "beta_diversity/"
if(!dir.exists(beta_diversity_out_dir)) dir.create(beta_diversity_out_dir)

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
ggsave(filename = glue("{beta_diversity_out_dir}/{output_prefix}{distance_method}_dendrogram{assay_suffix}.png"),
       plot = dendogram, width = 14, height = 10, dpi = 300, units = "in")

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

**Input Data:**

* **amplicon_runsheet.csv** (metadata table)
* **counts_GLAmpSeq.tsv** (count table)
* **taxonomy_GLAmpSeq.tsv** (taxonomy table)

**Output Data:**

* **beta_diversity/<output_prefix><distance_method>_dendrogram_GLAmpSeq.png** (Dendogram)
* **beta_diversity/<output_prefix><distance_method>_adonis_table_GLAmpSeq.csv** (Adonis Stats Table)
* **beta_diversity/<output_prefix><distance_method>_PCoA_without_labels_GLAmpSeq.png** (Unlabeled PCoA)
* **beta_diversity/<output_prefix><distance_method>_PCoA_w_labels_GLAmpSeq.png** (Labeled PCoA)

where distance_method is either bray or euclidean for Bray Curtis and Euclidean distance, respectively.

<br>

---


## 9. Taxonomy Plots


```bash
Rscript plot_taxonomy.R \
                  --metadata-table mapping/GLDS-487_amplicon_v1_runsheet.csv \
                  --feature-table data/counts_GLAmpSeq.tsv \
                  --taxonomy-table data/taxonomy_GLAmpSeq.tsv \
                  --group groups \
                  --samples-column 'Sample Name' \
                  --remove-rare FALSE \
                  --prevalence-cutoff 0.15 \
                  --library-cutoff 100

```
**Parameter Definitions:**

*	`--metadata-table` – specifies the path to a comma separated samples metadata file with the group/treatment to be analyzed
*	`--feature-table` – specifies the path to a tab separated samples feature table i.e. ASV or OTU table
*	`--taxonomy-table` – specifies the path to a feature taxonomy table i.e. ASV taxonomy table
*	`--group` – specifies the group column in metadata to be analyzed
* `--samples-column` – specifies the column in metadata containing the sample names in the feature table
* `--remove-rare` - should rare features be filtered out prior to analysis? If set, rare feature will be removed
* `--prevalence-cutoff` - If --remove-rare, a numerical fraction between 0 and 1. 
                         Taxa with prevalences(the proportion of samples in which the taxon is present) less than --prevalence-cutoff will be excluded in the analysis. Default is 0.15, i.e. exclude taxa / features that are not present in at least 15% of the samples.
* `--library-cutoff` - If --remove-rare, a numerical threshold for filtering samples based on library sizes. 
                       Samples with library sizes less than lib_cut will be excluded in the analysis. Default is 100. 
                       if you do not want to discard any sample then set to 0.
type `Rscript plot_taxonomy.R --help` at the commandline for a full list of available parameters

Content of `plot_taxonomy.R`  

```R
taxonomy_plots_out_dir <- "taxonomy_plots/"
if(!dir.exists(taxonomy_plots_out_dir)) dir.create(taxonomy_plots_out_dir)

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


x_lab <- "Samples"
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
                               facet_wrap(facet_by, scales = "free", nrow = 1) +
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
thresholds <- c(phylum=1,class=2, order=2, family=2, genus=2, species=2)

# Convert from wide to long format for every treatment group of interest
group_rare <- TRUE
maximum_number_of_taxa <- 500

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
plot_width <- 2 * number_of_groups
walk2(.x = group_relAbundace_tbs, .y = taxon_levels[-1], 
                           .f = function(relAbundace_tb=.x, taxon_level=.y){
                             
                             p <- ggplot(data =  relAbundace_tb, mapping = aes(x= !!sym(groups_colname), y = !!sym(y)   )) +
                               geom_col(aes(fill = !!sym(taxon_level))) + 
                               publication_format +
                               theme(axis.text.x=element_text(
                                 margin=margin(t=0.5,r=0,b=0,l=0,unit ="cm"),
                                 angle = 90, 
                                 hjust = 0.5, vjust = 0.5)) + 
                               labs(x = NULL , y = y_lab, fill = tools::toTitleCase(taxon_level)) + 
                               scale_fill_manual(values = custom_palette)
                             ggsave(filename = glue("{taxonomy_plots_out_dir}/{output_prefix}groups_{taxon_level}{assay_suffix}.png"),
                                    plot=p, width = plot_width, height = 10, dpi = 300)
                           })
```

**Input Data:**

* **amplicon_metdata.csv** (metadata table)
* **counts_GLAmpSeq.tsv** (count table)
* **taxonomy_GLAmpSeq.tsv** (taxonomy table)

**Output Data:**

* **taxonomy_plots/<output_prefix>samples_<taxon_level>_GLAmpSeq.png** (samples barplots)
* **taxonomy_plots/<output_prefix>groups_<taxon_level>_GLAmpSeq.png** (groups barplots)

where taxon_level is one of phylum, class, order, family, genus and species.

> please note that species plot should only be taken with a grain of salt as short amplicon sequences can't be used to accurately predict species.

<br>

---


## 10. Differential Abundance Testing


### 10a. ANCOMBC 1
```bash
Rscript pairwise_ancombc1.R \
                  --metadata-table amplicon_runsheet.csv \
                  --feature-table counts_GLAmpSeq.tsv \
                  --taxonomy-table taxonomy_GLAmpSeq.tsv \
                  --group groups \
                  --samples-column 'Sample Name' \
                  --target-region 16S \
                  --remove-rare FALSE \
                  --prevalence-cutoff 0.15 \
                  --library-cutoff 100 \
                  --cpus 5

```
**Parameter Definitions:**

*	`--metadata-table` – specifies the path to a comma separated samples metadata file with the group/treatment to be analyzed
*	`--feature-table` – specifies the path to a tab separated samples feature table i.e. ASV or OTU table
*	`--taxonomy-table` – specifies the path to a feature taxonomy table i.e. ASV taxonomy table
*	`--group` – specifies the group column in metadata to be analyzed
* `--samples-column` – specifies the column in metadata containing the sample names in the feature table
* `--target-region`  – specifies the amplicon target region. Options are either 16S, 18S or ITS
* `--remove-rare` - should rare features be filtered out prior to analysis? If set, rare feature will be removed
* `--prevalence-cutoff` - If --remove-rare, a numerical fraction between 0 and 1. 
                         Taxa with prevalences(the proportion of samples in which the taxon is present) less than --prevalence-cutoff will be excluded in the analysis. Default is 0.15, i.e. exclude taxa / features that are not present in at least 15% of the samples.
* `--library-cutoff` - If --remove-rare, a numerical threshold for filtering samples based on library sizes. 
                       Samples with library sizes less than lib_cut will be excluded in the analysis. Default is 100. 
                       if you do not want to discard any sample then set to 0.
* `--cpus ` - Specifies the number of cpus to use for parallel processing.

Type `Rscript pairwise_ancombc1.R --help` at the commandline for a full list of available parameters

Content of `pairwise_ancombc1.R`  

```R
# Create output directory if it doesn't already exist
diff_abund_out_dir <- "differential_abundance/ancombc1/"
if(!dir.exists(diff_abund_out_dir)) dir.create(diff_abund_out_dir, recursive = TRUE)


# Create phyloseq object
ps <- phyloseq(otu_table(feature_table, taxa_are_rows = TRUE),
               sample_data(metadata),
               tax_table(as.matrix(taxonomy_table)))

# Convert phyloseq to tree summarized experiment object
tse <-  mia::makeTreeSummarizedExperimentFromPhyloseq(ps)


# Get unique group comparison as a matrix
pairwise_comp.m <- utils::combn(metadata[,group] %>% unique, 2)
pairwise_comp_df <- pairwise_comp.m %>% as.data.frame 

colnames(pairwise_comp_df) <- map_chr(pairwise_comp_df,
                                      \(col) str_c(col, collapse = "v"))
comparisons <- colnames(pairwise_comp_df)
names(comparisons) <- comparisons


message("Running ANCOMBC1....")
set.seed(123)
final_results_bc1  <- map(pairwise_comp_df, function(col){
  
  group1 <- col[1]
  group2 <- col[2]
  
  tse_sub <-  tse[, tse[[group]] %in% c(group1, group2)]
  
  # Note that by default, levels of a categorical variable in R are sorted 
  # alphabetically. 
  # Changing the reference group by reordering the factor levels
  tse_sub[[group]] <- factor(tse_sub[[group]] , levels = c(group1, group2))
  
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
                  formula = group, 
                  p_adj_method = "fdr", prv_cut = prevalence_cutoff,
                  lib_cut = library_cutoff, 
                  group = group, struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
                  max_iter = 100, conserve = TRUE, alpha = 0.05, global = FALSE,
                  n_cl = threads, verbose = TRUE)
  
  # ------ Set data frame names ---------# 
  # LFC 
  lfc <- out$res$lfc %>%
    as.data.frame() %>% 
    select(-contains("Intercept")) %>% 
    set_names(
      c("taxon",
        glue("logFC_({group2})v({group1})"))
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
  
  
  res <-lfc %>%
    left_join(se) %>%
    left_join(W) %>% 
    left_join(p_val)  %>% 
    left_join(q_val) %>% 
    left_join(diff_abn)
  
  
  return(res)
  
})



# Create merged stats pairwise dataframe
# initialize the merged stats dataframe to contain the taxon column for joining
merged_stats_df <- final_results_bc1[[names(final_results_bc1)[1]]] %>%
  as.data.frame() %>% select(taxon)

walk(comparisons[names(final_results_bc1)], .f = function(comparison){
  
  df <-  final_results_bc1[[comparison]] %>% as.data.frame()
  
  merged_stats_df <<- merged_stats_df %>%
    dplyr::full_join(df, by = join_by("taxon"))
  
})

# Sort ASVs in ascending order
merged_stats_df <- merged_stats_df %>% 
  rename(!!feature := taxon) %>%
  mutate(!!feature := SortMixed(!!sym(feature)))



comp_names <- merged_stats_df %>% 
  select(starts_with("logFC")) %>%
  colnames() %>% str_remove_all("logFC_")
names(comp_names) <- comp_names

message("Making volcano plots...")
# -------------- Make volcano plots ------------------ #
volcano_plots <- map(comp_names, function(comparison){
  
  comp_col  <- c(
    glue("logFC_{comparison}"),
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
  
  p <- ggplot(sub_res_df, aes(x=logFC, y=-log10(pvalue), color=diff, label=!!sym(feature))) +
    geom_point(size=4) +
    scale_color_manual(values=c("TRUE"="red", "FALSE"="black")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    ggrepel::geom_text_repel() + 
    labs(x="logFC", y="-log10(Pvalue)", 
         title = comparison, color="Significant") + publication_format
  
  ggsave(filename = glue("{output_prefix}{comparison}_volcano{assay_suffix}.png"), plot = p, device = "png",
         width = 6, height = 8, units = "in", dpi = 300, path = diff_abund_out_dir)
  
  return(p)
})

number_of_columns <- 2
number_of_rows = ceiling(length(comp_names) / number_of_columns)
fig_height = 7.5 * number_of_rows

p <- wrap_plots(volcano_plots, ncol = 2)
#  Try to combine all the volcano plots in one figure
try(
ggsave(filename = glue("{output_prefix}{feature}_volcano{assay_suffix}.png"), plot = p, device = "png",
       width = 16, height = fig_height, units = "in", dpi = 300,
       path = diff_abund_out_dir, limitsize = FALSE)
)

# Add NCBI id to feature i.e. ASV
tax_names <- map_chr(str_replace_all(taxonomy_table$species, ";_","")  %>%
                       str_split(";"),
                     function(row) row[length(row)])

df <- data.frame(ASV=rownames(taxonomy_table), best_taxonomy=tax_names)

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


samples <- metadata[[samples_column]]
samplesdropped <- setdiff(x = samples, y = colnames(normalized_table)[-1])
missing_df <- data.frame(ASV=normalized_table[[feature]],
                         matrix(data = NA, 
                                nrow = nrow(normalized_table),
                                ncol = length(samplesdropped)
                         )
)
colnames(missing_df) <- c(feature,samplesdropped)


group_levels <- metadata[, group] %>% unique() %>% sort()
group_means_df <- normalized_table[feature]
walk(group_levels, function(group_level){
  
  
  mean_col <- glue("Group.Mean_({group_level})")
  std_col <- glue("Group.Stdev_({group_level})")
  
  # Samples that belong to the current group
  Samples <- metadata %>%
    filter(!!sym(group) == group_level) %>%
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


# Append Mean and standard deviation
normalized_table <- normalized_table %>%
  rowwise() %>%
  mutate(All.Mean=mean(c_across(where(is.numeric))),
         All.Stdev=sd(c_across(where(is.numeric))) )%>% 
  left_join(missing_df, by = feature) %>% 
  select(!!feature, all_of(samples), All.Mean, All.Stdev)


merged_df <- df  %>%
  left_join(taxonomy_table %>%
              as.data.frame() %>%
              rownames_to_column(feature)) %>% 
  select(!!feature, domain:species,everything()) # Try to generalize


merged_df <- merged_df %>%
  select(!!sym(feature):NCBI_id) %>%
  left_join(normalized_table, by = feature) %>%
  left_join(merged_df) %>% 
  left_join(group_means_df, by = feature) %>% 
  mutate(across(where(is.numeric), ~round(.x, digits=3))) %>% 
  mutate(across(where(is.matrix), as.numeric))

output_file <- glue("{diff_abund_out_dir}/{output_prefix}ancombc1_differential_abundance{assay_suffix}.csv")
message("Writing out results of differential abundance using ANCOMBC1...")
write_csv(merged_df,output_file)


#  --------------- Make log abundance box plots ------------------ #

df2 <- (metadata %>% select(!!samples_column, !!group)) %>% 
  left_join(feature_table %>%
              t %>%
              as.data.frame %>%
              rownames_to_column(samples_column))

message("Making abundance box plots...")
boxplots <- map( merged_stats_df[[feature]], function(feature){
  
  p <- ggplot(df2, aes(x=!!sym(group), y=log(!!sym(feature)+1), fill=!!sym(group) )) +
    geom_boxplot() + 
    labs(x=NULL, y="Log Abundance", fill=tools::toTitleCase(group), title = feature) +
    theme_light() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.title.y = element_text(face = "bold", size=12),
          legend.text = element_text(face = "bold", size=10), 
          legend.title = element_text(face = "bold", size=12))
  
  # Save feature boxplot as separate figures
  ggsave(plot = p, filename = glue("{output_prefix}{feature}_boxplot{assay_suffix}.png"), device = "png", 
         width = 8, height = 5, units = "in", dpi = 300, path = diff_abund_out_dir)
  
  return(p)
})

p <- wrap_plots(boxplots, ncol = 2, guides = 'collect')

number_of_features <- merged_stats_df[[feature]] %>% length
number_of_columns <- 2
number_of_rows = ceiling(number_of_features / number_of_columns)
fig_height = 5 * number_of_rows

# Try to Plot all features / ASVs in one figure
try(
ggsave(filename = glue("{output_prefix}{feature}_boxplots{assay_suffix}.png"), plot = p, device = "png",
       width = 14, height = fig_height, units = "in", dpi = 300,
       limitsize = FALSE, path = diff_abund_out_dir)  # There too many things to plot

)

```

**Input Data:**
* **amplicon_metdata.csv** (metadata table)
* **counts_GLAmpSeq.tsv** (count table)
* **taxonomy_GLAmpSeq.tsv** (taxonomy table)

**Output Data:**

* **differential_abundance/ancombc1/<output_prefix><comparison>_volcano_GLAmpSeq.png** (Comparion Volcano Plot)
* **differential_abundance/ancombc1/<output_prefix><feature>_volcano_GLAmpSeq.png** (optional - Combined Volcano Plots)
* **differential_abundance/ancombc1/<output_prefix>ancombc1_differential_abundance_GLAmpSeq.csv** (Statistics Table)
* **differential_abundance/ancombc1/<output_prefix><feature>_boxplot_GLAmpSeq.png** (ASV Boxplots)
* **differential_abundance/ancombc1/<output_prefix><feature>_boxplots_GLAmpSeq.png** (Combined Boxplots)

<br>

---

### 10b. ANCOMBC 2
```bash
Rscript pairwise_ancombc2.R \
                  --metadata-table amplicon_runsheet.csv \
                  --feature-table counts_GLAmpSeq.tsv \
                  --taxonomy-table taxonomy_GLAmpSeq.tsv \
                  --group groups \
                  --samples-column 'Sample Name' \
                  --target-region 16S \
                  --remove-rare FALSE \
                  --prevalence-cutoff 0.15 \
                  --library-cutoff 100 \
                  --cpus 5

```
**Parameter Definitions:**

*	`--metadata-table` – specifies the path to a comma separated samples metadata file with the group/treatment to be analyzed
*	`--feature-table` – specifies the path to a tab separated samples feature table i.e. ASV or OTU table
*	`--taxonomy-table` – specifies the path to a feature taxonomy table i.e. ASV taxonomy table
*	`--group` – specifies the group column in metadata to be analyzed
* `--samples-column` – specifies the column in metadata containing the sample names in the feature table
* `--target-region`  – specifies the amplicon target region. Options are either 16S, 18S or ITS
* `--remove-rare` - should rare features be filtered out prior to analysis? If set, rare feature will be removed
* `--prevalence-cutoff` - If --remove-rare, a numerical fraction between 0 and 1. 
                         Taxa with prevalences(the proportion of samples in which the taxon is present) less than --prevalence-cutoff will be excluded in the analysis. Default is 0.15, i.e. exclude taxa / features that are not present in at least 15% of the samples.
* `--library-cutoff` - If --remove-rare, a numerical threshold for filtering samples based on library sizes. 
                       Samples with library sizes less than lib_cut will be excluded in the analysis. Default is 100. 
                       if you do not want to discard any sample then set to 0.

* `--cpus ` - Specifies the number of cpus to use for parallel processing.

Type `Rscript pairwise_ancombc2.R --help` at the commandline for a full list of available parameters

Content of `pairwise_ancombc2.R`  

```R
diff_abund_out_dir <- "differential_abundance/ancombc2/"
if(!dir.exists(diff_abund_out_dir)) dir.create(diff_abund_out_dir, recursive = TRUE)

# Create phyloseq object
ps <- phyloseq(otu_table(feature_table, taxa_are_rows = TRUE),
               sample_data(metadata),
               tax_table(as.matrix(taxonomy_table)))

# Convert phyloseq to tree summarized experiment object
tse <-  mia::makeTreeSummarizedExperimentFromPhyloseq(ps)

# Getting the reference group and making sure that it is the reference 
# used in the analysis
group_levels <- metadata[, group] %>% unique() %>% sort()
ref_group <- group_levels[1]
tse[[group]] <- factor(tse[[group]] , levels = group_levels)

message("Running ANCOMBC2....")
# Run acombc2
output <- ancombc2(data = tse, assay_name = "counts", tax_level = NULL,
                   fix_formula = group, rand_formula = NULL,
                   p_adj_method = "fdr", pseudo_sens = TRUE,
                   prv_cut = prevalence_cutoff, 
                   lib_cut = library_cutoff, s0_perc = 0.05,
                   group = group, struc_zero = TRUE, neg_lb = FALSE,
                   alpha = 0.05, n_cl = threads, verbose = TRUE,
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
                          if(str_count(colname,group) == 1){
                            str_replace_all(string=colname, 
                                            pattern=glue("(.+)_{group}(.+)"),
                                            replacement=glue("\\1_(\\2)v({ref_group})")) %>% 
                            str_replace(pattern = "^lfc_", replacement = "logFC_") %>% 
                            str_replace(pattern = "^se_", replacement = "lfcSE_") %>% 
                            str_replace(pattern = "^W_", replacement = "Wstat_") %>%
                            str_replace(pattern = "^p_", replacement = "pvalue_") %>%
                            str_replace(pattern = "^q_", replacement = "qvalue_")
                            
                          # Columns with normal two groups comparison
                          } else if(str_count(colname,group) == 2){
                            
                            str_replace_all(string=colname, 
                                            pattern=glue("(.+)_{group}(.+)_{group}(.+)"),
                                            replacement=glue("\\1_(\\2)v(\\3)")) %>% 
                            str_replace(pattern = "^lfc_", replacement = "logFC_") %>% 
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


# ------ Sort columns by group comparisons --------#
# Create a data frame containing only the feature/ASV column
res_df <- paired_stats_df[1] 
walk(uniq_comps, function(comp){
  
  # Get the results for a comparison
  temp_df <- paired_stats_df %>% select(ASV, contains(comp))
  
  # Merge the current comparison to previous comparisons by feature/ASV id
  res_df <<- res_df %>% left_join(temp_df)
})



# --------- Add NCBI id to feature  ---------------#

# Get the best taxonomy assigned to each ASV
tax_names <- map_chr(str_replace_all(taxonomy_table$species, ";_","")  %>%
                       str_split(";"),
                     function(row) row[length(row)])

df <- data.frame(ASV=rownames(taxonomy_table), best_taxonomy=tax_names)

message("Querying NCBI...")
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


samples <- metadata[[samples_column]]
samplesdropped <- setdiff(x = samples, y = colnames(normalized_table)[-1])
missing_df <- data.frame(ASV=normalized_table[[feature]],
           matrix(data = NA, 
                  nrow = nrow(normalized_table),
                  ncol = length(samplesdropped)
                  )
           )
colnames(missing_df) <- c(feature,samplesdropped)


group_means_df <- normalized_table[feature]
walk(group_levels, function(group_level){
  
  
  mean_col <- glue("Group.Mean_({group_level})")
  std_col <- glue("Group.Stdev_({group_level})")
  
  # Samples that belong to the current group
  Samples <- metadata %>%
    filter(!!sym(group) == group_level) %>%
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


# Calculate global mean and standard deviation
normalized_table <- normalized_table %>%
  rowwise() %>%
  mutate(All.Mean=mean(c_across(where(is.numeric))),
         All.Stdev=sd(c_across(where(is.numeric))) ) %>% 
  left_join(missing_df, by = feature) %>% 
  select(!!feature, all_of(samples), All.Mean, All.Stdev)

# Append the taxonomy table to the ncbi and stats table
merged_df <- df  %>%
  left_join(taxonomy_table %>%
              as.data.frame() %>%
              rownames_to_column(feature)) %>% 
  select(!!feature,domain:species,everything())

# Add group means and normalized table
merged_df <- merged_df %>%
  select(!!sym(feature):NCBI_id) %>%
  left_join(normalized_table, by = feature) %>%
  left_join(merged_df) %>% 
  left_join(group_means_df, by = feature) %>% 
  mutate(across(where(is.numeric), ~round(.x, digits=3)))

message("Writing out results of differential abundance using ANCOMBC2...")
output_file <- glue("{diff_abund_out_dir}/{output_prefix}ancombc2_differential_abundance{assay_suffix}.csv")
write_csv(merged_df,output_file)


# ---------------------- Visualization --------------------------------------- #
message("Making volcano plots...")
# ------------ Make volcano ---------------- #
volcano_plots <- map(uniq_comps, function(comparison){
  
  comp_col  <- c(
    glue("logFC_{comparison}"),
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
  
  p <- ggplot(sub_res_df, aes(x=logFC, y=-log10(pvalue), color=diff, label=!!sym(feature))) +
    geom_point(size=4) + geom_point(size=4) +
    scale_color_manual(values=c("TRUE"="red", "FALSE"="black")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    ggrepel::geom_text_repel() + 
    labs(x="logFC", y="-log10(Pvalue)", 
         title = comparison, color="Significant") + publication_format


  ggsave(filename = glue("{output_prefix}{comparison}_volcano{assay_suffix}.png"), plot = p, device = "png",
         width = 6, height = 8, units = "in", dpi = 300, path = diff_abund_out_dir)
  
  return(p)
  
})

p <- wrap_plots(volcano_plots, ncol = 2)


number_of_columns <- 2
number_of_rows = ceiling(length(uniq_comps) / number_of_columns)
fig_height = 7.5 * number_of_rows

#  Try to combine all the volcano plots in one figure
try(
ggsave(filename = glue("{output_prefix}{feature}_volcano{assay_suffix}.png"), plot = p, device = "png", 
       width = 16, height = fig_height, units = "in",
       dpi = 300, limitsize = FALSE, path=diff_abund_out_dir)
)

# ------------- Box plots ---------------- #

df2 <- (metadata %>% select(!!samples_column, !!group)) %>% 
  left_join(feature_table %>%
              t %>%
              as.data.frame %>%
              rownames_to_column(samples_column))

message("Making abundance box plots...")
boxplots <- map(res_df[[feature]], function(feature){
  
  p <- ggplot(df2, aes(x=!!sym(group), y=log(!!sym(feature)+1), fill=!!sym(group) )) +
    geom_boxplot() + 
    labs(x=NULL, y="Log Abundance", fill=tools::toTitleCase(group), title = feature) +
    theme_light() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.title.y = element_text(face = "bold", size=12),
          legend.text = element_text(face = "bold", size=10), 
          legend.title = element_text(face = "bold", size=12))
  
  ggsave(filename = glue("{output_prefix}{feature}_boxplot{assay_suffix}.png"), plot = p, device = "png",
         width = 8, height = 5, units = "in", dpi = 300, path = diff_abund_out_dir)
  
  return(p)
})


p <- wrap_plots(boxplots, ncol = 2, guides = 'collect')

number_of_features <- res_df[[feature]] %>% length
number_of_columns <- 2
number_of_rows = ceiling(number_of_features / number_of_columns)
fig_height = 5 * number_of_rows

# Try to Plot all features / ASVs in one figure
try(
ggsave(filename = glue("{output_prefix}{feature}_boxplots{assay_suffix}.png"), plot = p, device = "png",
       width = 14, height = fig_height, units = "in", dpi = 300,
       path = diff_abund_out_dir, limitsize = FALSE)
)
```

**Input Data:**

* **amplicon_metdata.csv** (metadata table)
* **counts_GLAmpSeq.tsv** (count table)
* **taxonomy_GLAmpSeq.tsv** (taxonomy table)

**Output Data:**

* **differential_abundance/ancombc2/<output_prefix><comparison>_volcano_GLAmpSeq.png** (Comparion Volcano Plot)
* **differential_abundance/ancombc2/<output_prefix><feature>_volcano_GLAmpSeq.png** (optional - Combined Volcano Plots)
* **differential_abundance/ancombc2/<output_prefix>ancombc2_differential_abundance_GLAmpSeq.csv** (Statistics Table)
* **differential_abundance/ancombc2/<output_prefix><feature>_boxplot_GLAmpSeq.png** (ASV Boxplots)
* **differential_abundance/ancombc2/<output_prefix><feature>_boxplots_GLAmpSeq.png** (Combined Boxplots)


<br>

---

### 10c. DESeq2

```bash
Rscript run_deseq2.R \
                  --metadata-table amplicon_runsheet.csv \
                  --feature-table counts_GLAmpSeq.tsv \
                  --taxonomy-table taxonomy_GLAmpSeq.tsv \
                  --group groups \
                  --samples-column 'Sample Name' \
                  --target-region 16S \
                  --remove-rare FALSE \
                  --prevalence-cutoff 0.15 \
                  --library-cutoff 100

```
**Parameter Definitions:**

*	`--metadata-table` – specifies the path to a comma separated samples metadata file with the group/treatment to be analyzed
*	`--feature-table` – specifies the path to a tab separated samples feature table i.e. ASV or OTU table
*	`--taxonomy-table` – specifies the path to a feature taxonomy table i.e. ASV taxonomy table
*	`--group` – specifies the group column in metadata to be analyzed
* `--samples-column` – specifies the column in metadata containing the sample names in the feature table
* `--target-region`  – specifies the amplicon target region. Options are either 16S, 18S or ITS
* `--remove-rare` - should rare features be filtered out prior to analysis? If set, rare feature will be removed
* `--prevalence-cutoff` - If --remove-rare, a numerical fraction between 0 and 1. 
                         Taxa with prevalences(the proportion of samples in which the taxon is present) less than --prevalence-cutoff will be excluded in the analysis. Default is 0.15, i.e. exclude taxa / features that are not present in at least 15% of the samples.
* `--library-cutoff` - If --remove-rare, a numerical threshold for filtering samples based on library sizes. 
                       Samples with library sizes less than lib_cut will be excluded in the analysis. Default is 100. 
                       if you do not want to discard any sample then set to 0.

type `Rscript run_deseq2.R --help` at the commandline for a full list of available parameters

Content of `run_deseq2.R`  

```R
# Create output directory if it doesn't already exist
diff_abund_out_dir <- "differential_abundance/deseq2/"
if(!dir.exists(diff_abund_out_dir)) dir.create(diff_abund_out_dir, recursive = TRUE)

#### pairwise comparisons
unique_groups <- unique(metadata[[group]])

# Create phyloseq object
ASV_physeq <- phyloseq(otu_table(feature_table, taxa_are_rows = TRUE),
                       tax_table(as.matrix(taxonomy_table)),
                       sample_data(metadata))

deseq_obj <- phyloseq_to_deseq2(physeq = ASV_physeq,
                                design = reformulate(group))

# Add pseudocount if any 0 count samples are present
if (sum(colSums(counts(deseq_obj)) == 0) > 0) {
  count_data <- counts(deseq_obj) + 1 
  
  count_data <- as.matrix(apply(count_data, 2, as.integer))
  rownames(count_data) <- rownames(counts(deseq_obj))
  colnames(count_data) <- colnames(counts(deseq_obj))
  counts(deseq_obj) <- count_data
}

# Run Deseq
# https://rdrr.io/bioc/phyloseq/src/inst/doc/phyloseq-mixture-models.R 
deseq_modeled <- tryCatch({
  # Attempt to run DESeq
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
pairwise_comp.m <- utils::combn(metadata[,group] %>% unique, 2)
pairwise_comp_df <- pairwise_comp.m %>% as.data.frame 

colnames(pairwise_comp_df) <- map_chr(pairwise_comp_df,
                                      \(col) str_c(col, collapse = "v"))
comparisons <- colnames(pairwise_comp_df)
names(comparisons) <- comparisons

# Retrieve statistics table
merged_stats_df <-  data.frame(ASV=rownames(feature_table))
colnames(merged_stats_df) <- feature

walk(pairwise_comp_df, function(col){
  
  group1 <- col[1]
  group2 <- col[2]
  
df <- results(deseq_modeled, contrast = c(group, group1, group2)) %>%
  data.frame() %>%
  rownames_to_column(feature) %>% 
  set_names(c(feature ,
              glue("baseMean_({group1})v({group2})"),
              glue("log2FC_({group1})v({group2})"),
              glue("lfcSE_({group1})v({group2})"), 
              glue("stat_({group1})v({group2})"), 
              glue("pvalue_({group1})v({group2})"),
              glue("padj_({group1})v({group2})") 
            ))

            
  merged_stats_df <<- merged_stats_df %>% 
                          dplyr::left_join(df, join_by("ASV"))
})



# Add NCBI id to feature i.e. ASV
tax_names <- map_chr(str_replace_all(taxonomy_table$species, ";_","")  %>%
                       str_split(";"),
                     function(row) row[length(row)])

df <- data.frame(ASV=rownames(taxonomy_table), best_taxonomy=tax_names)

# Pull NCBI IDS for unique taxonomy names
df2 <- data.frame(best_taxonomy = df$best_taxonomy %>%
                    unique()) %>%
  mutate(NCBI_id=get_ncbi_ids(best_taxonomy, target_region),
         .after = best_taxonomy)

df <- df %>%
  left_join(df2, join_by("best_taxonomy")) %>% 
  right_join(merged_stats_df)




group_levels <- metadata[, group] %>% unique() %>% sort()
normalized_table <- counts(deseq_modeled, normalized=TRUE) %>% 
                        as.data.frame() %>%
                        rownames_to_column(feature)

# Creating a dataframe of samples that were dropped because they didn't
# meet are cut-off criteria
samples <- metadata[[samples_column]]
samplesdropped <- setdiff(x = samples, y = colnames(normalized_table)[-1])
missing_df <- data.frame(ASV=normalized_table[[feature]],
                         matrix(data = NA, 
                                nrow = nrow(normalized_table),
                                ncol = length(samplesdropped)
                         )
)
colnames(missing_df) <- c(feature,samplesdropped)


group_means_df <- normalized_table[feature]
walk(group_levels, function(group_level){
  
  
  mean_col <- glue("Group.Mean_({group_level})")
  std_col <- glue("Group.Stdev_({group_level})")
  
  # Samples that belong to the current group
  Samples <- metadata %>%
    filter(!!sym(group) == group_level) %>%
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


# Append Mean and standard deviation
normalized_table <- normalized_table %>%
  rowwise() %>%
  mutate(All.Mean=mean(c_across(where(is.numeric))),
         All.Stdev=sd(c_across(where(is.numeric))) )%>% 
  left_join(missing_df, by = feature) %>% 
  select(!!feature, all_of(samples), All.Mean, All.Stdev)


# Add taxonomy
merged_df <- df  %>%
  left_join(taxonomy_table %>%
              as.data.frame() %>%
              rownames_to_column(feature)) %>% 
  select(!!feature, domain:species,everything()) # Try to generalize

# Merge all prepared tables 
merged_df <- merged_df %>%
  select(!!sym(feature):NCBI_id) %>%
  left_join(normalized_table, by = feature) %>%
  left_join(merged_df) %>% 
  left_join(group_means_df, by = feature) %>% 
  mutate(across(where(is.numeric), ~round(.x, digits=3))) %>% 
  mutate(across(where(is.matrix), as.numeric))


output_file <- glue("{diff_abund_out_dir}/{output_prefix}deseq2_differential_abundance{assay_suffix}.csv")
message("Writing out results of differential abundance using DESeq2...")
write_csv(merged_df,output_file)



# Make volcano plots
walk(pairwise_comp_df, function(col){
  
  group1 <- col[1]
  group2 <- col[2]
  
  plot_width_inches <- 11.1
  plot_height_inches <- 8.33
  p_val <- 0.1 #also logfc cutoff?
  
  deseq_res <- results(deseq_modeled, contrast = c(group, group1, group2))
  volcano_data <- as.data.frame(deseq_res)
  
  
  volcano_data <- volcano_data[!is.na(volcano_data$padj), ]
  volcano_data$significant <- volcano_data$padj <= p_val #also logfc cutoff?
  
  ###### Long x-axis label adjustments ##########
  x_label <- paste("Log2 Fold Change\n(",group1," vs ",group2,")")
  label_length <- nchar(x_label)
  max_allowed_label_length <- plot_width_inches * 10
  
  # Construct x-axis label with new line breaks if was too long
  if (label_length > max_allowed_label_length){
    x_label <- paste("Log2 Fold Change\n\n(", group1, "\n vs \n", group2, ")", sep="")
  }
  #######################################
  
  # ASVs promoted in space on right, reduced on left
  p <- ggplot(volcano_data, aes(x=log2FoldChange, y=-log10(padj), 
                                color=significant)) +
    geom_point(alpha=0.7, size=2) +
    geom_hline(yintercept = -log10(p_val), linetype = "dashed") +
    scale_color_manual(values=c("black", "red"), 
                       labels=c(paste0("padj > ", p_val), 
                                paste0(" padj \u2264 ", p_val))) +
    theme_bw() +
    labs(title="Volcano Plot",
         x=x_label,
         y="-Log10 P-value",
         color=paste0("")) +
    theme(legend.position="top")
  
  # label points and plot
  top_points <- volcano_data %>%
    arrange(padj) %>%
    filter(significant) %>%
    head(10)
  
  volcano_plot <- p + geom_text_repel(data=top_points, 
                                      aes(label=row.names(top_points)),
                                      size=3)
  
  # Save volcano plot
  ggsave(filename=glue("{diff_abund_out_dir}{output_prefix}volcano_{group1}_vs_{group2}.png"),
         plot=volcano_plot,
         width = plot_width_inches, 
         height = plot_height_inches, 
         dpi = 300)
})
```

**Input Data:**

* **amplicon_metdata.csv** (metadata table)
* **counts_GLAmpSeq.tsv** (count table)
* **taxonomy_GLAmpSeq.tsv** (taxonomy table)

**Output Data:**

* **differential_abundance/deseq2/<output_prefix>volcano_<group1>_vs_<group2>.png** (Comparion Volcano Plot)
* **differential_abundance/deseq2/<output_prefix>deseq2_differential_abundance_GLAmpSeq.csv** (Statistics Table)


<br>

---