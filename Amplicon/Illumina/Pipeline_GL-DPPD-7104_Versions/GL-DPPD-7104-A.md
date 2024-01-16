# Bioinformatics pipeline for amplicon Illumina sequencing data  

> **This page holds an overview and instructions for how GeneLab processes Illumina amplicon datasets. Exact processing commands for specific datasets that have been released are available in the [GLDS_Processing_Scripts](../GLDS_Processing_Scripts) sub-directory and/or are provided with their processed data in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/).**  

---

**Date:** May 13, 2020  
**Revision:** A  
**Document Number:** GL-DPPD-7104-A  

**Submitted by:**  
Michael D. Lee (GeneLab Data Processing Team)

**Approved by:**  
Sylvain Costes (GeneLab Project Manager)  
Samrawit Gebre (GeneLab Deputy Project Manager)  
Homer Fogle (GeneLab Data Processing Representative)  
Jonathan Galazka (GeneLab Project Scientist)  
Anushree Sonic (Genelab Configuration Manager)  

---

# Table of contents  

- [**Software used**](#software-used)
- [**Reference databases used**](#reference-databases-used)
- [**General processing overview with example commands**](#general-processing-overview-with-example-commands)
  - [**1. Raw Data QC**](#1-raw-data-qc)
    - [Compile Raw Data QC](#compile-raw-data-qc)
  - [**2. Trim Primers**](#2-trim-primers)
  - [**3. Quality filtering**](#3-quality-filtering)
  - [**4. Filtered Data QC**](#4-filtered-data-qc)
    - [Compile Filtered Data QC](#compile-filtered-data-qc)
  - [**5. Calculate error model, apply DADA2 algorithm, assign taxonomy, and create output tables**](#5-calculate-error-model-apply-dada2-algorithm-assign-taxonomy-and-create-output-tables)
    - [Learning the error rates](#learning-the-error-rates)
    - [Inferring sequences](#inferring-sequences)
    - [Merging forward and reverse reads](#merging-forward-and-reverse-reads-not-needed-if-data-are-single-end)
    - [Generating sequence table with counts per sample](#generating-sequence-table-with-counts-per-sample)
    - [Removing putative chimeras](#removing-putative-chimeras)
    - [Assigning taxonomy](#assigning-taxonomy)
    - [Generating and writing standard outputs](#generating-and-writing-standard-outputs)

---

# Software used  

|Program|Version*|Relevant Links|
|:------|:------:|:-------------|
|FastQC|`0.11.9`|[https://www.bioinformatics.babraham.ac.uk/projects/fastqc/](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)|
|MultiQC|`1.9`|[https://multiqc.info/](https://multiqc.info/)|
|Cutadapt|`2.3`|[https://cutadapt.readthedocs.io/en/stable/](https://cutadapt.readthedocs.io/en/stable/)|
|DADA2|`1.20.0`|[https://www.bioconductor.org/packages/release/bioc/html/dada2.html](https://www.bioconductor.org/packages/release/bioc/html/dada2.html)|
|DECIPHER|`2.20.0`|[https://bioconductor.org/packages/release/bioc/html/DECIPHER.html](https://bioconductor.org/packages/release/bioc/html/DECIPHER.html)|
|biomformat|`1.20.0`|[https://github.com/joey711/biomformat](https://github.com/joey711/biomformat)|

>**\*** Exact versions are available along with the processing commands for each specific dataset.

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
* `*.fastq.gz` – the input reads are specified as a positional argument, and can be given all at once with wildcards like this, or as individual arguments with spaces in between them

**Input Data:**

* fastq, compressed or uncompressed

**Output Data:**

* fastqc.html (FastQC output html summary)
* fastqc.zip (FastQC output data)


<br>  

### Compile Raw Data QC  

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

## 3. Quality filtering
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

### Compile Filtered Data QC
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

## 5. Calculate error model, apply DADA2 algorithm, assign taxonomy, and create output tables
> The following is run in an R environment.  

These example commands as written assumes paired-end data, with notes included on what would be different if working with single-end data. The taxonomy reference database used below is as an example only, suitable for the example 16S dataset ([GLDS-200](https://osdr.nasa.gov/bio/repo/data/studies/OSD-200)) used here. But others designed for DECIPHER can be found here: [http://www2.decipher.codes/Downloads.html](http://www2.decipher.codes/Downloads.html)  

<br>

### Learning the error rates
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

### Inferring sequences
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

### Merging forward and reverse reads; not needed if data are single-end
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

### Generating sequence table with counts per sample
```R
seqtab <- makeSequenceTable(merged_contigs)
```

**Parameter Definitions:**  

*	If single-end data, instead of “merged_contigs”, the forward_seqs object would be provided to the `makeSequenceTable()` function here

<br>

### Removing putative chimeras
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

### Assigning taxonomy

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

### Generating and writing standard outputs

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

* **ASVs_GLAmpSeq.fasta** (inferred sequences)
* **counts_GLAmpSeq.tsv** (count table)
* **taxonomy_GLAmpSeq.tsv** (taxonomy table)
* **taxonomy-and-counts_GLAmpSeq.tsv** (combined taxonomy and count table)
* **taxonomy-and-counts_GLAmpSeq.biom** (count and taxonomy table in biom format)
* **read-count-tracking_GLAmpSeq.tsv** (read counts at each processing step)

