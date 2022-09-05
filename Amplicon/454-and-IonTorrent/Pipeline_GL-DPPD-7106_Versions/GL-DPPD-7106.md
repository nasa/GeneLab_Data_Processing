# Bioinformatics pipeline for 454 and IonTorrent amplicon sequencing data  

> **This page holds an overview and some example commands of how GeneLab processes 454 and IonTorrent amplicon datasets. Exact processing commands for specific datasets that have been released are provided with their processed data in the [GeneLab Data Systems (GLDS) repository](https://genelab-data.ndc.nasa.gov/genelab/projects).**  

---

**Date:** January 13, 2021  
**Revision:** -  
**Document Number:** GL-DPPD-7106  

**Submitted by:**  
Michael D. Lee (GeneLab Analysis Team)  

**Approved by:**  
Sylvain Costes (GeneLab Project Manager)  
Samrawit Gebre (GeneLab Deputy Project Manager and Interim GeneLab Configuration Manager)  
Amanda Saravia-Butler (GeneLab Data Processing Lead)  
Jonathan Galazka (GeneLab Project Scientist)  

---

# Table of contents  

- [**Software used**](#software-used)
- [**Reference databases used**](#reference-databases-used)
- [**General processing overview with example commands**](#general-processing-overview-with-example-commands)
  - [**1. Raw Data QC**](#1-raw-data-qc)
    - [1a. Compile Raw Data QC](#1a-compile-raw-data-qc)
  - [**2. Trim Primers**](#2-trim-primers)
  - [**3. Quality filtering**](#3-quality-filtering)
  - [**4. Filtered Data QC**](#4-filtered-data-qc)
    - [4a. Compile Filtered Data QC](#4a-compile-filtered-data-qc)
  - [**5. Generating OTUs and counts per sample**](#5-generating-otus-and-counts-per-sample)
    - [5a. Dereplicate individual samples](#5a-dereplicate-individual-samples)
    - [5b. Generate OTUs](#5b-generate-otus)
    - [5c. Map reads to OTUs](#5c-map-reads-to-otus)
  - [**6. Generating taxonomy and additional outputs**](#6-generating-taxonomy-and-additional-outputs)
    - [6a. Assigning taxonomy](#6a-assigning-taxonomy)
    - [6b. Generating and writing outputs](#6b-generating-and-writing-outputs)

---

# Software used

|Program|Version*|Relevant Links|
|:------|:-----:|-------------:|
|FastQC|`fastqc -v`|[https://www.bioinformatics.babraham.ac.uk/projects/fastqc/](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)|
|MultiQC|`multiqc -v`|[https://multiqc.info/](https://multiqc.info/)|
|Cutadapt|`cutadapt --version`|[https://cutadapt.readthedocs.io/en/stable/](https://cutadapt.readthedocs.io/en/stable/)|
|bbduk|`bbduk.sh --version`|[https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/)|
|vsearch|`vsearch --version`|[https://github.com/torognes/vsearch](https://github.com/torognes/vsearch)|
|R|`R --version` (at command line) | [https://www.r-project.org/](https://www.r-project.org/)|
|DECIPHER|`packageVersion("DECIPHER")`|[https://bioconductor.org/packages/release/bioc/html/DECIPHER.html](https://bioconductor.org/packages/release/bioc/html/DECIPHER.html)|
|biomformat|`packageVersion("biomformat")`|[https://github.com/joey711/biomformat](https://github.com/joey711/biomformat)|

>**\*** Exact versions utilized for a given dataset are available along with the processing commands for each specific dataset (this is due to how the system may need to be updated regularly).

# Reference databases used

|Program used| Database| Relevant Links|
|:-----|:-----:|--------:|
|DECIPHER| SILVA SSU r138 | [http://www2.decipher.codes/Classification/TrainingSets/SILVA_SSU_r138_2019.RData](http://www2.decipher.codes/Classification/TrainingSets/)|
|DECIPHER| UNITE v2020 | [http://www2.decipher.codes/Classification/TrainingSets/UNITE_v2020_February2020.RData](http://www2.decipher.codes/Classification/TrainingSets/)|

---

# General processing overview with example commands

> Exact processing commands for specific datasets are provided with their processed data in the [GeneLab Data Systems (GLDS) repository](https://genelab-data.ndc.nasa.gov/genelab/projects).  

---

## 1. Raw Data QC
**Uses [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)**

```
fastqc -o raw_fastqc_output *raw.fastq.gz
```

**Parameter Definitions:**

* `-o` – the output directory to store results
* `*raw.fastq.gz` – the input reads are specified as a positional argument, and can be given all at once with wildcards like this, or as individual arguments with spaces in between them

**Input data:**

* fastq, compressed or uncompressed

**Output data:**

* fastqc.html (FastQC output html summary)
* fastqc.zip (FastQC output data)


### 1a. Compile Raw Data QC
**Uses [MultiQC](https://multiqc.info/)**

```
multiqc -o raw_multiqc_output -n raw_multiqc -z raw_fastqc_output/
```

**Parameter Definitions:**

*	`-o` – the output directory to store results
* `-n` – the filename prefix of results
* `-z` – specifies to zip the output data directory
*	`raw_fastqc_output/` – the directory holding the output data from the fastqc run, provided as a positional argument

**Input data:**

* raw_fastqc_output/*fastqc.zip (FastQC output data)

**Output data:**

* raw_multiqc_output/raw_multiqc_report.html (multiqc output html summary)
* raw_multiqc_output/raw_multiqc_data.zip (zipped directory containing multiqc output data)

---

## 2. Trim Primers
**Uses [cutadapt](https://cutadapt.readthedocs.io/en/stable/)**

The location and orientation of primers in the data is important to understand in deciding how to do this step. `cutadapt` has many options for primer identification and removal. They are described in detail on their documentation page here: [https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types](https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types)  

The following example command shows how it was done for [GLDS-72](https://genelab-data.ndc.nasa.gov/genelab/accession/GLDS-72/), which was 454 sequencing of the 16S gene using these primers:  
* forward: 5’-AGAGTTTGATCCTGGCTCAG-3’  
* reverse: 5’- GCTGCCTCCCGTAGGAGT-3’  


```
cutadapt -g AGAGTTTGATCCTGGCTCAG -a GCTGCCTCCCGTAGGAGT \
         -o sample-1_trimmed.fastq.gz sample-1_raw.fastq.gz
```

**Parameter Definitions:**

*	`-g` – specifies the forward primer 

*	`-a` – specifies the reverse primer

*	`-o` – specifies output primer-trimmed reads

*	`sample-1_raw.fastq.gz` – positional argument specifying the input reads


**Input data:**

* sample-1_raw.fastq.gz (raw reads)

**Output data:**

* sample-1_trimmed.fastq.gz (primer-trimmed reads)

---

## 3. Quality filtering
**Uses [bbduk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/)**

```
bbduk.sh in=sample-1_trimmed.fastq.gz out=sample-1_filtered.fastq.gz \
         qtrim=r trimq=10 mlf=0.5 minavgquality=20 minlength=50
```

**Parameter Definitions:**

*	`in=` – the input file

*	`out=` – the output file

*	`qtrim=` – specifies the direction to trim low-quality bases from

*	`trimq=` – the minimum quality to keep while trimming low-quality bases

*	`mlf=` – allowed minimum length after trimming is half the initial read size (read filtered out if shorter)

*	`minavgquality=` – minimum average quality allowed after trimming (read filtered out if below)

*	`minlength=` – allowed minimum length (backup to `mlf` setting if starting read is shorter than 100 bps)

**Input data:**

* sample-1_trimmed.fastq.gz (primer-trimmed reads)

**Output data:**

* sample-1_filtered.fastq.gz (filtered reads)

---

## 4. Filtered Data QC
**Uses [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)**

```
fastqc -o filtered_fastqc_output *filtered.fastq.gz
```

**Parameter Definitions:**

*	`-o` – the output directory to store results  
*	`*filtered.fastq.gz` – the input reads are specified as a positional argument, and can be given all at once with wildcards like this, or as individual arguments with spaces in between them  

**Input data:**

* *filtered.fastq.gz (filtered reads)

**Output data:**

* *fastqc.html (FastQC output html summary)
* *fastqc.zip (FastQC output data)


### 4a. Compile Filtered Data QC
**Uses [MultiQC](https://multiqc.info/)**

```
multiqc -o filtered_multiqc_output -n filtered_multiqc -z filtered_fastqc_output
```

**Parameter Definitions:**

*	`-o` – the output directory to store results
* `-n` – the filename prefix of results
* `-z` – specifies to zip the output data directory
*	`filtered_fastqc_output` – the directory holding the output data from the fastqc run, provided as a positional argument

**Input data:**

* filtered_fastqc_output/*fastqc.zip (FastQC output data)

**Output data:**

* filtered_multiqc_output/filtered_multiqc_report.html (multiqc output html summary)
* filtered_multiqc_output/filtered_multiqc_data.zip (zipped directory containing multiqc output data)


---

## 5. Generating OTUs and counts per sample
**Uses [vsearch](https://github.com/torognes/vsearch)**

### 5a. Dereplicate individual samples
```
vsearch --derep_fulllength sample-1_filtered.fastq.gz --strand both \
        --output sample-1_derep.fasta --sizeout --relabel "sample=sample-1;seq_" 
```

**Parameter Definitions:**

*	`--derep_fulllength` – the vsearch command
*	`sample-1_filtered.fastq.gz` – input file, provided as a positional argument
*  `--strand both` - specifies to check both orientations
*  `--output` - designate the name of the output fasta file
*  `--sizeout` - incorporates abundance information in the sequence header
*  `--relabel` - relabel the headers of the sequences starting with this prefix

**Input data:**

* sample-1_filtered.fastq.gz (filtered reads)

**Output data:**

* sample-1_derep.fasta (dereplicated reads)


### 5b. Generate OTUs

#### Combining all individual sample dereplicated sequences for further processing
```bash
cat *_derep.fasta > all-samples-seqs.fasta
rm *_derep.fasta
```

#### Dereplicating combined sequences
```bash
vsearch --derep_fulllength all-samples-seqs.fasta --strand both \
        --output all-samples_derep.fasta --sizein --sizeout
```

**Parameter Definitions:**

*	`--derep_fulllength` – the vsearch command
*	`all-samples-seqs.fasta` – input file, provided as a positional argument
*  `--strand both` - specifies to check both orientations
*  `--output` - designate the name of the output fasta file
*  `--sizein` - take into account abundance information in header
*  `--sizeout` - incorporates abundance information in the sequence header


**Input data:**

* all-samples-seqs.fasta (combined, individual-sample-dereplicated reads)

**Output data:**

* all-samples_derep.fasta (combined dereplicated reads)


#### Clustering to get representative sequences
```bash
vsearch --cluster_size all-samples_derep.fasta --id 0.97 --strand both \
        --sizein --sizeout --relabel "OTU_" --centroids rep-seqs.fasta
```

**Parameter Definitions:**

*	`--cluster_size ` – the vsearch command
*	`all-samples_derep.fasta` – input file, provided as a positional argument
*  `--id` - specifies the percent identity to cluster at
*  `--strand both` - specifies to consider both orientations
*  `--sizein` - take into account abundance information in header
*  `--sizeout` - incorporates abundance information in the sequence header
*  `--relabel` - relabel the headers of the sequences starting with this prefix
*  `--centroids` - designate the name of the output fasta file holding representative sequences


**Input data:**

* all-samples_derep.fasta (combined dereplicated reads)

**Output data:**

* rep-seqs.fasta (representative sequences)


#### Removing singletons
```bash
vsearch --sortbysize rep-seqs.fasta --minsize 2 --output rep-seqs-no-singletons.fasta
```

**Parameter Definitions:**

*	`--sortbysize ` – the vsearch command
*	`rep-seqs.fasta` – input file, provided as a positional argument
*  `--minsize` - minimum cluster size to be retained
*  `--output` - designate the name of the output fasta file holding filtered representative sequences


**Input data:**

* rep-seqs.fasta (representative sequences)

**Output data:**

* rep-seqs-no-singletons.fasta (representative sequences, with singletons removed)

#### Chimera check and removal

```bash
vsearch --uchime_denovo rep-seqs-no-singletons.fasta --sizein --nonchimeras OTUs.fasta --relabel "OTU_"
```

**Parameter Definitions:**

*	`--uchime_denovo ` – the vsearch command
*	`rep-seqs-no-singletons.fasta` – input file, provided as a positional argument
*  `--sizein` - take into account abundance information in header
*  `--nonchimeras` - designate the name of the output fasta file holding filtered representative sequences
*  `--relabel` - relabel the headers of the sequences starting with this prefix

**Input data:**

* rep-seqs-no-singletons.fasta (representative sequences, with singletons removed)

**Output data:**

* OTUs.fasta (representative sequences, with singletons and chimeras removed)


### 5c. Map reads to OTUs
```bash
vsearch --usearch_global all-samples-seqs.fasta -db OTUs.fasta --sizein \
        --id 0.97 --otutabout - | sed 's/^#OTU ID/OTU_ID/' > counts.tsv
```

**Parameter Definitions:**

*	`--usearch_global ` – the vsearch command
*	`all-samples-seqs.fasta` – input file, provided as a positional argument (note this is our combined individual-dereplicated samples' fasta of reads generated in step 5b) 
*  `--db` - specifies the reference OTUs to map to
*  `--id` - specifies the minimum percent identity to enable mapping a read to an OTU
*  `--otutabout` - designates the output of the count table (here going to stout in order to be modified in the subsequent `sed` command
*  `| sed ... > counts.tsv` - renaming the first column header and writing to `counts.tsv`


**Input data:**

* all-samples-seqs.fasta (combined, individual-sample-dereplicated reads)
* OTUs.fasta (representative sequences, with singletons and chimeras removed)

**Output data:**

* counts.tsv (count table of representative sequences per sample)


---

## 6. Generating taxonomy and additional outputs

> The following is performed within [R](https://www.r-project.org/).  

**Uses [DECIPHER](https://bioconductor.org/packages/release/bioc/html/DECIPHER.html) and [biom-format](https://github.com/joey711/biomformat)**

### 6a. Assigning taxonomy

Reading in OTU sequences:
```R
dna <- readDNAStringSet("OTUs.fasta")
```

Downloading the reference R taxonomy object:
```R
download.file(url=“http://www2.decipher.codes/Classification/TrainingSets/SILVA_SSU_r138_2019.RData”,
              destfile=“SILVA_SSU_r138_2019.RData”)
```

**Parameter Definitions:**  

*	`download.file()` – the function we are calling, with the following parameters set within it

*	`url=` – specifying the url address of the file to download

*	`destfile=` – specifying the path/name of the file after downloading


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


### 6b. Generating and writing outputs
Creating table of taxonomy and setting any that are unclassified as "NA", and writing out:

```R
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")

tax_tab <- t(sapply(tax_info, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))

colnames(tax_tab) <- ranks
row.names(tax_tab) <- NULL
otu_ids <- names(tax_info)
tax_tab <- data.frame("OTU_ID"=otu_ids, tax_tab, check.names=FALSE)

write.table(tax_tab, "taxonomy.tsv", sep = "\t", quote=F, row.names=FALSE)
```

Reading in counts table and generating additional outputs:

```R
otu_tab <- read.table("counts.tsv", sep="\t", header=TRUE, check.names=FALSE)

    # generating and writing out biom file format
biom_object <- make_biom(data=otu_tab, observation_metadata=tax_tab)
write_biom(biom_object, "taxonomy-and-counts.biom")

    # making a tsv of combined tax and counts
tax_and_count_tab <- merge(tax_tab, otu_tab)
write.table(tax_and_count_tab, "taxonomy-and-counts.tsv", sep="\t", quote=FALSE, row.names=FALSE)
```

**Input data:**

* OTUs.fasta (representative sequences, with singletons and chimeras removed)
* counts.tsv (count table of representative sequences per sample)

**Output data:**

* taxonomy-and-counts.tsv (table with taxonomy of representative sequences and counts per sample)
* taxonomy-and-counts.biom (biom formatted file with taxonomy of representative sequences and counts per sample)

---

