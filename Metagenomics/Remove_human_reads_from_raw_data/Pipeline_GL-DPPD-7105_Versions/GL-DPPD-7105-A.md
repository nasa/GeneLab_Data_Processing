# GeneLab removal of human reads from metagenomics datasets

> **It is NASA's policy that any human reads are to be removed from metagenomics datasets prior to being hosted in [GeneLab's data repository](https://genelab-data.ndc.nasa.gov/genelab/projects). As such, all metagenomics datasets are screened against a human reference-genome [kraken2](https://github.com/DerrickWood/kraken2/wiki) database. This document holds an overview and some example commands of how GeneLab does performs this.**

---

**Date:**  January, 13, 2021  
**Revision:** A  
**Document Number:** GL-DPPD-7105-A  

**Submitted by:**  
Michael D. Lee (GeneLab Science Team)  

**Approved by:**  
Sylvain Costes (GeneLab Project Manager)  
Samrawit Gebre (GeneLab Deputy Project Manager and Interim GeneLab Configuration Manager)  
Amanda Saravia-Butler (GeneLab Data Processing Lead)  
Jonathan Galazka (GeneLab Project Scientist)  

## Updates from previous revision
* In [Step 1](#1-build-kraken2-database), added the `--no-masking` flag to the `--download-library` command to err on the side of being more conservative in that it will now filter out reads that match to anywhere in the human reference genome (including low-complexity regions). 
* Also in [Step 1](#1-build-kraken2-database), used the default `--minimizer-spaces` setting (7), which has higher sensitivity and virtually identical accuracy than the 6 that was set in prior version (e.g. see Fig. 1E [here](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1891-0)).
* In [Step 3](#3-generate-a-kraken2-summary-report), changed `total_fragments=$(wc -l sample-1-kraken2-output.txt | cut -f 1 -d " ")` to `total_fragments=$(wc -l sample-1-kraken2-output.txt | sed 's/^ *//' | cut -f 1 -d " ")` for better portability across Unix-like systems.

---

# Table of contents

- [**Software used**](#software-used)
- [**General processing overview with example commands**](#general-processing-overview-with-example-commands)
  - [**1. Build kraken2 database**](#1-build-kraken2-database)
  - [**2. Filter out human-classified reads**](#2-filter-out-human-classified-reads)
    - [Example if paired-end reads](#example-if-paired-end-reads)
    - [Example if single-end reads](#example-if-single-end-reads)
  - [**3. Generate a kraken2 summary report**](#3-generate-a-kraken2-summary-report)

---

# Software used

|Program|Version*|Relevant Links|
|:------|:-----:|------:|
|kraken2|`kraken2 -v`|[https://github.com/DerrickWood/kraken2/wiki](https://github.com/DerrickWood/kraken2/wiki)|

> \* Exact versions utilized for a given dataset are available along with the processing commands for each specific dataset (this is due to how the system may need to be updated regularly).

---

# General processing overview with example commands

## 1. Build kraken2 database

```bash
kraken2-build --download-library human --db kraken2-human-db --threads 30 --no-masking
kraken2-build --download-taxonomy --db kraken2-human-db/
kraken2-build --build --db kraken2-human-db/ --threads 30
kraken2-build --clean --db kraken2-human-db/
```

**Parameter Definitions:**

* `--download-library` - specifies the references to download (here just the human reference genome)
* `--db` - specifies the directory we are putting the database in
* `--threads` - specifies the number of threads to use
* `--no-masking` - prevents [masking](https://github.com/DerrickWood/kraken2/wiki/Manual#masking-of-low-complexity-sequences) of low-complexity sequences
* `--download-taxonomy` - downloads taxonomic mapping information
* `--build` - specifies to construct kraken2-formatted database
* `--clean` - specifies to remove unnecessarily intermediate files

**Input data:**

* None

**Output data:**

* kraken2 database files (hash.k2d, opts.k2d, and taxo.k2d)

---

## 2. Filter out human-classified reads

### Example if paired-end reads

```bash
kraken2 --db kraken2-human-db --gzip-compressed --threads 4 --use-names --paired \
        --output sample-1-kraken2-output.txt --report sample-1-kraken2-report.tsv \
        --unclassified-out sample-1_R#.fastq sample-1-R1.fq.gz sample-1-R2.fq.gz
        
# renaming and gzipping output files
mv sample-1_R_1.fastq sample-1-R1-human-reads-removed.fastq && gzip sample-1-R1-human-reads-removed.fastq
mv sample-1_R_2.fastq sample-1-R2-human-reads-removed.fastq && gzip sample-1-R2-human-reads-removed.fastq
```

**Parameter Definitions:**

* `--db` - specifies the directory holding the kraken2 database files created in step 1
* `--gzip-compressed` - specifies the input fastq files are gzip-compressed
* `--threads` - specifies the number of threads to use
* `--use-names` - specifies adding taxa names in addition to taxids
* `--paired` - specifies input reads are paired-end
* `--output` - specifies the name of the kraken2 read-based output file (one line per read)
* `--report` - specifies the name of the kraken2 report output file (one line per taxa, with number of reads assigned to it)
* `--unclassified-out` - name of output files of reads that were not classified (the `#` symbol gets replaced with "_1" and "_2" in the output file names)
* last two positional arguments are the input read files

**Input data:**

* sample-1-R1.fq.gz (gzipped forward-reads fastq file)
* sample-1-R2.fq.gz (gzipped reverse-reads fastq file)

**Output data:**

* sample-1-kraken2-output.txt (kraken2 read-based output file (one line per read))
* sample-1-kraken2-report.tsv (kraken2 report output file (one line per taxa, with number of reads assigned to it))
* sample-1-R1-human-reads-removed.fastq.gz (human-read removed, gzipped forward-reads fastq file)
* sample-1-R2-human-reads-removed.fastq.gz (human-read removed, gzipped reverse-reads fastq file)

### Example if single-end reads

```bash
kraken2 --db kraken2-human-db --gzip-compressed --threads 4 --use-names \
        --output sample-1-kraken2-output.txt --report sample-1-kraken2-report.tsv \
        --unclassified-out sample-1-human-reads-removed.fastq sample-1.fq.gz

# gzipping output file
gzip sample-1-human-reads-removed.fastq
```

**Parameter Definitions:**

* `--db` - specifies the directory holding the kraken2 database files created in step 1
* `--gzip-compressed` - specifies the input fastq files are gzip-compressed
* `--threads` - specifies the number of threads to use
* `--use-names` - specifies adding taxa names in addition to taxids
* `--output` - specifies the name of the kraken2 read-based output file (one line per read)
* `--report` - specifies the name of the kraken2 report output file (one line per taxa, with number of reads assigned to it)
* `--unclassified-out` - name of output files of reads that were not classified 
* last positional argument is the input read file

**Input data:**

* sample-1.fq.gz (gzipped reads fastq file)

**Output data:**

* sample-1-kraken2-output.txt (kraken2 read-based output file (one line per read))
* sample-1-kraken2-report.tsv (kraken2 report output file (one line per taxa, with number of reads assigned to it))
* sample-1-human-reads-removed.fastq.gz (human-read removed, gzipped reads fastq file)

---

## 3. Generate a kraken2 summary report
Utilizes a Unix-like command-line.

```bash
total_fragments=$(wc -l sample-1-kraken2-output.txt | sed 's/^ *//' | cut -f 1 -d " ")

fragments_retained=$(grep -w -m 1 "unclassified" sample-1-kraken2-report.tsv | cut -f 2)

perc_removed=$(printf "%.2f\n" $(echo "scale=4; 100 - ${fragments_retained} / ${total_fragments} * 100" | bc -l))

cat <( printf "Sample_ID\tTotal_fragments_before\tTotal_fragments_after\tPercent_human_reads_removed\n" ) \
    <( printf "Sample-1\t${total_fragments}\t${fragments_retained}\t${perc_removed}\n" ) > Human-read-removal-summary.tsv
```

**Input data:**

* sample-1-kraken2-output.txt (kraken2 read-based output file (one line per read))
* sample-1-kraken2-report.tsv (kraken2 report output file (one line per taxa, with number of reads assigned to it))

**Output data:**

* Human-read-removal-summary.tsv (a tab-separated file with 4 columns: "Sample_ID", "Total_fragments_before", "Total_fragments_after", "Percent_human_reads_removed")

---
