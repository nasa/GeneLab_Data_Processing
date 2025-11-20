# GeneLab removal of human reads from metagenomics datasets

> **It is NASA's policy that any human reads are to be removed from metagenomics datasets prior to being hosted in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/). As such, all metagenomics datasets are screened against a human reference-genome [kraken2](https://github.com/DerrickWood/kraken2/wiki) database. This document holds an overview and some example commands of how GeneLab does perform this.**

---

**Date:**  November X, 2025  
**Revision:** B 
**Document Number:** GL-DPPD-7105-B  

**Submitted by:**  
Jihan Yehia (GeneLab Data Processing Team)  

**Approved by:**  
Samrawit Gebre (OSDR Project Manager)  
Danielle Lopez (OSDR Deputy Project Manager)  
Jonathan Galazka (OSDR Project Scientist)  
Amanda Saravia-Butler (GeneLab Science Lead)  
Barbara Novak (GeneLab Data Processing Lead)  

## Updates from previous revision
* Updated kraken2 from version 2.1.1 to 2.1.6
* In [Step 1](#1-build-kraken2-database), used kraken2's `k2` wrapper script for `download-taxonomy` because the script supports HTTPS download as mentioned [here](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#introducing-k2)

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

> Output files listed in **bold** below are included with each Metagenomics dataset in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/).

## 1. Build kraken2 database
This depends on the appropriate host genome. This example is done with the human genome ([GRCh38.p13 | GCF_000001405.39](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39)).
> **Note:** It is recommended to use NCBI with kraken2 because sequences not downloaded from NCBI may require explicit assignment of taxonomy information before they can be used to build the database, as mentioned [here](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown).

```bash
# downloading and decompressing reference genome
wget -q https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz
gunzip GCF_000001405.39_GRCh38.p13_genomic.fna.gz


# building kraken2 database
k2 download-taxonomy --db kraken2-human-db/
kraken2-build --add-to-library GCF_000001405.39_GRCh38.p13_genomic.fna --no-masking --db kraken2-human-db/
kraken2-build --build --db kraken2-human-db/ --threads 30 --no-masking
kraken2-build --clean --db kraken2-human-db/
```

**Parameter Definitions:**

* `download-taxonomy` - downloads taxonomic mapping information via [k2 wrapper script](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#introducing-k2)
* `--add-to-library` - adds the fasta file to the library of sequences being included
* `--db` - specifies the directory we are putting the database in
* `--threads` - specifies the number of threads to use
* `--no-masking` - prevents [masking](https://github.com/DerrickWood/kraken2/wiki/Manual#masking-of-low-complexity-sequences) of low-complexity sequences
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
mv sample-1_R_1.fastq sample-1_R1_HRremoved_raw.fastq && gzip sample-1_R1_HRremoved_raw.fastq
mv sample-1_R_2.fastq sample-1_R2_HRremoved_raw.fastq && gzip sample-1_R2_HRremoved_raw.fastq
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
* **sample-1_R1_HRremoved_raw.fastq.gz** (host-read removed, gzipped forward-reads fastq file)
* **sample-1_R2_HRremoved_raw.fastq.gz** (host-read removed, gzipped reverse-reads fastq file)

### Example if single-end reads

```bash
kraken2 --db kraken2-human-db --gzip-compressed --threads 4 --use-names \
        --output sample-1-kraken2-output.txt --report sample-1-kraken2-report.tsv \
        --unclassified-out sample-1_HRremoved_raw.fastq sample-1.fq.gz

# gzipping output file
gzip sample-1_HRremoved_raw.fastq
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
* **sample-1_HRremoved_raw.fastq.gz** (host-read removed, gzipped reads fastq file)

---

## 3. Generate a kraken2 summary report
Utilizes a Unix-like command-line.

```bash
total_fragments=$(wc -l sample-1-kraken2-output.txt | sed 's/^ *//' | cut -f 1 -d " ")

fragments_retained=$(grep -w -m 1 "unclassified" sample-1-kraken2-report.tsv | cut -f 2)

perc_removed=$(printf "%.2f\n" $(echo "scale=4; 100 - ${fragments_retained} / ${total_fragments} * 100" | bc -l))

cat <( printf "Sample_ID\tTotal_fragments_before\tTotal_fragments_after\tPercent_host_reads_removed\n" ) \
    <( printf "Sample-1\t${total_fragments}\t${fragments_retained}\t${perc_removed}\n" ) > Human-read-removal-summary.tsv
```

**Input data:**

* sample-1-kraken2-output.txt (kraken2 read-based output file (one line per read))
* sample-1-kraken2-report.tsv (kraken2 report output file (one line per taxa, with number of reads assigned to it))

**Output data:**

* <Human>-read-removal-summary.tsv (a tab-separated file with 4 columns: "Sample_ID", "Total_fragments_before", "Total_fragments_after", "Percent_host_reads_removed")
* *Note: The percent human reads removed from each sample is provided in the assay table on the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/).*

---
