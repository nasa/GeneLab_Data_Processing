# GeneLab tracking of host reads in Illumina metagenomic datasets

> **In order to provide an estimate of host DNA in Illumina metagenomic datasets that are sequenced from host-derived samples, datasets are screened against an appropriate reference genome using [kraken2](https://github.com/DerrickWood/kraken2/wiki). Reads are not removed from the dataset, but the percentage of detected host reads is reported.**

---

**Date:**  February, 4, 2022  
**Revision:** -  
**Document Number:** GL-DPPD-7109  

**Submitted by:**  
Michael D. Lee (GeneLab Science Team)  

**Approved by:**  
Sylvain Costes (GeneLab Project Manager)  
Samrawit Gebre (GeneLab Deputy Project Manager and Interim GeneLab Configuration Manager)  
Amanda Saravia-Butler (GeneLab Data Processing Lead)  
Jonathan Galazka (GeneLab Project Scientist)  

---

# Table of contents

- [**Software used**](#software-used)
- [**General processing overview with example commands**](#general-processing-overview-with-example-commands)
  - [**1. Build kraken2 database**](#1-build-kraken2-database-of-host-genome)
  - [**2. Identify host-classified reads**](#2-identify-host-classified-reads)
    - [Example if paired-end reads](#example-if-paired-end-reads)
    - [Example if single-end reads](#example-if-single-end-reads)
  - [**3. Generate a summary report**](#3-generate-a-summary-report)

---

# Software used

|Program|Version|Relevant Links|
|:------|:-----:|------:|
|kraken2|2.1.1|[https://github.com/DerrickWood/kraken2/wiki](https://github.com/DerrickWood/kraken2/wiki)|

---

# General processing overview with example commands

## 1. Build kraken2 database of host genome
This depends on the appropriate host genome. This example is done with the mouse genome ([GRCm39 | GCF_000001635.27](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.27)).

```bash
# downloading and decompressing reference genome
curl -LO https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz
gunzip GCF_000001635.27_GRCm39_genomic.fna.gz

# building kraken2 database
kraken2-build --download-taxonomy --db kraken2-mouse-db/
kraken2-build --add-to-library GCF_000001635.27_GRCm39_genomic.fna --no-masking --db kraken2-mouse-db/
kraken2-build --build --db kraken2-mouse-db/ --threads 30 --no-masking
kraken2-build --clean --db kraken2-mouse-db/
```

**Parameter Definitions:**

* `--download-taxonomy` - downloads taxonomic mapping information
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
* reference genome used (*.fna)

---

## 2. Identify host-classified reads

### Example if paired-end reads

```bash
kraken2 --db kraken2-mouse-db --gzip-compressed --threads 4 --use-names --paired \
        --output sample-1-kraken2-output.txt --report sample-1-kraken2-report.tsv Sample-1_R1.fastq.gz Sample-1_R2.fastq.gz
```

**Parameter Definitions:**

* `--db` - specifies the directory holding the kraken2 database files created in step 1
* `--gzip-compressed` - specifies the input fastq files are gzip-compressed
* `--threads` - specifies the number of threads to use
* `--use-names` - specifies adding taxa names in addition to taxids
* `--paired` - specifies input reads are paired-end
* `--output` - specifies the name of the kraken2 read-based output file (one line per read)
* `--report` - specifies the name of the kraken2 report output file (one line per taxa, with number of reads assigned to it)
* last two positional arguments are the input read files

**Input data:**

* Sample-1_R1.fastq.gz (gzipped forward-reads fastq file)
* Sample-1_R2.fastq.gz (gzipped reverse-reads fastq file)

**Output data:**

* sample-1-kraken2-output.txt (kraken2 read-based output file (one line per read))
* sample-1-kraken2-report.tsv (kraken2 report output file (one line per taxa, with number of reads assigned to it))

### Example if single-end reads

```bash
kraken2 --db kraken2-mouse-db --gzip-compressed --threads 4 --use-names \
        --output sample-1-kraken2-output.txt --report sample-1-kraken2-report.tsv Sample-1.fastq.gz
```

**Parameter Definitions:**

* `--db` - specifies the directory holding the kraken2 database files created in step 1
* `--gzip-compressed` - specifies the input fastq files are gzip-compressed
* `--threads` - specifies the number of threads to use
* `--use-names` - specifies adding taxa names in addition to taxids
* `--output` - specifies the name of the kraken2 read-based output file (one line per read)
* `--report` - specifies the name of the kraken2 report output file (one line per taxa, with number of reads assigned to it)
* last positional argument is the input read file

**Input data:**

* Sample-1.fastq.gz (gzipped reads fastq file)

**Output data:**

* sample-1-kraken2-output.txt (kraken2 read-based output file (one line per read))
* sample-1-kraken2-report.tsv (kraken2 report output file (one line per taxa, with number of reads assigned to it))

---

## 3. Generate a summary report
Utilizes a Unix-like command-line.

```bash
total_fragments=$(wc -l sample-1-kraken2-output.txt | sed 's/^ *//' | cut -f 1 -d " ")

fragments_classified=$(grep -w -c "^C" sample-1-kraken2-output.txt)

perc_host=$(printf "%.2f\n" $(echo "scale=4; ${fragments_classified} / ${total_fragments} * 100" | bc -l))

cat <( printf "Sample_ID\tTotal_fragments\tTotal_host_fragments\tPercent_host\n" ) \
    <( printf "Sample-1\t${total_fragments}\t${fragments_classified}\t${perc_host}\n" ) > Host-read-count-summary.tsv
```

**Input data:**

* sample-1-kraken2-output.txt (kraken2 read-based output file (one line per read))
* sample-1-kraken2-report.tsv (kraken2 report output file (one line per taxa, with number of reads assigned to it))

**Output data:**

* Host-read-count-summary.tsv (a tab-separated file with 4 columns: "Sample\_ID", "Total\_fragments", "Total\_host\_fragments", "Percent\_host")

---
