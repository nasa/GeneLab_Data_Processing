# GeneLab removal of human reads from metagenomics datasets

> **It is NASA's policy that any human reads are to be removed from metagenomics datasets prior to being hosted in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/). As such, all metagenomics datasets are screened against a human reference-genome [kraken2](https://github.com/DerrickWood/kraken2/wiki) database. This document holds an overview and some example commands of how GeneLab does perform this.**

---

 **Date:**  March X, 2026  
**Revision:** B  
**Document Number:** GL-DPPD-7105-B  

**Submitted by:**  
Jihan Yehia (GeneLab Data Processing Team)  

**Approved by:**  
Jonathan Galazka (OSDR Project Manager)  
Danielle Lopez (OSDR Deputy Project Manager)  
Amanda Saravia-Butler (OSDR Subject Matter Expert)  
Barbara Novak (GeneLab Data Processing Lead)  

## Updates from previous revision
* Updated kraken2 from version 2.1.1 to 2.1.6
* In [Step 1](#1-build-kraken2-database), replaced direct `kraken2-build` calls with kraken2's `k2` wrapper script, which provides a higher-level interface for database construction including HTTPS download support, as described [here](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#introducing-k2)

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
|kraken2| 2.1.6 |[https://github.com/DerrickWood/kraken2/wiki](https://github.com/DerrickWood/kraken2/wiki)|

> \* Exact version utilized for a given dataset are available along with the processing commands for each specific dataset (this is due to how the system may need to be updated regularly).

---

# General processing overview with example commands

> Output files listed in **bold** below are included with each Metagenomics dataset in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/). 
> **Note:** the examples provided are for human read removal. For the more generic host read removal, please refer to the Workflow documentation.


### 1. Build kraken2 database
For human read removal, building the database relies on kraken2 downloading the sequences from NCBI, 
which may include multiple different versions of the genome (for example, both [GRCh38.p14](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.40) 
and [T2T-CHM13v2.0](https://www.ncbi.nlm.nih.gov/assembly/GCF_009914755.1)). To use this pipeline for 
read removal of hosts other than human, download the host fasta file and add it to the kraken2 database 
using the [`k2 add-to-library` function as described in the [kraken2 documentation](https://github.com/DerrickWood/kraken2/wiki/Manual#add-to-library) or 
in the NASA GeneLab [Estimate Host Reads pipeline](../../Estimate_host_reads_in_raw_data/Pipeline_GL-DPPD-7109_Versions/)

> **Note:** It is recommended to use NCBI genome files with kraken2 because sequences not downloaded from 
NCBI may require explicit assignment of taxonomy information before they can be used to build the 
database, as mentioned in the [Kraken2 Documentation](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown). 

```bash
# download human fasta sequences
k2 download-library --db kraken2-human-db/ --library human --threads 30 --no-masking

# Download NCBI taxonomic information 
k2 download-taxonomy --db kraken2-human-db/

# Build the database
k2 build --db kraken2-human-db/ --kmer-len 35 --minimizer-len 31 --threads 30

# Clean up intermediate files
k2 clean --db kraken2-human-db/
```

**Parameter Definitions:**

- `download-library` - downloads the references-containing library
  - `--library` - specifies the library to download (here the human reference genome)
  - `--no-masking` - disables masking of low-complexity sequences. For additional 
                   information see the [kraken documentation for masking](https://github.com/DerrickWood/kraken2/wiki/Manual#masking-of-low-complexity-sequences).
- `--db` - specifies the name of the directory for the kraken2 database
- `--threads` - specifies the number of threads to use
- `download-taxonomy` - downloads taxonomic mapping information
- `build` - builds a kraken2-formatted database from the library files
  - `--kmer-len` - k-mer length in bp (default: 35).
  - `--minimizer-len` - minimizer length in bp (default: 31)
- `clean` - removes unneeded intermediate files

**Input data:**

* None

**Output data:**

- kraken2_human_db/ (Kraken2 human database directory, containing hash.k2d, opts.k2d, and taxo.k2d files)

---

## 2. Filter out human-classified reads

>**Note:** HRrm is the file suffix NASA GeneLab uses for all files with human reads removed. If host reads are removed, the suffix is updated to "HostRm" instead.

### Example if paired-end reads

```bash
kraken2 --db kraken2-human-db --gzip-compressed --threads 4 --use-names --paired \
        --output sample-1-kraken2-output.txt --report sample-1-kraken2-report.tsv \
        --unclassified-out sample-1_R#.fastq sample-1-R1.fq.gz sample-1-R2.fq.gz
        
# renaming and gzipping output files
mv sample-1_R_1.fastq sample-1_R1_HRrm.fastq && gzip sample-1_R1_HRrm.fastq
mv sample-1_R_2.fastq sample-1_R2_HRrm.fastq && gzip sample-1_R2_HRrm.fastq
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

* kraken2-human-db (Kraken2 human database directory, from [Step 1](#1-build-kraken2-database))
* sample-1-R1.fq.gz (gzipped forward-reads fastq file)
* sample-1-R2.fq.gz (gzipped reverse-reads fastq file)

**Output data:**

* sample-1-kraken2-output.txt (kraken2 read-based output file (one line per read))
* sample-1-kraken2-report.tsv (kraken2 report output file (one line per taxa, with number of reads assigned to it))
* **sample-1_R1_HRrm.fastq.gz** (human-read removed, gzipped forward-reads fastq file)
* **sample-1_R2_HRrm.fastq.gz** (human-read removed, gzipped reverse-reads fastq file)

### Example if single-end reads

```bash
kraken2 --db kraken2-human-db --gzip-compressed --threads 4 --use-names \
        --output sample-1-kraken2-output.txt --report sample-1-kraken2-report.tsv \
        --unclassified-out sample-1_HRrm.fastq sample-1.fq.gz

# gzipping output file
gzip sample-1_HRrm.fastq
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
* **sample-1_HRrm.fastq.gz** (human-read removed, gzipped reads fastq file)

---

## 3. Generate a kraken2 summary report
Utilizes a Unix-like command-line.

>**Note:** In the case of host removal rather than human removal, Update the script below to change "human" to "host". 

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
* *Note: The percent human reads removed from each sample is provided in the assay table on the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/).*

---
