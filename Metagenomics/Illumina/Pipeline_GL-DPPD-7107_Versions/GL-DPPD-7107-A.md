# Bioinformatics pipeline for Illumina metagenomics data

> **This document holds an overview and some example commands of how GeneLab processes Illumina metagenomics datasets. Exact processing commands for specific datasets that have been released are provided with their processed data in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/).**  

---

**Date:** October 28, 2024  
**Revision:** -A  
**Document Number:** GL-DPPD-7107  

**Submitted by:**  
Olabiyi A. Obayomi (GeneLab Analysis Team)  

**Approved by:**  
Samrawit Gebre (OSDR Project Manager)  
Lauren Sanders (OSDR Project Scientist)  
Amanda Saravia-Butler (GeneLab Science Lead)  
Barbara Novak (GeneLab Data Processing Lead)  

---

## Updates from previous version  

- The following tool versions were updated:
  - FastQC
  - MultiQC
  - bowtie2
  - samtools
  - CAT
  - GTDB-Tk
  - HUMAnN
  - MetaPhlAn
- In [step 14d](#14d-mag-taxonomic-classification), MAG taxonomic classification, added the new `--skip_ani_screen` argument to `gtdbtk classify_wf` to continue classifying genomes as in previous versions of GTDB-Tk, using mash and skani.

---

# Table of contents

- [**Software used**](#software-used)
- [**General processing overview with example commands**](#general-processing-overview-with-example-commands)
  - [**Pre-processing**](#pre-processing)
    - [1. Raw Data QC](#1-raw-data-qc)
    - [2. Quality filtering/trimming](#2-quality-filteringtrimming)
    - [3. Filtered/Trimmed Data QC](#3-filteredtrimmed-data-qc)
  - [**Assembly-based processing**](#assembly-based-processing)
    - [4. Sample assembly](#4-sample-assembly)
    - [5. Renaming contigs and summarizing assemblies](#5-renaming-contigs-and-summarizing-assemblies)
    - [6. Gene prediction](#6-gene-prediction)
    - [7. Functional annotation](#7-functional-annotation)
    - [8. Taxonomic classification](#8-taxonomic-classification)
    - [9. Read-mapping](#9-read-mapping)
    - [10. Getting coverage information and filtering based on detection](#10-getting-coverage-information-and-filtering-based-on-detection)
    - [11. Combining gene-level coverage, taxonomy, and functional annotations into one table for each sample](#11-combining-gene-level-coverage-taxonomy-and-functional-annotations-into-one-table-for-each-sample)
    - [12. Combining contig-level coverage and taxonomy into one table for each sample](#12-combining-contig-level-coverage-and-taxonomy-into-one-table-for-each-sample)
    - [13. Generating normalized, gene- and contig-level coverage summary tables of KO-annotations and taxonomy across samples](#13-generating-normalized-gene--and-contig-level-coverage-summary-tables-of-ko-annotations-and-taxonomy-across-samples)
    - [14. **M**etagenome-**A**ssembled **G**enome (MAG) recovery](#14-metagenome-assembled-genome-mag-recovery)
    - [15. Generating MAG-level functional summary overview](#15-generating-mag-level-functional-summary-overview)
  - [**Read-based processing**](#read-based-processing)
    - [16. Taxonomic and functional profiling](#16-taxonomic-and-functional-profiling)

---

# Software used

|Program|Version|Relevant Links|
|:------|:-----:|------:|
|FastQC| 0.12.1 |[https://www.bioinformatics.babraham.ac.uk/projects/fastqc/](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)|
|MultiQC| 1.19 |[https://multiqc.info/](https://multiqc.info/)|
|bbduk| 38.86 |[https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/)|
|MEGAHIT| 1.2.9 |[https://github.com/voutcn/megahit#megahit](https://github.com/voutcn/megahit#megahit)|
|bit| 1.8.53 |[https://github.com/AstrobioMike/bioinf_tools#bioinformatics-tools-bit](https://github.com/AstrobioMike/bioinf_tools#bioinformatics-tools-bit)|
|bowtie2| 2.4.1 |[https://github.com/BenLangmead/bowtie2#overview](https://github.com/BenLangmead/bowtie2#overview)|
|samtools| 1.20 |[https://github.com/samtools/samtools#samtools](https://github.com/samtools/samtools#samtools)|
|Prodigal| 2.6.3 |[https://github.com/hyattpd/Prodigal#prodigal](https://github.com/hyattpd/Prodigal#prodigal)|
|KOFamScan| 1.3.0 |[https://github.com/takaram/kofam_scan#kofamscan](https://github.com/takaram/kofam_scan#kofamscan)|
|CAT| 5.2.3 |[https://github.com/dutilh/CAT#cat-and-bat](https://github.com/dutilh/CAT#cat-and-bat)|
|MetaBAT| 2.15 |[https://bitbucket.org/berkeleylab/metabat/src/master/](https://bitbucket.org/berkeleylab/metabat/src/master/)|
|CheckM| 1.1.3 |[https://github.com/Ecogenomics/CheckM](https://github.com/Ecogenomics/CheckM)|
|GTDB-Tk| 2.4.0 |[https://github.com/Ecogenomics/GTDBTk](https://github.com/Ecogenomics/GTDBTk)|
|KEGG-Decoder| 1.2.2 |[https://github.com/bjtully/BioData/tree/master/KEGGDecoder#kegg-decoder](https://github.com/bjtully/BioData/tree/master/KEGGDecoder#kegg-decoder)
|HUMAnN| 3.9 |[https://github.com/biobakery/humann](https://github.com/biobakery/humann)|
|MetaPhlAn| 4.1.0 |[https://github.com/biobakery/MetaPhlAn](https://github.com/biobakery/MetaPhlAn)|

---

# General processing overview with example commands

> Exact processing commands and output files listed in **bold** below are included with each Metagenomics Seq processed dataset in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/).  

## Pre-processing
### 1. Raw Data QC

```
fastqc -o raw_fastqc_output *raw.fastq.gz
```

**Parameter Definitions:**

* `-o` – the output directory to store results
* `*raw.fastq.gz` – the input reads are specified as a positional argument, and can be given all at once with wildcards like this, or as individual arguments with spaces in between them

**Input data:**

* *raw.fastq.gz (raw reads, after human read removal)

**Output data:**

* *fastqc.html (FastQC output html summary)
* *fastqc.zip (FastQC output data)


#### 1a. Compile Raw Data QC

```
multiqc -o raw_multiqc_output -n raw_multiqc raw_fastqc_output/
# this is how it's packaged with our workflow outputs
zip -r raw_multiqc_GLmetagenomics_report.zip raw_multiqc_output
```

**Parameter Definitions:**

*	`-o` – the output directory to store results
*	`-n` – the filename prefix of results
*	`raw_fastqc_output/` – the directory holding the output data from the fastqc run, provided as a positional argument

**Input data:**

* raw_fastqc_output/*fastqc.zip (FastQC output data)

**Output data:**

* **raw_multiqc_GLmetagenomics_report.zip** (zip containing the following)
  * **raw_multiqc.html** (multiqc output html summary)
  * **raw_multiqc_data** (directory containing multiqc output data)

<br>  

---

### 2. Quality filtering/trimming

```      
bbduk.sh in=sample-1-R1-raw.fastq.gz in2=sample-1-R2-raw.fastq.gz out1=sample-1_R1_filtered.fastq.gz \
         out2=sample-1_R2_filtered.fastq.gz ref=ref-adapters.fa ktrim=l k=17 ftm=5 qtrim=rl \
         trimq=10 mlf=0.5 maxns=0 > bbduk.log 2>&1
         
# if libraries were prepared with the Swift1S kit
# bbduk.sh in=sample-1-R1-raw.fastq.gz in2=sample-1-R2-raw.fastq.gz out1=sample-1_R1_filtered.fastq.gz \
         out2=sample-1_R2_filtered.fastq.gz ref=ref-adapters.fa ktrim=l k=17 ftm=5 qtrim=rl \
         trimq=10 mlf=0.5 maxns=0 swift=t > bbduk.log 2>&1

```

**Parameter Definitions:**

*	`in` and `in2` – specifies the forward and reverse input reads, respectively (no `in2` if working with single-end data)

*	`out1` and `out2` – specifies the forward and reverse output reads, respectively (no `out2` if working with single-end data)

*	`ref` – specifies a fasta file holding potential adapter sequences (comes with bbduk installation)

*	`ktrim` – specifies to trim adapters from the 5’ end (left) if found

*	`k` – sets minimum length of kmer match to identify adapter sequences (provided by the “ref” file above)

*	`ftm` – sets a multiple of expected length the sequence should be (handles poor additional bases that are sometimes present, see “Force-Trim Modulo” section on [this page](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/))

*	`qtrim` – sets quality-score-based trimming to be applied to left and right sides

*	`trimq` – sets the score to use for PHRED-algorithm trimming

*	`mlf` – sets the minimum length of reads retained based on their initial length

*	`maxns` – sets the maximum number of Ns allowed in a read before it will be filtered out

*	`swift` – tells the program to look for and trim low-complexity adaptase reminants from the Swift1S kit

*	`> bbduk.log 2>&1` – redirects the stderr and stdout to a log file for saving

**Input data:**

* *raw.fastq.gz (raw reads)

**Output data:**

* **\*_filtered.fastq.gz** (filtered/trimmed reads)
* **\*-bbduk.log** (log file of standard output and error from bbduk run)

<br>

---

### 3. Filtered/Trimmed Data QC
```
fastqc -o filtered_fastqc_output/ *filtered.fastq.gz
```

**Parameter Definitions:**

*	`-o` – the output directory to store results  
*	`*filtered.fastq.gz` – the input reads are specified as a positional argument, and can be given all at once with wildcards like this, or as individual arguments with spaces in between them  

**Input data:**

* *filtered.fastq.gz (filtered/trimmed reads)

**Output data:**

* *fastqc.html (FastQC output html summary)
* *fastqc.zip (FastQC output data)


#### 3a. Compile Filtered/Trimmed Data QC
```
multiqc -o filtered_multiqc_output -n filtered_multiqc filtered_fastqc_output/
# this is how it's packaged with our workflow outputs
zip -r filtered_multiqc_GLmetagenomics_report.zip filtered_multiqc_output
```

**Parameter Definitions:**

*	`-o` – the output directory to store results
*	`-n` – the filename prefix of results
*	`filtered_fastqc_output/` – the directory holding the output data from the fastqc run, provided as a positional argument

**Input data:**

* filtered_fastqc_output/*fastqc.zip (FastQC output data)

**Output data:**

* **filtered_multiqc_GLmetagenomics_report.zip** (zip containing the following)
  * **filtered_multiqc.html** (multiqc output html summary)
  * **filtered_multiqc_data** (directory containing multiqc output data)

<br>

---

## Assembly-based processing
### 4. Sample assembly
```
megahit -1 sample-1_R1_filtered.fastq.gz -2 sample-1_R2_filtered.fastq.gz \
        -o sample-1-assembly -t NumberOfThreads --min-contig-length 500 > sample-1-assembly.log 2>&1
```

**Parameter Definitions:**  

*	`-1 and -2` – specifies the input forward and reverse reads (if single-end data, then neither `-1` nor `-2` are used, instead single-end reads are passed to `-r`)

*	`-o` – specifies output directory

*	`-t` – specifies the number of threads to use

*	`--min-contig-length` – specifies the minimum contig length to write out

*	`> sample-1-assembly.log 2>&1` – sends stdout/stderr to log file


**Input data:**

* *fastq.gz (filtered/trimmed reads from [step 2](#2-quality-filteringtrimming) above)

**Output data:**

* sample-1-assembly/final.contigs.fa (assembly file)
* **sample-1-assembly.log** (log file)

<br>

---

### 5. Renaming contigs and summarizing assemblies

#### 5a. Renaming contig headers
```
bit-rename-fasta-headers -i sample-1-assembly/final.contigs.fa -w c_sample-1 -o sample-1-assembly.fasta
```

**Parameter Definitions:**  

*	`-i` – input fasta file

*	`-w` – wanted header prefix (a number will be appended for each contig), starts with a “c_” to ensure they won’t start with a number which can be problematic

*	`-o` – output fasta file


**Input data:**

* sample-1-assembly/final.contigs.fa (assembly file from [step 4](#4-sample-assembly))

**Output files:**

* **sample-1-assembly.fasta** (contig-renamed assembly file)


#### 5b. Summarizing assemblies

```
bit-summarize-assembly -o assembly-summaries_GLmetagenomics.tsv *assembly.fasta
```

**Parameter Definitions:**  

*	`-o` – output summary table

*	– multiple input assemblies can be provided as positional arguments


**Input data:**

* *-assembly.fasta (contig-renamed assembly files from [step 5a](#5a-renaming-contig-headers))

**Output files:**

* **assembly-summaries_GLmetagenomics.tsv** (table of assembly summary statistics)

<br>

---

### 6. Gene prediction
```
prodigal -a sample-1-genes.faa -d sample-1-genes.fasta -f gff -p meta -c -q \
         -o sample-1-genes.gff -i sample-1-assembly.fasta
```
**Parameter Definitions:**

*	`-a` – specifies the output amino acid sequences file

*	`-d` – specifies the output nucleotide sequences file

*	`-f` – specifies the output format gene-calls file

*	`-p` – specifies which mode to run the gene-caller in 

*	`-c` – no incomplete genes reported 

*	`-q` – run in quiet mode (don’t output process on each contig) 

*	`-o` – specifies the name of the output gene-calls file 

*	`-i` – specifies the input assembly

**Input data:**

* sample-1-assembly.fasta (contig-renamed assembly file from [step 5a](#5a-renaming-contig-headers))

**Output data:**

* **sample-1-genes.faa** (gene-calls amino-acid fasta file)
* **sample-1-genes.fasta** (gene-calls nucleotide fasta file)
* **sample-1-genes.gff** (gene-calls in general feature format)

<br>

---

### 7. Functional annotation
> **Notes**  
> The annotation process overwrites the same temporary directory by default. So if running multiple processses at a time, it is necessary to specify a specific temporary directory with the `--tmp-dir` argument as shown below.


#### 7a. Downloading reference database of HMM models (only needs to be done once)

```
curl -LO ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz
curl -LO ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz
tar -xzvf profiles.tar.gz
gunzip ko_list.gz 
```

#### 7b. Running KEGG annotation
```
exec_annotation -p profiles/ -k ko_list --cpu 15 -f detail-tsv -o sample-1-KO-tab.tmp \
                --tmp-dir sample-1-tmp-KO --report-unannotated sample-1-genes.faa 
```

**Parameter Definitions:**
*	`-p` – specifies the directory holding the downloaded reference HMMs

*	`-k` – specifies the downloaded reference KO  (Kegg Orthology) terms 

*	`--cpu` – specifies the number of searches to run in parallel

*	`-f` – specifies the output format

*	`-o` – specifies the output file name

*	`--tmp-dir` – specifies the temporary directory to write to (needed if running more than one process concurrently, see Notes above)

*	`--report-unannotated` – specifies to generate an output for each entry

*	`sample-1-genes.faa` – the input file is specified as a positional argument 


**Input data:**

* sample-1-genes.faa (amino-acid fasta file, from [step 6](#6-gene-prediction))
* profiles/ (reference directory holding the KO HMMs)
* ko_list (reference list of KOs to scan for)

**Output data:**

* sample-1-KO-tab.tmp (table of KO annotations assigned to gene IDs)


#### 7c. Filtering output to retain only those passing the KO-specific score and top hits
```
bit-filter-KOFamScan-results -i sample-1-KO-tab.tmp -o sample-1-annotations.tsv

  # removing temporary files
rm -rf sample-1-tmp-KO/ sample-1-KO-annots.tmp
```

**Parameter Definitions:**  

*	`-i` – specifies the input table

*	`-o` – specifies the output table


**Input data:**

* sample-1-KO-tab.tmp (table of KO annotations assigned to gene IDs from [step 7b](#7b-running-kegg-annotation))

**Output data:**

* sample-1-annotations.tsv (table of KO annotations assigned to gene IDs)

<br>

---

### 8. Taxonomic classification

#### 8a. Pulling and un-packing pre-built reference db (only needs to be done once)
```
wget tbb.bio.uu.nl/bastiaan/CAT_prepare/CAT_prepare_20200618.tar.gz
tar -xvzf CAT_prepare_20200618.tar.gz
```

#### 8b. Running taxonomic classification
```
CAT contigs -c sample-1-assembly.fasta -d CAT_prepare_20200618/2020-06-18_database/ \
            -t CAT_prepare_20200618/2020-06-18_taxonomy/ -p sample-1-genes.faa \
            -o sample-1-tax-out.tmp -n NumberOfThreads -r 3 --top 4 --I_know_what_Im_doing --no-stars
```

**Parameter Definitions:**  

*	`-c` – specifies the input assembly fasta file

*	`-d` – specifies the CAT reference sequence database

*	`-t` – specifies the CAT reference taxonomy database

*	`-p` – specifies the input protein fasta file

*	`-o` – specifies the output prefix

*	`-n` – specifies the number of CPU cores to use

*	`-r` – specifies the number of top protein hits to consider in assigning tax

*	`--top` – specifies the number of protein alignments to store

*	`--I_know_what_Im_doing` – allows us to alter the `--top` parameter

* `--no-stars` - suppress marking of suggestive taxonomic assignments


**Input data:**

* sample-1-assembly.fasta (assembly file from [step 5a](#5a-renaming-contig-headers))
* sample-1-genes.faa (gene-calls amino-acid fasta file from [step 6](#6-gene-prediction))

**Output data:**

* sample-1-tax-out.tmp.ORF2LCA.txt (gene-calls taxonomy file)
* sample-1-tax-out.tmp.contig2classification.txt (contig taxonomy file)

#### 8c. Adding taxonomy info from taxids to genes
```
CAT add_names -i sample-1-tax-out.tmp.ORF2LCA.txt -o sample-1-gene-tax-out.tmp \
              -t CAT_prepare_20200618/2020-06-18_taxonomy/ --only_official --exclude-scores
```

**Parameter Definitions:**  

*	`-i` – specifies the input taxonomy file

*	`-o` – specifies the output file 

*	`-t` – specifies the CAT reference taxonomy database

*	`--only_official` – specifies to add only standard taxonomic ranks

* `--exclude-scores` - specifies to exclude bit-score support scores in the lineage

**Input data:**

* sample-1-tax-out.tmp.ORF2LCA.txt (gene-calls taxonomy file from [step 8b](#8b-running-taxonomic-classification))

**Output data:**

* sample-1-gene-tax-out.tmp (gene-calls taxonomy file with lineage info added)



#### 8d. Adding taxonomy info from taxids to contigs
```
CAT add_names -i sample-1-tax-out.tmp.contig2classification.txt -o sample-1-contig-tax-out.tmp \
              -t CAT-ref/2020-06-18_taxonomy/ --only_official --exclude-scores
```

**Parameter Definitions:**  

*	`-i` – specifies the input taxonomy file

*	`-o` – specifies the output file 

*	`-t` – specifies the CAT reference taxonomy database

*	`--only_official` – specifies to add only standard taxonomic ranks

* `--exclude-scores` - specifies to exclude bit-score support scores in the lineage


**Input data:**

* sample-1-tax-out.tmp.contig2classification.txt (contig taxonomy file from [step 8b](#8b-running-taxonomic-classification))

**Output data:**

* sample-1-contig-tax-out.tmp (contig taxonomy file with lineage info added)


#### 8e. Formatting gene-level output with awk and sed
```
awk -F $'\t' ' BEGIN { OFS=FS } { if ( $3 == "lineage" ) { print $1,$3,$5,$6,$7,$8,$9,$10,$11 } \
    else if ( $2 == "ORF has no hit to database" || $2 ~ /^no taxid found/ ) \
    { print $1,"NA","NA","NA","NA","NA","NA","NA","NA" } else { n=split($3,lineage,";"); \
    print $1,lineage[n],$5,$6,$7,$8,$9,$10,$11 } } ' sample-1-gene-tax-out.tmp | \
    sed 's/no support/NA/g' | sed 's/superkingdom/domain/' | sed 's/# ORF/gene_ID/' | \
    sed 's/lineage/taxid/'  > sample-1-gene-tax-out.tsv
```

#### 8f. Formatting contig-level output with awk and sed
```
awk -F $'\t' ' BEGIN { OFS=FS } { if ( $2 == "classification" ) { print $1,$4,$6,$7,$8,$9,$10,$11,$12 } \
    else if ( $2 == "no taxid assigned" ) { print $1,"NA","NA","NA","NA","NA","NA","NA","NA" } \
    else { n=split($4,lineage,";"); print $1,lineage[n],$6,$7,$8,$9,$10,$11,$12 } } ' sample-1-contig-tax-out.tmp | \
    sed 's/no support/NA/g' | sed 's/superkingdom/domain/' | sed 's/^# contig/contig_ID/' | \
    sed 's/lineage/taxid/' > sample-1-contig-tax-out.tsv

  # clearing intermediate files
rm sample-1*.tmp*
```

**Input data:**

* sample-1-gene-tax-out.tmp (gene-calls taxonomy file with lineage info added from [step 8c](#8c-adding-taxonomy-info-from-taxids-to-genes))
* sample-1-contig-tax-out.tmp (contig taxonomy file with lineage info added from [step 8d](#8d-adding-taxonomy-info-from-taxids-to-contigs))


**Output data:**

* sample-1-gene-tax-out.tsv (gene-calls taxonomy file with lineage info added reformatted)
* sample-1-contig-tax-out.tsv (contig taxonomy file with lineage info added reformatted)

<br>

---

### 9. Read-mapping

#### 9a. Building reference index
```
bowtie2-build sample-1-assembly.fasta sample-1-assembly-bt-index
```

**Parameter Definitions:**  

*	`sample-1-assembly.fasta` - first positional argument specifies the input assembly

*	`sample-1-assembly-bt-index` - second positional argument specifies the prefix of the output index files


#### 9b. Performing mapping, conversion to bam, and sorting
```
bowtie2 --threads NumberOfThreads -x sample-1-assembly-bt-index -1 sample-1_R1_filtered.fastq.gz \
        -2 sample-1_R2_filtered.fastq.gz --no-unal 2> sample-1-mapping-info.txt | samtools view -b | samtools sort -@ NumberOfThreads > sample-1.bam
```

**Parameter Definitions:**  

*	`--threads` – specifies the number of threads to run in parallel

*	`-x` – specifies the prefix of the reference index files to map to (generated in [step 9a](#9a-building-reference-index))

*	`-1 and -2` – specifies the forward and reverse reads to map (if single-end data, neither `-1` nor `-2` are provided, and the single-end reads are passed to `-r`)

* `--no-unal` - suppress SAM records for unaligned reads

* `2> sample-1-mapping-info.txt` – capture the printed summary results in a log file

*	`samtools view -b` – convert the output directly to bam format (compressed)

*	`samtools sort -@` – sort the bam file using the specified number of threads

*	`>` – redirect the output to a file

#### 9c. Indexing
```
samtools index -@ NumberOfThreads sample-1.bam 
```

**Parameter Definitions:**  
*	`-@` – set number of threads to use 

*	`sample-1.bam` - input bam file is provided as a positional argument as generated in [step 9b](#9b-performing-mapping-conversion-to-bam-and-sorting)

**Input data:**

* sample-1-assembly.fasta (assembly file from [step 5](#5a-renaming-contig-headers))
* *.fastq.gz (filtered/trimmed reads from [step 2](#2-quality-filteringtrimming))

**Output data:**

* **sample-1.bam** (mapping file)
* sample-1.bam.bai (bam index file)
* **sample-1-mapping-info.txt** (read-mapping log file)

<br>

---

### 10. Getting coverage information and filtering based on detection
> **Notes**  
> “Detection” is a metric of what proportion of a reference sequence recruited reads (see [here](http://merenlab.org/2017/05/08/anvio-views/#detection)). Filtering based on detection is one way of helping to mitigate non-specific read-recruitment.

#### 10a. Filtering coverage levels based on detection

```
  # pileup.sh comes from the bbduk.sh package
pileup.sh -in sample-1.bam fastaorf=sample-1-genes.fasta outorf=sample-1-gene-cov-and-det.tmp \
          out=sample-1-contig-cov-and-det.tmp
```

**Parameter Definitions:**  

*	`-in` – the input bam file

*	`fastaorf=` – input gene-calls nucleotide fasta file

*	`outorf=` – the output gene-coverage tsv file

*	`out=` – the output contig-coverage tsv file


#### 10b. Filtering gene coverage based on requiring 50% detection and parsing down to just gene ID and coverage
```
grep -v "#" sample-1-gene-cov-and-det.tmp | awk -F $'\t' ' BEGIN { OFS=FS } { if ( $10 <= 0.5 ) $4 = 0 } \
     { print $1,$4 } ' > sample-1-gene-cov.tmp

cat <( printf "gene_ID\tcoverage\n" ) sample-1-gene-cov.tmp > sample-1-gene-coverages.tsv
```

Filtering contig coverage based on requiring 50% detection and parsing down to just contig ID and coverage:
```
grep -v "#" sample-1-contig-cov-and-det.tmp | awk -F $'\t' ' BEGIN { OFS=FS } { if ( $5 <= 50 ) $2 = 0 } \
     { print $1,$2 } ' > sample-1-contig-cov.tmp

cat <( printf "contig_ID\tcoverage\n" ) sample-1-contig-cov.tmp > sample-1-contig-coverages.tsv

  # removing intermediate files

rm sample-1-*.tmp
```

**Input data:**

* sample-1.bam (mapping file from [step 9b](#9b-performing-mapping-conversion-to-bam-and-sorting))
* sample-1-genes.fasta (gene-calls nucleotide fasta file from [step 6](#6-gene-prediction))

**Output data:**

* sample-1-gene-coverages.tsv (table with gene-level coverages)
* sample-1-contig-coverages.tsv (table with contig-level coverages)

<br>

---

### 11. Combining gene-level coverage, taxonomy, and functional annotations into one table for each sample
> **Notes**  
> Just uses `paste`, `sed`, and `awk`, all are standard in any Unix-like environment.  

```
paste <( tail -n +2 sample-1-gene-coverages.tsv | sort -V -k 1 ) <( tail -n +2 sample-1-annotations.tsv | sort -V -k 1 | cut -f 2- ) \
      <( tail -n +2 sample-1-gene-tax-out.tsv | sort -V -k 1 | cut -f 2- ) > sample-1-gene-tab.tmp

paste <( head -n 1 sample-1-gene-coverages.tsv ) <( head -n 1 sample-1-annotations.tsv | cut -f 2- ) \
      <( head -n 1 sample-1-gene-tax-out.tsv | cut -f 2- ) > sample-1-header.tmp

cat sample-1-header.tmp sample-1-gene-tab.tmp > sample-1-gene-coverage-annotation-and-tax.tsv

  # removing intermediate files
rm sample-1*tmp sample-1-gene-coverages.tsv sample-1-annotations.tsv sample-1-gene-tax-out.tsv
```

**Input data:**

* sample-1-gene-coverages.tsv (table with gene-level coverages from [step 10b](#10b-filtering-gene-coverage-based-on-requiring-50-detection-and-parsing-down-to-just-gene-id-and-coverage))
* sample-1-annotations.tsv (table of KO annotations assigned to gene IDs from [step 7c](#7c-filtering-output-to-retain-only-those-passing-the-ko-specific-score-and-top-hits))
* sample-1-gene-tax-out.tsv (gene-level taxonomic classifications from [step 8f](#8f-formatting-contig-level-output-with-awk-and-sed))


**Output data:**

* **sample-1-gene-coverage-annotation-and-tax.tsv** (table with combined gene coverage, annotation, and taxonomy info)

<br>

---

### 12. Combining contig-level coverage and taxonomy into one table for each sample
> **Notes**  
> Just uses `paste`, `sed`, and `awk`, all are standard in any Unix-like environment.  

```
paste <( tail -n +2 sample-1-contig-coverages.tsv | sort -V -k 1 ) \
      <( tail -n +2 sample-1-contig-tax-out.tsv | sort -V -k 1 | cut -f 2- ) > sample-1-contig.tmp

paste <( head -n 1 sample-1-contig-coverages.tsv ) <( head -n 1 sample-1-contig-tax-out.tsv | cut -f 2- ) \
      > sample-1-contig-header.tmp
      
cat sample-1-contig-header.tmp sample-1-contig.tmp > sample-1-contig-coverage-and-tax.tsv

  # removing intermediate files
rm sample-1*tmp sample-1-contig-coverages.tsv sample-1-contig-tax-out.tsv
```

**Input data:**

* sample-1-contig-coverages.tsv (table with contig-level coverages from [step 10b](#10b-filtering-gene-coverage-based-on-requiring-50-detection-and-parsing-down-to-just-gene-id-and-coverage))
* sample-1-contig-tax-out.tsv (contig-level taxonomic classifications from [step 8f](#8f-formatting-contig-level-output-with-awk-and-sed))


**Output data:**

* **sample-1-contig-coverage-and-tax.tsv** (table with combined contig coverage and taxonomy info)

<br>

---

### 13. Generating normalized, gene- and contig-level coverage summary tables of KO-annotations and taxonomy across samples

> **Notes**  
> * To combine across samples to generate these summary tables, we need the same "units". This is done for annotations based on the assigned KO terms, and all non-annotated functions are included together as "Not annotated". It is done for taxonomic classifications based on taxids (full lineages included in the table), and any not classified are included together as "Not classified". 
> * The values we are working with are coverage per gene (so they are number of bases recruited to the gene normalized by the length of the gene). These have been normalized by making the total coverage of a sample 1,000,000 and setting each individual gene-level coverage its proportion of that 1,000,000 total. So basically percent, but out of 1,000,000 instead of 100 to make the numbers more friendly. 

#### 13a. Generating gene-level coverage summary tables

```
bit-GL-combine-KO-and-tax-tables *-gene-coverage-annotation-and-tax.tsv -o Combined
```

**Parameter Definitions:**  

*	takes positional arguments specifying the input tsv files, can be provided as a space-delimited list of files, or with wildcards like above

-	`-o` – specifies the output prefix


**Input data:**

* *-gene-coverage-annotation-and-tax.tsv (tables with combined gene coverage, annotation, and taxonomy info generated for individual samples from [step 11](#11-combining-gene-level-coverage-taxonomy-and-functional-annotations-into-one-table-for-each-sample))

**Output data:**

* **Combined-gene-level-KO-function-coverages-CPM_GLmetagenomics.tsv** (table with all samples combined based on KO annotations; normalized to coverage per million genes covered)
* **Combined-gene-level-taxonomy-coverages-CPM_GLmetagenomics.tsv** (table with all samples combined based on gene-level taxonomic classifications; normalized to coverage per million genes covered)
* **Combined-gene-level-KO-function-coverages_GLmetagenomics.tsv** (table with all samples combined based on KO annotations)
* **Combined-gene-level-taxonomy-coverages_GLmetagenomics.tsv** (table with all samples combined based on gene-level taxonomic classifications)


#### 13b. Generating contig-level coverage summary tables

```
bit-GL-combine-contig-tax-tables *-contig-coverage-and-tax.tsv -o Combined
```
**Parameter Definitions:**  

*	takes positional arguments specifying the input tsv files, can be provided as a space-delimited list of files, or with wildcards like above

-	`-o` – specifies the output prefix


**Input data:**

* *-contig-coverage-annotation-and-tax.tsv (tables with combined contig coverage, annotation, and taxonomy info generated for individual samples from [step 12](#12-combining-contig-level-coverage-and-taxonomy-into-one-table-for-each-sample))

**Output data:**

* **Combined-contig-level-taxonomy-coverages-CPM_GLmetagenomics.tsv** (table with all samples combined based on contig-level taxonomic classifications; normalized to coverage per million genes covered)
* **Combined-contig-level-taxonomy-coverages_GLmetagenomics.tsv** (table with all samples combined based on contig-level taxonomic classifications)

<br>

---

### 14. **M**etagenome-**A**ssembled **G**enome (MAG) recovery

#### 14a. Binning contigs
```
jgi_summarize_bam_contig_depths --outputDepth sample-1-metabat-assembly-depth.tsv --percentIdentity 97 --minContigLength 1000 --minContigDepth 1.0  --referenceFasta sample-1-assembly.fasta sample-1.bam

metabat2  --inFile sample-1-assembly.fasta --outFile sample-1 --abdFile sample-1-metabat-assembly-depth.tsv -t NumberOfThreads

mkdir sample-1-bins
mv sample-1*bin*.fasta sample-1-bins
zip -r sample-1-bins.zip sample-1-bins
```

**Parameter Definitions:**  

*  `--outputDepth` – specifies the output depth file
*  `--percentIdentity` – minimum end-to-end percent identity of a mapped read to be included
*  `--minContigLength` – minimum contig length to include
*  `--minContigDepth` – minimum contig depth to include
*  `--referenceFasta` – the assembly fasta file generated in step 5a
*  `sample-1.bam` – final positional arguments are the bam files generated in step 9
*  `--inFile` - the assembly fasta file generated in step 5a
*  `--outFile` - the prefix of the identified bins output files
*  `--abdFile` - the depth file generated by the previous `jgi_summarize_bam_contig_depths` command
*  `-t` - specifies number of threads to use


**Input data:**

* sample-1-assembly.fasta (assembly fasta file created in [step 5a](#5a-renaming-contig-headers))
* sample-1.bam (bam file created in [step 9b](#9b-performing-mapping-conversion-to-bam-and-sorting))

**Output data:**

* **sample-1-metabat-assembly-depth.tsv** (tab-delimited summary of coverages)
* sample-1-bins/sample-1-bin\*.fasta (fasta files of recovered bins)
* **sample-1-bins.zip** (zip file containing fasta files of recovered bins)

#### 14b. Bin quality assessment
Utilizes the default `checkm` database available [here](https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz), `checkm_data_2015_01_16.tar.gz`.

```
checkm lineage_wf -f bins-overview_GLmetagenomics.tsv --tab_table -x fa ./ checkm-output-dir
```

**Parameter Definitions:**  

*  `lineage_wf` – specifies the workflow being utilized
*  `-f` – specifies the output summary file
*  `--tab_table` – specifies the output summary file should be a tab-delimited table
*  `-x` – specifies the extension that is on the bin fasta files that are being assessed
*  `./` – first positional argument at end specifies the directory holding the bins generated in step 14a
*  `checkm-output-dir` – second positional argument at end specifies the primary checkm output directory with detailed information

**Input data:**

* sample-1-bins/sample-1-bin\*.fasta (bin fasta files generated in [step 14a](#14a-binning-contigs))

**Output data:**

* **bins-overview_GLmetagenomics.tsv** (tab-delimited file with quality estimates per bin)
* checkm-output-dir (directory holding detailed checkm outputs)

#### 14c. Filtering MAGs

```
cat <( head -n 1 bins-overview_GLmetagenomics.tsv ) \
    <( awk -F $'\t' ' $12 >= 90 && $13 <= 10 && $14 == 0 ' bins-overview_GLmetagenomics.tsv | sed 's/bin./MAG-/' ) \
    > checkm-MAGs-overview.tsv
    
# copying bins into a MAGs directory in order to run tax classification
awk -F $'\t' ' $12 >= 90 && $13 <= 10 && $14 == 0 ' bins-overview_GLmetagenomics.tsv | cut -f 1 > MAG-bin-IDs.tmp

mkdir MAGs
for ID in MAG-bin-IDs.tmp
do
    MAG_ID=$(echo $ID | sed 's/bin./MAG-/')
    cp ${ID}.fasta MAGs/${MAG_ID}.fasta
done

for SAMPLE in $(cat MAG-bin-IDs.tmp | sed 's/-bin.*//' | sort -u);
do
  mkdir ${SAMPLE}-MAGs
  mv ${SAMPLE}-*MAG*.fasta ${SAMPLE}-MAGs
  zip -r ${SAMPLE}-MAGs.zip ${SAMPLE}-MAGs
done
```

**Input data:**

* bins-overview_GLmetagenomics.tsv (tab-delimited file with quality estimates per bin from [step 14b](#14b-bin-quality-assessment))

**Output data:**

* checkm-MAGs-overview.tsv (tab-delimited file with quality estimates per MAG)
* MAGs/\*.fasta (directory holding high-quality MAGs)
* **\*-MAGs.zip** (zip files containing directories of high-quality MAGs)


#### 14d. MAG taxonomic classification
Uses default `gtdbtk` database setup with program's `download.sh` command.

```
gtdbtk classify_wf --genome_dir MAGs/ -x fa --out_dir gtdbtk-output-dir  --skip_ani_screen
```

**Parameter Definitions:**  

*  `classify_wf` – specifies the workflow being utilized
*  `--genome_dir` – specifies the directory holding the MAGs generated in step 14c
*  `-x` – specifies the extension that is on the MAG fasta files that are being taxonomically classified
*  `--out_dir` – specifies the output directory
*  `--skip_ani_screen`  - specifies to skip ani_screening step to classify genomes using mash and skani

**Input data:**

* MAGs/\*.fasta (directory holding high-quality MAGs from [step 14c](#14c-filtering-mags))

**Output data:**

* gtdbtk-output-dir/gtdbtk.\*.summary.tsv (files with assigned taxonomy and info)

#### 14e. Generating overview table of all MAGs

```bash
# combine summaries
for MAG in $(cut -f 1 assembly-summaries_GLmetagenomics.tsv | tail -n +2); do

    grep -w -m 1 "^${MAG}" checkm-MAGs-overview.tsv | cut -f 12,13,14 \
        >> checkm-estimates.tmp

    grep -w "^${MAG}" gtdbtk-output-dir/gtdbtk.*.summary.tsv | \
    cut -f 2 | sed 's/^.__//' | \
    sed 's/;.__/\t/g' | \
    awk 'BEGIN{ OFS=FS="\t" } { for (i=1; i<=NF; i++) if ( $i ~ /^ *$/ ) $i = "NA" }; 1' \
        >> gtdb-taxonomies.tmp

done

# Add headers
cat <(printf "est. completeness\test. redundancy\test. strain heterogeneity\n") checkm-estimates.tmp \
    > checkm-estimates-with-headers.tmp

cat <(printf "domain\tphylum\tclass\torder\tfamily\\tgenus\tspecies\n") gtdb-taxonomies.tmp \
    > gtdb-taxonomies-with-headers.tmp

paste assembly-summaries_GLmetagenomics.tsv \
checkm-estimates-with-headers.tmp \
gtdb-taxonomies-with-headers.tmp \
    > MAGs-overview.tmp

# Ordering by taxonomy
head -n 1 MAGs-overview.tmp > MAGs-overview-header.tmp

tail -n +2 MAGs-overview.tmp | sort -t \$'\t' -k 14,20 > MAGs-overview-sorted.tmp

cat MAGs-overview-header.tmp MAGs-overview-sorted.tmp \
    > MAGs-overview_GLmetagenomics.tsv

```

**Input data:**

* assembly-summaries_GLmetagenomics.tsv (table of assembly summary statistics from [step 5b](#5b-summarizing-assemblies))
* MAGs/\*.fasta (directory holding high-quality MAGs from [step 14c](#14c-filtering-mags))
* checkm-MAGs-overview.tsv (tab-delimited file with quality estimates per MAG from [step 14c](#14c-filtering-mags))
* gtdbtk-output-dir/gtdbtk.\*.summary.tsv (directory of files with assigned taxonomy and info from [step 14d](#14d-mag-taxonomic-classification))

**Output data:**

* **MAGs-overview_GLmetagenomics.tsv** (a tab-delimited overview of all recovered MAGs)


<br>

---

### 15. Generating MAG-level functional summary overview

#### 15a. Getting KO annotations per MAG
This utilizes the helper script [`parse-MAG-annots.py`](../Workflow_Documentation/NF_MGIllumina/workflow_code/bin/parse-MAG-annots.py) 

```bash
for file in $( ls MAGs/*.fasta )
do

    MAG_ID=$( echo ${file} | cut -f 2 -d "/" | sed 's/.fasta//' )
    sample_ID=$( echo ${MAG_ID} | sed 's/-MAG-[0-9]*$//' )

    grep "^>" ${file} | tr -d ">" > ${MAG_ID}-contigs.tmp

    python parse-MAG-annots.py -i annotations-and-taxonomy/${sample_ID}-gene-coverage-annotation-and-tax.tsv \
                               -w ${MAG_ID}-contigs.tmp -M ${MAG_ID} \
                               -o MAG-level-KO-annotations_GLmetagenomics.tsv

    rm ${MAG_ID}-contigs.tmp

done
```

**Parameter Definitions:**  

*	`-i` – specifies the input sample gene-coverage-annotation-and-tax.tsv file generated in step 11

*  `-w` – specifies the appropriate temporary file holding all the contigs in the current MAG

*	`-M` – specifies the current MAG unique identifier

*	`-o` – specifies the output file

**Input data:**

* \*-gene-coverage-annotation-and-tax.tsv (sample gene-coverage-annotation-and-tax.tsv file generated in [step 11](#11-combining-gene-level-coverage-taxonomy-and-functional-annotations-into-one-table-for-each-sample))
* MAGs/\*.fasta (directory holding high-quality MAGs from [step 14c](#14c-filtering-mags))

**Output data:**

* **MAG-level-KO-annotations_GLmetagenomics.tsv** (tab-delimited table holding MAGs and their KO annotations)


#### 15b. Summarizing KO annotations with KEGG-Decoder

```bash
KEGG-decoder -v interactive -i MAG-level-KO-annotations_GLmetagenomics.tsv -o MAG-KEGG-Decoder-out_GLmetagenomics.tsv
```

**Parameter Definitions:**  

*  `-v interactive` – specifies to create an interactive html output
 
*	`-i` – specifies the input MAG-level-KO-annotations_GLmetagenomics.tsv file generated in [step 15a](#15a-getting-ko-annotations-per-mag)

*	`-o` – specifies the output table

**Input data:**

* MAG-level-KO-annotations_GLmetagenomics.tsv (tab-delimited table holding MAGs and their KO annotations, generated in [step 15a](#15a-getting-ko-annotations-per-mag))

**Output data:**

* **MAG-KEGG-Decoder-out_GLmetagenomics.tsv** (tab-delimited table holding MAGs and their proportions of genes held known to be required for specific pathways/metabolisms)

* **MAG-KEGG-Decoder-out_GLmetagenomics.html** (interactive heatmap html file of the above output table)

<br>
---

## Read-based processing
### 16. Taxonomic and functional profiling
The following uses the `humann` and `metaphlan` reference databases downloaded on 13-Jun-2024 as follows:

```bash
humann_databases --download chocophlan full
humann_databases --download uniref uniref90_diamond 
humann_databases --download utility_mapping full 
metaphlan --install
```

#### 16a. Running humann (which also runs metaphlan)
```bash
  # forward and reverse reads need to be provided combined if paired-end (if not paired-end, single-end reads are provided to the --input argument next)
cat sample-1_R1_filtered.fastq.gz sample-1_R2_filtered.fastq.gz > sample-1-combined.fastq.gz

humann --input sample-1-combined.fastq.gz --output sample-1-humann3-out-dir --threads NumberOfThreads \
       --output-basename sample-1 --metaphlan-options "--unknown_estimation --add_viruses \
       --sample_id sample-1"
```

**Parameter Definitions:**  

*	`--input` – specifies the input combined forward and reverse reads (if paired-end)

*	`--output` – specifies output directory

*	`--threads` – specifies the number of threads to use

*	`--output-basename` – specifies prefix of the output files

*	`--metaphlan-options` – options to be passed to metaphlan
	* `--unknown_estimation` – include unclassified in estimated relative abundances
	* `--add_viruses` – include viruses in the reference database
	* `--sample_id` – specifies the sample identifier we want in the table (rather than full filename)


#### 16b. Merging multiple sample functional profiles into one table
```bash
  # they need to be in their own directories
mkdir genefamily-results/ pathabundance-results/ pathcoverage-results/

  # copying results from previous running humann3 step (16a) to get them all together in their own directories (as is needed)
cp *-humann3-out-dir/*genefamilies.tsv genefamily-results/
cp *-humann3-out-dir/*abundance.tsv pathabundance-results/
cp *-humann3-out-dir/*coverage.tsv pathcoverage-results/

humann_join_tables -i genefamily-results/ -o gene-families.tsv
humann_join_tables -i pathabundance-results/ -o path-abundances.tsv
humann_join_tables -i pathcoverage-results/ -o path-coverages.tsv
```

**Parameter Definitions:**  

*	`-i` – the directory holding the input tables

*	`-o` – the name of the output combined table


#### 16c. Splitting results tables
The read-based functional annotation tables have taxonomic info and non-taxonomic info mixed together initially. `humann` comes with a helper script to split these. Here we are using that to generate both non-taxonomically grouped functional info files and taxonomically grouped ones.

```bash
humann_split_stratified_table -i gene-families.tsv -o ./
mv gene-families_stratified.tsv Gene-families-grouped-by-taxa_GLmetagenomics.tsv
mv gene-families_unstratified.tsv Gene-families_GLmetagenomics.tsv

humann_split_stratified_table -i path-abundances.tsv -o ./
mv path-abundances_stratified.tsv Path-abundances-grouped-by-taxa_GLmetagenomics.tsv
mv path-abundances_unstratified.tsv Path-abundances_GLmetagenomics.tsv

humann2_split_stratified_table -i path-coverages.tsv -o ./
mv path-coverages_stratified.tsv Path-coverages-grouped-by-taxa_GLmetagenomics.tsv
mv path-coverages_unstratified.tsv Path-coverages_GLmetagenomics.tsv
```

**Parameter Definitions:**  

*	`-i` – the input combined table

*	`-o` – output directory (here specifying current directory)


#### 16d. Normalizing gene families and pathway abundance tables
This generates some normalized tables of the read-based functional outputs from humann that are more readily suitable for across sample comparisons.

```bash
humann_renorm_table -i Gene-families_GLmetagenomics.tsv -o Gene-families-cpm_GLmetagenomics.tsv --update-snames
humann_renorm_table -i Path-abundances_GLmetagenomics.tsv -o Path-abundances-cpm_GLmetagenomics.tsv --update-snames
```

**Parameter Definitions:**  

*	`-i` – the input combined table

*	`-o` – name of the output normalized table

*	`--update-snames` – change suffix of column names in tables to "-CPM"


#### 16e. Generating a normalized gene-family table that is grouped by Kegg Orthologs (KOs)

```bash
humann_regroup_table -i Gene-families_GLmetagenomics.tsv -g uniref90_ko | humann_rename_table -n kegg-orthology | \
                     humann_renorm_table -o Gene-families-KO-cpm_GLmetagenomics.tsv --update-snames
```

**Parameter Definitions:**  

*	`-i` – the input table

*	`-g` – the map to use to group uniref IDs into Kegg Orthologs

*	`|` – sending that output into the next humann command to add human-readable Kegg Orthology names

*	`-n` – specifying we are converting Kegg orthology IDs into Kegg orthology human-readable names

*	`|` – sending that output into the next humann command to normalize to copies-per-million

*	`-o` – specifying the final output file name

*  `--update-snames` – change suffix of column names in tables to "-CPM"

#### 16f. Combining taxonomy tables

```bash
merge_metaphlan_tables.py *-humann3-out-dir/*_humann_temp/*_metaphlan_bugs_list.tsv > Metaphlan-taxonomy_GLmetagenomics.tsv
```

**Parameter Definitions:**  

*	input metaphlan tables are provided as position arguments (produced during humann3 run in [step 16a](#16a-running-humann-which-also-runs-metaphlan)

*  `>` – output is redirected from stdout to a file


**Input data:**

* *fastq.gz (filtered/trimmed reads from [step 2](#2-quality-filteringtrimming), forward and reverse reads concatenated if paired-end)

**Output data:**

* **Gene-families_GLmetagenomics.tsv** (gene-family abundances) 
* **Gene-families-grouped-by-taxa_GLmetagenomics.tsv** (gene-family abundances grouped by taxa)
* **Gene-families-cpm_GLmetagenomics.tsv** (gene-family abundances normalized to copies-per-million)
* **Gene-families-KO-cpm_GLmetagenomics.tsv** (KO term abundances normalized to copies-per-million)
* **Pathway-abundances_GLmetagenomics.tsv** (pathway abundances)
* **Pathway-abundances-grouped-by-taxa_GLmetagenomics.tsv** (pathway abundances grouped by taxa)
* **Pathway-abundances-cpm_GLmetagenomics.tsv** (pathway abundances normalized to copies-per-million)
* **Pathway-coverages_GLmetagenomics.tsv** (pathway coverages)
* **Pathway-coverages-grouped-by-taxa_GLmetagenomics.tsv** (pathway coverages grouped by taxa)
* **Metaphlan-taxonomy_GLmetagenomics.tsv** (metaphlan estimated taxonomic relative abundances)

---
