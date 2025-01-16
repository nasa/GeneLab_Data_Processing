# Bioinformatics pipeline for Methylation Sequencing (Methyl-Seq) data

> **This document holds an overview and some example commands of how GeneLab processes bisulfite sequencing (methylseq) datasets. Exact processing commands for specific datasets that have been released are provided with their processed data in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/).**

---

**Date:** February 13, 2023  
**Revision:** -  
**Document Number:** GL-DPPD-7113  

**Submitted by:**  
Michael D. Lee (GeneLab Analysis Team)  

**Approved by:**  
Sylvain Costes (GeneLab Project Manager)  
Samrawit Gebre (GeneLab Deputy Project Manager)  
Amanda Saravia-Butler (GeneLab Data Processing Lead)  
Jonathan Galazka (GeneLab Project Scientist)  

---

# Table of contents

- [Software used](#software-used)
- [General processing overview with example commands](#general-processing-overview-with-example-commands)
  - [1. Raw Data QC](#1-raw-data-qc)
    - [1a. Raw Data QC](#1a-raw-data-qc)
    - [1b. Compile Raw Data QC](#1b-compile-raw-data-qc)
  - [2. Adapter trimming/quality filtering](#2-adapter-trimmingquality-filtering)
    - [If not RRBS or if RRBS using MseI digestion](#if-not-rrbs-or-if-rrbs-using-msei-digestion)
    - [If using a random priming bisulfite method](#if-using-a-random-priming-post-bisulfite-method)
    - [If RRBS with MspI digestion](#if-rrbs-with-mspi-digestion)
    - [If RRBS with NuGEN ovation kit](#if-rrbs-with-nugen-ovation-kit)
      - [First adapter-trimming/quality-filtering with trimgalore](#first-adapter-trimmingquality-filtering-with-trimgalore)
      - [Now running NuGEN-specific script](#now-running-nugen-specific-script)
  - [3. Filtered/Trimmed Data QC](#3-filteredtrimmed-data-qc)
    - [3a. Trimmed Data QC](#3a-trimmed-data-qc)
    - [3b. Compile Trimmed Data QC](#3b-compile-trimmed-data-qc)
  - [4. Alignment](#4-alignment)
    - [4a. Generate reference](#4a-generate-reference)
    - [4b. Align](#4b-align)
    - [4c. Sort Alignment Files](#4c-sort-alignment-files)
  - [5. Alignment QC](#5-alignment-qc)
  - [6. Deduplicate (skip if data are RRBS)](#6-deduplicate-skip-if-data-are-rrbs)
    - [6a. Deduplicate](#6a-deduplicate)
    - [6b. Sort Deduplicated Alignment Files](#6b-sort-deduplicated-alignment-files)
  - [7. Extract methylation calls](#7-extract-methylation-calls)
  - [8. Generate individual sample report](#8-generate-individual-sample-report)
  - [9. Generate combined summary reports](#9-generate-combined-summary-reports)
    - [9a. Bismark summary report](#9a-bismark-summary-report)
    - [9b. Compile Alignment and Bismark QC](#9b-compile-alignment-and-bismark-qc)
  - [10. Generate reference genome annotation information](#10-generate-reference-genome-annotation-information)
    - [10a. GTF to BED conversion](#10a-gtf-to-bed-conversion)
    - [10b. Making a mapping file of genes to transcripts](#10b-making-a-mapping-file-of-genes-to-transcripts)
  - [11. Differential methylation analysis](#11-differential-methylation-analysis)
    - [11a. Set up R environment](#11a-set-up-r-environment)
    - [11b. Individual-base analysis](#11b-individual-base-analysis)
    - [11c. Tile-level analysis](#11c-tile-level-analysis)
    - [11d. Add feature information](#11d-add-feature-information)
    - [11e. Add functional annotations](#11e-add-functional-annotations)
    - [11f. Make and write out tables of percent methylation levels](#11f-make-and-write-out-tables-of-percent-methylation-levels)
    - [11g. Make overview figure of percent differential methylation across features](#11g-make-overview-figure-of-percent-differential-methylation-across-features)

---

# Software used

|Program|Version|Relevant Links|
|:------------|:------:|------:|
|FastQC       | 0.12.1 |[https://www.bioinformatics.babraham.ac.uk/projects/fastqc/](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)|
|MultiQC      | 1.22.1 |[https://multiqc.info/](https://multiqc.info/)|
|Cutadapt     | 4.8    |[https://cutadapt.readthedocs.io/en/stable/](https://cutadapt.readthedocs.io/en/stable/)|
|TrimGalore!  | 0.6.10 |[https://github.com/FelixKrueger/TrimGalore](https://github.com/FelixKrueger/TrimGalore)|
|Bismark      | 0.24.2 |[https://github.com/FelixKrueger/Bismark](https://github.com/FelixKrueger/Bismark)|
|bowtie2      | 2.5.4  |[https://github.com/BenLangmead/bowtie2#overview](https://github.com/BenLangmead/bowtie2#overview)|
|hisat2       | 2.2.1  |[https://github.com/DaehwanKimLab/hisat2](https://github.com/DaehwanKimLab/hisat2)|
|samtools     | 1.20   |[https://github.com/samtools/samtools#samtools](https://github.com/samtools/samtools#samtools)|
|qualimap     | 2.3    |[http://qualimap.conesalab.org/](http://qualimap.conesalab.org/)|
|gtfToGenePred| 447    |[http://hgdownload.cse.ucsc.edu/admin/exe/](http://hgdownload.cse.ucsc.edu/admin/exe/)|
|genePredToBed| 447    |[http://hgdownload.cse.ucsc.edu/admin/exe/](http://hgdownload.cse.ucsc.edu/admin/exe/)
|R            | 4.4.0  |[https://www.r-project.org](https://www.r-project.org)|
|tidyverse    | 2.0.0  |[https://tidyverse.tidyverse.org/](https://tidyverse.tidyverse.org/)|
|Bioconductor | 3.19   |[https://bioconductor.org](https://bioconductor.org)
|methylKit    | 1.30.0 |[https://bioconductor.org/packages/release/bioc/html/methylKit.html](https://bioconductor.org/packages/release/bioc/html/methylKit.html)|
|genomation   | 1.36.0 |[https://bioconductor.org/packages/release/bioc/html/genomation.html](https://bioconductor.org/packages/release/bioc/html/genomation.html)|

---

# General processing overview with example commands

> Exact processing commands and output files listed in **bold** below are included with each RNAseq processed dataset in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/).

---

## 1. Raw Data QC

<br>

### 1a. Raw Data QC

```bash
fastqc -o raw_fastqc_output/ *raw.fastq.gz
```

**Parameter Definitions:**

* `-o` – the output directory to store results
* `*raw.fastq.gz` – the input reads are specified as a positional argument, and can be given all at once with wildcards, or as individual arguments with spaces in between them

**Input data:**

* *raw.fastq.gz (raw reads)

**Output data:**

* *fastqc.html (FastQC output html summary)
* *fastqc.zip (FastQC output data)

<br>

### 1b. Compile Raw Data QC

```bash
multiqc --interactive \
  -o raw_multiqc_GLmethylSeq_data/ \
  -n raw_multiqc_GLmethylSeq \
  -z \
  raw_fastqc_output/
```

**Parameter Definitions:**

*	`--interactive` - force reports to use interactive plots
*	`-o` – the output directory to store results
*	`-n` – the filename prefix for output files
*	`-z` – specifies to zip the output data directory
*	`raw_fastqc_output/` – the directory holding the output data from the fastqc run, provided as a positional argument

**Input data:**

* raw_fastqc_output/*fastqc.zip (FastQC output data from [Step 1a](#1a-raw-data-qc))

**Output data:**

* **raw_multiqc_GLMethylSeq.html** (multiqc output html summary)
* **raw_multiqc_GLMethylSeq_data.zip** (zipped directory containing multiqc output data)

> If using RNA, file suffix will be **GLRNAMethylSeq** instead of **GLMethylSeq**

<br>  

---

## 2. Adapter trimming/quality filtering
See `trim_galore --help` or [TrimGalore User Guide](https://github.com/FelixKrueger/TrimGalore/blob/0.6.10/Docs/Trim_Galore_User_Guide.md) for more info on any of the below.

Additionally, the Bismark documentation also includes guidelines for specific MethylSeq library types: [Bismark library type guide](http://felixkrueger.github.io/Bismark/bismark/library_types/). Some library types will require additional 5' and/or 3' hard trimming to remove the signature of the oligos used for random priming. Leaving these bases may cause misalignments and methylation biases.

<br>

### If not RRBS or if RRBS using MseI digestion
Note that the `--rrbs` option is **not** appropriate when RRBS (reduced representation bisulfite sequencing) libraries were prepared with MseI digestion (see the TrimGalore User Guide [Note for RRBS using MseI](https://github.com/FelixKrueger/TrimGalore/blob/0.6.10/Docs/Trim_Galore_User_Guide.md#rrbs-specific-options-mspi-digested-material).

**Single-end example**  

```bash
trim_galore --gzip \
  --cores NumberOfThreads \
  --phred33 \
  --output_dir trimmed_reads_out_dir/ \
  sample-1_raw.fastq.gz

# renaming output to use GeneLab standard conventions
mv sample-1_raw_trimmed.fq.gz sample-1_trimmed.fastq.gz
```

**Paired-end example**  

```bash
trim_galore --gzip \
  --cores NumberOfThreads \
  --phred33 \
  --output_dir trimmed_reads_out_dir/ \
  --paired \
  sample-1_R1_raw.fastq.gz sample-1_R2_raw.fastq.gz

# renaming outputs to use GeneLab standard suffix
mv sample-1_R1_raw_val_1.fq.gz sample-1_R1_trimmed.fastq.gz
mv sample-1_R2_raw_val_2.fq.gz sample-1_R2_trimmed.fastq.gz
```

<br>

### If using a random priming post-bisulfite method
(such as TruSeq (formerly EpiGnome), PBAT, scBSSeq, Pico Methyl, Accel, etc.)
Random priming is not truly random and the signature left at the ends of the reads can introduce errors, indels, and methylation biases. Add the optional clipping parameters (`--clip_r1`, `--clip_r2`, `--three_prime_clip_r1`, and `--three_prime_clip_r2`) to trim off the random priming signature on the 5' ends of each read and next to the 3' end after adapter trimming. See [Bismark library type guide](http://felixkrueger.github.io/Bismark/bismark/library_types/) for more detailed information. 

**Paired-end example for TruSeq (EpiGnome) library prep**
```bash
trim_galore --gzip \
  --cores NumberOfThreads \
  --phred33 \
  --output_dir trimmed_reads_out_dir/ \
  --paired \
  --clip_R1 8 \
  --clip_R2 8 \
  --three_prime_clip_R1 8 \
  --three_prime_clip_R2 8 \
  sample-1_R1_raw.fastq.gz sample-1_R2_raw.fastq.gz

# renaming outputs to use GeneLab standard suffix
mv sample-1_R1_raw_val_1.fq.gz sample-1_R1_trimmed.fastq.gz
mv sample-1_R2_raw_val_2.fq.gz sample-1_R2_trimmed.fastq.gz
```

<br>

### If RRBS with MspI digestion
Note that if the library preparation was non-directional, the `--non_directional` flag needs to be added to this command (whether single-end or paired-end; see [TrimGalore User Guide](https://github.com/FelixKrueger/TrimGalore/blob/0.6.10/Docs/Trim_Galore_User_Guide.md#rrbs-specific-options-mspi-digested-material)). 

**Single-end example**  

```bash
trim_galore --gzip \
  --cores NumberOfThreads \
  --phred33 \
  --rrbs \
  --output_dir trimmed_reads_out_dir/ \
  sample-1_raw.fastq.gz

# renaming output to use GeneLab standard conventions
mv sample-1_raw_trimmed.fq.gz sample-1_trimmed.fastq.gz
```

**Paired-end example**  

```bash
trim_galore --gzip \
  --cores NumberOfThreads \
  --phred33 \
  --rrbs \
  --output_dir trimmed_reads_out_dir/ \
  --paired \
  sample-1_R1_raw.fastq.gz sample-1_R2_raw.fastq.gz

# renaming outputs to use GeneLab standard suffix
mv sample-1_R1_raw_val_1.fq.gz sample-1_R1_trimmed.fastq.gz
mv sample-1_R2_raw_val_2.fq.gz sample-1_R2_trimmed.fastq.gz
```

<br>

### If RRBS with NuGEN ovation kit
Libraries prepared with the NuGEN ovation kit need to be processed with an additional script provided by the company's [github](https://github.com/nugentechnologies/NuMetRRBS#analysis-guide-for-nugen-ovation-rrbs-methyl-seq). 

Following their instructions, we first run an adapter-trimming/quality-filtering step with trimgalore. Note that the `--rrbs` option is not appropriate to pass to trimgalore when this kit is used (see Bismark documentation for [RRBS NuGEN Ovation Methyl-Seq System](http://felixkrueger.github.io/Bismark/bismark/library_types/#rrbs-nugen-ovation-methyl-seq-system). Then we utilize the company's script to remove the random diversity sequences added by the kit. 

#### First adapter-trimming/quality-filtering with trimgalore

**Single-end example**  

```bash
trim_galore --gzip \
  --cores NumberOfThreads \
  --phred33 \
  -a AGATCGGAAGAGC \
  --output_dir trimmed_reads_out_dir/ \
  sample-1_raw.fastq.gz

# renaming output to use GeneLab standard conventions
mv sample-1_raw_trimmed.fq.gz sample-1_trimmed.fastq.gz
```

**Paired-end example**  

```bash
trim_galore --gzip \
  --cores NumberOfThreads \
  --phred33 \
  --paired \
  -a AGATCGGAAGAGC \
  -a2 AAATCAAAAAAAC \
  --output_dir trimmed_reads_out_dir/ \
  sample-1_R1_raw.fastq.gz sample-1_R2_raw.fastq.gz

# renaming output to use GeneLab standard conventions
mv sample-1_R1_raw_val_1.fq.gz sample-1_R1_trimmed.fastq.gz
mv sample-1_R2_raw_val_2.fq.gz sample-1_R2_trimmed.fastq.gz
```

#### Now running NuGEN-specific script

The NuGEN-specific script can be downloaded from their [github](https://github.com/nugentechnologies/NuMetRRBS#analysis-guide-for-nugen-ovation-rrbs-methyl-seq) with the following:

```bash
curl -LO \
     https://raw.githubusercontent.com/nugentechnologies/NuMetRRBS/master/trimRRBSdiversityAdaptCustomers.py
```

**Single-end example**  

```bash
python2 trimRRBSdiversityAdaptCustomers.py \
  -1 sample-1_trimmed.fastq.gz

# renaming outputs to use GeneLab standard suffix
mv sample-1_trimmed.fastq_trimmed.fq.gz sample-1_trimmed.fastq.gz
```

**Paired-end example**  

```bash
python2 trimRRBSdiversityAdaptCustomers.py \
  -1 sample-1_R1_trimmed.fastq.gz \
  -2 sample-1_R2_trimmed.fastq.gz

# renaming outputs to use GeneLab standard suffix
mv sample-1_R1_trimmed.fastq_trimmed.fq.gz sample-1_R1_trimmed.fastq.gz
mv sample-1_R2_trimmed.fastq_trimmed.fq.gz sample-1_R2_trimmed.fastq.gz
```

<br>

**Parameter Definitions for `trim_galore`:**  

* `--gzip` - gzip compress the output(s)
* `--cores` - number of cores to use (this value is dependent on the number of threads available on the system running trim galore)
* `--phred33` - instructs cutadapt to use ASCII+33 quality scores as Phred scores for quality trimming
* `--rrbs` - specific trimming suitable for RRBS data generated with MspI digestion only
* `-a` - specific adapter sequence to be trimmed off of forward or single reads (applicable for libraries prepared with the NuGEN ovation kit)
* `-a2` - specific adapter sequence to be trimmed off of reverse reads (applicable for libraries prepared with the NuGEN ovation kit)
* `--paired` - specifies data are paired-end
* `--output_dir` - the output directory to store results
* `--clip_R1` - number of bases to trim off the 5' end of each R1 read (optional, for use with library prep kits that use random priming, such as TruSeq(EpiGnome))
* `--clip_R2` - number of bases to trim off the 5' end of each R2 read (optional, for use with library prep kits that use random priming, such as TruSeq(EpiGnome))
* `--three_prime_clip_R1` - number of bases to trim off the 3' end of each R1 read AFTER adapter trimming. (optional, for use with library prep kits that use random priming, such as TruSeq(EpiGnome)) 
* `--three_prime_clip_R2` - number of bases to trim off the 3' end of each R2 read AFTER adapter trimming. (optional, for use with library prep kits that use random priming, such as TruSeq(EpiGnome)) 
* `sample-1_raw.fastq.gz` or `sample-1_R1_raw.fastq.gz sample-1_R2_raw.fastq.gz` – the input reads are specified as positional arguments, paired-end read files are listed pairwise such that the forward reads (\*R1_raw.fastq.gz) are immediately followed by the respective reverse reads (\*R2_raw.fastq.gz)


**Parameter Definitions for `trimRRBSdiversityAdaptCustomers.py `:**  

- `-1` - forward or single input read file
- `-2` - reverse read file if data is paired-end

**Input Data:**

* \*fastq.gz (raw reads)

**Output Data:**

* **\*trimmed.fastq.gz** (adapter-trimmed/quality-filtered reads)
* \*trimming_report.txt (trimming report)

<br>

---

## 3. Filtered/Trimmed Data QC

<br>

### 3a. Trimmed Data QC

```bash
fastqc -o trimmed_fastqc_output/ *trimmed.fastq.gz
```

**Parameter Definitions:**

*	`-o` – the output directory to store results  
*	`*trimmed.fastq.gz` – the input reads are specified as a positional argument, and can be given all at once with wildcards as in the example, or as individual arguments with spaces in between them  

**Input data:**

* *trimmed.fastq.gz (filtered/trimmed reads, output from [Step 2](#2-adapter-trimmingquality-filtering))

**Output data:**

* *fastqc.html (FastQC output html summary)
* *fastqc.zip (FastQC output data)

<br>

### 3b. Compile Trimmed Data QC

```bash
multiqc --interactive \
  -o trimmed_multiqc_GLmethylSeq_data/ \
  -n trimmed_multiqc_GLmethylSeq \
  -z \
  trimmed_fastqc_output/ trimgalore_output/
```

**Parameter Definitions:**

*	`--interactive` - force reports to use interactive plots
*	`-o` – the output directory to store results
*	`-n` – the filename prefix for output files
*	`-z` – specifies to zip the output data directory
*	`trimmed_fastqc_output/` – the directory holding the output data from the fastqc run, provided as a positional argument
* `trimgalore_output/` - the directory holding the trimgalore trimming reports, provided as a positional argument

**Input data:**

* trimmed_fastqc_output/*fastqc.zip (FastQC output data from [Step 3a](#3a-trimmed-data-qc))

**Output data:**

* **trimmed_multiqc_GLMethylSeq.html** (multiqc output html summary)
* **trimmed_multiqc_GLMethylSeq_data.zip** (directory containing multiqc output data)

> If using RNA, file suffix will be **GLRNAMethylSeq** instead of **GLMethylSeq**

<br>

---

## 4. Alignment

<br>

### 4a. Generate reference
The reference will need to be specific to the organism that was sequenced. Bismark operates on a directory holding the target reference genome in fasta format.

```bash
# creating directory to hold reference and moving it into there
mkdir bismark_reference_genome
mv ref-genome.fasta bismark_reference_genome/

bismark_genome_preparation --bowtie2 \
  --parallel NumberOfThreads \
  bismark_reference_genome/

bam2nuc --genomic_composition_only \
  --genome_folder bismark_reference_genome/
```

**Parameter Definitions:**

*bismark_genome_preparation*
*	`--bowtie2` - specify bismark to create bisulfite indexes for use with Bowtie2
*	`--parallel` – specifies how many threads to use (note these will be doubled as it operates on both strands simultaneously)
*  `bismark_reference_genome/` - positional argument specifying the directory holding the reference genome (should end in ".fa" or ".fasta", can be gzipped and including ".gz")

*bam2nuc*
* --genomic_composition_only - specifies creation of the (genome-specific) genomic_nucleotide_frequencies.txt report  
* --genome_folder - specifies the directory holding the reference genome (should end in ".fa" or ".fasta", can be gzipped and including ".gz")

 
**Input data:**

* a directory holding the reference genome in fasta format (this pipeline version uses the Ensembl fasta file indicated in the `fasta` column of the [GL-DPPD-7110-A_annotations.csv](../../GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv) GeneLab Annotations file))

**Output data:**

* the reference genome directory that was provided as input will now hold indexes for the bisulfite-converted reference genome (all `*.bt2` files are indexes, all `*.fa` files are converted versions of the reference genome)
* bismark_reference_genome/Bisulfite_Genome/
  * CT_converion/
    * BS_CT.1.bt2
    * BS_CT.2.bt2
    * BS_CT.3.bt2
    * BS_CT.4.bt2
    * BS_CT.rev.1.bt2
    * BS_CT.rev.2.bt2
    * genome_mfa.CT_conversion.fa
  * GA_conversion/
    * BS_GA.1.bt2
    * BS_GA.2.bt2
    * BS_GA.3.bt2
    * BS_GA.4.bt2
    * BS_GA.rev.1.bt2
    * BS_GA.rev.2.bt2
    * genome_mfa.GA_conversion.fa
  * \*.txt (captured standard output from the command)
* **bismark_reference_genome/genomic_nucleotide_frequencies.txt** (tab-delimited table of mono- and di-nucleotide frequencies in reference genome)



> **NOTE**  
> If using RNA, the `--hisat2` flag is added instead of `--bowtie2`, which specifies to use the splice-aware aligner [HISAT2](https://github.com/DaehwanKimLab/hisat2#hisat2), and the outputs will include 8 `*ht2` files in separate sub-directories along with each reference-genome conversion.

<br>

### 4b. Align

Note that if the library preparation was non-directional, the `--non_directional` flag needs to be added to this command (whether single-end or paired-end). For a full list of alignment option recommendations library type and/or commercially available kit, please see the library page in the [Bismark documentation](http://felixkrueger.github.io/Bismark/bismark/library_types/) 

**Single-end example**  

```bash
bismark --bowtie2 \
  --bam \
  --parallel NumberOfThreads \
  --non_bs_mm \
  --nucleotide_coverage \
  --gzip \
  --output_dir mapping_files_out_dir/ \
  --genome_folder bismark_reference_genome/ \
  sample-1_trimmed.fastq.gz

# renaming output files so they are cleaner and will work with sorted bam 
# file/auto-detection of bismark2summary later
mv sample-1_trimmed_bismark_bt2_SE_report.txt sample-1_bismark_bt2_sorted_SE_report.txt
mv sample-1_trimmed_bismark_bt2.nucleotide_stats.txt sample-1_bismark_bt2.nucleotide_stats.txt
mv sample-1_trimmed_bismark_bt2.bam sample-1_bismark_bt2.bam
```

**Paired-end example**  

```bash
bismark --bowtie2 \
  --bam \
  --parallel NumberOfThreads \
  --non_bs_mm \
  --nucleotide_coverage \
  --gzip \
  --output_dir mapping_files_out_dir/ \
  --genome_folder bismark_reference_genome/ \
  -1 sample-1_R1_trimmed.fastq.gz \
  -2 sample-1_R2_trimmed.fastq.gz

# renaming output files so they are cleaner and will work with sorted bam 
# file/auto-detection of bismark2summary later
mv sample-1_R1_trimmed_bismark_bt2_PE_report.txt sample-1_bismark_bt2_sorted_PE_report.txt
mv sample-1_R1_trimmed_bismark_bt2_pe.nucleotide_stats.txt sample-1_bismark_bt2_pe.nucleotide_stats.txt
mv sample-1_R1_trimmed_bismark_bt2_pe.bam sample-1_bismark_bt2_pe.bam
```

**Parameter Definitions:**

* `--bowtie2` - specifies to use bowtie2 for alignment (limited to end-to-end alignments)
* `--bam` - specifies to convert the default output sam format into compressed bam format
* `--parallel` - allows us to specify the number of threads to use (note: will consume 3-5X this value)
* `--non_bs_mm` - outputs an extra column in the bam file specifying the number of non-bisulfite mismatches each read has
* `--nucleotide_coverage` - outputs a table with mono- and di-nucleotide sequence compositions and coverage values compared to genomic compositions
* `--gzip` - write temporary bisulfite conversion files in gzip format to save disk space during alignment
* `--output_dir` - the output directory to store results
* `--genome_folder` - specifies the directory holding the reference genome indexes (the same that was provided in [Step 4a.](#4a-generate-reference) above)
* `-1` - for paired-end data, the forward trimmed reads (not used for single-end data)
* `-2` - for paired-end data, the reverse trimmed reads (not used for single-end data)
* `sample-1_trimmed.fastq.gz` - for single-end data, input trimmed-reads are provided as a positional argument


**Input data:**
* bismark_reference_genome/ (directory holding indexes of reference genome)
* gzip-compressed fastq files (adapter-trimmed/quality-filtered reads)

**Output data:**  

* sample-1_bismark_bt2*.bam (mapping file) 
* **\*_[SP]E_report.txt** (bismark mapping report file)
* **\*.nucleotide_stats.txt** (tab-delimited table with sample-specific mono- and di-nucleotide sequence compositions and coverage values compared to genomic compositions)


> **NOTE**  
> If using RNA, need to add the `--hisat2` and `--path_to_hisat2` flags. Output files will include "bismark_hisat2" instead of "bismark_bt2" in the name.

<br>

### 4c. Sort Alignment Files

```bash
samtools sort -@ NumberOfThreads \
  -o sample-1_bismark_bt2_sorted.bam \
  sample-1_bismark_bt2*.bam
```

**Parameter Definitions:**

* `sort` - specifies to use the `sort` sub-program of `samtools`
* `-@` - specifies the number of threads to use
* `-o` - specifies the output file name
* sample-1_bismark_bt2*.bam - the input bam file, provided as a positional argument

**Input data:**

* sample-1_bismark_bt2*.bam (bismark alignment bam file, output from [Step 4b.](#4b-align) above)
> **NOTE**  
> If using RNA, files will include "bismark_hisat2" instead of "bismark_bt2" in the name.

**Output data:**

* **sample-1_bismark_bt2_sorted.bam** (bismark bowtie2 alignment bam file sorted by chromosomal coordinates)

<br>

---

## 5. Alignment QC

```bash
qualimap bamqc -bam sample-1_bismark_bt2_sorted.bam \
  -gff reference.gtf \
  -outdir sample-1_bismark_bt2_qualimap/ \
  --collect-overlap-pairs \
  --java-mem-size=6G \
  -nt NumberOfThreads
```

**Parameter Definitions:**

* `bamqc` - specifies the `bamqc` sub-program of `qualimap`
* `-bam` - specifies the input bam file
* `-gff` - specifies the feature file containing regions of interest for the reference genome (can be gff, gtf, or bed format)
* `-outdir` - specifies the path to print the alignment QC output files to
* `--collect-overlap-pairs` - instructs the program to output statistics of overlapping paired-end reads (if data were paired-end, no effect if single-end)
* `--java-mem-size=6G` - specifies the amount of memory to use (here this is set to 6G; see [qualimap FAQ here](http://qualimap.conesalab.org/doc_html/faq.html?highlight=java-mem-size))
* `-nt` - specifies the number of threads to use

**Input data:**

* sample-1_bismark_bt2_sorted.bam (bismark bowtie2 alignment bam file sorted by chromosomal coordinates, output from [Step 4c](#4c-sort-alignment-files) above)
* a feature file containing regions of interest for the reference genome in gtf format (this pipeline version uses the Ensembl fasta file indicated in the `gtf` column of the [GL-DPPD-7110-A_annotations.csv](../../GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv) GeneLab Annotations file))

**Output data:** 

* **\*sample-1_bismark_bt2_qualimap/** (subdirectory of many alignment QC output files and formatting files for presenting in an html file (see [qualimap documentation](http://qualimap.conesalab.org/doc_html/analysis.html#output))

> **NOTE**  
> If using RNA, files will include "bismark_hisat2" instead of "bismark_bt2" in the name.

<br>

---

## 6. Deduplicate (skip if data are RRBS)
> **NOTE**  
> This step should **not** be done if the data are RRBS (reduced representation bisulfite sequencing; see e.g., [bismark documentation](https://felixkrueger.github.io/Bismark/bismark/deduplication/)). 

<br>

### 6a. Deduplicate 

```bash
deduplicate_bismark sample-1_bismark_bt2*.bam
```

**Parameter Definitions:**

* `sample-1_bismark_bt2*.bam` - the input bam file, provided as a positional argument

**Input data:**

* sample-1_bismark_bt2*.bam (unsorted bismark bowtie2 alignment bam file, output from [Step 4b](#4b-align) above)
> **NOTE**  
> If using RNA, files will include "bismark_hisat2" instead of "bismark_bt2" in the name.

**Output data:**

* **\*.deduplicated.bam** (unsorted bismark bowtie2 alignment bam file, with duplicates removed)
* **\*.deduplication_report.txt** (report file containing deduplication information) 


<br>

### 6b. Sort Deduplicated Alignment Files

```bash
samtools sort -@ NumberOfThreads \
  -o sample-1_bismark_bt2*.deduplicated_sorted.bam \
  sample-1_bismark_bt2*.deduplicated.bam
```

**Parameter Definitions:**

* `sort` - specifies to use the `sort` sub-program of `samtools`
* `-@` - specifies the number of threads to use
* `-o` - specifies the output file name
* `sample-1_bismark_bt2*.deduplicated.bam` - the input bam file, provided as a positional argument

**Input data:**

* sample-1_bismark_bt2*.deduplicated.bam (bismark bowtie2 alignment bam file, output from [Step 6a.](#6a-deduplicate-skip-if-data-are-rrbs) above)
> **NOTE**  
> If using RNA, files will include "bismark_hisat2" instead of "bismark_bt2" in the name.

**Output data:**

* **sample-1_bismark_bt2*.deduplicated_sorted.bam** (bismark bowtie2 alignment bam file sorted by chromosomal coordinates, with duplicates removed)

<br>

---

## 7. Extract methylation calls


**Single-end example**  

```bash
bismark_methylation_extractor --parallel NumberOfThreads \
  --bedGraph \
  --gzip \
  --comprehensive \
  --output_dir methylation_calls_out_dir/ \
  --cytosine_report \
  --genome_folder bismark_reference_genome/ \
  sample-1_bismark_bt2*.bam
```

**Paired-end example**  

```bash
bismark_methylation_extractor --parallel NumberOfThreads \
  --bedGraph \
  --gzip \
  --comprehensive \
  --output_dir methylation_calls_out_dir/ \
  --cytosine_report \
  --genome_folder bismark_reference_genome/ \
  --ignore_r2 2 \
  --ignore_3prime_r2 2 \
  sample-1_bismark_bt2*.bam
```


**Parameter Definitions:**

* `--parallel` - specifies the number of cores to use for methylation extraction, note: the program will utilize ~3X the number specified 
* `--bedGraph` - instructs the program to generate a sorted bedGraph file that reports the position of a given cytosine and its methlyation state (by default, only methylated CpGs are reported - see bedgraph options in [bismark documentation](https://felixkrueger.github.io/Bismark/options/methylation_extraction/#bedgraph-specific-options) for more info)
* `--gzip` - specifies to gzip-compress the methylation extractor output files
* `--comprehensive` - specifies to merge all four possible strand-specific methylation info into context-dependent output files
* `--output_dir` - the output directory to store results
* `--cytosine_report` - instructions the program to produce a genome-wide methylation report for all cytosines in the genome
* `--genome_folder` - a directory holding the reference genome in fasta format (this pipeline version uses the Ensembl fasta file indicated in the `fasta` column of the [GL-DPPD-7110-A_annotations.csv](../../GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv) GeneLab Annotations file))
* `--ignore_r2` - specifies how many bases to ignore from the 5' end of the reverse reads (bismark docs recommend 2, see [bismark documentation](https://felixkrueger.github.io/Bismark/options/methylation_extraction/#options))
  > **Note:** The first couple of bases in the reverse read of bisulfite sequence experiments show a severe bias towards non-methylation as a result of end-repairing sonicated fragments with unmethylated cytosines, so it is recommended to remove the first few basepairs
* `--ignore_3prime_r2` - specifies how many bases to ignore from the 3' end of the reverse reads to remove unwanted biases from the end of reads (For specific recommendations see Bismark documentation on [Library Types](https://felixkrueger.github.io/Bismark/bismark/library_types/)) 
* `sample-1_bismark_bt2*.bam` - the input bam file, provided as a positional argument

**Input data:**

* sample-1_bismark_bt2*.bam (bismark bowtie2 alignment bam file, output from [Step 4b](#4b-align) above if data are RRBS, or deduplicated bam file from [step 6a](#6a-deduplicate) if data are not RRBS and the bam file was deduplicated (e.g., sample-1_bismark_bt2.deduplicated.bam))
* a directory holding the reference genome in fasta format (this pipeline version uses the Ensembl fasta file indicated in the `fasta` column of the [GL-DPPD-7110-A_annotations.csv](../../GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv) GeneLab Annotations file))
> **NOTE**  
> If using RNA, files will include "bismark_hisat2" instead of "bismark_bt2" in the name.


**Output data:**

* **\*\_context\_\*.txt.gz** (bismark methylation-call files for CpG, CHG, and CHH contexts that were detected; see [bismark documentation](https://felixkrueger.github.io/Bismark/), namely [methylation call](http://felixkrueger.github.io/Bismark/bismark/alignment/#methylation-call) for symbols, and [methylation extraction output](http://felixkrueger.github.io/Bismark/bismark/methylation_extraction/#the-methylation-extractor-output-looks-like-this-tab-separated) for file format)
* **\*.bedGraph.gz** (gzip-compressed bedGraph-formatted file of methylation percentages of each CpG site; see [bismark documentation](https://felixkrueger.github.io/Bismark/options/methylation_extraction/#bedgraph-output))
* **\*.bismark.cov.gz** (gzip-compressed bedGraph-formatted file like above "\*.bedGraph.gz", but also including 2 more columns of methylated and unmethylated counts at the specified position; see [bismark documentation](https://felixkrueger.github.io/Bismark/options/methylation_extraction/#bedgraph-specific-options))
* **\*.M-bias.txt** (text file with methylation information in the context of the position in reads, helpful for investigating bias as a function of base position in the read; see [bismark documentation](http://felixkrueger.github.io/Bismark/bismark/methylation_extraction/#m-bias-plot))
* **\*_splitting_report.txt** (text file containing general methylation detection information)
* **\*.cytosine_context_summary.txt** (tsv file of detected cytosine-methylation information summed by nucleotide context)
* **\*.CpG_report.txt.gz** (a genome-wide methylation report for all CpG cytosines)

<br>

---

## 8. Generate individual sample report


```bash
bismark2report --dir sample-1_bismark_report_out_dir/ \
  --alignment_report sample-1_bismark_bt2*_report.txt \
  --nucleotide_report sample-1_bismark_bt2*.nucleotide_stats.txt \
  --splitting_report sample-1_bismark_bt2*_splitting_report.txt \
  --mbias_report sample-1_bismark_bt2*.M-bias.txt \
  --dedup_report sample-1_bismark_bt2*.deduplication_report.txt
```

**Parameter Definitions:**

* `--dir` - the output directory to store results
* `--alignment_report` - specifies the alignment report input file  
* `--nucleotide_report` - specifies the nucleotide stats report input file
* `--splitting_report` - specifies the splitting report input file  
* `--mbias_report` - specifies the methylation bias report input file 
* `--dedup_report` - specifies the deduplication report input file (optional, use only if deduplication was run)

**Input data:**

* sample-1_bismark_bt2*_report.txt (bismark mapping report file, output from [Step 4b.](#4b-align))
* sample-1_bismark_bt2*.nucleotide_stats.txt (bismark nucleotide stats report file, output from [Step 4b.](#4b-align))
* sample-1_bismark_bt2*_splitting_report.txt (splitting report file, output from [Step 7](#7-extract-methylation-calls) above)
* sample-1_bismark_bt2*.M-bias.txt (text file with methylation information in the context of the position in reads, output from [Step 7](#7-extract-methylation-calls) above)
* sample-1_bismark_bt2*.deduplication_report.txt (optional deduplication report, output from [Step 6a.](#6a-deduplicate) if deduplication was run)
> **NOTE**  
> If using RNA, files will include "bismark_hisat2" instead of "bismark_bt2" in the name.

**Output data:**

* **\*_report.html** (a bismark summary html file for the given sample)

<br>

---

## 9. Generate combined summary reports

### 9a. Bismark Summary Report

```bash
bismark2summary sample-1_bismark_bt2*.bam

#rename output files using standard GeneLab suffix
mv bismark_summary_report.txt bismark_summary_report_GLMethylSeq.txt
mv bismark_summary_report.html bismark_summary_report_GLMethylSeq.html

```

**Parameter Definitions:**

* `sample-1_bismark_bt2*.bam` - the input bam files are specified as a positional argument, and can be given all at once with wildcards, or as individual arguments with spaces in between them

**Input data:**  

* autodetects appropriate files based on the input bam files
  * the positional argument(s) can either be explicitly naming the bam file(s) as done above or can be the top-level directory holding the initial bam files and relevant report files
  * the autodetected files cannot be explicitly provided, but it looks for those named like these listed here and includes them if they exist for each individual starting bam file it is given or finds
    * sample-1_bismark_bt2*_report.txt generated from [Step 4b.](#4b-align) above
    * sample-1_bismark_bt2*_splitting_report.txt from [Step 7](#7-extract-methylation-calls) above
    * sample-1_bismark_bt2*.deduplication_report.txt if deduplication was performed in [Step 6a.](#6a-deduplicate)
> **NOTE**  
> If using RNA, files will include "bismark_hisat2" instead of "bismark_bt2" in the name.

**Output data:**  

* **bismark_summary_report_GLMethylSeq.txt** (summary table of general information on all included samples)
* **bismark_summary_report_GLMethylSeq.html** (html summary of general information on all included samples)

> **NOTE**
> If using RNA, file suffix will be **GLRNAMethylSeq** instead of **GLMethylSeq**

<br>

---

### 9b. Compile Alignment and Bismark QC

```bash
multiqc --interactive -o align_and_bismark_multiqc_data/ -n align_and_bismark_multiqc -z \
  qualimap_out_dir/ mapping_files_out_dir/ methylation_calls_out_dir/ deduplication_out_dir/
```

**Parameter Definitions:**

*	`--interactive` - force reports to use interactive plots
*	`-o` – the output directory to store results
*	`-n` – the filename prefix for output files
*	`-z` – specifies to zip the output data directory
*	`qualimap_out_dir/` – the directory holding the output data from the qualimap run, provided as a positional argument
*	`methylation_calls_out_dir/` – the directory holding the output data from the methylation extraction run, provided as a positional argument
*	`mapping_files_out_dir/` – the directory holding the output data from the alignment run, provided as a positional argument
*	`deduplication_out_dir/` – the directory holding the output data from the deduplication run, provided as a positional argument (omitted if RRBS data)

**Input data:**

* qualimap_out_dir/*.txt (qualimap output results from [Step 5](#5-alignment-qc))
* mapping_files_out_dir/*_[SP]E_report.txt (Bismark alignment results from [Step 4b](#4b-align))
* mapping_files_out_dir/*.nucleotide_stats.txt (Bismark nucleotide stats results from [Step 4b](#4b-align))
* deduplication_out_dir/*.deduplication_report.txt (Bismark deduplication results from [Step 6a](#6a-deduplicate) if not RRBS)
* methylation_calls_out_dir/*_splitting_report.txt (methylation extraction results from [Step 7](#7-extract-methylation-calls))
* methylation_calls_out_dir/*M-bias.txt (methylation bias results from [Step 7](#7-extract-methylation-calls))

**Output data:**

* **align_and_bismark_multiqc_GLmethylSeq.html** (multiqc output html summary)
* **align_and_bismark_multiqc_GLmethylSeq_data** (directory containing multiqc output data)

> If using RNA, file suffix will be **GLRNAMethylSeq** instead of **GLMethylSeq**

<br>


---

## 10. Generate reference genome annotation information

### 10a. GTF to BED conversion
A bed-formatted annotation file is needed for adding annotation information to the output from the differential methylation analysis. We utilize gtf files from [Ensembl](https://www.ensembl.org/) and convert them to the bed format as follows:

```bash
gtfToGenePred reference.gtf reference.genePred
genePredToBed reference.genePred reference.bed
```

**Parameter Definitions:**

* `reference.gtf` - genome annotation file in GTF format, provided as a positional argument
* `reference.genePred` - genome annotation file in genePred format, provided as a positional argument 
* `reference.bed` - genome annotation file in BED format, provided as a positional argument

**Input data:**

* reference.gtf (genome annotation, this pipeline version uses the Ensembl gtf file indicated in the `gtf` column of the [GL-DPPD-7110-A_annotations.csv](../../GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv) GeneLab Annotations file)

**Output data:**

* reference.genePred (intermediate genome annotation file in genePred format)
* reference.bed (genome annotation file in BED format)

<br>

---

### 10b. Making a mapping file of genes to transcripts
Making a mapping file of gene names to transcript names, which we need to link functional annotations in a primary output table. 

```bash
awk ' $3 == "transcript" ' reference.gtf | cut -f 9 | tr -s ";" "\t" | \
    cut -f 1,3 | tr -s " " "\t" | cut -f 2,4 | tr -d '"' \
    > reference-gene-to-transcript-map.tsv
```

**Parameter Definitions:**

* `reference.gtf` - genome annotation file in GTF format, provided as a positional argument
* `reference-gene-to-transcript-map.tsv` - output file, provided as a positional argument


**Input data:**

* reference.gtf (genome annotation, this pipeline version uses the Ensembl gtf file indicated in the `gtf` column of the [GL-DPPD-7110-A_annotations.csv](../../GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv) GeneLab Annotations file)

**Output data:**

* **reference-gene-to-transcript-map.tsv** (a gene-to-transcript mapping file with gene IDs in the first column and transcript IDs in the second)

<br>

---

## 11. Differential methylation analysis

<!-- Example data for the R code below can be downloaded and unpacked with the following:

```bash
curl -L -o test-methylkit-files.zip https://figshare.com/ndownloader/files/38765340
unzip test-methylkit-files.zip && rm test-methylkit-files.zip
``` -->

The remainder of this document is performed in R. 

<br>

---

### 11a. Set up R environment

```R
### Install R packages if not already installed ###

install.packages("remotes")
remotes::install_version("tidyverse", version = "2.0.0")

if (!require("BiocManager", quietly = TRUE)) 
    install.packages("BiocManager")
BiocManager::install(version = "3.19")
BiocManager::install("methylKit", version = "3.19")
BiocManager::install("genomation", version = "3.19")


### Import libraries (###

library(tidyverse)
library(methylKit)
library(genomation)

### Setting up variables and objects ###

# Define which organism is used in the study 
    # This should be consistent with the name in the "name" column of 
    # the GL-DPPD-7110-A_annotations.csv file 
organism <- "organism_that_samples_were_derived_from"

# Set path to the directory holding bismark coverage files (*.bismark.cov.gz)
coverage_files_dir_path <- file.path("bismark-coverage-files/")

# Set path to directory holding reference bed (reference.bed) and 
# transcript-to-gene-map (reference-gene-to-transcript-map.tsv) files
references_dir_path <- file.path("reference-files/")

# Create a list of unique sample names
sample.list <- list("sample-1", "sample-2", "sample-3",
                    "sample-4", "sample-5", "sample-6")

# Generate a vector holding input coverage files
input_cov_files_vec <- list.files(coverage_files_dir_path, 
                                  pattern = ".*.bismark.cov.gz$", 
                                  full.names = TRUE)

# Creating a list containing the detected input files and ensuring
# they are in the same order as the sample names 
file.list <- as.list(input_cov_files_vec[sapply(unlist(sample.list), 
                           function(x) grep(paste0(x,"_bismark.*.bismark.cov.gz$"), 
                                            input_cov_files_vec, value = FALSE))])

# Set variable for bed file
reference_bed_file <- list.files(references_dir_path, 
                                 pattern = ".*.bed$", full.names = TRUE)

# Set variable for transcript-to-gene-map file
gene_transcript_map_file <- list.files(references_dir_path, 
                                       pattern = ".*-gene-to-transcript-map.tsv$", 
                                       full.names = TRUE)

# Set variable for GL reference table
base_GL_github_link <- 
    "https://raw.githubusercontent.com/nasa/GeneLab_Data_Processing/master/"
ref_tab_location <-
    "GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv"
ref_tab_link <- paste0(base_GL_github_link, ref_tab_location)

# Read in reference table
ref_table <- read.csv(ref_tab_link)

# Set variable with organism and ensembl version number
ensembl_version <- ref_table %>% 
    filter(name == organism) %>% pull(ensemblVersion)

org_and_ensembl_version <- paste(organism, ensembl_version, sep = "_")

# Set variable with link to appropriate annotations table
annotations_tab_link <- ref_table %>% 
    filter(name == organism) %>% pull(genelab_annots_link)

# Read coverage files into R object and setting up contrast vector with the 'treatment' option
myobj <- methRead(location = file.list,
                  sample.id = sample.list,
                  assembly = org_and_ensembl_version,
                  pipeline = "bismarkCoverage",
                  header = FALSE,
                  treatment = c(1,1,1,0,0,0),
                  mincov = 10)

# Example of how to store as tables if memory requirements are too high
# myobj_storage <- methRead(location = file.list,
#                   sample.id = sample.list,
#                   assembly = org_and_ensembl_version,
#                   dbtype = "tabix",
#                   pipeline = "bismarkCoverage",
#                   header = FALSE,
#                   treatment = c(1,1,1,0,0,0),
#                   mincov = 10)
    # then 'myobj_storage' would be used in place of 'myobj' below
```

**Input data:**

* \*.bismark.cov.gz (gzip-compressed bedGraph-formatted files generated in [Step 7](#7-extract-methylation-calls))
* reference.bed (bed file generated in [Step 10a](#10a-gtf-to-bed-conversion))
* reference-gene-to-transcript-map.tsv (gene-to-transcript mapping file generated in [Step 10b](#10b-making-a-mapping-file-of-genes-to-transcripts))

**Output data:**

* `reference_bed_file` (variable holding the path to the reference.bed file)
* `gene_transcript_map_file` (variable holding the path to the \*gene-to-transcript-map.tsv file)
* `annotations_tab_link` (variable holding the link to the GeneLab reference annotations table)
* `myobj` (methylRawList object holding coverage information for all samples)

<br>

---

### 11b. Individual-base analysis

```R
# Merge samples
meth <- unite(myobj)

# Calculate differential methylation
myDiff <- calculateDiffMeth(meth, mc.cores = 4)

# Get all significantly different bases
myDiff.all_sig <- getMethylDiff(myDiff, difference = 25, qvalue = 0.01)

# Getting significantly hypermethylated bases
myDiff.hyper <- getMethylDiff(myDiff, difference = 25, qvalue = 0.01, type = "hyper")

# Get significantly hypomethylated bases
myDiff.hypo <- getMethylDiff(myDiff, difference = 25, qvalue = 0.01, type = "hypo")
```

**Input data:**

* `myobj` (methylRawList object created in [Step 11a](#11a-set-up-r-environment))

**Output data:**

* `meth` (methylBase object, at the base-level)
* `myDiff.all_sig` (object containing all significantly differentially methylated bases, based on q-value < 0.01)
* `myDiff.hyper` (object containing all significantly differentially hyper-methylated bases, based on q-value < 0.01)
* `myDiff.hypo` (object containing all significantly differentially hypo-methylated bases, based on q-value < 0.01)

<br>

---

### 11c. Tile-level analysis

```R
# Create tiles object
tiles_obj <- tileMethylCounts(myobj, win.size = 1000, step.size = 1000, cov.bases = 10)

# Merge tiled samples
tiles_meth <- unite(tiles_obj)

# Calculate differential methylation on tiles
tiles_diff <- calculateDiffMeth(tiles_meth, mc.cores = 4)

# Getting all significantly different tiles
tiles_myDiff.all_sig <- getMethylDiff(tiles_diff, difference = 25, qvalue = 0.01)

# Get significantly hypermethylated tiles
tiles_myDiff.hyper <- getMethylDiff(tiles_diff, difference = 25, qvalue = 0.01, type = "hyper")

# Get significantly hypomethylated tiles
tiles_myDiff.hypo <- getMethylDiff(tiles_diff, difference = 25, qvalue = 0.01, type = "hypo")
```

**Input data:**

* `myobj` (methylRawList object created in [Step 11a](#11a-set-up-r-environment))

**Output data:**

* `tiles_meth` (methylBase object, at the tile level)
* `tiles_myDiff.all_sig` (object containing all significantly differentially methylated titles, based on q-value < 0.01)
* `tiles_myDiff.hyper` (object containing all significantly differentially hyper-methylated tiles, based on q-value < 0.01)
* `tiles_myDiff.hypo` (object containing all significantly differentially hypo-methylated tiles, based on q-value < 0.01)

<br>

---

### 11d. Add feature information

```R
# Read in reference bed file
gene.obj <- readTranscriptFeatures(reference_bed_file, up.flank = 1000, 
                                   down.flank = 1000, remove.unusual = TRUE, 
                                   unique.prom = TRUE)

## Add features to individual-base objects
diffAnn.all_sig <- annotateWithGeneParts(as(myDiff.all_sig, "GRanges"), gene.obj)
diffAnn.hyper <- annotateWithGeneParts(as(myDiff.hyper, "GRanges"), gene.obj)
diffAnn.hypo <- annotateWithGeneParts(as(myDiff.hypo, "GRanges"), gene.obj)

# Make base-level significance table with features 
sig_all_bases_tab_with_features <- cbind(data.frame(myDiff.all_sig), 
                                         getAssociationWithTSS(diffAnn.all_sig), 
                                         as.data.frame(getMembers(diffAnn.all_sig))) %>% 
                                             .[, !names(.) %in% c("target.row")]

# Make base-level significantly hyper-methylated table with features
sig_bases_hyper_tab_with_features <- cbind(data.frame(myDiff.hyper), 
                                           getAssociationWithTSS(diffAnn.hyper), 
                                           as.data.frame(getMembers(diffAnn.hyper))) %>% 
                                               .[, !names(.) %in% c("target.row")]

# Make base-level significantly hypo-methylated table with features
sig_bases_hypo_tab_with_features <- cbind(data.frame(myDiff.hypo), 
                                          getAssociationWithTSS(diffAnn.hypo), 
                                          as.data.frame(getMembers(diffAnn.hypo))) %>% 
                                              .[, !names(.) %in% c("target.row")]


## Add features to tile objects
tiles_diffAnn.all_sig <- annotateWithGeneParts(as(tiles_myDiff.all_sig, "GRanges"), gene.obj)
tiles_diffAnn.hyper <- annotateWithGeneParts(as(tiles_myDiff.hyper, "GRanges"), gene.obj)
tiles_diffAnn.hypo <- annotateWithGeneParts(as(tiles_myDiff.hypo, "GRanges"), gene.obj)


# Make tiles significance table with features 
tiles_sig_all_out_tab_with_features <- cbind(data.frame(tiles_myDiff.all_sig), 
                                             getAssociationWithTSS(tiles_diffAnn.all_sig), 
                                             as.data.frame(getMembers(tiles_diffAnn.all_sig))) %>% 
                                                 .[, !names(.) %in% c("target.row")]
                                               
# Make tiles significantly hyper-methylated table with features
tiles_sig_hyper_tab_with_features <- cbind(data.frame(tiles_myDiff.hyper), 
                                           getAssociationWithTSS(tiles_diffAnn.hyper), 
                                           as.data.frame(getMembers(tiles_diffAnn.hyper))) %>% 
                                               .[, !names(.) %in% c("target.row")]

# Make tiles significantly hypo-methylated table with features
tiles_sig_hypo_tab_with_features <- cbind(data.frame(tiles_myDiff.hypo), 
                                          getAssociationWithTSS(tiles_diffAnn.hypo), 
                                          as.data.frame(getMembers(tiles_diffAnn.hypo))) %>% 
                                              .[, !names(.) %in% c("target.row")]
```

**Input data:**

* `reference_bed_file` (variable holding the path to the reference bed file created in [Step 11a](#11a-set-up-r-environment))
* `myDiff.all_sig`, `myDiff.hyper`, `myDiff.hypo` (base-level methylDiff objects created in [Step 11b](#11b-individual-base-analysis))
* `tiles_myDiff.all_sig`, `tiles_myDiff.hyper`, `tiles_myDiff.hypo` (tile-level methylDiff objects created in [Step 11c](#11c-tile-level-analysis))

**Output data:**

* `diffAnn.all_sig` (AnnotationByGeneParts object for all significantly differentially methylated bases)
* `tiles_diffAnn.all_sig` (AnnotationByGeneParts object for significantly differentially methylated tiles)
* `sig_all_bases_tab_with_features` (dataframe object holding all significantly differentially methylated bases and associated features)
* `sig_bases_hyper_tab_with_features` (dataframe object holding all significantly differentially hyper-methylated bases and associated features)
* `sig_bases_hypo_tab_with_features` (dataframe object holding all significantly differentially hypo-methylated bases and associated features)
* `tiles_sig_all_out_tab_with_features` (dataframe object holding all significantly differentially methylated tiles and associated features)
* `tiles_sig_hyper_tab_with_features` (dataframe object holding all significantly differentially hyper-methylated tiles and associated features)
* `tiles_sig_hypo_tab_with_features` (dataframe object holding all significantly differentially hypo-methylated tiles and associated features)

<br>

---

### 11e. Add functional annotations

```R
# Adjust time-out to enable downloading reference annotation table
options(timeout = 600)

# Read in annotations table
functional_annots_tab <- read.table(annotations_tab_link, 
                                    sep = "\t", quote = "", 
                                    header = TRUE)

# Read in gene-to-transcript mapping file
gene_transcript_map <- read.table(gene_transcript_map_file, sep = "\t", 
                                  col.names = c("gene_ID", "feature.name"))

## Add annotations to individual-base objects ##
    # Retrieving the corresponding gene ID for each transcript ID in the
    # *_bases*tab_with_features tables, and adding it to the table

sig_all_bases_tab_with_features_and_gene_IDs <- 
    left_join(sig_all_bases_tab_with_features, gene_transcript_map)

sig_hyper_bases_tab_with_features_and_gene_IDs <-
    left_join(sig_bases_hyper_tab_with_features, gene_transcript_map)

sig_hypo_bases_tab_with_features_and_gene_IDs <-
    left_join(sig_bases_hypo_tab_with_features, gene_transcript_map)

    # now adding functional annotations and renaming column of output tables 
    # to be primary keytype instead of "gene_ID"
sig_all_bases_tab_with_features_and_annots <- 
    left_join(sig_all_bases_tab_with_features_and_gene_IDs, 
              functional_annots_tab, by = c("gene_ID" = "ENSEMBL"))
              
sig_all_bases_tab_with_features_and_annots <- 
    rename(sig_all_bases_tab_with_features_and_annots, "gene_ID" = "ENSEMBL")


sig_hyper_bases_tab_with_features_and_annots <- 
    left_join(sig_hyper_bases_tab_with_features_and_gene_IDs, 
              functional_annots_tab, by = c("gene_ID" = "ENSEMBL"))

sig_hyper_bases_tab_with_features_and_annots <- 
    rename(sig_hyper_bases_tab_with_features_and_annots, "gene_ID" = "ENSEMBL")


sig_hypo_bases_tab_with_features_and_annots <- 
    left_join(sig_hypo_bases_tab_with_features_and_gene_IDs, 
              functional_annots_tab, by = c("gene_ID" = "ENSEMBL"))

sig_hypo_bases_tab_with_features_and_annots <- 
    rename(sig_hypo_bases_tab_with_features_and_annots, "gene_ID" = "ENSEMBL")

# Write out tables
write.table(sig_all_bases_tab_with_features_and_annots, 
            "sig-diff-methylated-bases.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(sig_hyper_bases_tab_with_features_and_annots, 
            "sig-diff-hypermethylated-bases.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(sig_hypo_bases_tab_with_features_and_annots, 
            "sig-diff-hypomethylated-bases.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)

## Add annotations to tiles objects ##
    # Retrieving the corresponding gene ID for each transcript ID in the
    # *tiles*tab_with_features tables, and adding it to the table

sig_all_tiles_tab_with_features_and_gene_IDs <- 
    left_join(tiles_sig_all_out_tab_with_features, gene_transcript_map)

sig_hyper_tiles_tab_with_features_and_gene_IDs <-
    left_join(tiles_sig_hyper_tab_with_features, gene_transcript_map)

sig_hypo_tiles_tab_with_features_and_gene_IDs <-
    left_join(tiles_sig_hypo_tab_with_features, gene_transcript_map)

    # now adding functional annotations and renaming column of output tables 
    # to be primary keytype instead of "gene_ID"
sig_all_tiles_tab_with_features_and_annots <- 
    left_join(sig_all_tiles_tab_with_features_and_gene_IDs, 
              functional_annots_tab, by = c("gene_ID" = "ENSEMBL"))

sig_all_tiles_tab_with_features_and_annots <- 
    rename(sig_all_tiles_tab_with_features_and_annots, "gene_ID" = "ENSEMBL")


sig_hyper_tiles_tab_with_features_and_annots <- 
    left_join(sig_hyper_tiles_tab_with_features_and_gene_IDs, 
              functional_annots_tab, by = c("gene_ID" = "ENSEMBL"))

sig_hyper_tiles_tab_with_features_and_annots <- 
    rename(sig_hyper_tiles_tab_with_features_and_annots, "gene_ID" = "ENSEMBL")


sig_hypo_tiles_tab_with_features_and_annots <- 
    left_join(sig_hypo_tiles_tab_with_features_and_gene_IDs, 
              functional_annots_tab, by = c("gene_ID" = "ENSEMBL"))

sig_hypo_tiles_tab_with_features_and_annots <- 
    rename(sig_hypo_tiles_tab_with_features_and_annots, "gene_ID" = "ENSEMBL")


# Write out tables
write.table(sig_all_tiles_tab_with_features_and_annots, 
            "sig-diff-methylated-tiles.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(sig_hyper_tiles_tab_with_features_and_annots, 
            "sig-diff-hypermethylated-tiles.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(sig_hypo_tiles_tab_with_features_and_annots, 
            "sig-diff-hypomethylated-tiles.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)
```
**Input data:**


* `annotations_tab_link` (variable holding the link to the GeneLab reference annotations table created in [Step 11a](#11a-set-up-r-environment))
* `gene_transcript_map_file` (variable holding the path to the \*gene-to-transcript-map.tsv file created in [Step 11a](#11a-set-up-r-environment))
* `sig_all_bases_tab_with_features`, `sig_bases_hyper_tab_with_features`, `sig_bases_hypo_tab_with_features` (dataframe objects holding significantly differentially methylated bases and associated features created in [Step 11d](#11d-add-feature-information))
* `tiles_sig_all_out_tab_with_features`, `tiles_sig_hyper_tab_with_features`, `tiles_sig_hypo_tab_with_features` (dataframe objects holding significantly differentially methylated tiles and associated features created in [Step 11d](#11d-add-feature-information))

**Output data:**

* **\*sig-diff-methylated-bases.tsv** (table containing all significantly differentially methylated cytosines)
* **\*sig-diff-hypermethylated-bases.tsv** (table containing all cytosines with significantly elevated methylation levels)
* **\*sig-diff-hypomethylated-bases.tsv** (table containing all cytosines with significantly reduced methylation levels)
* **\*sig-diff-methylated-tiles.tsv** (table containing all significantly differentially methylated tiles)
* **\*sig-diff-hypermethylated-tiles.tsv** (table containing all tiles with significantly elevated methylation levels)
* **\*sig-diff-hypomethylated-tiles.tsv** (table containing all tiles with significantly reduced methylation levels)

> NOTE:  
> In practice these will be based on a contrast listed in the filename, e.g., they will be something like "Flight_vs_Ground-sig-diff-methylated-bases.tsv".

<br>

---

### 11f. Make and write out tables of percent methylation levels

```R
# Make and write out a table of base-level percent methylated
perc.meth <- percMethylation(meth, rowids = TRUE)
perc.meth <- perc.meth %>% data.frame(check.names = FALSE) %>% rownames_to_column("location")

write.table(perc.meth, "base-level-percent-methylated.tsv", sep = "\t", 
            quote = FALSE, row.names = FALSE)

# Make and write out a table of tile percent methylated
tiles_perc.meth <- percMethylation(tiles_meth, rowids = TRUE)
tiles_perc.meth <- tiles_perc.meth %>% data.frame(check.names = FALSE) %>% rownames_to_column("location")

write.table(tiles_perc.meth, "tile-level-percent-methylated.tsv", sep = "\t", 
            quote = FALSE, row.names = FALSE)
```

**Input data:**
* `meth` (base-level methylBase object created in [Step 11b](#11b-individual-base-analysis))
* `tiles_meth` (tile-level methylBase object created in [Step 11c](#11c-tile-level-analysis))

**Output data:**

* **base-level-percent-methylated.tsv** (table containing the percent methylation levels of all cytosines across all samples)
* **tile-level-percent-methylated.tsv** (table containing the percent methylation levels of all tiles across all samples)

<br>

---

### 11g. Make overview figure of percent differential methylation across features

```R
# Based on individual bases
pdf("sig-diff-methylated-bases-across-features.pdf")
plotTargetAnnotation(diffAnn.all_sig, precedence = TRUE, 
                     main = "% of sig. diff. methylated sites across features")
dev.off()

# Based on tiles
pdf("sig-diff-methylated-tiles-across-features.pdf")
plotTargetAnnotation(tiles_diffAnn.all_sig, precedence = TRUE, 
                     main = "% of sig. diff. methylated tiles across features")
dev.off()
```

**Input data:**

* `diffAnn.all_sig`, `tiles_diffAnn.all_sig` (AnnotationByGeneParts objects for all differentially methylated bases and tiles created in [Step 11d](#11d-add-feature-information))

**Output data:**

* **\*sig-diff-methylated-bases-across-features.pdf** (overview figure of the percent of significantly differentially methylated cytosines identified in specific features (promoter, exon, intron))
* **\*sig-diff-methylated-tiles-across-features.pdf** (overview figure of the percent of significantly differentially methylated tiles identified in specific features (promoter, exon, intron))

> NOTE:  
> In practice these will be based on a contrast listed in the filename, e.g., they will be something like "Flight_vs_Ground-sig-diff-methylated-bases-across-features.pdf".

<br>

---
---
