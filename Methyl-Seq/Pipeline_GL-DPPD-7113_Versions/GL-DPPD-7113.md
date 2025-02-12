# Bioinformatics pipeline for Methylation Sequencing (Methyl-Seq) data

> **This document holds an overview and some example commands of how GeneLab processes bisulfite sequencing (methylseq) 
  datasets. Exact processing commands for specific datasets that have been released are provided with their processed 
  data in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/).**

---

**Date:** February 13, 2025  
**Revision:** -  
**Document Number:** GL-DPPD-7113  

**Submitted by:**  
Michael D. Lee and Barbara Novak (GeneLab Data Processing Team)  

**Approved by:**  
Barbara Novak (GeneLab Data Processing Lead)  
Amanda Saravia-Butler (GeneLab Science Lead)  
Lauren Sanders (OSDR Project Scientist)  
Danielle Lopez (OSDR Deputy Project Manager)  
Samrawit Gebre (OSDR Project Manager)  

---

# Table of contents

- [**Software used**](#software-used)
- [**General processing overview with example commands**](#general-processing-overview-with-example-commands)
  - [**1. Raw Data QC**](#1-raw-data-qc)
    - [1a. Raw Data QC](#1a-raw-data-qc)
    - [1b. Compile Raw Data QC](#1b-compile-raw-data-qc)
  - [**2. Adapter trimming/quality filtering**](#2-adapter-trimmingquality-filtering)
    - [If not RRBS or if RRBS using MseI digestion](#if-not-rrbs-or-if-rrbs-using-msei-digestion)
    - [If using a random priming bisulfite method](#if-using-a-random-priming-post-bisulfite-method)
    - [If RRBS with MspI digestion](#if-rrbs-with-mspi-digestion)
    - [If RRBS with NuGEN ovation kit](#if-rrbs-with-nugen-ovation-kit)
      - [First adapter-trimming/quality-filtering with trimgalore](#first-adapter-trimmingquality-filtering-with-trimgalore)
      - [Now running NuGEN-specific script](#now-running-nugen-specific-script)
  - [**3. Filtered/Trimmed Data QC**](#3-filteredtrimmed-data-qc)
    - [3a. Trimmed Data QC](#3a-trimmed-data-qc)
    - [3b. Compile Trimmed Data QC](#3b-compile-trimmed-data-qc)
  - [**4. Alignment**](#4-alignment)
    - [4a. Generate reference](#4a-generate-reference)
    - [4b. Align](#4b-align)
    - [4c. Sort Alignment Files](#4c-sort-alignment-files)
  - [**5. Alignment QC**](#5-alignment-qc)
  - [**6. Deduplicate (skip if data are RRBS)**](#6-deduplicate-skip-if-data-are-rrbs)
    - [6a. Deduplicate](#6a-deduplicate)
    - [6b. Sort Deduplicated Alignment Files](#6b-sort-deduplicated-alignment-files)
  - [**7. Extract methylation calls**](#7-extract-methylation-calls)
  - [**8. Generate individual sample report**](#8-generate-individual-sample-report)
  - [**9. Generate combined summary reports**](#9-generate-combined-summary-reports)
    - [9a. Bismark summary report](#9a-bismark-summary-report)
    - [9b. Compile Alignment and Bismark QC](#9b-compile-alignment-and-bismark-qc)
  - [**10. Generate reference genome annotation information**](#10-generate-reference-genome-annotation-information)
    - [10a. GTF to BED conversion](#10a-gtf-to-bed-conversion)
    - [10b. Making a mapping file of genes to transcripts](#10b-making-a-mapping-file-of-genes-to-transcripts)
  - [**11. Differential methylation analysis**](#11-differential-methylation-analysis)
    - [11a. Set up R environment](#11a-set-up-r-environment)
    - [11b. Configure Metadata, Sample Grouping, and Group Comparisons](#11b-configure-metadata-sample-grouping-and-group-comparisons)
    - [11c. Import Methylation Calls](#11c-import-methylation-calls)
    - [11d. Individual-base analysis](#11d-individual-base-analysis)
    - [11e. Tile analysis](#11e-tile-analysis)
    - [11f. Export Tables](#11f-export-tables)


---

# Software used

|Program|Version|Relevant Links|
|:------------|:------:|------:|
|FastQC       | 0.12.1 |[https://www.bioinformatics.babraham.ac.uk/projects/fastqc/](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)|
|MultiQC      | 1.26   |[https://multiqc.info/](https://multiqc.info/)|
|Cutadapt     | 4.8    |[https://cutadapt.readthedocs.io/en/stable/](https://cutadapt.readthedocs.io/en/stable/)|
|TrimGalore!  | 0.6.10 |[https://github.com/FelixKrueger/TrimGalore](https://github.com/FelixKrueger/TrimGalore)|
|Bismark      | 0.24.2 |[https://github.com/FelixKrueger/Bismark](https://github.com/FelixKrueger/Bismark)|
|bowtie2      | 2.5.4  |[https://github.com/BenLangmead/bowtie2#overview](https://github.com/BenLangmead/bowtie2#overview)|
|hisat2       | 2.2.1  |[https://github.com/DaehwanKimLab/hisat2](https://github.com/DaehwanKimLab/hisat2)|
|samtools     | 1.21   |[https://github.com/samtools/samtools#samtools](https://github.com/samtools/samtools#samtools)|
|qualimap     | 2.3    |[http://qualimap.conesalab.org/](http://qualimap.conesalab.org/)|
|gtfToGenePred| 469    |[http://hgdownload.cse.ucsc.edu/admin/exe/](http://hgdownload.cse.ucsc.edu/admin/exe/)|
|genePredToBed| 469    |[http://hgdownload.cse.ucsc.edu/admin/exe/](http://hgdownload.cse.ucsc.edu/admin/exe/)
|R            | 4.4.2  |[https://www.r-project.org](https://www.r-project.org)|
|tidyverse    | 2.0.0  |[https://tidyverse.tidyverse.org/](https://tidyverse.tidyverse.org/)|
|Bioconductor | 3.20   |[https://bioconductor.org](https://bioconductor.org)
|methylKit    | 1.32.0 |[https://bioconductor.org/packages/release/bioc/html/methylKit.html](https://bioconductor.org/packages/release/bioc/html/methylKit.html)|
|genomation   | 1.38.0 |[https://bioconductor.org/packages/release/bioc/html/genomation.html](https://bioconductor.org/packages/release/bioc/html/genomation.html)|

---

# General processing overview with example commands

> Exact processing commands and output files listed in **bold** below are included with each Methyl-Seq processed dataset in 
  the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/).

---

## 1. Raw Data QC

<br>

### 1a. Raw Data QC

```bash
fastqc -o raw_fastqc_output/ *raw.fastq.gz
```

**Parameter Definitions:**

* `-o` – the output directory to store results
* `*raw.fastq.gz` – the input reads are specified as a positional argument, and can be given all at once with wildcards, 
    or as individual arguments with spaces in between them

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

> If processing RNA Methylation Sequencing data, file suffix will be **GLRNAMethylSeq** instead of **GLMethylSeq**

<br>  

---

## 2. Adapter trimming/quality filtering
See `trim_galore --help` or [TrimGalore User Guide](https://github.com/FelixKrueger/TrimGalore/blob/0.6.10/Docs/Trim_Galore_User_Guide.md) 
for more info on any of the below.

Additionally, the Bismark documentation also includes guidelines for specific MethylSeq library types: 
[Bismark library type guide](http://felixkrueger.github.io/Bismark/bismark/library_types/). Some library types will 
require additional 5' and/or 3' hard trimming to remove the signature of the oligos used for random priming. Leaving 
these bases may cause misalignments and methylation biases.

<br>

### If not RRBS or if RRBS using MseI digestion
Note that the `--rrbs` option is **not** appropriate when RRBS (reduced representation bisulfite sequencing) libraries 
were prepared with MseI digestion (see the TrimGalore User Guide 
[Note for RRBS using MseI](https://github.com/FelixKrueger/TrimGalore/blob/0.6.10/Docs/Trim_Galore_User_Guide.md#rrbs-specific-options-mspi-digested-material).

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
Random priming is not truly random and the signature left at the ends of the reads can introduce errors, indels, and 
methylation biases. Add the optional clipping parameters (`--clip_r1`, `--clip_r2`, `--three_prime_clip_r1`, and 
`--three_prime_clip_r2`) to trim off the random priming signature on the 5' ends of each read and next to the 3' end 
after adapter trimming. See [Bismark library type guide](http://felixkrueger.github.io/Bismark/bismark/library_types/) 
for more detailed information. 

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
Note that if the library preparation was non-directional, the `--non_directional` flag needs to be added to this command 
(whether single-end or paired-end; see 
[TrimGalore User Guide](https://github.com/FelixKrueger/TrimGalore/blob/0.6.10/Docs/Trim_Galore_User_Guide.md#rrbs-specific-options-mspi-digested-material)). 

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
Libraries prepared with the NuGEN ovation kit need to be processed with an additional script provided by the company's 
[github](https://github.com/nugentechnologies/NuMetRRBS#analysis-guide-for-nugen-ovation-rrbs-methyl-seq). 

Following their instructions, we first run an adapter-trimming/quality-filtering step with trimgalore. Note that the 
`--rrbs` option is not appropriate to pass to trimgalore when this kit is used (see Bismark documentation for 
[RRBS NuGEN Ovation Methyl-Seq System](http://felixkrueger.github.io/Bismark/bismark/library_types/#rrbs-nugen-ovation-methyl-seq-system). 
Then we utilize the company's script to remove the random diversity sequences added by the kit. 

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

The NuGEN-specific script can be downloaded from their 
[github](https://github.com/nugentechnologies/NuMetRRBS#analysis-guide-for-nugen-ovation-rrbs-methyl-seq) 
with the following command:

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
* `--cores` - number of cores to use (this value is dependent on the number of threads available on the system 
  running trim galore)
* `--phred33` - instructs cutadapt to use ASCII+33 quality scores as Phred scores for quality trimming
* `--rrbs` - specific trimming suitable for RRBS data generated with MspI digestion only
* `-a` - specific adapter sequence to be trimmed off of forward or single reads (applicable for libraries prepared with 
  the NuGEN ovation kit)
* `-a2` - specific adapter sequence to be trimmed off of reverse reads (applicable for libraries prepared with the 
  NuGEN ovation kit)
* `--paired` - specifies data are paired-end
* `--output_dir` - the output directory to store results
* `--clip_R1` - number of bases to trim off the 5' end of each R1 read (optional, for use with library prep kits that 
  use random priming, such as TruSeq(EpiGnome))
* `--clip_R2` - number of bases to trim off the 5' end of each R2 read (optional, for use with library prep kits that use 
  random priming, such as TruSeq(EpiGnome))
* `--three_prime_clip_R1` - number of bases to trim off the 3' end of each R1 read AFTER adapter trimming. (optional, 
  for use with library prep kits that use random priming, such as TruSeq(EpiGnome)) 
* `--three_prime_clip_R2` - number of bases to trim off the 3' end of each R2 read AFTER adapter trimming. (optional, 
  for use with library prep kits that use random priming, such as TruSeq(EpiGnome)) 
* `sample-1_raw.fastq.gz` or `sample-1_R1_raw.fastq.gz sample-1_R2_raw.fastq.gz` – the input reads are specified as 
  positional arguments, paired-end read files are listed pairwise such that the forward reads (\*R1_raw.fastq.gz) are 
  immediately followed by the respective reverse reads (\*R2_raw.fastq.gz)


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
*	`*trimmed.fastq.gz` – the input reads are specified as a positional argument, and can be given all at once with 
  wildcards as in the example, or as individual arguments with spaces in between them  

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
  trimmed_fastqc_output/ trimmed_reads_out_dir/
```

**Parameter Definitions:**

*	`--interactive` - force reports to use interactive plots
*	`-o` – the output directory to store results
*	`-n` – the filename prefix for output files
*	`-z` – specifies to zip the output data directory
*	`trimmed_fastqc_output/` – the directory holding the output data from the fastqc run, provided as a positional argument
* `trimmed_reads_out_dir/` - the directory holding the trimgalore trimming reports, provided as a positional argument

**Input data:**

* trimmed_fastqc_output/*fastqc.zip (FastQC output data from [Step 3a](#3a-trimmed-data-qc))
* trimmed_reads_out_dir/*trimming_report.txt (TrimGalore! trimming report, output from [Step 2](#2-adapter-trimmingquality-filtering))

**Output data:**

* **trimmed_multiqc_GLMethylSeq.html** (multiqc output html summary)
* **trimmed_multiqc_GLMethylSeq_data.zip** (directory containing multiqc output data)

> If processing RNA Methylation Sequencing data, file suffix will be **GLRNAMethylSeq** instead of **GLMethylSeq**

<br>

---

## 4. Alignment

<br>

### 4a. Generate reference
The reference will need to be specific to the organism that was sequenced. Bismark operates on a directory holding the 
target reference genome in fasta format.

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
*  `bismark_reference_genome/` - positional argument specifying the directory holding the reference genome (should end 
  in ".fa" or ".fasta", can be gzipped and including ".gz")

*bam2nuc*
* --genomic_composition_only - specifies creation of the (genome-specific) genomic_nucleotide_frequencies.txt report  
* --genome_folder - specifies the directory holding the reference genome (should end in ".fa" or ".fasta", can be 
  gzipped and including ".gz")

 
**Input data:**

* a directory holding the reference genome in fasta format (this pipeline version uses the Ensembl fasta file indicated 
  in the `fasta` column of the 
  [GL-DPPD-7110-A_annotations.csv](../../GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv) 
  GeneLab Annotations file)

**Output data:**

* the reference genome directory that was provided as input will now hold indexes for the bisulfite-converted reference 
  genome (all `*.bt2` files are indexes, all `*.fa` files are converted versions of the reference genome)
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
* **bismark_reference_genome/genomic_nucleotide_frequencies.txt** (tab-delimited table of mono- and di-nucleotide 
  frequencies in reference genome)



> **NOTE**  
> If processing RNA Methylation Sequencing data, the `--hisat2` flag is added instead of `--bowtie2`, which specifies using the splice-aware aligner 
  [HISAT2](https://github.com/DaehwanKimLab/hisat2#hisat2), and the outputs will include 8 `*ht2` files in separate 
  sub-directories along with each reference-genome conversion.

<br>

### 4b. Align

> **NOTE**  
> If the library preparation was non-directional, the `--non_directional` flag needs to be added to this command 
  (whether single-end or paired-end). For a full list of alignment option recommendations library type and/or commercially 
  available kit, please see the library page in the [Bismark documentation](http://felixkrueger.github.io/Bismark/bismark/library_types/) 

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
* `--nucleotide_coverage` - outputs a table with mono- and di-nucleotide sequence compositions and coverage values 
    compared to genomic compositions
* `--gzip` - write temporary bisulfite conversion files in gzip format to save disk space during alignment
* `--output_dir` - the output directory to store results
* `--genome_folder` - specifies the directory holding the reference genome indexes (the same that was provided in 
    [Step 4a.](#4a-generate-reference) above)
* `-1` - for paired-end data, the forward trimmed reads (not used for single-end data)
* `-2` - for paired-end data, the reverse trimmed reads (not used for single-end data)
* `sample-1_trimmed.fastq.gz` - for single-end data, input trimmed-reads are provided as a positional argument


**Input data:**
* bismark_reference_genome/ (directory holding indexes of reference genome, output from [Step 4a](#4a-generate-reference))
* gzip-compressed fastq files (adapter-trimmed/quality-filtered reads, output from [Step 2](#2-adapter-trimmingquality-filtering))

**Output data:**  

* sample-1_bismark_bt2*.bam (mapping file) 
* **\*_[SP]E_report.txt** (bismark mapping report file)
* **\*.nucleotide_stats.txt** (tab-delimited table with sample-specific mono- and di-nucleotide sequence compositions 
  and coverage values compared to genomic compositions)


> **NOTE**  
> If processing RNA Methylation Sequencing data, command-line parameters `--hisat2` and `--path_to_hisat2` need to be added. Output file names will include 
  "bismark_hisat2" instead of "bismark_bt2".

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
> If processing RNA Methylation Sequencing data, file names will include "bismark_hisat2" instead of "bismark_bt2".

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
* `--collect-overlap-pairs` - instructs the program to output statistics of overlapping paired-end reads (if data were 
    paired-end, no effect if single-end)
* `--java-mem-size=6G` - specifies the amount of memory to use (here this is set to 6G; see 
    [qualimap FAQ here](http://qualimap.conesalab.org/doc_html/faq.html?highlight=java-mem-size))
* `-nt` - specifies the number of threads to use

**Input data:**

* sample-1_bismark_bt2_sorted.bam (bismark bowtie2 alignment bam file sorted by chromosomal coordinates, output from 
  [Step 4c](#4c-sort-alignment-files) above)
* a feature file containing regions of interest for the reference genome in gtf format (this pipeline version uses the 
  Ensembl fasta file indicated in the `gtf` column of the 
  [GL-DPPD-7110-A_annotations.csv](../../GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv) 
  GeneLab Annotations file)

**Output data:** 

* **\*sample-1_bismark_bt2_qualimap/** (subdirectory of many alignment QC output files and formatting files for presenting in an html file (see [qualimap documentation](http://qualimap.conesalab.org/doc_html/analysis.html#output))

> **NOTE**  
> If processing RNA Methylation Sequencing data, file names will include "bismark_hisat2" instead of "bismark_bt2".

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
> If processing RNA Methylation Sequencing data, file names will include "bismark_hisat2" instead of "bismark_bt2".

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
> If processing RNA Methylation Sequencing data, file names will include "bismark_hisat2" instead of "bismark_bt2".

**Output data:**

* **sample-1_bismark_bt2\*.deduplicated_sorted.bam** (bismark bowtie2 alignment bam file sorted by chromosomal coordinates, with duplicates removed)

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

* `--parallel` - specifies the number of cores to use for methylation extraction, note: the program will utilize ~3X the 
  number specified 
* `--bedGraph` - instructs the program to generate a sorted bedGraph file that reports the position of a given cytosine 
  and its methlyation state (by default, only methylated CpGs are reported - see bedgraph options in 
  [bismark documentation](https://felixkrueger.github.io/Bismark/options/methylation_extraction/#bedgraph-specific-options) 
  for more info)
* `--gzip` - specifies to gzip-compress the methylation extractor output files
* `--comprehensive` - specifies to merge all four possible strand-specific methylation info into context-dependent output files
* `--output_dir` - the output directory to store results
* `--cytosine_report` - instructions the program to produce a genome-wide methylation report for all cytosines in the genome
* `--genome_folder` - a directory holding the reference genome in fasta format (this pipeline version uses the Ensembl 
  fasta file indicated in the `fasta` column of the 
  [GL-DPPD-7110-A_annotations.csv](../../GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv) 
  GeneLab Annotations file)
* `--ignore_r2` - specifies how many bases to ignore from the 5' end of the reverse reads (bismark docs recommend 2, see 
  [bismark documentation](https://felixkrueger.github.io/Bismark/options/methylation_extraction/#options))
  > **Note:** The first several bases in the reverse read of bisulfite sequence experiments show a severe bias towards 
    non-methylation as a result of end-repairing sonicated fragments with unmethylated cytosines. It is 
    recommended to remove the first few basepairs.
* `--ignore_3prime_r2` - specifies how many bases to ignore from the 3' end of the reverse reads to remove unwanted biases 
  from the end of reads (For specific recommendations see Bismark documentation on 
  [Library Types](https://felixkrueger.github.io/Bismark/bismark/library_types/)) 
* `sample-1_bismark_bt2*.bam` - the input bam file, provided as a positional argument

**Input data:**

* sample-1_bismark_bt2*.bam (bismark bowtie2 alignment bam file, output from [Step 4b](#4b-align) above if data are RRBS, 
  or deduplicated bam file from [step 6a](#6a-deduplicate) if data are not RRBS and the bam file was deduplicated 
  (e.g., sample-1_bismark_bt2.deduplicated.bam))
* a directory holding the reference genome in fasta format (this pipeline version uses the Ensembl fasta file indicated 
  in the `fasta` column of the 
  [GL-DPPD-7110-A_annotations.csv](../../GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv) 
  GeneLab Annotations file)
> **NOTE**  
> If processing RNA Methylation Sequencing data, file names will include "bismark_hisat2" instead of "bismark_bt2".


**Output data:**

* **\*\_context\_\*.txt.gz** (bismark methylation-call files for CpG, CHG, and CHH contexts that were detected; see 
  [bismark documentation](https://felixkrueger.github.io/Bismark/), namely 
  [methylation call](http://felixkrueger.github.io/Bismark/bismark/alignment/#methylation-call) for symbols, and 
  [methylation extraction output](http://felixkrueger.github.io/Bismark/bismark/methylation_extraction/#the-methylation-extractor-output-looks-like-this-tab-separated) 
  for file format)
* **\*.bedGraph.gz** (gzip-compressed bedGraph-formatted file of methylation percentages of each CpG site; see 
  [bismark documentation](https://felixkrueger.github.io/Bismark/options/methylation_extraction/#bedgraph-output))
* **\*.bismark.cov.gz** (gzip-compressed bedGraph-formatted file like above "\*.bedGraph.gz", but also including 2 more 
  columns of methylated and unmethylated counts at the specified position; see 
  [bismark documentation](https://felixkrueger.github.io/Bismark/options/methylation_extraction/#bedgraph-specific-options))
* **\*.M-bias.txt** (text file with methylation information in the context of the position in reads, helpful for 
  investigating bias as a function of base position in the read; see 
  [bismark documentation](http://felixkrueger.github.io/Bismark/bismark/methylation_extraction/#m-bias-plot))
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
* sample-1_bismark_bt2*.M-bias.txt (text file with methylation information in the context of the position in reads, 
  output from [Step 7](#7-extract-methylation-calls) above)
* sample-1_bismark_bt2*.deduplication_report.txt (optional deduplication report, output from [Step 6a.](#6a-deduplicate) if deduplication was run)
> **NOTE**  
> If processing RNA Methylation Sequencing data, file names will include "bismark_hisat2" instead of "bismark_bt2".

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
> If processing RNA Methylation Sequencing data, file names will include "bismark_hisat2" instead of "bismark_bt2".

**Output data:**  

* **bismark_summary_report_GLMethylSeq.txt** (summary table of general information on all included samples)
* **bismark_summary_report_GLMethylSeq.html** (html summary of general information on all included samples)

> **NOTE**
> If processing RNA Methylation Sequencing data, file suffix will be **GLRNAMethylSeq** instead of **GLMethylSeq**

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

> If processing RNA Methylation Sequencing data, file suffix will be **GLRNAMethylSeq** instead of **GLMethylSeq**

<br>


---

## 10. Generate reference genome annotation information

### 10a. GTF to BED conversion
A bed-formatted annotation file is needed for adding annotation information to the output from the differential 
methylation analysis. We utilize gtf files from [Ensembl](https://www.ensembl.org/) and convert them to the bed format 
as follows:

```bash
gtfToGenePred reference.gtf reference.genePred
genePredToBed reference.genePred reference.bed
```

**Parameter Definitions:**

* `reference.gtf` - genome annotation file in GTF format, provided as a positional argument
* `reference.genePred` - genome annotation file in genePred format, provided as a positional argument 
* `reference.bed` - genome annotation file in BED format, provided as a positional argument

**Input data:**

* reference.gtf (genome annotation, this pipeline version uses the Ensembl gtf file indicated in the `gtf` column of the 
[GL-DPPD-7110-A_annotations.csv](../../GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv) 
GeneLab Annotations file)

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
* `> reference-gene-to-transcript-map.tsv` - re-directs output to the provided file name


**Input data:**

* reference.gtf (genome annotation, this pipeline version uses the Ensembl gtf file indicated in the `gtf` column of the 
[GL-DPPD-7110-A_annotations.csv](../../GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv) 
GeneLab Annotations file)

**Output data:**

* reference-gene-to-transcript-map.tsv (a gene-to-transcript mapping file with gene IDs in the first column and 
  transcript IDs in the second)

<br>

---

## 11. Differential methylation analysis

The remainder of this document is performed in R. 

<br>


### 11a. Set up R environment

```R
### Install and load required packages ###

# List of required packages
cran_packages <- c("tidyverse")
bioc_packages <- c("methylKit", "genomation")

# Install missing packages
for (pkg in cran_packages) {
    if (!require(pkg, character.only = TRUE)) {
        install.packages(pkg)
    }
}
for (pkg in bioc_packages) {
    if (!require(pkg, character.only = TRUE)) {
        BiocManager::install(pkg)
    }
}

# Load all packages
library(tidyverse)
library(methylKit)
library(genomation)


### Set up working environment, variables, and object ###

# Define which organism is used in the study - this should be consistent with the species name in 
# the "species" column of the # GL-DPPD-7110-A_annotations.csv file
organism <- "organism_that_samples_were_derived_from"

### Define the location of the input data and where the output data will be printed to ###

runsheet_path <- file.path("/path/to/directory/containing/runsheet.csv/file")
work_dir <- file.path("/path/to/working/directory/where/script/is/executed/from")

# define paths to annotation files produced in Step 10a (reference.bed) and Step 10b (reference-gene-to-transcript-map.tsv)
ref_bed_path <- file.path("/path/to/reference-files/my_organism.bed")
ref_gene_transcript_map_path <- file.path("path/to/reference-files/my_organism-gene-to-transcript-map.tsv")

# path to bismark.cov.gz files from step 7 above
coverage_files_dir_path <- file.path("/path/to/bismark-coverage-files/")  

# Reference annotation table link (GL-DPPD-7110-A_annotations.csv) file
ref_table_link <- "https://raw.githubusercontent.com/nasa/GeneLab_Data_Processing/master/GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv"


### Pull in the GeneLab annotation table ###

# read in annotations.csv table
ref_table <- read.csv(ref_table_link)

# Set variable with link to appropriate annotations table
annotations_tab_link <- ref_table %>% 
    filter(species == organism) %>% pull(genelab_annots_link)


### Pull in annotation data from Step 10 ###

# Read in gene-to-transcript mapping file
gene_transcript_map <- read.table(ref_gene_transcript_map_path, sep = "\t", 
                                  col.names = c("gene_ID", "feature.name"))

# Read in reference bed file
gene.obj <- readTranscriptFeatures(ref_bed_path, up.flank = 1000, 
                                   down.flank = 1000, remove.unusual = TRUE, 
                                   unique.prom = TRUE)


##### Set variable with organism and ensembl version number for methylKit #####
ensembl_version <- ref_table %>% 
    filter(species == organism) %>% pull(ensemblVersion)

org_and_ensembl_version <- paste(organism, ensembl_version, sep = "_")

setwd(work_dir)

```

**Input data:**

* runsheet.csv (CSV formatted runsheet file containing unique sample IDs and factors as defined in the 
  [file specification](../Workflow_Documentation/examples/runsheet/README.md))
* `organism` (name of organism samples were derived from, found in the `species` column of 
  [GL-DPPD-7110-A_annotations.csv](../../GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv) file)
* /path/to/bismark-coverage-files/ (directory containing `*.bismark.cov.gz` (gzip-compressed bedGraph formatted files 
  with two additional coverage columns generated in [Step 7](#7-extract-methylation-calls)))
* reference.bed (genome annotation file in BED format generated in [Step 10a](#10a-gtf-to-bed-conversion))
* reference-gene-to-transcript-map.tsv (gene-to-transcript mapping file generated in 
  [Step 10b](#10b-making-a-mapping-file-of-genes-to-transcripts))

**Output data:**

* `runsheet_path` (variable containing path to runsheet file as defined in the 
  [file specification](../Workflow_Documentation/examples/runsheet/README.md)) 
* `gene.obj` (a GRangesList object containing locations of exon/intron/promoter/TSS)
* `gene_transcript_map` (DataFrame holding the gene-to-transcript mappings)
* `annotations_tab_link` (variable containing URL to GeneLab gene annotation table for the organism of interest)
* `org_and_ensembl_version` (variable containing the organism and ensembl version to use in methylKit input object generation)
* `coverage_files_dir_path` (variable containing path to bismark coverage files (`*.bismark.cov.gz`) generated in [Step7](#7-extract-methylation-calls))

<br>

### 11b. Configure Metadata, Sample Grouping, and Group Comparisons

```R
### Pull all factors for each sample in the study from the runsheet provided ###

compare_csv_from_runsheet <- function(runsheet_path) {
    df <- read.csv(runsheet_path)
    factors <- df %>%
        dplyr::select(starts_with("Factor.Value", ignore.case = TRUE)) %>%
        rename_with(~ paste0("factor_", seq_along(.)))
    result <- df %>%
        dplyr::select(sample_id = "Sample.Name") %>%
        bind_cols(factors)
    return(result)
}


### Load metadata from runsheet csv file and create data frame containing all samples and respective factors ###
study <- compare_csv_from_runsheet(runsheet_path) %>% column_to_rownames(var = "sample_id")

##### Format groups and indicate the group that each sample belongs to #####
group <- apply(study, 1, paste, collapse = " & ")
group_names <- paste0("(", group, ")") ## human readable group names
# generate group naming compatible with R models, maintains the default behaviour of make.names but ensures 'X' is never prepended
group <- sub("^BLOCKER_", "", make.names(paste0("BLOCKER_", group)))  
names(group) <- group_names
rm(group_names)

# create a lookup table to map back from R naming to human readable group names
group_name_lookup <- data.frame(group = group, group_names = names(group)) %>%
    distinct() %>%
    column_to_rownames(var = "group")


##### Format contrasts table, defining pairwise comparisons for all groups, one-way contrasts #####

# generate matrix of pairwise group combinations for comparison
contrast.names <- combn(levels(factor(names(group))), 2)

# create computationally friendly contrast values
contrasts <- apply(
    contrast.names,
    MARGIN = 2,
    function(col) {
        # limited make.names call for each group (also removes leading parentheses)
        sub("^BLOCKER_", "", make.names(paste0("BLOCKER_", stringr::str_sub(col, 2, -2)))) 
    }
)

# format contrast combinations for output tables
contrast.names <- paste(contrast.names[1, ], contrast.names[2, ], sep = "v") 
colnames(contrasts) <- contrast.names
rm(contrast.names)

# Create sampleTable
sampleTable <- data.frame("sample_id" = rownames(study), "condition" = group)

```

**Input Data:**

* `runsheet_path` (variable containing path to runsheet file as defined in the 
  [file specification](../Workflow_Documentation/examples/runsheet/README.md))

**Output Data:**

* `study` (data frame specifying factor levels assigned to each sample)
* `group` (named vector specifying the group or set of factor levels for each sample)
* `contrasts` (matrix defining pairwise comparisons between groups)
* `sampleTable` (data frame mapping sample IDs to groups (conditions))
* `group_name_lookup` (data frame mapping groups with R naming scheme to human-readable group names)

<br>

### 11c. Import Methylation Calls

```R
### Import methylation call data ###

order_input_files <- function(sample_names, paths) {
    # this function takes in a vector of sample_names and
    # a vector of paths to coverage files, and returns a vector
    # of coverage files in the same order as the sample names

    ordered_paths <- c()

    for (sample in sample_names) {
        search_pattern <- paste0(sample, "_bismark")
        hits <- paths[grep(search_pattern, paths)]

        # making sure there is exactly one match
        if (length(hits) != 1) {
            stop("\n  There was a problem matching up sample names with their coverage files.\nCannot proceed.\n\n", call. = FALSE)
        }

        ordered_paths <- c(ordered_paths, hits)
    }

    return(ordered_paths)
}

# Get list of *.bismark.cov.gz files
bismark_cov_paths <- list.files(
  path = coverage_files_dir_path,
  pattern = ".*.bismark.cov.gz",
  full.names = TRUE,
  recursive = TRUE
)

# Make sure files list matches number of samples in runsheet
if (dim(sampleTable)[1] != length(bismark_cov_paths)) {
    error_message <- paste0(
        "\n  The number of '*.bismark.cov.gz' files found in the ",
        coverage_files_dir_path,
        " directory\n  does not match the number of samples specified in the runsheet.\nCannot proceed.\n\n"
    )

    stop(error_message, call. = FALSE)
}

# Make sure coverage-file-paths vector is in the same order as the sample names
bismark_cov_paths <- order_input_files(sampleTable$sample_id, bismark_cov_paths)

# Make a single table with info needed for methylkit
sample_meth_info_df <- full_join(
    data.frame("sample_id" = sampleTable$sample_id, "coverage_file_path" = bismark_cov_paths),
    sampleTable,
    "sample_id")

# Read in methylation counts
# NOTE: treatment vector is composed of integers corresponding to the groups. 
#       It is used during the min.per.group filtering in the unite command
meth_obj <- methRead(
    location = as.list(sample_meth_info_df %>% pull(coverage_file_path)),
    sample.id = as.list(sample_meth_info_df %>% pull(sample_id)),
    treatment = as.list(sample_meth_info_df %>% mutate(group_no = as.integer(factor(condition))) %>% pull(group_no)),
    pipeline = "bismarkCoverage",
    assembly = org_and_ensembl_version,
    header = FALSE,
    mincov = 10, # set minimum read coverage needed to include a base in the methylKit object.
    dbtype = "NA"  # indicates where the object should be stored ("NA" == in memory, "tabix" == on the file system)
)

```

**Input Data:**

* `coverage_files_dir_path` (variable containing path to bismark coverage files (`*.bismark.cov.gz`) generated in 
  [Step7](#7-extract-methylation-calls))
* `org_and_ensembl_version` (variable containing the organism name and ensembl version from [Step 11a](#11a-set-up-r-environment))
* `sampleTable` (data frame containing sample condition values, output from 
  [Step 11b](#11b-configure-metadata-sample-grouping-and-group-comparisons))

**Output Data:**

* `meth_obj` (methylRawList object holding methylation and coverage information for all samples)
* `sample_meth_info_df` (data frame containing sample IDs mapped to conditions and coverage file paths)

<br>

### 11d. Individual-base analysis

```R
### Function for computing differential methylation per base for each contrasts ###
compute_contrast_bases <- function(i) {
    # Get current contrast
    curr_contrasts_vec <- contrasts[, i]

    # Get subset sample info table relevant to current contrast
    curr_sample_info_df <- sample_meth_info_df %>% filter(condition %in% curr_contrasts_vec)

    # Get which samples are relevant to current contrast
    curr_samples_vec <- curr_sample_info_df %>% pull(sample_id)

    # Make binary vector for treatment argument to methRead()
    curr_treatment_vec <- c()
    for (value in curr_sample_info_df$condition) {
        if (value == curr_sample_info_df$condition[1]) {
            curr_treatment_vec <- c(curr_treatment_vec, 1)
        } else {
            curr_treatment_vec <- c(curr_treatment_vec, 0)
        }
    }

    # return methylation statistics for bases
    return(
        as.data.frame(getData(
            calculateDiffMeth(
                methylKit::unite(reorganize(
                    meth_obj,
                    sample.ids = as.list(curr_samples_vec),
                    treatment = as.integer(factor(curr_sample_info_df$condition)) - 1
                ),
                min.per.group = 1L), # Keep only bases with coverage in at least one sample per group
                mc.cores = myargs$mc_cores)
            )
        ) %>%
        relocate("meth.diff", .before = "pvalue") %>% # move methyl diff value before stats values
        rename_with(~ paste0(.x, "_", colnames(contrasts)[i]), any_of(c("pvalue", "qvalue", "meth.diff")))
    )
}


### Normalize the data ###

# First, filter samples by coverage to account for PCR bias or over-amplification
meth_obj <- filterByCoverage(meth_obj, lo.count = 10, lo.perc = NULL, hi.count = NULL, hi.perc = 99.9)
# Normalize coverage between samples using a scaling factor derived from the median coverage distributions
norm_meth_obj <- normalizeCoverage(meth_obj, method = "median")


### Base analysis ###

# create initial output object at base resolution
#   - compute pair-wise differential methylation for each contrast
#   - add percent methylation values for all samples
#   - add All.mean and All.stdev stats
# NOTE: when retrieving percent methylation, keep all bases (min.per.group = 0L) to
#       include bases that have coverage in some groups (and data for only some
#       contrasts) but no coverage in others
output_bases_df <- lapply(seq_len(dim(contrasts)[2]), compute_contrast_bases) %>%
    purrr::reduce(full_join, by = c("chr", "start", "end", "strand")) %>%
    left_join(as.data.frame(percMethylation(methylKit::unite(norm_meth_obj, min.per.group = 0L), rowids = TRUE)) %>%
        rownames_to_column() %>%
        separate_wider_regex(rowname, c(chr = ".*", "\\.", start = ".*", "\\.", end = ".*")) %>%
        mutate_at(vars(start, end), as.integer) %>%
        rowwise() %>%
        mutate(
            All.mean = mean(c_across(all_of(sample_meth_info_df$sample_id)), na.rm = TRUE),
            All.stdev = sd(c_across(all_of(sample_meth_info_df$sample_id)), na.rm = TRUE)
        ), by = c("chr", "start", "end"))

# Add group mean and stdev
for (g in rownames(group_name_lookup)) {
    colnames_to_process <- sample_meth_info_df %>% filter(condition == g) %>% pull(sample_id)
    output_bases_df <- output_bases_df %>%
        rowwise() %>%
        mutate(
            gmean = mean(c_across(all_of(colnames_to_process)), na.rm = TRUE),
            gstdev = sd(c_across(all_of(colnames_to_process)), na.rm = TRUE)
        ) %>%
        rename("gmean" = paste0("Group.Mean_", group_name_lookup[g, ]),
               "gstdev" = paste0("Group.Stdev_", group_name_lookup[g, ]))
}

# Annotate with gene parts (promoter/exon/intron)
suppressWarnings(output_ann_df <- annotateWithGeneParts(as(output_bases_df, "GRanges"), gene.obj))

# Associate the nearest TSS annotations with gene parts
# Add annotations to output bases dataframe
df_list <- list(
        as.tibble(getAssociationWithTSS(output_ann_df) %>% dplyr::rename(rowid = target.row)),
        tibble::rowid_to_column(as.data.frame(genomation::getMembers(output_ann_df))),
        tibble::rowid_to_column(output_bases_df))

# Create one combined table with all annotations, methylation values, and statistics
# reorder the columns to final expected format
bases_tab_with_features_and_annots <- df_list %>%
        purrr::reduce(full_join, by = "rowid") %>%
        mutate(rowid = NULL) %>%
    left_join(gene_transcript_map) %>%
    left_join(functional_annots_tab, by = c("gene_ID" = myargs$primary_keytype)) %>%
    rename("gene_ID" = myargs$primary_keytype) %>%
    relocate(c("ENSEMBL", "SYMBOL", "GENENAME", "REFSEQ", "ENTREZID", "STRING_id", "GOSLIM_IDS", "feature.name", "chr", "start", "end", "strand")) %>%
    relocate(sample_meth_info_df$sample_id, .after = "intron")
rm(df_list)

```

**Input data:**

* `contrasts` (matrix defining pairwise comparisons between groups from [Step 11b](#11b-configure-metadata-sample-grouping-and-group-comparisons))
* `group_name_lookup` (data frame mapping groups with R naming scheme to human-readable group names from [Step 11b](#11b-configure-metadata-sample-grouping-and-group-comparisons))
* `meth_obj` (methylRawList object created in [Step 11c](#11c-import-methylation-calls)
* `sample_meth_info_df` (data frame containing sample IDs mapped to conditions and coverage file paths from [Step 11c](#11c-import-methylation-calls))

**Output data:**

* `bases_tab_with_features_and_annots` (data frame containing the following columns:
  - Gene identifier column (ENSEMBL for non-plant or TAIR for plant studies) 
  - Additional organism-specific gene annotation columns
  - Feature.name (ID of feature that the annotation derives from)
  - Chromosomal coordinate (of the methylated base)
    - Chr
    - Start
    - End
    - Strand (`*` if unstranded)
  - dist.to.feature (distance between annotated feature and methylated base)
  - feature.strand (strand of the feature)
  - Gene part annotation (0 if false, 1 if true)
    - prom (promoter)
    - exon
    - intron
  - Normalized percent methylation in this chromosomal coordinate for each sample
  - For each pair-wise comparison:
    - meth.diff (differential methylation value)
    - pvalue
    - qvalue
  - Statistics for percent methylation values:
    - All.mean (mean across all samples)
    - All.stdev (standard deviation across all samples) 
    - For each experimental group:
      - Group.Mean_(group) (mean within group)
      - Group.Stdev_(group) (standard deviation within group))
* `norm_meth_obj` (median-normalized methylRawList object)

<br>

### 11e. Tile analysis

```R
### Tile analysis ###

# Group bases into tiles
meth_tiles_obj <- tileMethylCounts(norm_meth_obj,
                                   win.size = 1000,
                                   step.size = 1000,
                                   cov.bases = 10)

# create initial output object at tile resolution
#   - compute pair-wise differential methylation for each contrast
#   - add percent methylation values for all samples
#   - add All.mean and All.stdev stats
# NOTE: when retrieving percent methylation, keep all tiles (min.per.group = 0L) to
#       include tiles that have coverage in some groups (and data for only some
#       contrasts) but no coverage in others
output_tiles_df <- lapply(seq_len(dim(contrasts)[2]), compute_contrast_tiles) %>%
    purrr::reduce(full_join, by = c("chr", "start", "end", "strand")) %>%
    left_join(as.data.frame(percMethylation(methylKit::unite(meth_tiles_obj, min.per.group = 0L), rowids = TRUE)) %>%
        rownames_to_column() %>%
        separate_wider_regex(rowname, c(chr = ".*", "\\.", start = ".*", "\\.", end = ".*")) %>%
        mutate_at(vars(start, end), as.integer) %>%
        rowwise() %>%
        mutate(
            All.mean = mean(c_across(all_of(sample_meth_info_df$sample_id)), na.rm = TRUE),
            All.stdev = sd(c_across(all_of(sample_meth_info_df$sample_id)), na.rm = TRUE)
        ), by = c("chr", "start", "end"))

# Add group means and stdev
for (g in rownames(group_name_lookup)) {
    colnames_to_process <- sample_meth_info_df %>% filter(condition == g) %>% pull(sample_id)
    output_tiles_df <- output_tiles_df %>%
        rowwise() %>%
        mutate(
            gmean = mean(c_across(all_of(colnames_to_process)), na.rm = TRUE),
            gstdev = sd(c_across(all_of(colnames_to_process)), na.rm = TRUE)
        ) %>%
        rename("gmean" = paste0("Group.Mean_", group_name_lookup[g, ]),
               "gstdev" = paste0("Group.Stdev_", group_name_lookup[g, ]))
}

# Annotate with gene parts (promoter/exon/intron)
suppressWarnings(output_tiles_ann_df <- annotateWithGeneParts(as(output_tiles_df, "GRanges"), gene.obj))

# Associate the nearest TSS annotations with gene parts
# Add annotations to output bases dataframe
df_list <- list(
        as.tibble(getAssociationWithTSS(output_tiles_ann_df)) %>% dplyr::rename(rowid = target.row),
        tibble::rowid_to_column(as.data.frame(genomation::getMembers(output_tiles_ann_df))),
        tibble::rowid_to_column(output_tiles_df))

# Create one combined table with all annotations, methylation values, and statistics
# reorder the columns to final expected format
tiles_tab_with_features_and_annots <- df_list %>%
        purrr::reduce(full_join, by = "rowid") %>%
        mutate(rowid = NULL) %>%
    left_join(gene_transcript_map) %>%
    left_join(functional_annots_tab, by = c("gene_ID" = myargs$primary_keytype)) %>%
    rename("gene_ID" = myargs$primary_keytype) %>%
    relocate(c("ENSEMBL", "SYMBOL", "GENENAME", "REFSEQ", "ENTREZID", "STRING_id", "GOSLIM_IDS", "feature.name", "chr", "start", "end", "strand")) %>%
    relocate(sample_meth_info_df$sample_id, .after = "intron")
rm(df_list)
```

**Input data:**

* `contrasts` (matrix defining pairwise comparisons between groups from [Step 11b](#11b-configure-metadata-sample-grouping-and-group-comparisons))
* `group_name_lookup` (data frame mapping groups with R naming scheme to human-readable group names from [Step 11b](#11b-configure-metadata-sample-grouping-and-group-comparisons))
* `sample_meth_info_df` (data frame containing sample IDs mapped to conditions and coverage file paths from [Step 11c](#11c-import-methylation-calls))
* `norm_meth_obj` (median-normalized methylRawList object created in [Step 11d](#11d-individual-base-analysis)

**Output data:**

* `tiles_tab_with_features_and_annots` (data frame containing the following columns:
  - Gene identifier column (ENSEMBL for non-plant or TAIR for plant studies) 
  - Additional organism-specific gene annotation columns
  - Feature.name (ID of feature that the annotation derives from)
  - Chromosomal coordinate (of the methylated region)
    - Chr
    - Start
    - End
    - Strand (`*` if unstranded)
  - dist.to.feature (distance between annotated feature and methylated base)
  - feature.strand (strand of the feature)
  - Gene part annotation (0 if false, 1 if true)
    - prom (promoter)
    - exon
    - intron
  - Normalized percent methylation in this chromosomal coordinate for each sample
  - For each pair-wise comparison:
    - meth.diff (differential methylation value)
    - pvalue
    - qvalue
  - Statistics for percent methylation values:
    - All.mean (mean across all samples)
    - All.stdev (standard deviation across all samples) 
    - For each experimental group:
      - Group.Mean_(group) (mean within group)
      - Group.Stdev_(group) (standard deviation within group))

<br>

### 11f. Export Tables

```R
### Output contrasts table ###
write.csv(
    contrasts,
    file = file.path(myargs$methylkit_output_dir,
                     paste0("contrasts", myargs$file_suffix, ".csv"))
)

### Output Sample table ###
write.csv(
    sampleTable,
    file = file.path(myargs$methylkit_output_dir,
                     paste0("SampleTable", myargs$file_suffix, ".csv")),
    row.names = FALSE
)

### Export differentially methylated bases ###
data.table::fwrite(
    bases_tab_with_features_and_annots,
    row.names = FALSE,
    file = file.path(myargs$methylkit_output_dir,
        paste0("differential_methylation_bases", myargs$file_suffix, ".csv")
    ),
    quote = TRUE, na = "NA"
)

### Export differentially methylated tiles ###
data.table::fwrite(
    tiles_tab_with_features_and_annots,
    row.names = FALSE,
    file = file.path(myargs$methylkit_output_dir,
        paste0("differential_methylation_tiles", myargs$file_suffix, ".csv")
    ),
    quote = TRUE, na = "NA"
)

### print session info ###
print("Session Info: ")
sessionInfo()
print(paste0("BioC_version_associated_with_R_version: ",BiocManager::version()))

```

**Input data:**

* `contrasts` (matrix defining pairwise comparisons between groups from [Step 11b](#11b-configure-metadata-sample-grouping-and-group-comparisons))
* `sampleTable` (data frame mapping sample IDs to groups (conditions) from [Step 11b](#11b-configure-metadata-sample-grouping-and-group-comparisons))
* `bases_tab_with_features_and_annots` (Methylated base output table from [Step 11d](#11d-individual-base-analysis))
* `tiles_tab_with_features_and_annots` (Methylated tile output tables from [Step 11e](#11e-tile-analysis))

**Output data:**
* **SampleTable_GLMethylSeq.csv** (table specifying the group or set of factor levels for each sample)
* **contrasts_GLMethylSeq.csv** (table listing all pairwise group comparisons)
* **differential_methylation_bases_GLMethylSeq.csv** (Methylated bases table containing the following columns:
  - Gene identifier column (ENSEMBL for non-plant or TAIR for plant studies) 
  - Additional organism-specific gene annotation columns
  - Feature.name (ID of feature that the annotation derives from)
  - Chromosomal coordinate (of the methylated base)
    - Chr
    - Start
    - End
    - Strand (`*` if unstranded)
  - dist.to.feature (distance between annotated feature and methylated base)
  - feature.strand (strand of the feature)
  - Gene part annotation (0 if false, 1 if true)
    - prom (promoter)
    - exon
    - intron
  - Normalized percent methylation in this chromosomal coordinate for each sample
  - For each pair-wise comparison:
    - meth.diff (differential methylation value)
    - pvalue
    - qvalue
  - Statistics for percent methylation values:
    - All.mean (mean across all samples)
    - All.stdev (standard deviation across all samples) 
    - For each experimental group:
      - Group.Mean_(group) (mean within group)
      - Group.Stdev_(group) (standard deviation within group))
* **differential_methylation_tiles_GLMethylSeq.csv** (Methylated tiles table containing the following columns:
  - Gene identifier column (ENSEMBL for non-plant or TAIR for plant studies) 
  - Additional organism-specific gene annotation columns
  - Feature.name (ID of feature that the annotation derives from)
  - Chromosomal coordinate (of the methylated region)
    - Chr
    - Start
    - End
    - Strand (`*` if unstranded)
  - dist.to.feature (distance between annotated feature and methylated base)
  - feature.strand (strand of the feature)
  - Gene part annotation (0 if false, 1 if true)
    - prom (promoter)
    - exon
    - intron
  - Normalized percent methylation in this chromosomal coordinate for each sample
  - For each pair-wise comparison:
    - meth.diff (differential methylation value)
    - pvalue
    - qvalue
  - Statistics for percent methylation values:
    - All.mean (mean across all samples)
    - All.stdev (standard deviation across all samples) 
    - For each experimental group:
      - Group.Mean_(group) (mean within group)
      - Group.Stdev_(group) (standard deviation within group))

<br>

---
