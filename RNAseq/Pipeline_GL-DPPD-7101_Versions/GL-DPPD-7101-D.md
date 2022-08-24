# GeneLab bioinformatics processing pipeline for Illumina RNA-sequencing data

> **This page holds an overview and instructions for how GeneLab processes RNAseq datasets. Exact processing commands and GL-DPPD-7101 revision used for specific datasets are available in the [GLDS_Processing_Scripts](../GLDS_Processing_Scripts) sub-directory and are also provided with their processed data in the [GeneLab Data Systems (GLDS) repository](https://genelab-data.ndc.nasa.gov/genelab/projects).**  

---

**Date:** September 13, 2021  
**Revision:** D  
**Document Number:** GL-DPPD-7101-D  

**Submitted by:**  
Amanda Saravia-Butler (GeneLab Data Processing Team)

**Approved by:**  
Sylvain Costes (GeneLab Project Manager)  
Samrawit Gebre (GeneLab Deputy Project Manager and Interim GeneLab Configuration Manager)   
Jonathan Galazka (GeneLab Project Scientist)

---

## Updates from previous revision

Two additional sub-steps were added to step 4:
- A step to compile the alignment log files using multiQC, [step 4b](#4b-compile-alignment-logs)
- A step to index the alignment files, [step 4c](#4c-index-aligned-reads), which is required to assess read strandedness

Two additional steps were added prior to aligned read quantitation:
- Step 5, [5a](#5a-convert-gtf-to-genepred-file) and [5b](#5b-convert-genepred-to-bed-file), was added to create a reference annotation BED file required to assess read strandedness
- Step 6 was added to [determine read strandedness](#6a-determine-read-strandedness) and [compile read strandedness reports](#6b-compile-strandedness-reports), to determine which RSEM `--strandedness` setting to use during aligned read quantitation

The aligned read quantitation step, now [step 8](#8-count-aligned-reads-with-rsem), was modified to use the results of the read strandedness step to inform the correct RSEM `--strandedness` setting
> Note: A subset of samples from all datasets previously processed were evaluated for strandedness, and those datasets identified to have been processed with the incorrect RSEM `--strandedness` setting were reprocessed with the correct setting

The DESeq2 Normalization and DGE step for datasets with ERCC spike-in, now [step 9a](#9a-for-datasets-with-ercc-spike-in), was modified as follows:
- Perform ERCC normalization using only ERCC group B genes, since the concentration of these genes are the same in ERCC mix 1 and mix 2
- Remove any samples that do not contain detectible ERCC group B spike-ins prior to generation and subsequent analysis of ERCC-normalized count data
- Account for the edge case in which rescaling using ERCC size factors fails due to zero gene counts

Added "Stat_" column containing the Wald Statistic (similar to a Z-score) to the DGE output tables for datasets both with and without ERCC spike-in, now [step 9a](#9a-for-datasets-with-ercc-spike-in) and [step 9b](#9b-for-datasets-without-ercc-spike-in), respectively, which will be used for GSEA visualizations

---

# Table of contents  

- [**Software used**](#software-used)
- [**General processing overview with example commands**](#general-processing-overview-with-example-commands)
  - **1. Raw Data QC**
    - [**1a. Raw Data QC**](#1a-raw-data-qc)
    - [**1b. Compile Raw Data QC**](#1b-compile-raw-data-qc)
  - **2. Trim/Filter Raw Data and Trimmed Data QC**
    - [**2a. Trim/Filter Raw Data**](#2a-trimfilter-raw-data)
    - [**2b. Trimmed Data QC**](#2b-trimmed-data-qc)
    - [**2c. Compile Trimmed Data QC**](#2c-compile-trimmed-data-qc)
  - [**3. Build STAR Reference**](#3-build-star-reference)
  - **4. Align Reads to Reference Genome then Index**
    - [**4a. Align Reads to Reference Genome with STAR**](#4a-align-reads-to-reference-genome-with-star)
    - [**4b. Compile Alignment Logs**](#4b-compile-alignment-logs)
    - [**4c. Index Aligned Reads**](#4c-index-aligned-reads)
  - **5. Create Reference BED File**
    - [**5a. Convert GTF to genePred File**](#5a-convert-gtf-to-genepred-file)
    - [**5b. Convert genePred to BED File**](#5b-convert-genepred-to-bed-file)
  - **6. Assess Strandedness**
    - [**6a. Determine Read Strandedness**](#6a-determine-read-strandedness)
    - [**6b. Compile Strandedness Reports**](#6b-compile-strandedness-reports)
  - [**7. Build RSEM Reference**](#7-build-rsem-reference)
  - [**8. Count Aligned Reads with RSEM**](#8-count-aligned-reads-with-rsem)
  - [**9. Normalize Read Counts, Perform Differential Gene Expression Analysis, and Add Gene Annotations in R**](#9-normalize-read-counts-perform-differential-gene-expression-analysis-and-add-gene-annotations-in-r)
    - [**9a. For Datasets with ERCC Spike-In**](#9a-for-datasets-with-ercc-spike-in)
    - [**9b. For Datasets without ERCC Spike-In**](#9b-for-datasets-without-ercc-spike-in)
  
---

# Software used  

|Program|Version*|Relevant Links|
|:------|:------:|:-------------|
|FastQC|`fastqc -v`|[https://www.bioinformatics.babraham.ac.uk/projects/fastqc/](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)|
|MultiQC|`multiqc -v`|[https://multiqc.info/](https://multiqc.info/)|
|Cutadapt|`cutadapt --version`|[https://cutadapt.readthedocs.io/en/stable/](https://cutadapt.readthedocs.io/en/stable/)|
|TrimGalore!|`trim_galore -v`|[https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)|
|STAR|`STAR --version`|[https://github.com/alexdobin/STAR](https://github.com/alexdobin/STAR)|
|RSEM|`rsem-calculate-expression --version`|[https://github.com/deweylab/RSEM](https://github.com/deweylab/RSEM)|
|Samtools|`samtools --version`|[http://www.htslib.org/](http://www.htslib.org/)|
|gtfToGenePred|N/A|[http://hgdownload.cse.ucsc.edu/admin/exe/](http://hgdownload.cse.ucsc.edu/admin/exe/)|
|genePredToBed|N/A|[http://hgdownload.cse.ucsc.edu/admin/exe/](http://hgdownload.cse.ucsc.edu/admin/exe/)|
|infer_experiment|`infer_experiment.py --version`|[http://rseqc.sourceforge.net/#infer-experiment-py](http://rseqc.sourceforge.net/#infer-experiment-py)|
|Bioconductor|`BiocManager::version()`|[https://bioconductor.org](https://bioconductor.org)|
|DESeq2|`packageVersion("DESeq2")`|[https://bioconductor.org/packages/release/bioc/html/DESeq2.html](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)|
|tximport|`packageVersion("tximport")`|[https://bioconductor.org/packages/release/bioc/html/tximport.html](https://bioconductor.org/packages/release/bioc/html/tximport.html)|
|tidyverse|`packageVersion("tidyverse")`|[https://www.tidyverse.org](https://www.tidyverse.org)|
|Risa|`packageVersion("Risa")`|[https://www.bioconductor.org/packages/release/bioc/html/Risa.html](https://www.bioconductor.org/packages/release/bioc/html/Risa.html)|
|STRINGdb|`packageVersion("STRINGdb")`|[https://www.bioconductor.org/packages/release/bioc/html/STRINGdb.html](https://www.bioconductor.org/packages/release/bioc/html/STRINGdb.html)|
|PANTHER.db|`packageVersion("PANTHER.db")`|[https://bioconductor.org/packages/release/data/annotation/html/PANTHER.db.html](https://bioconductor.org/packages/release/data/annotation/html/PANTHER.db.html)|
|org.Hs.eg.db|`packageVersion("org.Hs.eg.db")`|[https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html](https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html)|
|org.Mm.eg.db|`packageVersion("org.Mm.eg.db")`|[https://bioconductor.org/packages/release/data/annotation/html/org.Mm.eg.db.html](https://bioconductor.org/packages/release/data/annotation/html/org.Mm.eg.db.html)|
|org.Dm.eg.db|`packageVersion("org.Dm.eg.db")`|[https://bioconductor.org/packages/release/data/annotation/html/org.Dm.eg.db.html](https://bioconductor.org/packages/release/data/annotation/html/org.Dm.eg.db.html)|
|org.Ce.eg.db|`packageVersion("org.Ce.eg.db")`|[https://bioconductor.org/packages/release/data/annotation/html/org.Ce.eg.db.html](https://bioconductor.org/packages/release/data/annotation/html/org.Ce.eg.db.html)|
|org.At.tair.db|`packageVersion("org.At.tair.db")`|[https://bioconductor.org/packages/release/data/annotation/html/org.At.tair.db.html](https://bioconductor.org/packages/release/data/annotation/html/org.At.tair.db.html)|
|org.EcK12.eg.db|`packageVersion("org.EcK12.eg.db")`|[https://bioconductor.org/packages/release/data/annotation/html/org.EcK12.eg.db.html](https://bioconductor.org/packages/release/data/annotation/html/org.EcK12.eg.db.html)|
|org.Sc.sgd.db|`packageVersion("org.Sc.sgd.db")`|[https://bioconductor.org/packages/release/data/annotation/html/org.Sc.sgd.db.html](https://bioconductor.org/packages/release/data/annotation/html/org.Sc.sgd.db.html)|

>**\*** Exact versions are available along with the processing commands for each specific dataset in the [GLDS_Processing_Scripts](../GLDS_Processing_Scripts) sub-directory. 

---

# General processing overview with example commands  

> Exact processing commands for specific datasets are provided in the [GLDS_Processing_Scripts](../GLDS_Processing_Scripts) sub-directory and are also provided with their processed data in the [GeneLab Data Systems (GLDS) repository](https://genelab-data.ndc.nasa.gov/genelab/projects).  

---

## 1a. Raw Data QC  

```
fastqc -o /path/to/raw_fastqc/output/directory *.fastq.gz
```

**Parameter Definitions:**

* `-o` – the output directory to store results
* `*.fastq.gz` – the input reads are specified as a positional argument, and can be given all at once with wildcards like this, or as individual arguments with spaces inbetween them

**Input Data:**
- *fastq.gz (raw reads)

**Output Data:**
- *fastqc.html (FastQC report)
- *fastqc.zip (FastQC data)

<br>

## 1b. Compile Raw Data QC  

```
multiqc -n raw_multiqc -o /path/to/raw_multiqc/output/directory /path/to/directory/containing/raw_fastqc/files
```

**Parameter Definitions:**

* `-n` - prefix name for output files
* `-o` – the output directory to store results
* `/path/to/directory/containing/raw_fastqc/files` – the directory holding the output data from the fastqc run, provided as a positional argument

**Input Data:**
- *fastqc.html (FastQC report)
- *fastqc.zip (FastQC data)

**Output Data:**
- raw_multiqc.html (multiqc report)
- raw_multiqc_data (directory containing multiqc data)

<br>

---

## 2a. Trim/Filter Raw Data  

```
trim_galore --gzip \
  --path_to_cutadapt /path/to/cutadapt \
  --phred33 \
  --illumina \ # if adapters are not illumina, replace with adapters used
  --output_dir /path/to/TrimGalore/output/directory \
  --paired \ # only for PE studies, remove this paramater if raw data are SE
  sample1_R1_raw.fastq.gz sample1_R2_raw.fastq.gz sample2_R1_raw.fastq.gz sample2_R2_raw.fastq.gz
# if SE, replace the last line with only the forward reads (R1) of each sample

```

**Parameter Definitions:**

* `--gzip` – compress the output files with `gzip`
* `--path_to_cutadapt` - specify path to cutadapt software if it is not in your `$PATH`
* `--phred33` - instructs cutadapt to use ASCII+33 quality scores as Phred scores for quality trimming
* `--illumina` - defines the adapter sequence to be trimmed as the first 13bp of the Illumina universal adapter `AGATCGGAAGAGC`
* `--output_dir` - the output directory to store results
* `--paired` - indicates paired-end reads - both reads, forward (R1) and reverse (R2) must pass length threshold or else both reads are removed
* `sample1_R1_raw.fastq.gz sample1_R2_raw.fastq.gz sample2_R1_raw.fastq.gz sample2_R2_raw.fastq.gz` – the input reads are specified as a positional argument, paired-end read files are listed pairwise such that the forward reads (*R1_raw.fastq.gz) are immediately followed by the respective reverse reads (*R2_raw.fastq.gz) for each sample

**Input Data:**
- *fastq.gz (raw reads)

**Output Data:**
- *fastq.gz (trimmed reads)
- *trimming_report.txt (trimming report)

<br>

## 2b. Trimmed Data QC  

```
fastqc -o /path/to/trimmed_fastqc/output/directory *.fastq.gz
```

**Parameter Definitions:**

* `-o` – the output directory to store results
* `*.fastq.gz` – the input reads are specified as a positional argument, and can be given all at once with wildcards like this, or as individual arguments with spaces inbetween them

**Input Data:**
- *fastq.gz (trimmed reads)

**Output Data:**
- *fastqc.html (FastQC report)
- *fastqc.zip (FastQC data)

<br>

## 2c. Compile Trimmed Data QC  

```
multiqc -n trimmed_multiqc -o /path/to/trimmed_multiqc/output/directory /path/to/directory/containing/trimmed_fastqc/files
```

**Parameter Definitions:**

* `-n` - prefix name for output files
* `-o` – the output directory to store results
* `/path/to/directory/containing/trimmed_fastqc/files` – the directory holding the output data from the fastqc run, provided as a positional argument

**Input Data:**
- *fastqc.html (FastQC report)
- *fastqc.zip (FastQC data)

**Output Data:**
- trimmed_multiqc.html (multiqc report)
- trimmed_multiqc_data (directory containing multiqc data)

<br>

---

## 3. Build STAR Reference  

```
STAR --runThreadN NumberOfThreads \
  --runMode genomeGenerate \
  --limitGenomeGenerateRAM 55000000000 \
  --genomeSAindexNbases 14 \
  --genomeDir /path/to/STAR/genome/directory \
  --genomeFastaFiles /path/to/genome/fasta/file \
  --sjdbGTFfile /path/to/annotation/gtf/file \
  --sjdbOverhang ReadLength-1

```

**Parameter Definitions:**

* `--runThreadN` – number of threads available on server node to create STAR reference
* `--runMode` - instructs STAR to run genome indices generation job
* `--limitGenomeGenerateRAM` - maximum RAM available (in bytes) to generate STAR reference, at least 35GB are needed for mouse and the example above shows 55GB
* `--genomeSAindexNbases` - length (in bases) of the SA pre-indexing string, usually between 10 and 15. Longer strings require more memory but allow for faster searches. This value should be scaled down for smaller genomes (like bacteria) to min(14, log2(GenomeLength)/2 - 1). For example, for a 1 megaBase genome this value would be 9.
* `--genomeDir` - specifies the path to the directory where the STAR reference will be stored. At least 100GB of available disk space is required for mammalian genomes.
* `--genomeFastaFiles` - specifies one or more fasta file(s) containing the genome reference sequences
* `--sjdbGTFfile` – specifies the file(s) containing annotated transcripts in the standard gtf format
* `--sjdbOverhang` - indicates the length of the genomic sequence around the annotated junction to be used in constructing the splice junctions database. The length should be one less than the maximum length of the reads.

**Input Data:**
- *.fasta (genome sequence\#)
- *.gtf (genome annotation\#)

\#See document(s) in the [GeneLab_Reference_and_Annotation_Files](../GeneLab_Reference_and_Annotation_Files) sub-directory for a list of the ensembl fasta genome sequences and associated gtf annotation files used to generate the RNAseq processed data available in the [GLDS repository](https://genelab-data.ndc.nasa.gov/genelab/projects).

**Output Data:**

STAR genome reference, which consists of the following files:
- chrLength.txt
- chrNameLength.txt
- chrName.txt
- chrStart.txt
- exonGeTrInfo.tab
- exonInfo.tab
- geneInfo.tab
- Genome
- genomeParameters.txt
- SA
- SAindex
- sjdbInfo.txt
- sjdbList.fromGTF.out.tab
- sjdbList.out.tab
- transcriptInfo.tab

<br>

---

## 4a. Align Reads to Reference Genome with STAR

```
STAR --twopassMode Basic \
	--limitBAMsortRAM 65000000000 \
	--genomeDir /path/to/STAR/genome/directory \
	--outSAMunmapped Within \
	--outFilterType BySJout \
	--outSAMattributes NH HI AS NM MD MC \
	--outFilterMultimapNmax 20 \
	--outFilterMismatchNmax 999 \
	--outFilterMismatchNoverReadLmax 0.04 \
	--alignIntronMin 20 \
	--alignIntronMax 1000000 \
	--alignMatesGapMax 1000000 \ # for PE only
	--alignSJoverhangMin 8 \
	--alignSJDBoverhangMin 1 \
	--sjdbScore 1 \
	--readFilesCommand zcat \
	--runThreadN NumberOfThreads \
	--outSAMtype BAM SortedByCoordinate \
	--quantMode TranscriptomeSAM \
	--outSAMheaderHD @HD VN:1.4 SO:coordinate \
	--outFileNamePrefix /path/to/STAR/output/directory/<sample_id> \
	--readFilesIn /path/to/trimmed_forward_reads \
	/path/to/trimmed_reverse_reads # only needed for PE studies

```

**Parameter Definitions:**

* `--twopassMode` – specifies 2-pass mapping mode; the `Basic` option instructs STAR to perform the 1st pass mapping, then automatically extract junctions, insert them into the genome index, and re-map all reads in the 2nd mapping pass
* `--limitBAMsortRAM` - maximum RAM available (in bytes) to sort the bam files, the example above indicates 65GB
* `--genomeDir` - specifies the path to the directory where the STAR reference is stored
* `--outSAMunmapped` - specifies ouput of unmapped reads in the sam format; the `Within` option instructs STAR to output the unmapped reads within the main sam file
* `--outFilterType` - specifies the type of filtering; the `BySJout` option instructs STAR to keep only those reads that contain junctions that passed filtering in the SJ.out.tab output file
* `--outSAMattributes` - list of desired sam attributes in the order desired for the output sam file; sam attribute descriptions can be found [here](https://samtools.github.io/hts-specs/SAMtags.pdf)
* `--outFilterMultimapNmax` – specifies the maximum number of loci the read is allowed to map to; all alignments will be output only if the read maps to no more loci than this value
* `--outFilterMismatchNmax` - maximum number of mismatches allowed to be included in the alignment output
* `--outFilterMismatchNoverReadLmax` - ratio of mismatches to read length allowed to be included in the alignment output; the `0.04` value indicates that up to 4 mismatches are allowed per 100 bases
* `--alignIntronMin` - minimum intron size; a genomic gap is considered an intron if its length is equal to or greater than this value, otherwise it is considered a deletion
* `--alignIntronMax` - maximum intron size
* `--alignMatesGapMax` - maximum genomic distance (in bases) between two mates of paired-end reads; this option should be removed for single-end reads
* `--alignSJoverhangMin` - minimum overhang (i.e. block size) for unannotated spliced alignments
* `--alignSJDBoverhangMin` - minimum overhang (i.e. block size) for annotated spliced alignments
* `--sjdbScore` - additional alignment score for alignments that cross database junctions
* `--readFilesCommand` - specifies command needed to interpret input files; the `zcat` option indicates input files are compressed with gzip and zcat will be used to uncompress the gzipped input files
* `--runThreadN` - indicates the number of threads to be used for STAR alignment and should be set to the number of available cores on the server node
* `--outSAMtype` - specifies desired output format; the `BAM SortedByCoordinate` options specify that the output file will be sorted by coordinate and be in the bam format
* `--quantMode` - specifies the type(s) of quantification desired; the `TranscriptomeSAM` option instructs STAR to output a separate sam/bam file containing alignments to the transcriptome
* `--outSAMheaderHD` - indicates a header line for the sam/bam file
* `--outFileNamePrefix` - specifies the path to and prefix for the output file names; for GeneLab the prefix is the sample id
* `--readFilesIn` - path to input read 1 (forward read) and read 2 (reverse read); for paired-end reads, read 1 and read 2 should be separated by a space; for single-end reads only read 1 should be indicated 


**Input Data:**
- STAR genome reference (output from Step 3)
- *fastq.gz (trimmed reads)

**Output Data:**
- *Aligned.sortedByCoord.out.bam# (sorted mapping to genome)
- *Aligned.toTranscriptome.out.bam# (sorted mapping to transcriptome)
- *Log.final.out# (log file containing alignment info/stats such as reads mapped, etc)
- *Log.out
- *Log.progress.out
- *SJ.out.tab\#
- *_STARgenome (directory containing the following:)
  - sjdbInfo.txt
  - sjdbList.out.tab
- *_STARpass1 (directory containing the following:)
  - Log.final.out
  - SJ.out.tab
- *_STARtmp (directory containing the following:)
  - BAMsort (directory containing subdirectories that are empty – this was the location for temp files that were automatically removed after successful completion)

\#Output files available with RNAseq processed data in the [GLDS repository](https://genelab-data.ndc.nasa.gov/genelab/projects). 

<br>

## 4b. Compile Alignment Logs 

```
multiqc -n align_multiqc -o /path/to/aligned_multiqc/output/directory /path/to/*Log.final.out/files
```

**Parameter Definitions:**

* `-n` - prefix name for output files
* `-o` – the output directory to store results
* `/path/to/*Log.final.out/files` – the directory holding the *Log.final.out output files from the [STAR alignment step](#4a-align-reads-to-reference-genome-with-star), provided as a positional argument

**Input Data:**
- *Log.final.out (log file conting alignment info/stats such as reads mapped, etc. from step 4a)

**Output Data:**
- align_multiqc.html (multiqc report)
- align_multiqc_data (directory containing multiqc data)

<br>

## 4c. Index Aligned Reads 

```
samtools index -@ NumberOfThreads /path/to/*Aligned.sortedByCoord.out.bam/files
```

**Parameter Definitions:**

* `-@` - number of threads available on server node to index alignment files
* `/path/to/*Aligned.sortedByCoord.out.bam/files` – the directory holding the *Aligned.sortedByCoord.out.bam output files from the [STAR alignment step](#4a-align-reads-to-reference-genome-with-star), provided as a positional argument

**Input Data:**
- *Aligned.sortedByCoord.out.bam (sorted mapping to genome file, from [step 4a](#4a-align-reads-to-reference-genome-with-star))

**Output Data:**
- *Aligned.sortedByCoord.out.bam.bai (index of sorted mapping to genome file)

<br>

---

## 5a. Convert GTF to genePred File  

```
gtfToGenePred -geneNameAsName2 \
  /path/to/annotation/gtf/file \
  /path/to/output/genePred/file

```

**Parameter Definitions:**

* `-geneNameAsName2` – instructs gtfToGenePred to use gene_name for the name2 field, instead of the default gene_id
* `/path/to/annotation/gtf/file` – specifies the file(s) containing annotated reference transcripts in the standard gtf format, provided as a positional argument
* `/path/to/output/genePred/file` – specifies the location and name of the output genePred file(s), provided as a positional argument

**Input Data:**
- *.gtf (genome annotation\#)

\#See document(s) in the [GeneLab_Reference_and_Annotation_Files](../GeneLab_Reference_and_Annotation_Files) sub-directory for a list of the ensembl fasta genome sequences and associated gtf annotation files used to generate the RNAseq processed data available in the [GLDS repository](https://genelab-data.ndc.nasa.gov/genelab/projects).

**Output Data:**
- *.genePred (genome annotation in genePred format)

<br>

## 5b. Convert genePred to BED File  

```
genePredToBed /path/to/annotation/genePred/file \
  /path/to/output/BED/file

```

**Parameter Definitions:**

* `/path/to/annotation/genePred/file` – specifies the file(s) containing annotated reference transcripts in the genePred format, provided as a positional argument
* `/path/to/output/BED/file` – specifies the location and name of the output BED file(s), provided as a positional argument

**Input Data:**
- *.genePred (genome annotation in genePred format, from step 5a)

**Output Data:**
- *.bed (genome annotation in BED format)

<br>

---

## 6a. Determine Read Strandedness

```
infer_experiment.py -r /path/to/annotation/BED/file \
	-i /path/to/*Aligned.sortedByCoord.out.bam \
	-s 15000000 > /path/to/*infer_expt.out
```

**Parameter Definitions:**

* `-r` – specifies the path to the reference annotation BED file
* `-i` - specifies the path to the input bam file(s)
* `-s` - specifies the number of reads to be sampled from the input bam file(s), 15M reads are sampled
* `>` - redirects standard output to specified file
* `/path/to/*infer_expt.out` - specifies the location and name of the file containing the infer_experiment standard output

**Input Data:**
- *.bed (genome annotation in BED format, output from step 5b)
- *Aligned.sortedByCoord.out.bam (sorted mapping to genome file, output from [step 4a](#4a-align-reads-to-reference-genome-with-star))
- *Aligned.sortedByCoord.out.bam.bai (index of sorted mapping to genome file, output from step [4c](#4c-index-aligned-read), although not indicated in the command, this file must be present in the same directory as the respective *Aligned.sortedByCoord.out.bam file)

**Output Data:**
- *infer_expt.out (file containing the infer_experiment standard output)

<br>

## 6b. Compile Strandedness Reports

```
multiqc -n infer_exp_multiqc -o /path/to/infer_exp_multiqc/output/directory /path/to/*infer_expt.out/files
```

**Parameter Definitions:**

* `-n` - prefix name for output files
* `-o` – the output directory to store results
* `/path/to/*infer_expt.out/files` – the directory holding the *infer_expt.out output files from the [read strandedness step](#6a-determine-read-strandedness), provided as a positional argument

**Input Data:**
- *infer_expt.out (file containing the infer_experiment standard output, from step 6a)

**Output Data:**
- infer_exp_multiqc.html (multiqc report)
- infer_exp_multiqc_data (directory containing multiqc data)

<br>

---

## 7. Build RSEM Reference

```
rsem-prepare-reference --gtf /path/to/annotation/gtf/file \
	/path/to/genome/fasta/file \
	/path/to/RSEM/genome/directory/RSEM_ref_prefix

```

**Parameter Definitions:**

* `--gtf` – specifies the file(s) containing annotated transcripts in the standard gtf format
* `/path/to/genome/fasta/file` – specifies one or more fasta file(s) containing the genome reference sequences, provided as a positional argument
* `/path/to/RSEM/genome/directory/RSEM_ref_prefix` - specifies the path to the directory where the RSEM reference will be stored and the prefix desired for the RSEM reference files, provided as a positional argument 

**Input Data:**
- *.fasta (genome sequence\#)
- *.gtf (genome annotation\#)

\#See document(s) in the [GeneLab_Reference_and_Annotation_Files](../GeneLab_Reference_and_Annotation_Files) sub-directory for a list of the ensembl fasta genome sequences and associated gtf annotation files used to generate the RNAseq processed data available in the [GLDS repository](https://genelab-data.ndc.nasa.gov/genelab/projects).

**Output Data:**

RSEM genome reference, which consists of the following files:
- RSEM_ref_prefix.chrlist
- RSEM_ref_prefix.grp
- RSEM_ref_prefix.idx.fa
- RSEM_ref_prefix.n2g.idx.fa
- RSEM_ref_prefix.seq
- RSEM_ref_prefix.ti
- RSEM_ref_prefix.transcripts.fa

<br>

---

## 8. Count Aligned Reads with RSEM

```
rsem-calculate-expression --num-threads NumberOfThreads \
	--alignments \
	--bam \
	--paired-end \
	--seed 12345 \
	--estimate-rspd \
	--no-bam-output \
	--strandedness reverse|forward|none \
	/path/to/*Aligned.toTranscriptome.out.bam \
	/path/to/RSEM/genome/directory/RSEM_ref_prefix \
	/path/to/RSEM/counts/output/directory/<sample_id>

```

**Parameter Definitions:**

* `--num-threads` – specifies the number of threads to use
* `--alignments` - indicates that the input file contains alignments in sam, bam, or cram format
* `--bam` - specifies that the input alignments are in bam format
* `--paired-end` - indicates that the input reads are paired-end reads; this option should be removed if the input reads are single-end
* `--seed` - the seed for the random number generators used in calculating posterior mean estimates and credibility intervals; must be a non-negative 32-bit integer
* `--estimate-rspd` - instructs RSEM to estimate the read start position distribution (rspd) from the data 
* `--no-bam-output` - instructs RSEM not to output any bam file
* `--strandedness` - defines the strandedness of the RNAseq reads; the `reverse` option is used if read strandedness (output from [step 6](#6a-determine-read-strandedness)) is antisense, `forward` is used with sense strandedness, and `none` is used if strandedness is half sense half antisense 
* `/path/to/*Aligned.toTranscriptome.out.bam` - specifies path to input bam files, provided as a positional argument
* `/path/to/RSEM/genome/directory/RSEM_ref_prefix` - specifies the path to the directory where the RSEM reference is stored and its prefix, provided as a positional argument
* `/path/to/RSEM/counts/output/directory` – specifies the path to and prefix for the output file names; for GeneLab the prefix is the sample id

**Input Data:**
- RSEM genome reference (output from Step 7)
- *Aligned.toTranscriptome.out.bam (sorted mapping to transcriptome, output from [step 4a](#4a-align-reads-to-reference-genome-with-star))

**Output Data:**
- *genes.results (counts per gene)
- *isoforms.results (counts per isoform)
- *stat (directory containing the following stats files)
	- *cnt
	- *model
	- *theta

<br>

---

## 9. Normalize Read Counts, Perform Differential Gene Expression Analysis, and Add Gene Annotations in R

<br>


### 9a. For Datasets with ERCC Spike-In


```R
## Install R packages if not already installed

install.packages("tidyverse")
source("https://bioconductor.org/biocLite.R")
if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("tximport")
BiocManager::install("DESeq2")
BiocManager::install("Risa")
BiocManager::install("STRINGdb")
BiocManager::install("PANTHER.db")


## Install annotation R packages if not already installed - only the annotation package for the organism that the data were derived from is required

BiocManager::install("org.Hs.eg.db")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("org.Rn.eg.db")
BiocManager::install("org.Dr.eg.db")
BiocManager::install("org.Dm.eg.db")
BiocManager::install("org.Ce.eg.db")
BiocManager::install("org.Sc.sgd.db")
BiocManager::install("org.At.tair.db")
BiocManager::install("org.EcK12.eg.db")
BiocManager::install("MeSH.Bsu.168.eg.db")


##### Set up your environment #####

## Import libraries (tximport, DESeq2, tidyverse, Risa)

library(tximport)
library(DESeq2)
library(tidyverse)
library(Risa)


## Define which organism is used in the study - this should be consistent with the name in the organisms.csv file, which matches the abbreviations used in the Panther database for each organism

organism <- "organism_that_samples_were_derived_from"

## Create a directory for the metadata and in that directory, download the *ISA.zip file for the study being analyzed, which is located in the [GLDS repository](https://genelab-data.ndc.nasa.gov/genelab/projects) under 'STUDY FILES' -> 'Study Metadata Files'

## Define the location of the input data and where the ouput data will be printed to
metadata_dir="/path/to/directory/containing/*ISA.zip/file"
work_dir="/path/to/working/directory/where/script/is/executed/from" ## Must contain organisms.csv file
counts_dir="/path/to/directory/containing/RSEM/counts/files"
norm_output="/path/to/normalized/counts/output/directory"
DGE_output="/path/to/DGE/output/directory"
DGE_output_ERCC="/path/to/ERCC-normalized/DGE/output/directory"


## Set your working directory to the directory containing the *ISA.zip file for the GLDS dataset being processed
setwd(file.path(metadata_dir))


##### Pull all factors for each sample in the study from the metadata in the ISA files #####

td = tempdir()
unzip(Sys.glob(file.path(metadata_dir,"*ISA.zip")), exdir = td)
isa <- readISAtab(path = td)
n = as.numeric(which(isa@assay.technology.types == "RNA Sequencing (RNA-Seq)"))
isa_tabs<-isa@assay.tabs[[n]]@assay.file
factors <- as.data.frame(isa@factors[[1]], stringsAsFactors = FALSE)
colnames(factors)<-paste("factor",1:dim(factors)[2], sep = "_")
compare_csv <- data.frame(sample_id = isa_tabs$`Sample Name`, factors) # note: both the s_* table and a_* table in ISA MUST have rows in the same order, otherwise samples may be assigned to the incorrect factor group

## Create data frame containing all samples and respective factors

study <- as.data.frame(compare_csv[,2:dim(compare_csv)[2]])
colnames(study) <- colnames(compare_csv)[2:dim(compare_csv)[2]]
rownames(study) <- compare_csv[,1]

## Set your working directory to the directory containing the organisms.csv file

setwd(file.path(work_dir))


##### Format groups and indicate the group that each sample belongs to #####

if (dim(study) >= 2){
	group<-apply(study,1,paste,collapse = " & ") # concatenate multiple factors into one condition per sample
} else{
	group<-study[,1]
}
group_names <- paste0("(",group,")",sep = "") # human readable group names
group <- make.names(group) # group naming compatible with R models
names(group) <- group_names
rm(group_names)


##### Format contrasts table, defining pairwise comparisons for all groups #####

contrasts <- combn(levels(factor(group)),2) # generate matrix of pairwise group combinations for comparison
contrast.names <- combn(levels(factor(names(group))),2)
contrast.names <- c(paste(contrast.names[1,],contrast.names[2,],sep = "v"),paste(contrast.names[2,],contrast.names[1,],sep = "v")) # format combinations for output table files names
contrasts <- cbind(contrasts,contrasts[c(2,1),])
colnames(contrasts) <- contrast.names
rm(contrast.names) 


##### Import RSEM raw (gene) count data #####

files <- list.files(file.path(counts_dir),pattern = ".genes.results", full.names = TRUE)

## Reorder the *genes.results files to match the ordering of the ISA samples
files <- files[sapply(rownames(study), function(x)grep(x, files, value=FALSE, fixed=TRUE))]

names(files) <- rownames(study)
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)

## Add 1 to genes with lengths of zero - needed to make DESeqDataSet object 
txi.rsem$length[txi.rsem$length == 0] <- 1


##### Make DESeqDataSet object #####

## Create data frame defining which group each sample belongs to
sampleTable <- data.frame(condition=factor(group))
rownames(sampleTable) <- colnames(txi.rsem$counts)

dds <- DESeqDataSetFromTximport(txi.rsem, sampleTable, ~condition)
summary(dds)


##### Create ERCC unfiltered raw counts table #####

## Note: These data are used internally at GeneLab for QC

## Make a DESeqDataSet object using only unfiltered ERCC genes
ercc_rows_all <- grep("ERCC-",rownames(dds))
ercc_dds_all <- dds[ercc_rows_all,]

## Print ERCC unfiltered raw counts table
ERCC_rawCounts_all = as.data.frame(counts(ercc_dds_all))
write.csv(ERCC_rawCounts_all,file='ERCC_rawCounts_unfiltered.csv')


##### Filter out genes with counts of less than 10 in all samples #####

keepGenes <- rowSums(counts(dds)) > 10
dds <- dds[keepGenes,]
summary(dds)
dim(dds)


##### Prepare filtered data to be normalized with and without considering ERCC genes #####

## Make a DESeqDataSet object using only filtered ERCC genes, which will be used to generate ERCC counts table
ercc_rows <- grep("ERCC-",rownames(dds))
ercc_dds <- dds[ercc_rows,]

## Create list of rows containing ERCC group B genes to use for ERCC-normalization
## Note: ERCC group B genes should be the same concentration in all samples

ercc_rows_gpB <- grep("ERCC-00096|ERCC-00171|ERCC-00009|ERCC-00042|ERCC-00060|ERCC-00035|ERCC-00025|ERCC-00051|ERCC-00053|ERCC-00148|ERCC-00126|ERCC-00034|ERCC-00150|ERCC-00067|ERCC-00031|ERCC-00109|ERCC-00073|ERCC-00158|ERCC-00104|ERCC-00142|ERCC-00138|ERCC-00117|ERCC-00075",rownames(dds))
ercc_dds_gpB <- dds[ercc_rows_gpB,]
summary(ercc_dds_gpB)
dim(ercc_dds_gpB)

## Identify and list samples that do not contain counts for ERCC genes in general and specifically group B ERCC genes
## This is to make sure all samples indeed contain ERCC spike-in
cat("Samples that do not have detectable ERCC spike-ins: ", colnames(ercc_dds[,colSums(counts(ercc_dds))==0]), sep="\n")

cat("Samples that do not have detectable ERCC group B spike-ins: ", colnames(dds[,colSums(counts(ercc_dds_gpB))==0]), sep="\n")

# Create a new group object WITHOUT the samples that don't have detectable ERCC group B spike-ins
remove <- colnames(dds[,colSums(counts(ercc_dds_gpB))==0])

study_sub <- subset(study,!rownames(study) %in% remove) # new study object with non-ERCC-gpB samples removed

## New group object with non-ERCC-gpB samples removed
if (dim(study_sub) >= 2){
  group_sub<-apply(study_sub,1,paste,collapse = " & ") # concatenate multiple factors into one condition per sample
} else{
  group_sub<-study_sub[,1]
}
group_names <- paste0("(",group_sub,")",sep = "") # human readable group names
group_sub <- make.names(group_sub) # group naming compatible with R models
names(group_sub) <- group_names
rm(group_names)

## Create new contrasts object that only reflects the groups in the subset group object
contrasts_sub <- combn(levels(factor(group_sub)),2) # generate matrix of pairwise group combinations for comparison
contrasts_sub.names <- combn(levels(factor(names(group_sub))),2)
contrasts_sub.names <- c(paste(contrasts_sub.names[1,],contrasts_sub.names[2,],sep = "v"),paste(contrasts_sub.names[2,],contrasts_sub.names[1,],sep = "v")) # format combinations for output table files names
contrasts_sub <- cbind(contrasts_sub,contrasts_sub[c(2,1),])
colnames(contrasts_sub) <- contrasts_sub.names
rm(contrasts_sub.names)

## If all samples contain ERCC spike-in (i.e. there are no samples to remove), reassign group_sub, study_sub and contrasts_sub back to the original variable contents
if (length(remove) == 0) {
  group_sub <- group
  study_sub <- study
  contrasts_sub <- contrasts
}

## Print ERCC filtered raw counts table
## Note: These data are used internally at GeneLab for QC
ERCC_rawCounts = as.data.frame(counts(ercc_dds))
write.csv(ERCC_rawCounts,file='ERCC_rawCounts_filtered.csv')




##### Generate a DESeqDataSet object using only non-ERCC genes #####

## dds_1 will be used to generate data without considering ERCC genes
dds_1 <- dds[-c(ercc_rows),] # remove ERCCs from full counts table

## dds_2 will be used to generate data with considering ERCC genes
## First remove samples that do not contain ERCC counts
## Note: All samples should contain ERCC spike-in and thus ERCC counts, if some samples do not contain ERCC counts, those samples should be removed and not used for downstream analysis

dds_2 <- dds[,colSums(counts(ercc_dds_gpB)) > 0] 
sampleTable_sub <- data.frame(condition=factor(group_sub)) # create a new sampleTable with only the kept samples
dds_2$condition <- sampleTable_sub$condition # reassign the condition to the subset condition to keep the groups and levels consistent with the "keep" samples
summary(dds_2)
dim(dds_2)



##### Perform DESeq analysis with and without considering ERCC genes #####

## Run DESeq analysis with ERCC-normalization by replacing size factor object with ERCC size factors for rescaling
## Try first to use the default type="median", but if there is an error (usually due to zeros in genes), use type="poscounts"
## From DESeq2 manual: "The "poscounts" estimator deals with a gene with some zeros, by calculating a modified geometric mean by taking the n-th root of the product of the non-zero counts."

dds_2 <- tryCatch(
      expr = { estimateSizeFactors(dds_2, controlGenes=ercc_rows_gpB) },
      error = function(e) { estimateSizeFactors(dds_2, type="poscounts", controlGenes=ercc_rows_gpB)}
)

dds_2 <- dds_2[-c(ercc_rows),] # remove ERCCs from counts table after normalization
dds_2 <- estimateDispersions(dds_2)
dds_2 <- nbinomWaldTest(dds_2)

## Run DESeq analysis without considering ERCC genes
dds_1 <- DESeq(dds_1)


##### Export unnormalized, normalized, and ERCC-normalized counts as well as the sample table #####

normCounts = as.data.frame(counts(dds_1, normalized=TRUE))
ERCCnormCounts = as.data.frame(counts(dds_2, normalized=TRUE))
setwd(file.path(norm_output))
write.csv(txi.rsem$counts,file='Unnormalized_Counts.csv')
write.csv(normCounts,file='Normalized_Counts.csv')
write.csv(ERCCnormCounts,file='ERCC_Normalized_Counts.csv')
write.csv(sampleTable,file='SampleTable.csv')

setwd(file.path(work_dir))


##### Generate F statistic p-value (similar to ANOVA p-value) using DESeq2 likelihood ratio test (LRT) design #####

## Add 1 to all counts to avoid issues with log transformation 
normCounts <- normCounts +1
ERCCnormCounts <- ERCCnormCounts +1

## For non-ERCC normalized data
dds_1_lrt <- DESeq(dds_1, test = "LRT", reduced = ~ 1)
res_1_lrt <- results(dds_1_lrt)

## For ERCC-normalized data
dds_2_lrt <- DESeq(dds_2, test = "LRT", reduced = ~ 1)
res_2_lrt <- results(dds_2_lrt)


##### Generate annotated DGE tables #####

## Import table with organism db objects for annotation
organism_table <- read.csv(file.path(work_dir,"organisms.csv"))

## Load annotation libraries
## If working through VPN, internet via proxy may need to be defined
#assign("has_internet_via_proxy", TRUE, environment(curl::has_internet))
library(STRINGdb) # for String database annotations
library(PANTHER.db) # for GOSLIM annotations

## Begin building anotation database
ann.dbi <- organism_table$annotations[organism_table$name == organism] # Organism specific gene annotation database
ann.dbi=as.character(ann.dbi)
if(!require(ann.dbi, character.only=TRUE)) {
	BiocManager::install(ann.dbi, ask = FALSE)
	library(ann.dbi, character.only=TRUE)
}


#### Generate annotated DGE table for non-ERCC normalized counts ####

## Start by creating output tables with (non-ERCC) normalized sample expression values

## reduced output table 1 will be used to generate human-readable DGE table
reduced_output_table_1 <- normCounts

## output table 1 will be used to generate computer-readable DGE table, which is used to create GeneLab visualization plots
output_table_1 <- normCounts

## Iterate through Wald Tests to generate pairwise comparisons of all groups
for (i in 1:dim(contrasts)[2]){
	res_1 <- results(dds_1, contrast=c("condition",contrasts[1,i],contrasts[2,i]))
	res_1 <- as.data.frame(res_1@listData)[,c(2,4,5,6)]
	colnames(res_1) <-c(paste0("Log2fc_",colnames(contrasts)[i]),paste0("Stat_",colnames(contrasts)[i]),paste0("P.value_",colnames(contrasts)[i]),paste0("Adj.p.value_",colnames(contrasts)[i]))
	output_table_1<-cbind(output_table_1,res_1)
	reduced_output_table_1 <- cbind(reduced_output_table_1,res_1)
	rm(res_1)
}

## Create annotation table and add gene annotation columns
keytype = "ENSEMBL" # will be different if primary annotations are not ENSEMBL
annot <- data.frame(rownames(output_table_1), stringsAsFactors = FALSE)
colnames(annot)[1]<-keytype
if ("SYMBOL" %in% columns(eval(parse(text = ann.dbi),env=.GlobalEnv))){
	annot$SYMBOL<-mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = rownames(output_table_1),keytype = keytype, column = "SYMBOL", multiVals = "first")
}
if ("GENENAME" %in% columns(eval(parse(text = ann.dbi),env=.GlobalEnv))){
        annot$GENENAME<-mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = rownames(output_table_1),keytype = keytype, column = "GENENAME", multiVals = "first")
}
if ("ENSEMBL" %in% columns(eval(parse(text = ann.dbi),env=.GlobalEnv))){
        annot$ENSEMBL<-mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = rownames(output_table_1),keytype = keytype, column = "ENSEMBL", multiVals = "first")
}
if ("REFSEQ" %in% columns(eval(parse(text = ann.dbi),env=.GlobalEnv))){
        annot$REFSEQ<-mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = rownames(output_table_1),keytype = keytype, column = "REFSEQ", multiVals = "first")
}
if ("ENTREZID" %in% columns(eval(parse(text = ann.dbi),env=.GlobalEnv))){
        annot$ENTREZID<-mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = rownames(output_table_1),keytype = keytype, column = "ENTREZID", multiVals = "first")
}

## Create and add string annotation columns to the annotation table
string_db <- STRINGdb$new( version="11", species=organism_table$taxon[organism_table$name == organism],score_threshold=0)
string_map <- string_db$map(annot,"SYMBOL",removeUnmappedRows = FALSE, takeFirst = TRUE)[,c(1,6)]
string_map <- string_map[!duplicated(string_map$SYMBOL),]
annot <- dplyr::left_join(annot,string_map, by = "SYMBOL")

## Create and add columns containing GOSLIM ids using the panther annotation database to the annotation table
pthOrganisms(PANTHER.db) <- organism
panther <- mapIds(PANTHER.db,keys = annot$ENTREZID,keytype = "ENTREZ",column = "GOSLIM_ID", multiVals = "list")
panther <- na.omit(panther)
annot$GOSLIM_IDS <- panther
rm(string_db,string_map,panther,keytype)

## Generate and add all sample mean column to the (non-ERCC) normalized counts table
output_table_1$All.mean <- rowMeans(normCounts, na.rm = TRUE, dims = 1)
reduced_output_table_1$All.mean <- rowMeans(normCounts, na.rm = TRUE, dims = 1)

## Generate and add all sample stdev column to the (non-ERCC) normalized counts table
output_table_1$All.stdev <- rowSds(as.matrix(normCounts), na.rm = TRUE, dims = 1)
reduced_output_table_1$All.stdev <- rowSds(as.matrix(normCounts), na.rm = TRUE, dims = 1)

## Add F statistic p-value (similar to ANOVA p-value) column to the (non-ERCC) normalized counts table
output_table_1$LRT.p.value <- res_1_lrt@listData$padj
reduced_output_table_1$LRT.p.value <- res_1_lrt@listData$padj

## Generate and add group mean and stdev columns to the (non-ERCC) normalized counts table
tcounts <- as.data.frame(t(normCounts))
tcounts$group <- group
group_means <- as.data.frame(t(aggregate(. ~ group,data = tcounts,mean)))
group_means <- group_means[-c(1),]
colnames(group_means) <- paste0("Group.Mean_",levels(factor(names(group))))
group_stdev <- as.data.frame(t(aggregate(. ~ group,data = tcounts,sd)))
group_stdev <- group_stdev[-c(1),]
colnames(group_stdev) <- paste0("Group.Stdev_",levels(factor(names(group))))

output_table_1 <- cbind(output_table_1,group_means)
reduced_output_table_1 <- cbind(reduced_output_table_1,group_means)

output_table_1 <- cbind(output_table_1,group_stdev)
reduced_output_table_1 <- cbind(reduced_output_table_1,group_stdev)

rm(group_stdev,group_means,tcounts)

### Add columns needed to generate GeneLab visulaization plots to the (non-ERCC) normalized counts table

## Add column to indicate the sign (positive/negative) of log2fc for each pairwise comparison
updown_table <- sign(output_table_1[,grep("Log2fc_",colnames(output_table_1))])
colnames(updown_table) <- gsub("Log2fc","Updown",grep("Log2fc_",colnames(output_table_1),value = TRUE))
output_table_1 <- cbind(output_table_1,updown_table)
rm(updown_table)

## Add column to indicate contrast significance with p <= 0.1
sig.1_table <- output_table_1[,grep("P.value_",colnames(output_table_1))]<=.1
colnames(sig.1_table) <- gsub("P.value","Sig.1",grep("P.value_",colnames(output_table_1),value = TRUE))
output_table_1 <- cbind(output_table_1,sig.1_table)
rm(sig.1_table)

## Add column to indicate contrast significance with p <= 0.05
sig.05_table <- output_table_1[,grep("P.value_",colnames(output_table_1))]<=.05
colnames(sig.05_table) <- gsub("P.value","Sig.05",grep("P.value_",colnames(output_table_1),value = TRUE))
output_table_1 <- cbind(output_table_1,sig.05_table)
rm(sig.05_table)

## Add columns for the volcano plot with p-value and adjusted p-value
log_pval_table <- log2(output_table_1[,grep("P.value_",colnames(output_table_1))])
colnames(log_pval_table) <- paste0("Log2_",colnames(log_pval_table))
output_table_1 <- cbind(output_table_1,log_pval_table)
rm(log_pval_table)
log_adj_pval_table <- log2(output_table_1[,grep("Adj.p.value_",colnames(output_table_1))])
colnames(log_adj_pval_table) <- paste0("Log2_",colnames(log_adj_pval_table))
output_table_1 <- cbind(output_table_1,log_adj_pval_table)
rm(log_adj_pval_table)

### Combine annotations table and the (non-ERCC) normalized counts table

output_table_1 <- cbind(annot,output_table_1)
reduced_output_table_1 <- cbind(annot,reduced_output_table_1)
rownames(output_table_1) <- NULL
rownames(reduced_output_table_1) <- NULL

output_table_1$GOSLIM_IDS <- vapply(output_table_1$GOSLIM_IDS, paste, collapse = ", ", character(1L))
reduced_output_table_1$GOSLIM_IDS <- vapply(reduced_output_table_1$GOSLIM_IDS, paste, collapse = ", ", character(1L))

### Export human- and computer/visualization- readable DGE tables

write.csv(output_table_1,file.path(DGE_output, "visualization_output_table.csv"), row.names = FALSE)
write.csv(contrasts,file.path(DGE_output, "contrasts.csv"))
write.csv(reduced_output_table_1,file.path(DGE_output, "differential_expression.csv"), row.names = FALSE)

### Generate and export PCA table for GeneLab visualization plots

exp_raw <- log2(normCounts)
PCA_raw <- prcomp(t(exp_raw), scale = FALSE)
write.csv(PCA_raw$x,file.path(DGE_output, "visualization_PCA_table.csv"), row.names = TRUE)
rm(exp_raw,PCA_raw)


#### Generate annotated DGE table for ERCC-normalized counts ####

## Start by creating output tables with ERCC-normalized sample expression values

## reduced output table 2 will be used to generate human-readable DGE table
reduced_output_table_2 <- ERCCnormCounts

## output table 2 will be used to generate computer-readable DGE table, which is used to create GeneLab visualization plots
output_table_2 <- ERCCnormCounts

## Iterate through Wald Tests to generate pairwise comparisons of all groups
for (i in 1:dim(contrasts_sub)[2]){
  res_2 <- results(dds_2, contrast=c("condition",contrasts_sub[1,i],contrasts_sub[2,i]))
  res_2 <- as.data.frame(res_2@listData)[,c(2,4,5,6)]
  colnames(res_2)<-c(paste0("Log2fc_",colnames(contrasts_sub)[i]),paste0("Stat_",colnames(contrasts_sub)[i]),paste0("P.value_",colnames(contrasts_sub)[i]),paste0("Adj.p.value_",colnames(contrasts_sub)[i]))
  output_table_2<-cbind(output_table_2,res_2)
  reduced_output_table_2 <- cbind(reduced_output_table_2,res_2)
  rm(res_2)
}


## Create annotation table and add gene annotation columns
keytype = "ENSEMBL" # will be different if primary annotations are not ENSEMBL
annot <- data.frame(rownames(output_table_2), stringsAsFactors = FALSE)
colnames(annot)[1]<-keytype
if ("SYMBOL" %in% columns(eval(parse(text = ann.dbi),env=.GlobalEnv))){
        annot$SYMBOL<-mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = rownames(output_table_2),keytype = keytype, column = "SYMBOL", multiVals = "first")
}
if ("GENENAME" %in% columns(eval(parse(text = ann.dbi),env=.GlobalEnv))){
        annot$GENENAME<-mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = rownames(output_table_2),keytype = keytype, column = "GENENAME", multiVals = "first")
}
if ("ENSEMBL" %in% columns(eval(parse(text = ann.dbi),env=.GlobalEnv))){
        annot$ENSEMBL<-mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = rownames(output_table_2),keytype = keytype, column = "ENSEMBL", multiVals = "first")
}
if ("REFSEQ" %in% columns(eval(parse(text = ann.dbi),env=.GlobalEnv))){
        annot$REFSEQ<-mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = rownames(output_table_2),keytype = keytype, column = "REFSEQ", multiVals = "first")
}
if ("ENTREZID" %in% columns(eval(parse(text = ann.dbi),env=.GlobalEnv))){
        annot$ENTREZID<-mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = rownames(output_table_2),keytype = keytype, column = "ENTREZID", multiVals = "first")
}

## Create and add string annotation columns to the annotation table
string_db <- STRINGdb$new( version="11", species=organism_table$taxon[organism_table$name == organism],score_threshold=0)
string_map <- string_db$map(annot,"SYMBOL",removeUnmappedRows = FALSE, takeFirst = TRUE)[,c(1,6)]
string_map <- string_map[!duplicated(string_map$SYMBOL),]
annot <- dplyr::left_join(annot,string_map, by = "SYMBOL")

## Create and add columns containing GOSLIM ids using the panther annotation database to the annotation table
pthOrganisms(PANTHER.db) <- organism
panther <- mapIds(PANTHER.db,keys = annot$ENTREZID,keytype = "ENTREZ",column = "GOSLIM_ID", multiVals = "list")
panther <- na.omit(panther)
annot$GOSLIM_IDS <- panther
rm(string_db,string_map,panther,keytype)

## Generate and add all sample mean column to the ERCC-normalized counts table
output_table_2$All.mean <- rowMeans(ERCCnormCounts, na.rm = TRUE, dims = 1)
reduced_output_table_2$All.mean <- rowMeans(ERCCnormCounts, na.rm = TRUE, dims = 1)

## Generate and add all sample stdev column to the ERCC-normalized counts table
output_table_2$All.stdev <- rowSds(as.matrix(ERCCnormCounts), na.rm = TRUE, dims = 1)
reduced_output_table_2$All.stdev <- rowSds(as.matrix(ERCCnormCounts), na.rm = TRUE, dims = 1)

## Add F statistic p-value (similar to ANOVA p-value) column to the ERCC-normalized counts table
output_table_2$LRT.p.value <- res_2_lrt@listData$padj
reduced_output_table_2$LRT.p.value <- res_2_lrt@listData$padj

## Generate and add group mean and stdev columns to the ERCC-normalized counts table
tcounts <- as.data.frame(t(ERCCnormCounts))
tcounts$group_sub <- group_sub
group_means <- as.data.frame(t(aggregate(. ~ group_sub,data = tcounts,mean)))
group_means <- group_means[-c(1),]
colnames(group_means) <- paste0("Group.Mean_",levels(factor(names(group_sub))))
group_stdev <- as.data.frame(t(aggregate(. ~ group_sub,data = tcounts,sd)))
group_stdev <- group_stdev[-c(1),]
colnames(group_stdev) <- paste0("Group.Stdev_",levels(factor(names(group_sub))))

output_table_2 <- cbind(output_table_2,group_means)
reduced_output_table_2 <- cbind(reduced_output_table_2,group_means)

output_table_2 <- cbind(output_table_2,group_stdev)
reduced_output_table_2 <- cbind(reduced_output_table_2,group_stdev)

rm(group_stdev,group_means,tcounts)

### Add columns needed to generate GeneLab visulaization plots to the ERCC-normalized counts table

## Add column to indicate the sign (positive/negative) of log2fc for each pairwise comparison
updown_table <- sign(output_table_2[,grep("Log2fc_",colnames(output_table_2))])
colnames(updown_table) <- gsub("Log2fc","Updown",grep("Log2fc_",colnames(output_table_2),value = TRUE))
output_table_2 <- cbind(output_table_2,updown_table)
rm(updown_table)

## Add column to indicate contrast significance with p <= 0.1
sig.1_table <- output_table_2[,grep("P.value_",colnames(output_table_2))]<=.1
colnames(sig.1_table) <- gsub("P.value","Sig.1",grep("P.value_",colnames(output_table_2),value = TRUE))
output_table_2 <- cbind(output_table_2,sig.1_table)
rm(sig.1_table)

## Add column to indicate contrast significance with p <= 0.05
sig.05_table <- output_table_2[,grep("P.value_",colnames(output_table_2))]<=.05
colnames(sig.05_table) <- gsub("P.value","Sig.05",grep("P.value_",colnames(output_table_2),value = TRUE))
output_table_2 <- cbind(output_table_2,sig.05_table)
rm(sig.05_table)

## Add columns for the volcano plot with p-value and adjusted p-value
log_pval_table <- log2(output_table_2[,grep("P.value_",colnames(output_table_2))])
colnames(log_pval_table) <- paste0("Log2_",colnames(log_pval_table))
output_table_2 <- cbind(output_table_2,log_pval_table)
rm(log_pval_table)
log_adj_pval_table <- log2(output_table_2[,grep("Adj.p.value_",colnames(output_table_2))])
colnames(log_adj_pval_table) <- paste0("Log2_",colnames(log_adj_pval_table))
output_table_2 <- cbind(output_table_2,log_adj_pval_table)
rm(log_adj_pval_table)

### Combine annotations table and the ERCC-normalized counts table

output_table_2 <- cbind(annot,output_table_2)
reduced_output_table_2 <- cbind(annot,reduced_output_table_2)
rownames(output_table_2) <- NULL
rownames(reduced_output_table_2) <- NULL

output_table_2$GOSLIM_IDS <- vapply(output_table_2$GOSLIM_IDS, paste, collapse = ", ", character(1L))
reduced_output_table_2$GOSLIM_IDS <- vapply(reduced_output_table_2$GOSLIM_IDS, paste, collapse = ", ", character(1L))

### Export human- and computer/visualization- readable DGE tables

write.csv(output_table_2,file.path(DGE_output_ERCC, "visualization_output_table_ERCCnorm.csv"), row.names = FALSE)
write.csv(contrasts,file.path(DGE_output_ERCC, "ERCCnorm_contrasts.csv"))
write.csv(reduced_output_table_2,file.path(DGE_output_ERCC, "ERCCnorm_differential_expression.csv"), row.names = FALSE)

### Generate and export PCA table for GeneLab visualization plots

exp_raw <- log2(ERCCnormCounts)
PCA_raw <- prcomp(t(exp_raw), scale = FALSE)
write.csv(PCA_raw$x,file.path(DGE_output_ERCC, "visualization_PCA_table_ERCCnorm.csv"), row.names = TRUE)
rm(exp_raw,PCA_raw)


## print session info ##
print("Session Info below: ")
sessionInfo()

```

**Input Data:**
- *ISA.zip (compressed ISA directory containing Investigation, Study, and Assay (ISA) metadata files for the respective GLDS dataset, used to define sample groups - the *ISA.zip file is located in the [GLDS repository](https://genelab-data.ndc.nasa.gov/genelab/projects) under 'STUDY FILES' -> 'Study Metadata Files')
- [organisms.csv](../organisms.csv) (csv file containing short name, species name, taxon ID, and annotation db object of model organisms hosted on GeneLab)
- *genes.results (RSEM counts per gene, output from step 6)

**Output Data:**

Output data without considering ERCC spike-in genes: 
- Unnormalized_Counts.csv\#
- Normalized_Counts.csv\#
- SampleTable.csv\#
- visualization_output_table.csv\#\#
- visualization_PCA_table.csv\#\#
- differential_expression.csv\#
- contrasts.csv\#

Output data with considering ERCC spike-in genes:
- ERCC_rawCounts_unfiltered.csv (ERCC unfiltered counts table)
- ERCC_rawCounts_filtered.csv (ERCC counts table after removing ERCC genes with low counts)
- ERCC_Normalized_Counts.csv\#
- visualization_output_table_ERCCnorm.csv\#\#
- visualization_PCA_table_ERCCnorm.csv\#\#
- ERCCnorm_differential_expression.csv\#
- ERCCnorm_contrasts.csv\#

\#Output files available with RNAseq processed data in the [GLDS repository](https://genelab-data.ndc.nasa.gov/genelab/projects).

\#\#Output files used to generate interactive plots from RNAseq processed data in the [GLDS visualization portal](https://genelab-data.ndc.nasa.gov/genelab/projects?page=1&paginate_by=25&viz=true).

<br>

---


### 9b. For Datasets without ERCC Spike-In


```R
## Install R packages if not already installed

install.packages("tidyverse")
source("https://bioconductor.org/biocLite.R")
if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("tximport")
BiocManager::install("DESeq2")
BiocManager::install("Risa")
BiocManager::install("STRINGdb")
BiocManager::install("PANTHER.db")


## Install annotation R packages if not already installed - only the annotation package for the organism that the data were derived from is required

BiocManager::install("org.Hs.eg.db")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("org.Rn.eg.db")
BiocManager::install("org.Dr.eg.db")
BiocManager::install("org.Dm.eg.db")
BiocManager::install("org.Ce.eg.db")
BiocManager::install("org.Sc.sgd.db")
BiocManager::install("org.At.tair.db")
BiocManager::install("org.EcK12.eg.db")
BiocManager::install("MeSH.Bsu.168.eg.db")


##### Set up your environment #####

## Import libraries (tximport, DESeq2, tidyverse, Risa)

library(tximport)
library(DESeq2)
library(tidyverse)
library(Risa)


## Define which organism is used in the study - this should be consistent with the name in the organisms.csv file, which matches the abbreviations used in the Panther database for each organism

organism <- "organism_that_samples_were_derived_from"

## Create a directory for the metadata and in that directory, download the *ISA.zip file for the study being analyzed, which is located in the [GLDS repository](https://genelab-data.ndc.nasa.gov/genelab/projects) under 'STUDY FILES' -> 'Study Metadata Files'

## Define the location of the input data and where the ouput data will be printed to
metadata_dir="/path/to/directory/containing/*ISA.zip/file"
work_dir="/path/to/working/directory/where/script/is/executed/from" ## Must contain organisms.csv file
counts_dir="/path/to/directory/containing/RSEM/counts/files"
norm_output="/path/to/normalized/counts/output/directory"
DGE_output="/path/to/DGE/output/directory"


## Set your working directory to the directory containing the *ISA.zip file for the GLDS being processed
setwd(file.path(metadata_dir))


##### Pull all factors for each sample in the study from the metadata in the ISA files #####

td = tempdir()
unzip(Sys.glob(file.path(metadata_dir,"*ISA.zip")), exdir = td)
isa <- readISAtab(path = td)
n = as.numeric(which(isa@assay.technology.types == "RNA Sequencing (RNA-Seq)"))
isa_tabs<-isa@assay.tabs[[n]]@assay.file
factors <- as.data.frame(isa@factors[[1]], stringsAsFactors = FALSE)
colnames(factors)<-paste("factor",1:dim(factors)[2], sep = "_")
compare_csv <- data.frame(sample_id = isa_tabs$`Sample Name`, factors) # note: both the s_* table and a_* table in ISA MUST have rows in the same order, otherwise samples may be assigned to the incorrect factor group

## Create data frame containing all samples and respective factors

study <- as.data.frame(compare_csv[,2:dim(compare_csv)[2]])
colnames(study) <- colnames(compare_csv)[2:dim(compare_csv)[2]]
rownames(study) <- compare_csv[,1]

## Set your working directory to the directory containing the organisms.csv file

setwd(file.path(work_dir))


##### Format groups and indicate the group that each sample belongs to #####

if (dim(study) >= 2){
	group<-apply(study,1,paste,collapse = " & ") # concatenate multiple factors into one condition per sample
} else{
	group<-study[,1]
}
group_names <- paste0("(",group,")",sep = "") # human readable group names
group <- make.names(group) # group naming compatible with R models
names(group) <- group_names
rm(group_names)


##### Format contrasts table, defining pairwise comparisons for all groups #####

contrasts <- combn(levels(factor(group)),2) # generate matrix of pairwise group combinations for comparison
contrast.names <- combn(levels(factor(names(group))),2)
contrast.names <- c(paste(contrast.names[1,],contrast.names[2,],sep = "v"),paste(contrast.names[2,],contrast.names[1,],sep = "v")) # format combinations for output table files names
contrasts <- cbind(contrasts,contrasts[c(2,1),])
colnames(contrasts) <- contrast.names
rm(contrast.names) 


##### Import RSEM raw (gene) count data #####

files <- list.files(file.path(counts_dir),pattern = ".genes.results", full.names = TRUE)

## Reorder the *genes.results files to match the ordering of the ISA samples
files <- files[sapply(rownames(study), function(x)grep(x, files, value=FALSE, fixed=TRUE))]

names(files) <- rownames(study)
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)

## Add 1 to genes with lengths of zero - needed to make DESeqDataSet object 
txi.rsem$length[txi.rsem$length == 0] <- 1


##### Make DESeqDataSet object #####

## Create data frame defining which group each sample belongs to
sampleTable <- data.frame(condition=factor(group))
rownames(sampleTable) <- colnames(txi.rsem$counts)

dds <- DESeqDataSetFromTximport(txi.rsem, sampleTable, ~condition)
summary(dds)


##### Filter out genes with counts of less than 10 in all samples #####

keepGenes <- rowSums(counts(dds)) > 10
dds <- dds[keepGenes,]
summary(dds)
dim(dds)


##### Perform DESeq analysis #####

dds_1 <- DESeq(dds)


##### Export unnormalized and normalized counts as well as the sample table #####

normCounts = as.data.frame(counts(dds_1, normalized=TRUE))
setwd(file.path(norm_output))
write.csv(txi.rsem$counts,file='Unnormalized_Counts.csv')
write.csv(normCounts,file='Normalized_Counts.csv')
write.csv(sampleTable,file='SampleTable.csv')

setwd(file.path(work_dir))


##### Generate F statistic p-value (similar to ANOVA p-value) using DESeq2 likelihood ratio test (LRT) design #####

## Add 1 to all counts to avoid issues with log transformation 
normCounts <- normCounts +1

dds_1_lrt <- DESeq(dds_1, test = "LRT", reduced = ~ 1)
res_1_lrt <- results(dds_1_lrt)


##### Generate annotated DGE tables #####

## Import table with organism db objects for annotation
organism_table <- read.csv(file.path(work_dir,"organisms.csv"))

## Load annotation libraries
## If working through VPN, internet via proxy may need to be defined
#assign("has_internet_via_proxy", TRUE, environment(curl::has_internet))
library(STRINGdb) # for String database annotations
library(PANTHER.db) # for GOSLIM annotations

## Begin building anotation database
ann.dbi <- organism_table$annotations[organism_table$name == organism] # Organism specific gene annotation database
ann.dbi=as.character(ann.dbi)
if(!require(ann.dbi, character.only=TRUE)) {
	BiocManager::install(ann.dbi, ask = FALSE)
	library(ann.dbi, character.only=TRUE)
}


#### Generate annotated DGE table containing normalized counts ####

## Start by creating output tables with normalized sample expression values

## reduced output table 1 will be used to generate human-readable DGE table
reduced_output_table_1 <- normCounts

## output table 1 will be used to generate computer-readable DGE table, which is used to create GeneLab visualization plots
output_table_1 <- normCounts

## Iterate through Wald Tests to generate pairwise comparisons of all groups
for (i in 1:dim(contrasts)[2]){
  res_1 <- results(dds_1, contrast=c("condition",contrasts[1,i],contrasts[2,i]))
  res_1 <- as.data.frame(res_1@listData)[,c(2,4,5,6)]
  colnames(res_1)<-c(paste0("Log2fc_",colnames(contrasts)[i]),paste0("Stat_",colnames(contrasts)[i]),paste0("P.value_",colnames(contrasts)[i]),paste0("Adj.p.value_",colnames(contrasts)[i]))
  output_table_1<-cbind(output_table_1,res_1)
  reduced_output_table_1 <- cbind(reduced_output_table_1,res_1)
  rm(res_1)
}

## Create annotation table and add gene annotation columns
keytype = "ENSEMBL" # will be different if primary annotations are not ENSEMBL
annot <- data.frame(rownames(output_table_1), stringsAsFactors = FALSE)
colnames(annot)[1]<-keytype
if ("SYMBOL" %in% columns(eval(parse(text = ann.dbi),env=.GlobalEnv))){
	annot$SYMBOL<-mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = rownames(output_table_1),keytype = keytype, column = "SYMBOL", multiVals = "first")
}
if ("GENENAME" %in% columns(eval(parse(text = ann.dbi),env=.GlobalEnv))){
        annot$GENENAME<-mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = rownames(output_table_1),keytype = keytype, column = "GENENAME", multiVals = "first")
}
if ("ENSEMBL" %in% columns(eval(parse(text = ann.dbi),env=.GlobalEnv))){
        annot$ENSEMBL<-mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = rownames(output_table_1),keytype = keytype, column = "ENSEMBL", multiVals = "first")
}
if ("REFSEQ" %in% columns(eval(parse(text = ann.dbi),env=.GlobalEnv))){
        annot$REFSEQ<-mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = rownames(output_table_1),keytype = keytype, column = "REFSEQ", multiVals = "first")
}
if ("ENTREZID" %in% columns(eval(parse(text = ann.dbi),env=.GlobalEnv))){
        annot$ENTREZID<-mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = rownames(output_table_1),keytype = keytype, column = "ENTREZID", multiVals = "first")
}

## Create and add string annotation columns to the annotation table
string_db <- STRINGdb$new( version="11", species=organism_table$taxon[organism_table$name == organism],score_threshold=0)
string_map <- string_db$map(annot,"SYMBOL",removeUnmappedRows = FALSE, takeFirst = TRUE)[,c(1,6)]
string_map <- string_map[!duplicated(string_map$SYMBOL),]
annot <- dplyr::left_join(annot,string_map, by = "SYMBOL")

## Create and add columns containing GOSLIM ids using the panther annotation database to the annotation table
pthOrganisms(PANTHER.db) <- organism
panther <- mapIds(PANTHER.db,keys = annot$ENTREZID,keytype = "ENTREZ",column = "GOSLIM_ID", multiVals = "list")
panther <- na.omit(panther)
annot$GOSLIM_IDS <- panther
rm(string_db,string_map,panther,keytype)

## Generate and add all sample mean column to the normalized counts table
output_table_1$All.mean <- rowMeans(normCounts, na.rm = TRUE, dims = 1)
reduced_output_table_1$All.mean <- rowMeans(normCounts, na.rm = TRUE, dims = 1)

## Generate and add all sample stdev column to the normalized counts table
output_table_1$All.stdev <- rowSds(as.matrix(normCounts), na.rm = TRUE, dims = 1)
reduced_output_table_1$All.stdev <- rowSds(as.matrix(normCounts), na.rm = TRUE, dims = 1)

## Add F statistic p-value (similar to ANOVA p-value) column to the normalized counts table
output_table_1$LRT.p.value <- res_1_lrt@listData$padj
reduced_output_table_1$LRT.p.value <- res_1_lrt@listData$padj

## Generate and add group mean and stdev columns to the normalized counts table
tcounts <- as.data.frame(t(normCounts))
tcounts$group <- group
group_means <- as.data.frame(t(aggregate(. ~ group,data = tcounts,mean)))
group_means <- group_means[-c(1),]
colnames(group_means) <- paste0("Group.Mean_",levels(factor(names(group))))
group_stdev <- as.data.frame(t(aggregate(. ~ group,data = tcounts,sd)))
group_stdev <- group_stdev[-c(1),]
colnames(group_stdev) <- paste0("Group.Stdev_",levels(factor(names(group))))

output_table_1 <- cbind(output_table_1,group_means)
reduced_output_table_1 <- cbind(reduced_output_table_1,group_means)

output_table_1 <- cbind(output_table_1,group_stdev)
reduced_output_table_1 <- cbind(reduced_output_table_1,group_stdev)

rm(group_stdev,group_means,tcounts)

### Add columns needed to generate GeneLab visulaization plots to the normalized counts table

## Add column to indicate the sign (positive/negative) of log2fc for each pairwise comparison
updown_table <- sign(output_table_1[,grep("Log2fc_",colnames(output_table_1))])
colnames(updown_table) <- gsub("Log2fc","Updown",grep("Log2fc_",colnames(output_table_1),value = TRUE))
output_table_1 <- cbind(output_table_1,updown_table)
rm(updown_table)

## Add column to indicate contrast significance with p <= 0.1
sig.1_table <- output_table_1[,grep("P.value_",colnames(output_table_1))]<=.1
colnames(sig.1_table) <- gsub("P.value","Sig.1",grep("P.value_",colnames(output_table_1),value = TRUE))
output_table_1 <- cbind(output_table_1,sig.1_table)
rm(sig.1_table)

## Add column to indicate contrast significance with p <= 0.05
sig.05_table <- output_table_1[,grep("P.value_",colnames(output_table_1))]<=.05
colnames(sig.05_table) <- gsub("P.value","Sig.05",grep("P.value_",colnames(output_table_1),value = TRUE))
output_table_1 <- cbind(output_table_1,sig.05_table)
rm(sig.05_table)

## Add columns for the volcano plot with p-value and adjusted p-value
log_pval_table <- log2(output_table_1[,grep("P.value_",colnames(output_table_1))])
colnames(log_pval_table) <- paste0("Log2_",colnames(log_pval_table))
output_table_1 <- cbind(output_table_1,log_pval_table)
rm(log_pval_table)
log_adj_pval_table <- log2(output_table_1[,grep("Adj.p.value_",colnames(output_table_1))])
colnames(log_adj_pval_table) <- paste0("Log2_",colnames(log_adj_pval_table))
output_table_1 <- cbind(output_table_1,log_adj_pval_table)
rm(log_adj_pval_table)

### Combine annotations table and the normalized counts table

output_table_1 <- cbind(annot,output_table_1)
reduced_output_table_1 <- cbind(annot,reduced_output_table_1)
rownames(output_table_1) <- NULL
rownames(reduced_output_table_1) <- NULL

output_table_1$GOSLIM_IDS <- vapply(output_table_1$GOSLIM_IDS, paste, collapse = ", ", character(1L))
reduced_output_table_1$GOSLIM_IDS <- vapply(reduced_output_table_1$GOSLIM_IDS, paste, collapse = ", ", character(1L))

### Export human- and computer/visualization- readable DGE tables

write.csv(output_table_1,file.path(DGE_output, "visualization_output_table.csv"), row.names = FALSE)
write.csv(contrasts,file.path(DGE_output, "contrasts.csv"))
write.csv(reduced_output_table_1,file.path(DGE_output, "differential_expression.csv"), row.names = FALSE)

### Generate and export PCA table for GeneLab visualization plots

exp_raw <- log2(normCounts)
PCA_raw <- prcomp(t(exp_raw), scale = FALSE)
write.csv(PCA_raw$x,file.path(DGE_output, "visualization_PCA_table.csv"), row.names = TRUE)
rm(exp_raw,PCA_raw)


## print session info ##
print("Session Info below: ")
sessionInfo()

```

**Input Data:**
- *ISA.zip (compressed ISA directory containing Investigation, Study, and Assay (ISA) metadata files for the respective GLDS dataset, used to define sample groups - the *ISA.zip file is located in the [GLDS repository](https://genelab-data.ndc.nasa.gov/genelab/projects) under 'STUDY FILES' -> 'Study Metadata Files')
- [organisms.csv](../organisms.csv) (csv file containing short name, species name, taxon ID, and annotation db object of model organisms hosted on GeneLab)
- *genes.results (RSEM counts per gene, output from step 6)

**Output Data:**

- Unnormalized_Counts.csv\#
- Normalized_Counts.csv\#
- SampleTable.csv\#
- visualization_output_table.csv\#\#
- visualization_PCA_table.csv\#\#
- differential_expression.csv\#
- contrasts.csv\#

\#Output files available with RNAseq processed data in the [GLDS repository](https://genelab-data.ndc.nasa.gov/genelab/projects).

\#\#Output files used to generate interactive plots from RNAseq processed data in the [GLDS visualization portal](https://genelab-data.ndc.nasa.gov/genelab/projects?page=1&paginate_by=25&viz=true).

<br>
