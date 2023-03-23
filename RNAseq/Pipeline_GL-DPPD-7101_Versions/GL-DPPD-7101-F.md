# GeneLab bioinformatics processing pipeline for Illumina RNA-sequencing data

> **This page holds an overview and instructions for how GeneLab processes RNAseq datasets. Exact processing commands, GL-DPPD-7101 version used, and processed data output files for specific datasets are provided in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/).**  

---

**Date:** August 18, 2022  
**Revision:** F  
**Document Number:** GL-DPPD-7101-F  

**Submitted by:**  
Jonathan Oribello (GeneLab Data Processing Team)

**Approved by:**  
Amanda Saravia-Butler (GeneLab Data Processing Lead)  
Sylvain Costes (GeneLab Project Manager)  
Samrawit Gebre (GeneLab Deputy Project Manager and Interim GeneLab Configuration Manager)  
Jonathan Galazka (GeneLab Project Scientist)

---

## Updates from previous version

Updated [Ensembl Reference Files](../../GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110/GL-DPPD-7110_annotations.csv) now use:
- Animals: Ensembl release 107
- Plants: Ensembl plants release 54
- Bacteria: Ensembl bacteria release 54

The DESeq2 Normalization and DGE step, [step 9](#9-normalize-read-counts-perform-differential-gene-expression-analysis-and-add-gene-annotations-in-r), was modified as follows:

- A separate sub-step, [step 9a](#9a-create-sample-runsheet), was added to use the [dp_tools](https://github.com/J-81/dp_tools) program to create a runsheet containing all the metadata needed for running DESeq2, including ERCC spike-in status and sample grouping. This runsheet is imported in the DESeq2 script in place of parsing the ISA.zip file associated with the GLDS dataset. 

- GeneLab now creates a custom reference annotation table as detailed in the [GeneLab_Reference_Annotations](../../GeneLab_Reference_Annotations) directory. The GeneLab Reference Annotation tables for each model organism created with [version GL-DPPD-7110](../../GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110) is now imported in the DESeq2 script to add gene annotations in [step 9f](#9f-prepare-genelab-dge-tables-with-annotations-on-datasets-with-ercc-spike-in) and [step 9i](#9i-prepare-genelab-dge-tables-with-annotations-on-datasets-without-ercc-spike-in). 

- Added the `ERCCnorm_SampleTable.csv` output file in [step 9g](#9g-export-genelab-dge-tables-with-annotations-for-datasets-with-ercc-spike-in) to indicate the samples used in the DESeq2 Normalization and DGE step for datasets with ERCC spike-in.
  > Note: In most cases, the ERCCnorm_SampleTable.csv and SampleTable.csv files are the same. They will only differ when, for the ERCC-based analysis, samples are removed due to a lack of detectable Group B ERCC spike-in genes.

- Fixed edge case where `contrasts.csv` and `ERCCnorm_contrasts.csv` table header and rows could become out of sync with each other in [step 9c](#9c-configure-metadata-sample-grouping-and-group-comparisons) and [step 9e](#9e-perform-dge-on-datasets-with-ercc-spike-in) by generating rows from header rather than generating both separately.

- Updated R version from 4.1.2 to 4.1.3.

- Fixed edge case where certain scripts would crash if sample names were prefixes of other sample names. This had affected [step 4c](#4c-tablulate-star-counts-in-r), [step 8c](#8c-calculate-total-number-of-genes-expressed-per-sample-in-r), and [step 9d](#9d-import-rsem-genecounts).

- Fixed rare edge case where groupwise mean and standard deviations could become misassociated to incorrect groups. This had affected [step 9f](#9f-prepare-genelab-dge-tables-with-annotations-on-datasets-with-ercc-spike-in) and [step 9i](#9i-prepare-genelab-dge-tables-with-annotations-on-datasets-without-ercc-spike-in).

---

# Table of contents  

- [**Software used**](#software-used)
- [**General processing overview with example commands**](#general-processing-overview-with-example-commands)
  - [**1. Raw Data QC**](#1-raw-data-qc)
    - [1a. Raw Data QC](#1a-raw-data-qc)
    - [1b. Compile Raw Data QC](#1b-compile-raw-data-qc)
  - [**2. Trim/Filter Raw Data and Trimmed Data QC**](#2-trimfilter-raw-data-and-trimmed-data-qc)
    - [2a. Trim/Filter Raw Data](#2a-trimfilter-raw-data)
    - [2b. Trimmed Data QC](#2b-trimmed-data-qc)
    - [2c. Compile Trimmed Data QC](#2c-compile-trimmed-data-qc)
  - [**3. Build STAR Reference**](#3-build-star-reference)
  - [**4. Align Reads to Reference Genome then Sort and Index**](#4-align-reads-to-reference-genome-then-sort-and-index)
    - [4a. Align Reads to Reference Genome with STAR](#4a-align-reads-to-reference-genome-with-star)
    - [4b. Compile Alignment Logs](#4b-compile-alignment-logs)
    - [4c. Tablulate STAR Counts in R](#4c-tablulate-star-counts-in-r)
    - [4d. Sort Aligned Reads](#4d-sort-aligned-reads)
    - [4e. Index Sorted Aligned Reads](#4e-index-sorted-aligned-reads)
  - [**5. Create Reference BED File**](#5-create-reference-bed-file)
    - [5a. Convert GTF to genePred File](#5a-convert-gtf-to-genepred-file)
    - [5b. Convert genePred to BED File](#5b-convert-genepred-to-bed-file)
  - [**6. Assess Strandedness, GeneBody Coverage, Inner Distance, and Read Distribution with RSeQC**](#6-assess-strandedness-genebody-coverage-inner-distance-and-read-distribution-with-rseqc)
    - [6a. Determine Read Strandedness](#6a-determine-read-strandedness)
    - [6b. Compile Strandedness Reports](#6b-compile-strandedness-reports)
    - [6c. Evaluate GeneBody Coverage](#6c-evaluate-genebody-coverage)
    - [6d. Compile GeneBody Coverage Reports](#6d-compile-genebody-coverage-reports)
    - [6e. Determine Inner Distance (For Paired End Datasets)](#6e-determine-inner-distance-for-paired-end-datasets-only)
    - [6f. Compile Inner Distance Reports](#6f-compile-inner-distance-reports)
    - [6g. Assess Read Distribution](#6g-assess-read-distribution)
    - [6h. Compile Read Distribution Reports](#6h-compile-read-distribution-reports)
  - [**7. Build RSEM Reference**](#7-build-rsem-reference)
  - [**8. Quantitate Aligned Reads**](#8-quantitate-aligned-reads)
    - [8a. Count Aligned Reads with RSEM](#8a-count-aligned-reads-with-rsem)
    - [8b. Compile RSEM Count Logs](#8b-compile-rsem-count-logs)
    - [8c. Calculate Total Number of Genes Expressed Per Sample in R](#8c-calculate-total-number-of-genes-expressed-per-sample-in-r)
  - [**9. Normalize Read Counts, Perform Differential Gene Expression Analysis, and Add Gene Annotations in R**](#9-normalize-read-counts-perform-differential-gene-expression-analysis-and-add-gene-annotations-in-r)
    - [9a. Create Sample RunSheet](#9a-create-sample-runsheet)
    - [9b. Environment Set Up](#9b-environment-set-up)
    - [9c. Configure Metadata, Sample Grouping, and Group Comparisons](#9c-configure-metadata-sample-grouping-and-group-comparisons)
    - [9d. Import RSEM GeneCounts](#9d-import-rsem-genecounts)
    - [9e. Perform DGE on Datasets With ERCC Spike-In](#9e-perform-dge-on-datasets-with-ercc-spike-in)
    - [9f. Prepare GeneLab DGE Tables with Annotations on Datasets With ERCC Spike-In](#9f-prepare-genelab-dge-tables-with-annotations-on-datasets-with-ercc-spike-in)
    - [9g. Export GeneLab DGE Tables with Annotations for Datasets With ERCC Spike-In](#9g-export-genelab-dge-tables-with-annotations-for-datasets-with-ercc-spike-in)
    - [9h. Perform DGE on Datasets Without ERCC Spike-In](#9h-perform-dge-on-datasets-without-ercc-spike-in)
    - [9i. Prepare GeneLab DGE Tables with Annotations on Datasets Without ERCC Spike-In](#9i-prepare-genelab-dge-tables-with-annotations-on-datasets-without-ercc-spike-in)
    - [9j. Export GeneLab DGE Tables with Annotations for Datasets Without ERCC Spike-In](#9j-export-genelab-dge-tables-with-annotations-for-datasets-without-ercc-spike-in)

  - [**10. Evaluate ERCC Spike-In Data**](#10-evaluate-ercc-spike-in-data)
    - [10a. Evaluate ERCC Count Data in Python](#10a-evaluate-ercc-count-data-in-python)
    - [10b. Perform DESeq2 Analysis of ERCC Counts in R](#10b-perform-deseq2-analysis-of-ercc-counts-in-r)
    - [10c. Analyze ERCC DESeq2 Results in Python](#10c-analyze-ercc-deseq2-results-in-python)

---

# Software used  

|Program|Version|Relevant Links|
|:------|:------:|:-------------|
|FastQC|0.11.9|[https://www.bioinformatics.babraham.ac.uk/projects/fastqc/](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)|
|MultiQC|1.12|[https://multiqc.info/](https://multiqc.info/)|
|Cutadapt|3.7|[https://cutadapt.readthedocs.io/en/stable/](https://cutadapt.readthedocs.io/en/stable/)|
|TrimGalore!|0.6.7|[https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)|
|STAR|2.7.10a|[https://github.com/alexdobin/STAR](https://github.com/alexdobin/STAR)|
|RSEM|1.3.1|[https://github.com/deweylab/RSEM](https://github.com/deweylab/RSEM)|
|Samtools|1.15|[http://www.htslib.org/](http://www.htslib.org/)|
|gtfToGenePred|377|[http://hgdownload.cse.ucsc.edu/admin/exe/](http://hgdownload.cse.ucsc.edu/admin/exe/)|
|genePredToBed|377|[http://hgdownload.cse.ucsc.edu/admin/exe/](http://hgdownload.cse.ucsc.edu/admin/exe/)|
|infer_experiment|4.0.0|[http://rseqc.sourceforge.net/#infer-experiment-py](http://rseqc.sourceforge.net/#infer-experiment-py)|
|geneBody_coverage|4.0.0|[http://rseqc.sourceforge.net/#genebody-coverage-py](http://rseqc.sourceforge.net/#genebody-coverage-py)|
|inner_distance|4.0.0|[http://rseqc.sourceforge.net/#inner-distance-py](http://rseqc.sourceforge.net/#inner-distance-py)|
|read_distribution|4.0.0|[http://rseqc.sourceforge.net/#read-distribution-py](http://rseqc.sourceforge.net/#read-distribution-py)|
|R|4.1.3|[https://www.r-project.org/](https://www.r-project.org/)|
|Bioconductor|3.14.0|[https://bioconductor.org](https://bioconductor.org)|
|DESeq2|1.34|[https://bioconductor.org/packages/release/bioc/html/DESeq2.html](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)|
|tximport|1.27.1|[https://github.com/mikelove/tximport](https://github.com/mikelove/tximport)|
|tidyverse|1.3.1|[https://www.tidyverse.org](https://www.tidyverse.org)|
|stringr|1.4.1|[https://github.com/tidyverse/stringr](https://github.com/tidyverse/stringr)|
|dp_tools|1.1.8|[https://github.com/J-81/dp_tools](https://github.com/J-81/dp_tools)|
|pandas|1.5.0|[https://github.com/pandas-dev/pandas](https://github.com/pandas-dev/pandas)|
|seaborn|0.12.0|[https://seaborn.pydata.org/](https://seaborn.pydata.org/)|
|matplotlib|3.6.0|[https://matplotlib.org/stable](https://matplotlib.org/stable)|
|jupyter notebook|6.4.12|[https://jupyter-notebook.readthedocs.io/](https://jupyter-notebook.readthedocs.io/)|
|numpy|1.23.3|[https://numpy.org/](https://numpy.org/)|
|scipy|1.9.1|[https://scipy.org/](https://scipy.org/)|
|singularity|3.9|[https://sylabs.io/](https://sylabs.io/)|

---

# General processing overview with example commands  

> Exact processing commands for specific datasets are provided in the [GLDS_Processing_Scripts](../GLDS_Processing_Scripts) directory.
> 
> All output files marked with a \# are published for each RNAseq processed dataset in the [GLDS repository](https://genelab-data.ndc.nasa.gov/genelab/projects). 

---

## 1. Raw Data QC

<br>

### 1a. Raw Data QC  

```bash
fastqc -o /path/to/raw_fastqc/output/directory *.fastq.gz
```

**Parameter Definitions:**

- `-o` – the output directory to store results
- `*.fastq.gz` – the input reads are specified as a positional argument, and can be given all at once with wildcards like this, or as individual arguments with spaces inbetween them

**Input Data:**

- *fastq.gz (raw reads)

**Output Data:**

- *fastqc.html (FastQC report)
- *fastqc.zip (FastQC data)

<br>

### 1b. Compile Raw Data QC  

```bash
multiqc --interactive -n raw_multiqc -o /path/to/raw_multiqc/output/directory /path/to/directory/containing/raw_fastqc/files
```

**Parameter Definitions:**

- `--interactive` - force reports to use interactive plots
- `-n` - prefix name for output files
- `-o` – the output directory to store results
- `/path/to/directory/containing/raw_fastqc/files` – the directory holding the output data from the fastqc run, provided as a positional argument

**Input Data:**

- *fastqc.zip (FastQC data, output from [Step 1a](#1a-raw-data-qc))

**Output Data:**

- raw_multiqc.html\# (multiqc report)
- /raw_multiqc_data\# (directory containing multiqc data)

<br>

---

## 2. Trim/Filter Raw Data and Trimmed Data QC

<br>

### 2a. Trim/Filter Raw Data  

```bash
trim_galore --gzip \
  --path_to_cutadapt /path/to/cutadapt \
  --cores NumberOfThreads \
  --phred33 \
  --illumina \ # if adapters are not illumina, replace with adapters used
  --output_dir /path/to/TrimGalore/output/directory \
  --paired \ # only for PE studies, remove this parameter if raw data are SE
  sample1_R1_raw.fastq.gz sample1_R2_raw.fastq.gz sample2_R1_raw.fastq.gz sample2_R2_raw.fastq.gz
# if SE, replace the last line with only the forward reads (R1) of each sample

```

**Parameter Definitions:**

- `--gzip` – compress the output files with `gzip`
- `--path_to_cutadapt` - specify path to cutadapt software if it is not in your `$PATH`
- `--cores` - specify the number of threads available on the server node to perform trimming
- `--phred33` - instructs cutadapt to use ASCII+33 quality scores as Phred scores for quality trimming
- `--illumina` - defines the adapter sequence to be trimmed as the first 13bp of the Illumina universal adapter `AGATCGGAAGAGC`
- `--output_dir` - the output directory to store results
- `--paired` - indicates paired-end reads - both reads, forward (R1) and reverse (R2) must pass length threshold or else both reads are removed
- `sample1_R1_raw.fastq.gz sample1_R2_raw.fastq.gz sample2_R1_raw.fastq.gz sample2_R2_raw.fastq.gz` – the input reads are specified as a positional argument, paired-end read files are listed pairwise such that the forward reads (*R1_raw.fastq.gz) are immediately followed by the respective reverse reads (*R2_raw.fastq.gz) for each sample

**Input Data:**

- *fastq.gz (raw reads)

**Output Data:**

- *fastq.gz\# (trimmed reads)
- *trimming_report.txt\# (trimming report)

<br>

### 2b. Trimmed Data QC  

```bash
fastqc -o /path/to/trimmed_fastqc/output/directory *.fastq.gz
```

**Parameter Definitions:**

- `-o` – the output directory to store results
- `*.fastq.gz` – the input reads are specified as a positional argument, and can be given all at once with wildcards like this, or as individual arguments with spaces inbetween them

**Input Data:**

- *fastq.gz (trimmed reads, output from [Step 2a](#2a-trimfilter-raw-data))

**Output Data:**

- *fastqc.html (FastQC report)
- *fastqc.zip (FastQC data)

<br>

### 2c. Compile Trimmed Data QC  

```bash
multiqc --interactive -n trimmed_multiqc -o /path/to/trimmed_multiqc/output/directory /path/to/directory/containing/trimmed_fastqc/files
```

**Parameter Definitions:**

- `--interactive` - force reports to use interactive plots
- `-n` - prefix name for output files
- `-o` – the output directory to store results
- `/path/to/directory/containing/trimmed_fastqc/files` – the directory holding the output data from the fastqc run, provided as a positional argument

**Input Data:**

- *fastqc.zip (FastQC data, output from [Step 2b](#2b-trimmed-data-qc))

**Output Data:**

- trimmed_multiqc.html\# (multiqc report)
- /trimmed_multiqc_data\# (directory containing multiqc data)

<br>

---

## 3. Build STAR Reference  

```bash
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

- `--runThreadN` – number of threads available on server node to create STAR reference
- `--runMode` - instructs STAR to run genome indices generation job
- `--limitGenomeGenerateRAM` - maximum RAM available (in bytes) to generate STAR reference, at least 35GB are needed for mouse and the example above shows 55GB
- `--genomeSAindexNbases` - length (in bases) of the SA pre-indexing string, usually between 10 and 15. Longer strings require more memory but allow for faster searches. This value should be scaled down for smaller genomes (like bacteria) to min(14, log2(GenomeLength)/2 - 1). For example, for a 1 megaBase genome this value would be 9.
- `--genomeDir` - specifies the path to the directory where the STAR reference will be stored. At least 100GB of available disk space is required for mammalian genomes.
- `--genomeFastaFiles` - specifies one or more fasta file(s) containing the genome reference sequences
- `--sjdbGTFfile` – specifies the file(s) containing annotated transcripts in the standard gtf format
- `--sjdbOverhang` - indicates the length of the genomic sequence around the annotated junction to be used in constructing the splice junctions database. The length should be one less than the maximum length of the reads.

**Input Data:**

- *.fasta (genome sequence, this scRCP version uses the Ensembl fasta file indicated in the `fasta` column of the [GL-DPPD-7110_annotations.csv](../../GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110/GL-DPPD-7110_annotations.csv) GeneLab Annotations file)
- *.gtf (genome annotation, this scRCP version uses the Ensembl gtf file indicated in the `gtf` column of the [GL-DPPD-7110_annotations.csv](../../GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110/GL-DPPD-7110_annotations.csv) GeneLab Annotations file)

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

## 4. Align Reads to Reference Genome then Sort and Index

<br>

### 4a. Align Reads to Reference Genome with STAR

```bash
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
 --quantMode TranscriptomeSAM GeneCounts \
 --outSAMheaderHD @HD VN:1.4 SO:coordinate \
 --outFileNamePrefix /path/to/STAR/output/directory/<sample_id> \
 --readFilesIn /path/to/trimmed_forward_reads \
 /path/to/trimmed_reverse_reads # only needed for PE studies

```

**Parameter Definitions:**

- `--twopassMode` – specifies 2-pass mapping mode; the `Basic` option instructs STAR to perform the 1st pass mapping, then automatically extract junctions, insert them into the genome index, and re-map all reads in the 2nd mapping pass
- `--limitBAMsortRAM` - maximum RAM available (in bytes) to sort the bam files, the example above indicates 65GB
- `--genomeDir` - specifies the path to the directory where the STAR reference is stored
- `--outSAMunmapped` - specifies output of unmapped reads in the sam format; the `Within` option instructs STAR to output the unmapped reads within the main sam file
- `--outFilterType` - specifies the type of filtering; the `BySJout` option instructs STAR to keep only those reads that contain junctions that passed filtering in the SJ.out.tab output file
- `--outSAMattributes` - list of desired sam attributes in the order desired for the output sam file; sam attribute descriptions can be found [here](https://samtools.github.io/hts-specs/SAMtags.pdf)
- `--outFilterMultimapNmax` – specifies the maximum number of loci the read is allowed to map to; all alignments will be output only if the read maps to no more loci than this value
- `--outFilterMismatchNmax` - maximum number of mismatches allowed to be included in the alignment output
- `--outFilterMismatchNoverReadLmax` - ratio of mismatches to read length allowed to be included in the alignment output; the `0.04` value indicates that up to 4 mismatches are allowed per 100 bases
- `--alignIntronMin` - minimum intron size; a genomic gap is considered an intron if its length is equal to or greater than this value, otherwise it is considered a deletion
- `--alignIntronMax` - maximum intron size
- `--alignMatesGapMax` - maximum genomic distance (in bases) between two mates of paired-end reads; this option should be removed for single-end reads
- `--alignSJoverhangMin` - minimum overhang (i.e. block size) for unannotated spliced alignments
- `--alignSJDBoverhangMin` - minimum overhang (i.e. block size) for annotated spliced alignments
- `--sjdbScore` - additional alignment score for alignments that cross database junctions
- `--readFilesCommand` - specifies command needed to interpret input files; the `zcat` option indicates input files are compressed with gzip and zcat will be used to uncompress the gzipped input files
- `--runThreadN` - indicates the number of threads to be used for STAR alignment and should be set to the number of available cores on the server node
- `--outSAMtype` - specifies desired output format; the `BAM SortedByCoordinate` options specify that the output file will be sorted by coordinate and be in the bam format
- `--quantMode` - specifies the type(s) of quantification desired; the `TranscriptomeSAM` option instructs STAR to output a separate sam/bam file containing alignments to the transcriptome and the `GeneCounts` option instructs STAR to output a tab delimited file containing the number of reads per gene
- `--outSAMheaderHD` - indicates a header line for the sam/bam file
- `--outFileNamePrefix` - specifies the path to and prefix for the output file names; for GeneLab the prefix is the sample id
- `--readFilesIn` - path to input read 1 (forward read) and read 2 (reverse read); for paired-end reads, read 1 and read 2 should be separated by a space; for single-end reads only read 1 should be indicated

**Input Data:**

- STAR genome reference (output from [Step 3](#3-build-star-reference))
- *fastq.gz (trimmed reads, output from [Step 2a](#2a-trimfilter-raw-data))

**Output Data:**

- *Aligned.sortedByCoord.out.bam (sorted mapping to genome)
- *Aligned.toTranscriptome.out.bam\# (sorted mapping to transcriptome)
- *Log.final.out\# (log file containing alignment info/stats such as reads mapped, etc)
- *ReadsPerGene.out.tab (tab delimitated file containing STAR read counts per gene with 4 columns that correspond to different strandedness options: column 1 = gene ID, column 2 = counts for unstranded RNAseq, column 3 = counts for 1st read strand aligned with RNA, column 4 = counts for 2nd read strand aligned with RNA)
- *Log.out (main log file containing detailed info about the STAR run)
- *Log.progress.out (minute-by-minute report containing job progress statistics, such as the number of processed reads, % of mapped reads etc.)
- *SJ.out.tab\# (high confidence collapsed splice junctions in tab-delimited format)
- *_STARgenome (directory containing the following:)
  - sjdbInfo.txt
  - sjdbList.out.tab
- *_STARpass1 (directory containing the following:)
  - Log.final.out
  - SJ.out.tab
- *_STARtmp (directory containing the following:)
  - BAMsort (directory containing subdirectories that are empty – this was the location for temp files that were automatically removed after successful completion)

<br>

### 4b. Compile Alignment Logs

```bash
multiqc --interactive -n align_multiqc -o /path/to/aligned_multiqc/output/directory /path/to/*Log.final.out/files
```

**Parameter Definitions:**

- `--interactive` - force reports to use interactive plots
- `-n` - prefix name for output files
- `-o` – the output directory to store results
- `/path/to/*Log.final.out/files` – the directory holding the *Log.final.out output files from the [STAR alignment step](#4a-align-reads-to-reference-genome-with-star), provided as a positional argument

**Input Data:**

- *Log.final.out (log file containing alignment info/stats such as reads mapped, etc., output from [Step 4a](#4a-align-reads-to-reference-genome-with-star))

**Output Data:**

- align_multiqc.html\# (multiqc report)
- /align_multiqc_data\# (directory containing multiqc data)

<br>

### 4c. Tablulate STAR Counts in R

```R
print("Make STAR counts table")
print("")

work_dir="/path/to/working/directory/where/script/is/executed/from" ## Must contain samples.txt file
align_dir="/path/to/directory/containing/STAR/counts/files"

setwd(file.path(work_dir))

### Pull in sample names where the "samples.txt" file is a single column list of sample names ###
study <- read.csv(Sys.glob(file.path(work_dir,"samples.txt")), header = FALSE, row.names = 1, stringsAsFactors = TRUE)

##### Import Data
ff <- list.files(file.path(align_dir), pattern = "ReadsPerGene.out.tab", recursive=TRUE, full.names = TRUE)

## Reorder the *genes.results files to match the ordering of the ISA samples
ff <- ff[sapply(rownames(study), function(x)grep(paste0(align_dir, '/', x,'_ReadsPerGene.out.tab$'), ff, value=FALSE))]

# Remove the first 4 lines
counts.files <- lapply( ff, read.table, skip = 4 )

# Get counts aligned to either strand for unstranded data by selecting col 2, to the first (forward) strand by selecting col 3 or to the second (reverse) strand by selecting col 4
counts <- as.data.frame( sapply( counts.files, function(x) x[ , 3 ] ) )

# Add column and row names
colnames(counts) <- rownames(study)
row.names(counts) <- counts.files[[1]]$V1


##### Export unnormalized counts table
setwd(file.path(align_dir))
write.csv(counts,file='STAR_Unnormalized_Counts.csv')


## print session info ##
print("Session Info below: ")
print("")
sessionInfo()
```

**Input Data:**

- samples.txt (A newline delimited list of sample IDs)
- *ReadsPerGene.out.tab (STAR counts per gene, output from [Step 4a](#4a-align-reads-to-reference-genome-with-star))

**Output Data:**

- STAR_Unnormalized_Counts.csv\# (Table containing raw STAR counts for each sample)

<br>

### 4d. Sort Aligned Reads

```bash
samtools sort -m 3G \
	--threads NumberOfThreads \
	-o /path/to/*Aligned.sortedByCoord_sorted.out.bam \
  /path/to/*Aligned.sortedByCoord.out.bam
```

**Parameter Definitions:**

- `-m` - memory available per thread, `3G` indicates 3 gigabytes, this can be changed based on user resources
- `--threads` - number of threads available on server node to sort genome alignment files
- `/path/to/*Aligned.sortedByCoord.out.bam` – path to the *Aligned.sortedByCoord.out.bam output files from the [STAR alignment step](#4a-align-reads-to-reference-genome-with-star), provided as a positional argument

**Input Data:**

- *Aligned.sortedByCoord.out.bam (sorted mapping to genome file, output from [Step 4a](#4a-align-reads-to-reference-genome-with-star))

**Output Data:**

- *Aligned.sortedByCoord_sorted.out.bam\# (samtools sorted genome aligned bam file)

<br>

### 4e. Index Sorted Aligned Reads

```bash
samtools index -@ NumberOfThreads /path/to/*Aligned.sortedByCoord_sorted.out.bam
```

**Parameter Definitions:**

- `-@` - number of threads available on server node to index the sorted alignment files
- `/path/to/*Aligned.sortedByCoord_sorted.out.bam` – the path to the sorted *Aligned.sortedByCoord_sorted.out.bam output files from the [step 4d](#4d-sort-aligned-reads), provided as a positional argument

**Input Data:**

- *Aligned.sortedByCoord_sorted.out.bam (sorted mapping to genome file, output from [Step 4d](#4d-sort-aligned-reads))

**Output Data:**

- *Aligned.sortedByCoord_sorted.out.bam.bai\# (index of sorted mapping to genome file)

<br>

---

## 5. Create Reference BED File

<br>

### 5a. Convert GTF to genePred File  

```bash
gtfToGenePred /path/to/annotation/gtf/file \
  /path/to/output/genePred/file

```

**Parameter Definitions:**

- `/path/to/annotation/gtf/file` – specifies the file(s) containing annotated reference transcripts in the standard gtf format, provided as a positional argument
- `/path/to/output/genePred/file` – specifies the location and name of the output genePred file(s), provided as a positional argument

**Input Data:**

- *.gtf (genome annotation, this scRCP version uses the Ensembl gtf file indicated in the `gtf` column of the [GL-DPPD-7110_annotations.csv](../../GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110/GL-DPPD-7110_annotations.csv) GeneLab Annotations file)

**Output Data:**

- *.genePred (genome annotation in genePred format)

<br>

### 5b. Convert genePred to BED File  

```bash
genePredToBed /path/to/annotation/genePred/file \
  /path/to/output/BED/file
```

**Parameter Definitions:**

- `/path/to/annotation/genePred/file` – specifies the file(s) containing annotated reference transcripts in the genePred format, provided as a positional argument
- `/path/to/output/BED/file` – specifies the location and name of the output BED file(s), provided as a positional argument

**Input Data:**

- *.genePred (genome annotation in genePred format, output from [Step 5a](#5a-convert-gtf-to-genepred-file))

**Output Data:**

- *.bed (genome annotation in BED format)

<br>

---

## 6. Assess Strandedness, GeneBody Coverage, Inner Distance, and Read Distribution with RSeQC

<br>

### 6a. Determine Read Strandedness

```bash
infer_experiment.py -r /path/to/annotation/BED/file \
 -i /path/to/*Aligned.sortedByCoord_sorted.out.bam \
 -s 15000000 > /path/to/*infer_expt.out
```

**Parameter Definitions:**

- `-r` – specifies the path to the reference annotation BED file
- `-i` - specifies the path to the input bam file(s)
- `-s` - specifies the number of reads to be sampled from the input bam file(s), 15M reads are sampled
- `>` - redirects standard output to specified file
- `/path/to/*infer_expt.out` - specifies the location and name of the file containing the infer_experiment standard output

**Input Data:**

- *.bed (genome annotation in BED format, output from [Step 5b](#5b-convert-genepred-to-bed-file))
- *Aligned.sortedByCoord_sorted.out.bam (sorted mapping to genome file, output from [Step 4d](#4d-sort-aligned-reads))
- *Aligned.sortedByCoord_sorted.out.bam.bai (index of sorted mapping to genome file, output from [Step 4e](#4e-index-sorted-aligned-reads), although not indicated in the command, this file must be present in the same directory as the respective \*Aligned.sortedByCoord_sorted.out.bam file)

**Output Data:**

- *infer_expt.out (file containing the infer_experiment standard output)

<br>

### 6b. Compile Strandedness Reports

```bash
multiqc --interactive -n infer_exp_multiqc -o /path/to/infer_exp_multiqc/output/directory /path/to/*infer_expt.out/files
```

**Parameter Definitions:**

- `--interactive` - force reports to use interactive plots
- `-n` - prefix name for output files
- `-o` – the output directory to store results
- `/path/to/*infer_expt.out/files` – the directory holding the *infer_expt.out output files from the [read strandedness step](#6a-determine-read-strandedness), provided as a positional argument

**Input Data:**

- *infer_expt.out (file containing the infer_experiment standard output, output from [Step 6a](#6a-determine-read-strandedness))

**Output Data:**

- infer_exp_multiqc.html\# (multiqc report)
- /infer_exp_multiqc_data\# (directory containing multiqc data)

<br>

### 6c. Evaluate GeneBody Coverage

```bash
geneBody_coverage.py -r /path/to/annotation/BED/file \
 -i /path/to/*Aligned.sortedByCoord_sorted.out.bam \
 -o /path/to/geneBody_coverage/output/directory/<sample_id>
```

**Parameter Definitions:**

- `-r` – specifies the path to the reference annotation BED file
- `-i` - specifies the path to the input bam file(s)
- `-o` - specifies the path to the output directory
- `/path/to/geneBody_coverage/output/directory/<sample_id>` - specifies the location and name of the directory containing the geneBody_coverage output files

**Input Data:**

- *.bed (genome annotation in BED format, output from [Step 5b](#5b-convert-genepred-to-bed-file))
- *Aligned.sortedByCoord_sorted.out.bam (sorted mapping to genome file, output from [Step 4d](#4d-sort-aligned-reads))
- *Aligned.sortedByCoord_sorted.out.bam.bai (index of sorted mapping to genome file, output from [Step 4e](#4e-index-sorted-aligned-reads), although not indicated in the command, this file must be present in the same directory as the respective \*Aligned.sortedByCoord_sorted.out.bam file)

**Output Data:**

- *.geneBodyCoverage.curves.pdf (genebody coverage line plot)
- *.geneBodyCoverage.r (R script that generates the genebody coverage line plot)
- *.geneBodyCoverage.txt (tab delimited file containing genebody coverage values used to generate the line plot)

<br>

### 6d. Compile GeneBody Coverage Reports

```bash
multiqc --interactive -n genebody_cov_multiqc -o /path/to/geneBody_coverage_multiqc/output/directory /path/to/geneBody_coverage/output/files
```

**Parameter Definitions:**

- `--interactive` - force reports to use interactive plots
- `-n` - prefix name for output files
- `-o` – the output directory to store results
- `/path/to/geneBody_coverage/output/files` – the directory holding the geneBody_coverage output files from [step 6c](#6c-evaluate-genebody-coverage), provided as a positional argument

**Input Data:**

- *.geneBodyCoverage.txt (tab delimited file containing genebody coverage values, output from [Step 6c](#6c-evaluate-genebody-coverage))

**Output Data:**

- geneBody_cov_multiqc.html\# (multiqc report)
- /geneBody_cov_multiqc_data\# (directory containing multiqc data)

<br>

### 6e. Determine Inner Distance (For Paired End Datasets ONLY)

```bash
inner_distance.py -r /path/to/annotation/BED/file \
 -i /path/to/*Aligned.sortedByCoord_sorted.out.bam \
 -k 15000000 \
 -l -150 \
 -u 350 \
 -o  /path/to/inner_distance/output/directory
```

**Parameter Definitions:**

- `-r` – specifies the path to the reference annotation BED file
- `-i` - specifies the path to the input bam file(s)
- `-k` - specifies the number of reads to be sampled from the input bam file(s), 15M reads are sampled
- `-l` - specifies the lower bound of inner distance (bp).
- `-u` - specifies the upper bound of inner distance (bp)
- `/path/to/inner_distance/output/directory` - specifies the location and name of the directory containing the inner_distance output files

**Input Data:**

- *.bed (genome annotation in BED format, output from [Step 5b](#5b-convert-genepred-to-bed-file))
- *Aligned.sortedByCoord_sorted.out.bam (sorted mapping to genome file, output from [Step 4d](#4d-sort-aligned-reads))
- *Aligned.sortedByCoord_sorted.out.bam.bai (index of sorted mapping to genome file, output from [Step 4e](#4e-index-sorted-aligned-reads), although not indicated in the command, this file must be present in the same directory as the respective \*Aligned.sortedByCoord_sorted.out.bam file)

**Output Data:**

- *.inner_distance.txt (log of read-wise inner distance results)
- *.inner_distance_freq.txt (tab delimited table of inner distances mapped to number of reads with that distance)
- *.inner_distance_plot.pdf (histogram plot of inner distance distribution)
- *.inner_distance_plot.r (R script that generates the histogram plot)

<br>

### 6f. Compile Inner Distance Reports

```bash
multiqc --interactive -n inner_dist_multiqc /path/to/inner_dist_multiqc/output/directory /path/to/inner_dist/output/files
```

**Parameter Definitions:**

- `--interactive` - force reports to use interactive plots
- `-n` - prefix name for output files
- `-o` – the output directory to store results
- `/path/to/inner_dist/output/files` – the directory holding the inner_distance output files from [Step 6e](#6e-determine-inner-distance-for-paired-end-datasets-only), provided as a positional argument

**Input Data:**

- *.inner_distance_freq.txt (tab delimited table of inner distances from [step 6e](#6e-determine-inner-distance-for-paired-end-datasets-only))

**Output Data:**

- inner_distance_multiqc.html\# (multiqc report)
- /inner_distance_multiqc_data\# (directory containing multiqc data)

<br>

### 6g. Assess Read Distribution

```bash
read_distribution.py -r /path/to/annotation/BED/file \
 -i /path/to/*Aligned.sortedByCoord_sorted.out.bam > /path/to/*read_dist.out
```

**Parameter Definitions:**

- `-r` – specifies the path to the reference annotation BED file
- `-i` - specifies the path to the input bam file(s)
- `>` - redirects standard output to specified file
- `/path/to/*read_dist.out` - specifies the location and name of the file containing the read_distribution standard output

**Input Data:**

- *.bed (genome annotation in BED format, output from [Step 5b](#5b-convert-genepred-to-bed-file))
- *Aligned.sortedByCoord_sorted.out.bam (sorted mapping to genome file, output from [Step 4d](#4d-sort-aligned-reads))
- *Aligned.sortedByCoord_sorted.out.bam.bai (index of sorted mapping to genome file, output from [Step 4e](#4e-index-sorted-aligned-reads), although not indicated in the command, this file must be present in the same directory as the respective \*Aligned.sortedByCoord_sorted.out.bam file)

**Output Data:**

- *read_dist.out (file containing the read distribution standard output)

<br>

### 6h. Compile Read Distribution Reports

```bash
multiqc --interactive -n read_dist_multiqc -o /path/to/read_dist_multiqc/output/directory /path/to/*read_dist.out/files
```

**Parameter Definitions:**

- `--interactive` - force reports to use interactive plots
- `-n` - prefix name for output files
- `-o` – the output directory to store results
- `/path/to/*read_dist.out/files` – the directory holding the *read_dist.out output files from [Step 6g](#6g-assess-read-distribution) provided as a positional argument

**Input Data:**

- *read_dist.out (files containing the read_distribution standard output, output from [Step 6g](#6g-assess-read-distribution))

**Output Data:**

- read_dist_multiqc.html\# (multiqc report)
- /read_dist_multiqc_data\# (directory containing multiqc data)

<br>

---

## 7. Build RSEM Reference

```bash
rsem-prepare-reference --gtf /path/to/annotation/gtf/file \
 /path/to/genome/fasta/file \
 /path/to/RSEM/genome/directory/RSEM_ref_prefix

```

**Parameter Definitions:**

- `--gtf` – specifies the file(s) containing annotated transcripts in the standard gtf format
- `/path/to/genome/fasta/file` – specifies one or more fasta file(s) containing the genome reference sequences, provided as a positional argument
- `/path/to/RSEM/genome/directory/RSEM_ref_prefix` - specifies the path to the directory where the RSEM reference will be stored and the prefix desired for the RSEM reference files, provided as a positional argument

**Input Data:**

- *.fasta (genome sequence, this scRCP version uses the Ensembl fasta file indicated in the `fasta` column of the [GL-DPPD-7110_annotations.csv](../../GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110/GL-DPPD-7110_annotations.csv) GeneLab Annotations file)
- *.gtf (genome annotation, this scRCP version uses the Ensembl gtf file indicated in the `gtf` column of the [GL-DPPD-7110_annotations.csv](../../GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110/GL-DPPD-7110_annotations.csv) GeneLab Annotations file)

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

## 8. Quantitate Aligned Reads

<br>

### 8a. Count Aligned Reads with RSEM

```bash
rsem-calculate-expression --num-threads NumberOfThreads \
 --alignments \
 --bam \
 --paired-end \
 --seed 12345 \
 --seed-length 20 \
 --estimate-rspd \
 --no-bam-output \
 --strandedness reverse|forward|none \
 /path/to/*Aligned.toTranscriptome.out.bam \
 /path/to/RSEM/genome/directory/RSEM_ref_prefix \
 /path/to/RSEM/counts/output/directory/<sample_id>
```

**Parameter Definitions:**

- `--num-threads` – specifies the number of threads to use
- `--alignments` - indicates that the input file contains alignments in sam, bam, or cram format
- `--bam` - specifies that the input alignments are in bam format
- `--paired-end` - indicates that the input reads are paired-end reads; this option should be removed if the input reads are single-end
- `--seed` - the seed for the random number generators used in calculating posterior mean estimates and credibility intervals; must be a non-negative 32-bit integer
- `--seed-length 20` - instructs RSEM to ignore any aligned read if it or its mates' (for paired-end reads) length is less than 20bp
- `--estimate-rspd` - instructs RSEM to estimate the read start position distribution (rspd) from the data
- `--no-bam-output` - instructs RSEM not to output any bam file
- `--strandedness` - defines the strandedness of the RNAseq reads; the `reverse` option is used if read strandedness (output from [step 6](#6a-determine-read-strandedness)) is antisense, `forward` is used with sense strandedness, and `none` is used if strandedness is half sense half antisense
- `/path/to/*Aligned.toTranscriptome.out.bam` - specifies path to input bam files, provided as a positional argument
- `/path/to/RSEM/genome/directory/RSEM_ref_prefix` - specifies the path to the directory where the RSEM reference is stored and its prefix, provided as a positional argument
- `/path/to/RSEM/counts/output/directory` – specifies the path to and prefix for the output file names; for GeneLab the prefix is the sample id

**Input Data:**

- RSEM genome reference (output from [Step 7](#7-build-rsem-reference))
- *Aligned.toTranscriptome.out.bam (sorted mapping to transcriptome, output from [Step 4a](#4a-align-reads-to-reference-genome-with-star))

**Output Data:**

- *genes.results\# (counts per gene)
- *isoforms.results\# (counts per isoform)
- *stat (directory containing the following stats files)
  - *cnt
  - *model
  - *theta

<br>

### 8b. Compile RSEM Count Logs

```bash
multiqc --interactive -n RSEM_count_multiqc -o /path/to/RSEM_count_multiqc/output/directory /path/to/*stat/files
```

**Parameter Definitions:**

- `--interactive` - force reports to use interactive plots
- `-n` - prefix name for output files
- `-o` – the output directory to store results
- `/path/to/*stat/files` – the directories holding the *stat output files from the [RSEM Counts step](#8a-count-aligned-reads-with-rsem), provided as a positional argument

**Input Data:**

- *stat (directory containing the following stats files, output from [Step 8a](#8a-count-aligned-reads-with-rsem))
  - *cnt
  - *model
  - *theta

**Output Data:**

- RSEM_count_multiqc.html\# (multiqc report)
- /RSEM_count_multiqc_data\# (directory containing multiqc data)

<br>

### 8c. Calculate Total Number of Genes Expressed Per Sample in R

```R
library(tximport)
library(tidyverse)

work_dir="/path/to/working/directory/where/script/is/executed/from" ## Must contain samples.txt file
counts_dir="/path/to/directory/containing/RSEM/counts/files"

setwd(file.path(work_dir))

### Pull in sample names where the "samples.txt" file is a single column list of sample names ###
samples <- read.csv(Sys.glob(file.path(work_dir,"samples.txt")), header = FALSE, row.names = 1, stringsAsFactors = TRUE)

##### Import RSEM Gene Count Data
files <- list.files(file.path(counts_dir),pattern = ".genes.results", full.names = TRUE)

### reorder the genes.results files to match the ordering of the samples in the metadata file
files <- files[sapply(rownames(samples), function(x)grep(paste0(counts_dir, '/',  x,'.genes.results$'), files, value=FALSE))]

names(files) <- rownames(samples)
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)

##### Count the number of genes with non-zero counts for each sample 
rawCounts <- txi.rsem$counts
NumNonZeroGenes <- (as.matrix(colSums(rawCounts > 0), row.names = 1))
colnames(NumNonZeroGenes) <- c("Number of genes with non-zero counts")

##### Export the number of genes with non-zero counts for each sample
setwd(file.path(counts_dir))
write.csv(NumNonZeroGenes,file='NumNonZeroGenes.csv')

## print session info ##
print("Session Info below: ")
print("")
sessionInfo()
```

**Input Data:**

- samples.txt (A newline delimited list of sample IDs)
- *genes.results (RSEM counts per gene, output from [Step 8a](#8a-count-aligned-reads-with-rsem))

**Output Data:**

- NumNonZeroGenes.csv (A samplewise table of the number of genes expressed)

<br>

---

## 9. Normalize Read Counts, Perform Differential Gene Expression Analysis, and Add Gene Annotations in R

<br>

### 9a. Create Sample RunSheet

> Note: Rather than running the command below to create the runsheet needed for processing, the runsheet may also be created manually by following the [file specification](../Workflow_Documentation/NF_RCP-F/examples/runsheet/README.md).

```bash
### Download the *ISA.zip file from the GeneLab Repository ###

dpt-get-isa-archive \
 --accession GLDS-###

### Parse the metadata from the *ISA.zip file to create a sample runsheet ###

dpt-isa-to-runsheet --accession GLDS-### \
 --config-type bulkRNASeq \
 --config-version Latest \
 --isa-archive *ISA.zip
```

**Parameter Definitions:**

- `--accession GLDS-###` - GLDS accession ID (replace ### with the GLDS number being processed), used to retrieve the urls for the ISA archive and raw reads hosted on the GeneLab Repository
- `--config-type` - Instructs the script to extract the metadata required for `bulkRNAseq` processing from the ISA archive
- `--config-version` - Specifies the `dp-tools` configuration version to use, a value of `Latest` will specify the most recent version
- `--isa-archive` - Specifies the *ISA.zip file for the respective GLDS dataset, downloaded in the `dpt-get-isa-archive` command


**Input Data:**

- No input data required but the GLDS accession ID needs to be indicated, which is used to download the respective ISA archive 

**Output Data:**

- *ISA.zip (compressed ISA directory containing Investigation, Study, and Assay (ISA) metadata files for the respective GLDS dataset, used to define sample groups - the *ISA.zip file is located in the [GLDS repository](https://genelab-data.ndc.nasa.gov/genelab/projects) under 'Study Files' -> 'metadata')

- {GLDS-Accession-ID}_bulkRNASeq_v{version}_runsheet.csv\# (table containing metadata required for processing, version denotes the dp_tools schema used to specify the metadata to extract from the ISA archive)

<br>

### 9b. Environment Set Up

```R
### Install R packages if not already installed ###


install.packages("tidyverse")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.14")
BiocManager::install("tximport")
BiocManager::install("DESeq2")


### Import libraries (tximport, DESeq2, tidyverse, stringr) ###

library(tximport)
library(DESeq2)
library(tidyverse)
library(stringr)


### Define which organism is used in the study - this should be consistent with the name in the "name" column of the GL-DPPD-7110_annotations.csv file, which matches the abbreviations used in the Panther database for each organism ###

organism <- "organism_that_samples_were_derived_from"


### Define the location of the input data and where the output data will be printed to ###

runsheet_path="/path/to/directory/containing/runsheet.csv/file" ## This is the runsheet created in Step 9a above
work_dir="/path/to/working/directory/where/script/is/executed/from" 
counts_dir="/path/to/directory/containing/RSEM/counts/files"
norm_output="/path/to/normalized/counts/output/directory"
DGE_output="/path/to/DGE/output/directory"
DGE_output_ERCC="/path/to/ERCC-normalized/DGE/output/directory" ## Only needed for datasets with ERCC spike-in


### Pull in the GeneLab annotation table (GL-DPPD-7110_annotations.csv) file ###

org_table_link <- "https://raw.githubusercontent.com/nasa/GeneLab_Data_Processing/master/GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110/GL-DPPD-7110_annotations.csv"

org_table <- read.table(org_table_link, sep = ",", header = TRUE)


### Define the link to the GeneLab annotation table for the organism of interest ###

annotations_link <- org_table[org_table$name == organism, "genelab_annots_link"]


### Set your working directory to the directory where you will execute your DESeq2 script from ###

setwd(file.path(work_dir))

```

<br>

### 9c. Configure Metadata, Sample Grouping, and Group Comparisons

```R
### Pull all factors for each sample in the study from the runsheet created in Step 9a ###

compare_csv_from_runsheet <- function(runsheet_path) {
    df = read.csv(runsheet_path)
    # get only Factor Value columns
    factors = as.data.frame(df[,grep("Factor.Value", colnames(df), ignore.case=TRUE)])
    colnames(factors) = paste("factor",1:dim(factors)[2], sep= "_")
    result = data.frame(sample_id = df[,c("Sample.Name")], factors)	
    return(result)
}


### Load metadata from runsheet csv file ###

compare_csv <- compare_csv_from_runsheet(runsheet_path)


### Create data frame containing all samples and respective factors ###

study <- as.data.frame(compare_csv[,2:dim(compare_csv)[2]])
colnames(study) <- colnames(compare_csv)[2:dim(compare_csv)[2]]
rownames(study) <- compare_csv[,1]


### Format groups and indicate the group that each sample belongs to ###

if (dim(study) >= 2){
	group<-apply(study,1,paste,collapse = " & ") ## concatenate multiple factors into one condition per sample
} else{
	group<-study[,1]
}
group_names <- paste0("(",group,")",sep = "") ## human readable group names
group <- sub("^BLOCKER_", "",  make.names(paste0("BLOCKER_", group))) # group naming compatible with R models, this maintains the default behaviour of make.names with the exception that 'X' is never prepended to group names
names(group) <- group_names
rm(group_names)


### Format contrasts table, defining pairwise comparisons for all groups ###

contrast.names <- combn(levels(factor(names(group))),2) ## generate matrix of pairwise group combinations for comparison
contrasts <- apply(contrast.names, MARGIN=2, function(col) sub("^BLOCKER_", "",  make.names(paste0("BLOCKER_", stringr::str_sub(col, 2, -2))))) # limited make.names call for each group (also removes leading parentheses)
contrast.names <- c(paste(contrast.names[1,],contrast.names[2,],sep = "v"),paste(contrast.names[2,],contrast.names[1,],sep = "v")) ## format combinations for output table files names
contrasts <- cbind(contrasts,contrasts[c(2,1),])
colnames(contrasts) <- contrast.names
rm(contrast.names) 

```

<br>

### 9d. Import RSEM GeneCounts

```R
### Import RSEM raw (gene) count data ###

files <- list.files(file.path(counts_dir),pattern = ".genes.results", full.names = TRUE)


### Reorder the *genes.results files to match the ordering of the ISA samples ###

files <- files[sapply(rownames(study), function(x)grep(paste0(counts_dir, '/', x,".genes.results$"), files, value=FALSE))]

names(files) <- rownames(study)
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)


### Add 1 to genes with lengths of zero - needed to make DESeqDataSet object ###

txi.rsem$length[txi.rsem$length == 0] <- 1

```

<br>

### 9e. Perform DGE on Datasets With ERCC Spike-In
> Note: For datasets that do not contain ERCC spike-in, skip to [Step 9h](#9h-perform-dge-on-datasets-without-ercc-spike-in)

```R
### Create data frame defining which group each sample belongs to ###

sampleTable <- data.frame(condition=factor(group))
rownames(sampleTable) <- colnames(txi.rsem$counts)


### Make DESeqDataSet object ###

dds <- DESeqDataSetFromTximport(txi.rsem, sampleTable, ~condition)
summary(dds)


############################################################
######### Create ERCC unfiltered raw counts table ##########
############################################################
## Note: These data are used internally at GeneLab for QC ##
############################################################


### Make a DESeqDataSet object using only unfiltered ERCC genes ###

ercc_rows_all <- grep("ERCC-",rownames(dds))
ercc_dds_all <- dds[ercc_rows_all,]

### Print ERCC unfiltered raw counts table ###

ERCC_rawCounts_all = as.data.frame(counts(ercc_dds_all))
write.csv(ERCC_rawCounts_all,file='ERCC_rawCounts_unfiltered.csv')


#############################################################################
### Prepare data to be normalized with and without considering ERCC genes ###
#############################################################################

### Filter out genes with counts of less than 10 in all samples ###

keepGenes <- rowSums(counts(dds)) > 10
dds <- dds[keepGenes,]
summary(dds)
dim(dds)


### Make a DESeqDataSet object using only filtered ERCC genes, which will be used to generate ERCC counts table ###

ercc_rows <- grep("ERCC-",rownames(dds))
ercc_dds <- dds[ercc_rows,]


### Print ERCC filtered raw counts table ###
## Note: These data are used internally at GeneLab for QC

ERCC_rawCounts = as.data.frame(counts(ercc_dds))
write.csv(ERCC_rawCounts,file='ERCC_rawCounts_filtered.csv')


### Create a list of rows containing ERCC group B genes to use for ERCC-normalization ###
## Note: ERCC group B genes should be the same concentration in all samples

ercc_rows_gpB <- grep("ERCC-00096|ERCC-00171|ERCC-00009|ERCC-00042|ERCC-00060|ERCC-00035|ERCC-00025|ERCC-00051|ERCC-00053|ERCC-00148|ERCC-00126|ERCC-00034|ERCC-00150|ERCC-00067|ERCC-00031|ERCC-00109|ERCC-00073|ERCC-00158|ERCC-00104|ERCC-00142|ERCC-00138|ERCC-00117|ERCC-00075",rownames(dds))
ercc_dds_gpB <- dds[ercc_rows_gpB,]
summary(ercc_dds_gpB)
dim(ercc_dds_gpB)


### Identify and list samples that do not contain any counts for ERCC genes and specifically group B ERCC genes ###
## Note: All samples should contain ERCC spike-in and thus ERCC counts, if some samples do not contain ERCC counts, those samples should be removed and not used for downstream ERCC-associated analysis

cat("Samples that do not have detectable ERCC spike-ins: ", colnames(ercc_dds[,colSums(counts(ercc_dds))==0]), sep="\n")

cat("Samples that do not have detectable ERCC group B spike-ins: ", colnames(dds[,colSums(counts(ercc_dds_gpB))==0]), sep="\n")


### Create a new study object WITHOUT the samples that don't have detectable ERCC group B spike-ins ###

remove <- colnames(dds[,colSums(counts(ercc_dds_gpB))==0])
study_sub <- subset(study,!rownames(study) %in% remove) # new study object with non-ERCC-gpB samples removed


### Create a new group object WITHOUT the samples that don't have detectable ERCC group B spike-ins ###

if (dim(study_sub) >= 2){
  group_sub<-apply(study_sub,1,paste,collapse = " & ") # concatenate multiple factors into one condition per sample
} else{
  group_sub<-study_sub[,1]
}
group_names <- paste0("(",group_sub,")",sep = "") # human readable group names
group_sub <- sub("^BLOCKER_", "",  make.names(paste0("BLOCKER_", group_sub))) # group naming compatible with R models, this maintains the default behaviour of make.names with the exception that 'X' is never prepended to group names
names(group_sub) <- group_names
rm(group_names)


### Create new contrasts object that only contains the groups in the subset group object ###

contrasts_sub.names <- combn(levels(factor(names(group_sub))),2) # generate matrix of pairwise group combinations for comparison
contrasts_sub <- apply(contrasts_sub.names, MARGIN=2, function(col) sub("^BLOCKER_", "",  make.names(paste0("BLOCKER_", stringr::str_sub(col, 2, -2))))) # limited make.names call for each group (also removes leading parentheses)
contrasts_sub.names <- c(paste(contrasts_sub.names[1,],contrasts_sub.names[2,],sep = "v"),paste(contrasts_sub.names[2,],contrasts_sub.names[1,],sep = "v")) # format combinations for output table files names
contrasts_sub <- cbind(contrasts_sub,contrasts_sub[c(2,1),])
colnames(contrasts_sub) <- contrasts_sub.names
rm(contrasts_sub.names)

### If all samples contain ERCC spike-in (i.e. there are no samples to remove), reassign group_sub, study_sub and contrasts_sub back to the original variable contents ###

if (length(remove) == 0) {
  group_sub <- group
  study_sub <- study
  contrasts_sub <- contrasts
}


########################################################################################################
### Prepare DESeqDataSet objects to be used for DGE analysis with and without considering ERCC genes ###
########################################################################################################

### Generate a DESeqDataSet object using only non-ERCC genes ###
## dds_1 will be used to generate data without considering ERCC genes

dds_1 <- dds[-c(ercc_rows),] ## remove ERCCs from full counts table


### Generate a DESeqDataSet object using only samples that contain ERCC group B genes ###
## dds_2 will be used to generate data with considering ERCC genes

dds_2 <- dds[,colSums(counts(ercc_dds_gpB)) > 0] ## samples that do not contain ERCC group B counts are removed
sampleTable_sub <- data.frame(condition=factor(group_sub)) ## create a new sampleTable only with samples that contain ERCC group B counts
rownames(sampleTable_sub) <- rownames(study_sub)
dds_2$condition <- sampleTable_sub$condition ## reassign the dds_2 condition to the subset condition containing only samples with ERCC group B counts
summary(dds_2)
dim(dds_2)


#######################################################################
### Perform DESeq2 analysis with and without considering ERCC genes ###
#######################################################################

### Run DESeq analysis with ERCC-normalization by replacing size factor object with ERCC size factors for rescaling ###
## Try first to use the default type="median", but if there is an error (usually due to zeros in genes), use type="poscounts"
## From DESeq2 manual: "The "poscounts" estimator deals with a gene with some zeros, by calculating a modified geometric mean by taking the n-th root of the product of the non-zero counts."

dds_2 <- tryCatch(
      expr = { estimateSizeFactors(dds_2, controlGenes=ercc_rows_gpB) },
      error = function(e) { estimateSizeFactors(dds_2, type="poscounts", controlGenes=ercc_rows_gpB)}
)

dds_2 <- dds_2[-c(ercc_rows),] # remove ERCCs from counts table after normalization
dds_2 <- estimateDispersions(dds_2)
dds_2 <- nbinomWaldTest(dds_2)


### Run DESeq analysis without considering ERCC genes ###

dds_1 <- DESeq(dds_1)


### Generate F statistic p-value (similar to ANOVA p-value) using DESeq2 likelihood ratio test (LRT) design ###

## For ERCC-normalized data

dds_2_lrt <- DESeq(dds_2, test = "LRT", reduced = ~ 1)
res_2_lrt <- results(dds_2_lrt)

## For non-ERCC normalized data

dds_1_lrt <- DESeq(dds_1, test = "LRT", reduced = ~ 1)
res_1_lrt <- results(dds_1_lrt)

```

<br>

### 9f. Prepare GeneLab DGE Tables with Annotations on Datasets With ERCC Spike-In
> Note: For datasets that do not contain ERCC spike-in, skip to [Step 9h](#9h-perform-dge-on-datasets-without-ercc-spike-in)

```R
### Create two data frames, one containing (non-ERCC) normalized counts and the other containing ERCC-normalized counts ###

normCounts <- as.data.frame(counts(dds_1, normalized=TRUE))
ERCCnormCounts <- as.data.frame(counts(dds_2, normalized=TRUE))


### Add 1 to all (non-ERCC) normalized counts and to all ERCC-normalized counts to avoid issues with downstream calculations ###
normCounts <- normCounts +1
ERCCnormCounts <- ERCCnormCounts +1


########################################################
### Prepare DGE table without considering ERCC genes ###
########################################################

### Start the DGE output table with the (non-ERCC) normalized counts for all samples ###
## reduced output table 1 will be used to generate human-readable DGE table

reduced_output_table_1 <- normCounts

## output tables 1 will be used to generate computer-readable DGE table, which is used to create GeneLab visualization plots

output_table_1 <- normCounts


### Iterate through Wald Tests to generate pairwise comparisons of all groups ###

for (i in 1:dim(contrasts)[2]){
	res_1 <- results(dds_1, contrast=c("condition",contrasts[1,i],contrasts[2,i]))
	res_1 <- as.data.frame(res_1@listData)[,c(2,4,5,6)]
	colnames(res_1) <- c(paste0("Log2fc_", colnames(contrasts)[i]), paste0("Stat_",colnames(contrasts)[i]), paste0("P.value_",colnames(contrasts)[i]), paste0("Adj.p.value_",colnames(contrasts)[i]))
	output_table_1 <- cbind(output_table_1,res_1)
	reduced_output_table_1 <- cbind(reduced_output_table_1,res_1)
	rm(res_1)
}


### Generate and add all sample mean column to the (non-ERCC) DGE table ###

output_table_1$All.mean <- rowMeans(normCounts, na.rm = TRUE, dims = 1)
reduced_output_table_1$All.mean <- rowMeans(normCounts, na.rm = TRUE, dims = 1)


### Generate and add all sample stdev column to the (non-ERCC) DGE table ###

output_table_1$All.stdev <- rowSds(as.matrix(normCounts), na.rm = TRUE, dims = 1)
reduced_output_table_1$All.stdev <- rowSds(as.matrix(normCounts), na.rm = TRUE, dims = 1)


### Add F statistic p-value (similar to ANOVA p-value) column to the (non-ERCC) DGE table ###

output_table_1$LRT.p.value <- res_1_lrt@listData$padj
reduced_output_table_1$LRT.p.value <- res_1_lrt@listData$padj


### Generate and add group mean and stdev columns to the (non-ERCC) DGE table ###

tcounts <- as.data.frame(t(normCounts))
tcounts$group <- names(group) # Used final table group name formatting (e.g. '( Space Flight & Blue Light )' )

group_means <- as.data.frame(t(aggregate(. ~ group,data = tcounts,mean))) # Compute group name group-wise means
colnames(group_means) <- paste0("Group.Mean_", group_means['group',]) # assign group name as column names

group_stdev <- as.data.frame(t(aggregate(. ~ group,data = tcounts,sd))) # Compute group name group-wise standard deviation
colnames(group_stdev) <- paste0("Group.Stdev_", group_stdev['group',]) # assign group name as column names

group_means <- group_means[-c(1),] # Drop group name row from data rows (now present as column names)
group_stdev <- group_stdev[-c(1),] # Drop group name row from data rows (now present as column names)
output_table_1 <- cbind(output_table_1,group_means, group_stdev) # Column bind the group-wise data
reduced_output_table_1 <- cbind(reduced_output_table_1,group_means, group_stdev) # Column bind the group-wise data

rm(group_stdev,group_means,tcounts)


### Add columns needed to generate GeneLab visualization plots to the (non-ERCC) DGE table ###

## Add column to indicate the sign (positive/negative) of log2fc for each pairwise comparison ##

updown_table <- sign(output_table_1[,grep("Log2fc_",colnames(output_table_1))])
colnames(updown_table) <- gsub("Log2fc","Updown",grep("Log2fc_",colnames(output_table_1),value = TRUE))
output_table_1 <- cbind(output_table_1,updown_table)
rm(updown_table)


## Add column to indicate contrast significance with p <= 0.1 ##

sig.1_table <- output_table_1[,grep("P.value_",colnames(output_table_1))]<=.1
colnames(sig.1_table) <- gsub("P.value","Sig.1",grep("P.value_",colnames(output_table_1),value = TRUE))
output_table_1 <- cbind(output_table_1,sig.1_table)
rm(sig.1_table)


## Add column to indicate contrast significance with p <= 0.05 ##

sig.05_table <- output_table_1[,grep("P.value_",colnames(output_table_1))]<=.05
colnames(sig.05_table) <- gsub("P.value","Sig.05",grep("P.value_",colnames(output_table_1),value = TRUE))
output_table_1 <- cbind(output_table_1,sig.05_table)
rm(sig.05_table)


## Add columns for the volcano plot with p-value and adjusted p-value ##

log_pval_table <- log2(output_table_1[,grep("P.value_",colnames(output_table_1))])
colnames(log_pval_table) <- paste0("Log2_",colnames(log_pval_table))
output_table_1 <- cbind(output_table_1,log_pval_table)
rm(log_pval_table)
log_adj_pval_table <- log2(output_table_1[,grep("Adj.p.value_",colnames(output_table_1))])
colnames(log_adj_pval_table) <- paste0("Log2_",colnames(log_adj_pval_table))
output_table_1 <- cbind(output_table_1,log_adj_pval_table)
rm(log_adj_pval_table)


## Prepare PCA table for GeneLab visualization plots ##

exp_raw <- log2(normCounts)
PCA_raw <- prcomp(t(exp_raw), scale = FALSE)


### Read in GeneLab annotation table for the organism of interest ###

annot <- read.table(annotations_link, sep = "\t", header = TRUE, quote = "", comment.char = "", row.names = 1)


### Combine annotations table and the (non-ERCC) DGE table ###

output_table_1 <- merge(annot, output_table_1, by='row.names', all.y=TRUE)
output_table_1 <- output_table_1 %>% 
  rename(
    ENSEMBL = Row.names ## Change ENSEMBL to TAIR for plant studies ##
  )


reduced_output_table_1 <- merge(annot, reduced_output_table_1, by='row.names', all.y=TRUE)
reduced_output_table_1 <- reduced_output_table_1 %>% 
  rename(
    ENSEMBL = Row.names ## Change ENSEMBL to TAIR for plant studies ##
  )


#####################################################
### Prepare DGE table with considering ERCC genes ###
#####################################################

### Start the DGE output table with the ERCC-normalized counts for all samples ###
## reduced output table 2 will be used to generate human-readable DGE table

reduced_output_table_2 <- ERCCnormCounts

## output table 2 will be used to generate computer-readable DGE table, which is used to create GeneLab visualization plots

output_table_2 <- ERCCnormCounts


### Iterate through Wald Tests to generate pairwise comparisons of all groups ###

for (i in 1:dim(contrasts_sub)[2]){
  res_2 <- results(dds_2, contrast=c("condition",contrasts_sub[1,i],contrasts_sub[2,i]))
  res_2 <- as.data.frame(res_2@listData)[,c(2,4,5,6)]
  colnames(res_2)<-c(paste0("Log2fc_",colnames(contrasts_sub)[i]), paste0("Stat_",colnames(contrasts_sub)[i]), paste0("P.value_",colnames(contrasts_sub)[i]), paste0("Adj.p.value_",colnames(contrasts_sub)[i]))
  output_table_2<-cbind(output_table_2,res_2)
  reduced_output_table_2 <- cbind(reduced_output_table_2,res_2)
  rm(res_2)
}


### Generate and add all sample mean column to the ERCC-normalized DGE table ###

output_table_2$All.mean <- rowMeans(ERCCnormCounts, na.rm = TRUE, dims = 1)
reduced_output_table_2$All.mean <- rowMeans(ERCCnormCounts, na.rm = TRUE, dims = 1)


### Generate and add all sample stdev column to the ERCC-normalized DGE table ###

output_table_2$All.stdev <- rowSds(as.matrix(ERCCnormCounts), na.rm = TRUE, dims = 1)
reduced_output_table_2$All.stdev <- rowSds(as.matrix(ERCCnormCounts), na.rm = TRUE, dims = 1)


### Add F statistic p-value (similar to ANOVA p-value) column to the ERCC-normalized DGE table ###

output_table_2$LRT.p.value <- res_2_lrt@listData$padj
reduced_output_table_2$LRT.p.value <- res_2_lrt@listData$padj


### Generate and add group mean and stdev columns to the ERCC-normalized DGE table ###

tcounts <- as.data.frame(t(ERCCnormCounts))
tcounts$group_sub <- names(group_sub) # Used final table group name formatting (e.g. '( Space Flight & Blue Light )' )

group_means <- as.data.frame(t(aggregate(. ~ group_sub,data = tcounts,mean))) # Compute group name group-wise means
colnames(group_means) <- paste0("Group.Mean_", group_means['group_sub',]) # assign group name as column names

group_stdev <- as.data.frame(t(aggregate(. ~ group_sub,data = tcounts,sd))) # Compute group name group-wise standard deviation
colnames(group_stdev) <- paste0("Group.Stdev_", group_stdev['group_sub',]) # assign group name as column names

group_means <- group_means[-c(1),] # Drop group name row from data rows (now present as column names)
group_stdev <- group_stdev[-c(1),] # Drop group name row from data rows (now present as column names)
output_table_2 <- cbind(output_table_2,group_means, group_stdev) # Column bind the group-wise data
reduced_output_table_2 <- cbind(reduced_output_table_2,group_means, group_stdev) # Column bind the group-wise data

rm(group_stdev,group_means,tcounts)


### Add columns needed to generate GeneLab visualization plots to the ERCC-normalized DGE table ###

## Add column to indicate the sign (positive/negative) of log2fc for each pairwise comparison ##

updown_table <- sign(output_table_2[,grep("Log2fc_",colnames(output_table_2))])
colnames(updown_table) <- gsub("Log2fc","Updown",grep("Log2fc_",colnames(output_table_2),value = TRUE))
output_table_2 <- cbind(output_table_2,updown_table)
rm(updown_table)


## Add column to indicate contrast significance with p <= 0.1 ##

sig.1_table <- output_table_2[,grep("P.value_",colnames(output_table_2))]<=.1
colnames(sig.1_table) <- gsub("P.value","Sig.1",grep("P.value_",colnames(output_table_2),value = TRUE))
output_table_2 <- cbind(output_table_2,sig.1_table)
rm(sig.1_table)


## Add column to indicate contrast significance with p <= 0.05 ##

sig.05_table <- output_table_2[,grep("P.value_",colnames(output_table_2))]<=.05
colnames(sig.05_table) <- gsub("P.value","Sig.05",grep("P.value_",colnames(output_table_2),value = TRUE))
output_table_2 <- cbind(output_table_2,sig.05_table)
rm(sig.05_table)


## Add columns for the volcano plot with p-value and adjusted p-value ##

log_pval_table <- log2(output_table_2[,grep("P.value_",colnames(output_table_2))])
colnames(log_pval_table) <- paste0("Log2_",colnames(log_pval_table))
output_table_2 <- cbind(output_table_2,log_pval_table)
rm(log_pval_table)
log_adj_pval_table <- log2(output_table_2[,grep("Adj.p.value_",colnames(output_table_2))])
colnames(log_adj_pval_table) <- paste0("Log2_",colnames(log_adj_pval_table))
output_table_2 <- cbind(output_table_2,log_adj_pval_table)
rm(log_adj_pval_table)


## Prepare PCA table for GeneLab visualization plots ##

exp_raw_ERCCnorm <- log2(ERCCnormCounts)
PCA_raw_ERCCnorm <- prcomp(t(exp_raw_ERCCnorm), scale = FALSE)


### Combine annotations table and the ERCC-normalized DGE table ###

output_table_2 <- merge(annot, output_table_2, by='row.names', all.y=TRUE)
output_table_2 <- output_table_2 %>% 
  rename(
    ENSEMBL = Row.names ## Change ENSEMBL to TAIR for plant studies ##
  )


reduced_output_table_2 <- merge(annot, reduced_output_table_2, by='row.names', all.y=TRUE)
reduced_output_table_2 <- reduced_output_table_2 %>% 
  rename(
    ENSEMBL = Row.names ## Change ENSEMBL to TAIR for plant studies ##
  )

```

<br>

### 9g. Export GeneLab DGE Tables with Annotations for Datasets With ERCC Spike-In
> Note: For datasets that do not contain ERCC spike-in, skip to [Step 9h](#9h-perform-dge-on-datasets-without-ercc-spike-in)

```R
### Export unnormalized, (non-ERCC) normalized, and ERCC-normalized counts tables ###

normCounts_exp <- as.data.frame(counts(dds_1, normalized=TRUE))
ERCCnormCounts_exp <- as.data.frame(counts(dds_2, normalized=TRUE))

write.csv(txi.rsem$counts,file.path(norm_output, "RSEM_Unnormalized_Counts.csv"))
write.csv(normCounts_exp,file.path(norm_output, "Normalized_Counts.csv"))
write.csv(ERCCnormCounts_exp,file.path(norm_output, "ERCC_Normalized_Counts.csv"))


### Export sample grouping and contrasts tables for (non-ERCC) normalized and ERCC-normalized data ###

write.csv(sampleTable,file.path(DGE_output, "SampleTable.csv"))
write.csv(sampleTable_sub,file.path(DGE_output_ERCC, "ERCCnorm_SampleTable.csv"))

write.csv(contrasts,file.path(DGE_output, "contrasts.csv"))
write.csv(contrasts_sub,file.path(DGE_output_ERCC, "ERCCnorm_contrasts.csv"))


### Export human-readable (non-ERCC) normalized and ERCC-normalized DGE tables ###

write.csv(reduced_output_table_1,file.path(DGE_output, "differential_expression.csv"), row.names = FALSE)
write.csv(reduced_output_table_2,file.path(DGE_output_ERCC, "ERCCnorm_differential_expression.csv"), row.names = FALSE)


### Export computer-readable DGE and PCA tables used for GeneLab visualization ###

write.csv(output_table_1,file.path(DGE_output, "visualization_output_table.csv"), row.names = FALSE)
write.csv(output_table_2,file.path(DGE_output_ERCC, "visualization_output_table_ERCCnorm.csv"), row.names = FALSE)

write.csv(PCA_raw$x,file.path(DGE_output, "visualization_PCA_table.csv"), row.names = TRUE)
write.csv(PCA_raw_ERCCnorm$x,file.path(DGE_output_ERCC, "visualization_PCA_table_ERCCnorm.csv"), row.names = TRUE)


### print session info ###

print("Session Info below: ")
sessionInfo()

```

<br>

### 9h. Perform DGE on Datasets Without ERCC Spike-In

```R
### Create data frame defining which group each sample belongs to ###

sampleTable <- data.frame(condition=factor(group))
rownames(sampleTable) <- colnames(txi.rsem$counts)


### Make DESeqDataSet object ###

dds <- DESeqDataSetFromTximport(txi.rsem, sampleTable, ~condition)
summary(dds)


### Filter out genes with counts of less than 10 in all samples ###

keepGenes <- rowSums(counts(dds)) > 10
dds_1 <- dds[keepGenes,]
summary(dds_1)
dim(dds_1)


### Run DESeq analysis ###

dds_1 <- DESeq(dds_1)


### Generate F statistic p-value (similar to ANOVA p-value) using DESeq2 likelihood ratio test (LRT) design ###

dds_1_lrt <- DESeq(dds_1, test = "LRT", reduced = ~ 1)
res_1_lrt <- results(dds_1_lrt)

```

<br>

### 9i. Prepare GeneLab DGE Tables with Annotations on Datasets Without ERCC Spike-In

```R
### Create a data frame containing normalized counts ###

normCounts <- as.data.frame(counts(dds_1, normalized=TRUE))


### Add 1 to all normalized counts to avoid issues with downstream calculations ###
normCounts <- normCounts +1


### Start the DGE output table with the normalized counts for all samples ###
## reduced output table 1 will be used to generate human-readable DGE table

reduced_output_table_1 <- normCounts

## output tables 1 will be used to generate computer-readable DGE table, which is used to create GeneLab visualization plots

output_table_1 <- normCounts


### Iterate through Wald Tests to generate pairwise comparisons of all groups ###

for (i in 1:dim(contrasts)[2]){
	res_1 <- results(dds_1, contrast=c("condition",contrasts[1,i],contrasts[2,i]))
	res_1 <- as.data.frame(res_1@listData)[,c(2,4,5,6)]
	colnames(res_1) <- c(paste0("Log2fc_", colnames(contrasts)[i]), paste0("Stat_",colnames(contrasts)[i]), paste0("P.value_",colnames(contrasts)[i]), paste0("Adj.p.value_",colnames(contrasts)[i]))
	output_table_1 <- cbind(output_table_1,res_1)
	reduced_output_table_1 <- cbind(reduced_output_table_1,res_1)
	rm(res_1)
}


### Generate and add all sample mean column to the DGE table ###

output_table_1$All.mean <- rowMeans(normCounts, na.rm = TRUE, dims = 1)
reduced_output_table_1$All.mean <- rowMeans(normCounts, na.rm = TRUE, dims = 1)


### Generate and add all sample stdev column to the DGE table ###

output_table_1$All.stdev <- rowSds(as.matrix(normCounts), na.rm = TRUE, dims = 1)
reduced_output_table_1$All.stdev <- rowSds(as.matrix(normCounts), na.rm = TRUE, dims = 1)


### Add F statistic p-value (similar to ANOVA p-value) column to the DGE table ###

output_table_1$LRT.p.value <- res_1_lrt@listData$padj
reduced_output_table_1$LRT.p.value <- res_1_lrt@listData$padj


### Generate and add group mean and stdev columns to the DGE table ###

tcounts <- as.data.frame(t(normCounts))
tcounts$group <- names(group) # Used final table group name formatting (e.g. '( Space Flight & Blue Light )' )

group_means <- as.data.frame(t(aggregate(. ~ group,data = tcounts,mean))) # Compute group name group-wise means
colnames(group_means) <- paste0("Group.Mean_", group_means['group',]) # assign group name as column names

group_stdev <- as.data.frame(t(aggregate(. ~ group,data = tcounts,sd))) # Compute group name group-wise standard deviation
colnames(group_stdev) <- paste0("Group.Stdev_", group_stdev['group',]) # assign group name as column names

group_means <- group_means[-c(1),] # Drop group name row from data rows (now present as column names)
group_stdev <- group_stdev[-c(1),] # Drop group name row from data rows (now present as column names)
output_table_1 <- cbind(output_table_1,group_means, group_stdev) # Column bind the group-wise data
reduced_output_table_1 <- cbind(reduced_output_table_1,group_means, group_stdev) # Column bind the group-wise data

rm(group_stdev,group_means,tcounts)


### Add columns needed to generate GeneLab visualization plots to the DGE table ###

## Add column to indicate the sign (positive/negative) of log2fc for each pairwise comparison ##

updown_table <- sign(output_table_1[,grep("Log2fc_",colnames(output_table_1))])
colnames(updown_table) <- gsub("Log2fc","Updown",grep("Log2fc_",colnames(output_table_1),value = TRUE))
output_table_1 <- cbind(output_table_1,updown_table)
rm(updown_table)


## Add column to indicate contrast significance with p <= 0.1 ##

sig.1_table <- output_table_1[,grep("P.value_",colnames(output_table_1))]<=.1
colnames(sig.1_table) <- gsub("P.value","Sig.1",grep("P.value_",colnames(output_table_1),value = TRUE))
output_table_1 <- cbind(output_table_1,sig.1_table)
rm(sig.1_table)


## Add column to indicate contrast significance with p <= 0.05 ##

sig.05_table <- output_table_1[,grep("P.value_",colnames(output_table_1))]<=.05
colnames(sig.05_table) <- gsub("P.value","Sig.05",grep("P.value_",colnames(output_table_1),value = TRUE))
output_table_1 <- cbind(output_table_1,sig.05_table)
rm(sig.05_table)


## Add columns for the volcano plot with p-value and adjusted p-value ##

log_pval_table <- log2(output_table_1[,grep("P.value_",colnames(output_table_1))])
colnames(log_pval_table) <- paste0("Log2_",colnames(log_pval_table))
output_table_1 <- cbind(output_table_1,log_pval_table)
rm(log_pval_table)
log_adj_pval_table <- log2(output_table_1[,grep("Adj.p.value_",colnames(output_table_1))])
colnames(log_adj_pval_table) <- paste0("Log2_",colnames(log_adj_pval_table))
output_table_1 <- cbind(output_table_1,log_adj_pval_table)
rm(log_adj_pval_table)


## Prepare PCA table for GeneLab visualization plots ##

exp_raw <- log2(normCounts)
PCA_raw <- prcomp(t(exp_raw), scale = FALSE)


### Read in GeneLab annotation table for the organism of interest ###

annot <- read.table(annotations_link, sep = "\t", header = TRUE, quote = "", comment.char = "", row.names = 1)


### Combine annotations table and the DGE table ###

output_table_1 <- merge(annot, output_table_1, by='row.names', all.y=TRUE)
output_table_1 <- output_table_1 %>% 
  rename(
    ENSEMBL = Row.names ## Change ENSEMBL to TAIR for plant studies ##
  )


reduced_output_table_1 <- merge(annot, reduced_output_table_1, by='row.names', all.y=TRUE)
reduced_output_table_1 <- reduced_output_table_1 %>% 
  rename(
    ENSEMBL = Row.names ## Change ENSEMBL to TAIR for plant studies ##
  )

```

<br>

### 9j. Export GeneLab DGE Tables with Annotations for Datasets Without ERCC Spike-In

```R
### Export unnormalized and normalized counts tables ###

normCounts_exp <- as.data.frame(counts(dds_1, normalized=TRUE))

write.csv(txi.rsem$counts,file.path(norm_output, "RSEM_Unnormalized_Counts.csv"))
write.csv(normCounts_exp,file.path(norm_output, "Normalized_Counts.csv"))


### Export sample grouping and contrasts tables ###

write.csv(sampleTable,file.path(DGE_output, "SampleTable.csv"))

write.csv(contrasts,file.path(DGE_output, "contrasts.csv"))


### Export human-readable DGE table ###

write.csv(reduced_output_table_1,file.path(DGE_output, "differential_expression.csv"), row.names = FALSE)


### Export computer-readable DGE and PCA tables used for GeneLab visualization ###

write.csv(output_table_1,file.path(DGE_output, "visualization_output_table.csv"), row.names = FALSE)

write.csv(PCA_raw$x,file.path(DGE_output, "visualization_PCA_table.csv"), row.names = TRUE)


### print session info ###

print("Session Info below: ")
sessionInfo()

```


**Input Data:**

- *runsheet.csv file (table containing metadata required for analysis, output from [step 9a](#9a-create-sample-runsheet))
- [GL-DPPD-7110_annotations.csv](../../GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110/GL-DPPD-7110_annotations.csv) (csv file containing link to GeneLab annotations) 
- *genes.results (RSEM counts per gene, output from [Step 8a](#8a-count-aligned-reads-with-rsem))

<a id=ERCCspikeOut></a>
**Output Data for Datasets with ERCC Spike-In:**

Output data without considering ERCC spike-in genes:

- RSEM_Unnormalized_Counts.csv\# (table containing raw RSEM gene counts for each sample)
- Normalized_Counts.csv\# (table containing normalized gene counts for each sample)
- SampleTable.csv\# (table containing samples and their respective groups)
- visualization_output_table.csv (file used to generate GeneLab DGE visualizations)
- visualization_PCA_table.csv (file used to generate GeneLab PCA plots)
- differential_expression.csv\# (table containing normalized counts for each sample, group statistics, DESeq2 DGE results for each pairwise comparison, and gene annotations) 
- contrasts.csv\# (table containing all pairwise comparisons)

Output data with considering ERCC spike-in genes:

- ERCC_rawCounts_unfiltered.csv (table containing raw ERCC unfiltered counts)
- ERCC_rawCounts_filtered.csv (ERCC counts table after removing ERCC genes with low counts)
- ERCC_Normalized_Counts.csv\# (table containing ERCC-normalized gene counts for each sample)
- ERCCnorm_SampleTable.csv\# (table containing samples with detectable ERCC group B genes and their respective groups)
- visualization_output_table_ERCCnorm.csv (file used to generate GeneLab DGE visualizations for ERCC-normalized data)
- visualization_PCA_table_ERCCnorm.csv (file used to generate GeneLab PCA plots for ERCC-normalized data)
- ERCCnorm_differential_expression.csv\# (table containing ERCC-normalized counts for each sample, group statistics, DESeq2 DGE results for each pairwise comparison, and gene annotations)
- ERCCnorm_contrasts.csv\# (table containing all pairwise comparisons for samples containing ERCC spike-in)


**Output Data for Datasets without ERCC Spike-In:**

- RSEM_Unnormalized_Counts.csv\# (table containing raw RSEM gene counts for each sample)
- Normalized_Counts.csv\# (table containing normalized gene counts for each sample)
- SampleTable.csv\# (table containing samples and their respective groups)
- visualization_output_table.csv (file used to generate GeneLab DGE visualizations)
- visualization_PCA_table.csv (file used to generate GeneLab PCA plots)
- differential_expression.csv\# (table containing normalized counts for each sample, group statistics, DESeq2 DGE results for each pairwise comparison, and gene annotations) 
- contrasts.csv\# (table containing all pairwise comparisons)

> Note: RNAseq processed data interactive tables and plots are found in the [GLDS visualization portal](https://visualization.genelab.nasa.gov/data/studies).

<br>

---

## 10. Evaluate ERCC Spike-In Data 
> Note: This is only applicable for datasets with ERCC spike-in

<br>

### 10a. Evaluate ERCC Count Data in Python

```python
### Setting up the notebook

# import python packages
import pandas as pd
pd.set_option('mode.chained_assignment', None) # suppress chained indexing warnings
import numpy as np
from json import loads
from re import search
import zipfile
import seaborn as sns
from scipy.stats import linregress
import matplotlib.pyplot as plt


### Get and parse data and metadata

# Get and unzip ISA.zip to extract metadata.

accession = 'GLDS-NNN' # Replace Ns with GLDS number
isaPath = '/path/to/GLDS-NNN_metadata_GLDS-NNN-ISA.zip' # Replace with path to ISA archive file
zip_file_object =  zipfile.ZipFile(isaPath, "r")
list_of_ISA_files = zip_file_object.namelist()
UnnormalizedCountsPath = '/path/to/GLDS-NNN_rna_seq_RSEM_Unnormalized_Counts.csv'

# Print contents of ISA zip file to view file order
list_of_ISA_files

# There are datasets that have multiple assays (including microarray), so the RNAseq ISA files from the above output must be selected. 
# Txt files outputted above are indexed as 0, 1, 2, etc. Fill in the indexed number corresponding to the sample (s_*txt) and assay files for RNAseq (a_*_(RNA-Seq).txt) in the code block below.

# Extract metadata from the sample file (s_*txt)

sample_file = list_of_ISA_files[1] # replace [1] with index corresponding to the (s_*txt) file
file = zip_file_object.open(sample_file)
sample_table = pd.read_csv(zip_file_object.open(sample_file), sep='\t')


# Extract metadata from the assay (a_*_(RNA-Seq).txt) file

assay_file = list_of_ISA_files[0] # replace [0] with index corresponding to the (a_*_(RNA-Seq).txt) file
file = zip_file_object.open(assay_file)
assay_table = pd.read_csv(zip_file_object.open(assay_file), sep='\t')


# Check the sample table

pd.set_option('display.max_columns', None)
print(sample_table.head(n=3))


# Check the assay table

pd.set_option('display.max_columns', None)
assay_table.head(n=3)



# Get raw counts table

raw_counts_table = pd.read_csv(UnnormalizedCountsPath, index_col=0) 
raw_counts_table.index.rename('Gene_ID', inplace=True)
print(raw_counts_table.head(n=3))

raw_counts_transcripts = raw_counts_table[raw_counts_table.index.str.contains('^ENSMUSG')] # Change according to organism of interest
raw_counts_transcripts = raw_counts_transcripts.sort_values(by=list(raw_counts_transcripts), ascending=False)
print(raw_counts_transcripts)


# Get ERCC counts

ercc_counts = raw_counts_table[raw_counts_table.index.str.contains('^ERCC-')] 
ercc_counts.reset_index(inplace=True)
ercc_counts = ercc_counts.rename(columns={'Gene_ID':'ERCC ID'})
ercc_counts = ercc_counts.sort_values(by=list(ercc_counts), ascending=False)
print(ercc_counts.head())

# Get files containing ERCC gene concentrations and metadata

ercc_url = 'https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_095046.txt'
ercc_table = pd.read_csv(ercc_url, '\t')
print(ercc_table.head(n=3))



### Calculate the number of ERCC genes detected in each of the 4 (A, B, C and D) groups for each sample

# Extract ERCC counts and calculate the log(2)

meltERCC = ercc_counts.melt(id_vars=['ERCC ID'])
meltERCC['log2 Count'] = meltERCC['value']+1
meltERCC['log2 Count'] = np.log2(meltERCC['log2 Count'])
meltERCC = meltERCC.rename(columns={'variable':'Sample Name', 'value':'Count'})
print(meltERCC.head(n=3))

# Build Mix dictionary to link sample name to mix added and read depth using the assay table

mix_dict = assay_table.filter(['Sample Name','Parameter Value[Spike-in Mix Number]', 
                       'Parameter Value[Read Depth]'])
mix_dict = mix_dict.rename(columns={'Parameter Value[Spike-in Mix Number]':'Mix',
                                    'Parameter Value[Read Depth]':
                                    'Total Reads'})
print(mix_dict.head(n=3))


# Make combined ercc counts and assay table

merged_ercc = meltERCC.merge(mix_dict, on='Sample Name')
print(merged_ercc)

# Read ERCC info including concentrations from merged_ercc table

groupA = ercc_table.loc[ercc_table['subgroup'] == 'A']['ERCC ID']
groupB = ercc_table.loc[ercc_table['subgroup'] == 'B']['ERCC ID']
groupC = ercc_table.loc[ercc_table['subgroup'] == 'C']['ERCC ID']
groupD = ercc_table.loc[ercc_table['subgroup'] == 'D']['ERCC ID']


# Make a dictionary for ERCC groups

group_dict = dict(zip(ercc_table['ERCC ID'], ercc_table['subgroup']))


# Calculate ERCC counts per million and log(2) counts per million

merged_ercc['Count per million'] = merged_ercc['Count'] / (merged_ercc['Total Reads'] / 1000000.0)
merged_ercc['log2 Count per million'] = np.log2(merged_ercc['Count per million']+1)


# Add ERCC group column

merged_ercc['ERCC group'] = merged_ercc['ERCC ID'].map(group_dict)
merged_ercc = merged_ercc.sort_values(by=['Mix'], ascending=True)
print(merged_ercc)


### Filter and calculate mean counts per million of Mix1 and Mix2 spiked samples in each of the 4 groups
## Check the 'Parameter Value[Spike-in Mix Number]' column of the assay table
## If the values do not have a space, change 'Mix 1' and 'Mix 2' to 'Mix1' and 'Mix2', respectively

# Filter Mix1 CPM and Mix2 CPM in group A 

Adf = merged_ercc.loc[merged_ercc['ERCC group'] == 'A']
Amix1df = Adf.loc[Adf['Mix']=='Mix 1']
Amix1df['Mix1 CPM'] = Amix1df[Amix1df['Count per million'] > 0]['Count per million'].dropna()
Amix1df = Amix1df.groupby('ERCC ID')['Mix1 CPM'].agg(np.mean).rename('Avg Mix1 CPM')
Amix1df = Amix1df.to_frame()
Amix2df = Adf.loc[Adf['Mix']=='Mix 2']
Amix2df['Mix2 CPM'] = Amix2df[Amix2df['Count per million'] > 0]['Count per million'].dropna()
Amix2df = Amix2df.groupby('ERCC ID')['Mix2 CPM'].agg(np.mean).rename('Avg Mix2 CPM')
Amix2df = Amix2df.to_frame()

adf = Amix1df.merge(Amix2df, on='ERCC ID', suffixes=('', '_2'))
adf = adf.reset_index()
adf['Avg Mix1 CPM/ Avg Mix2 CPM'] = (adf['Avg Mix1 CPM'] / adf['Avg Mix2 CPM'])


# Filter Mix1 CPM and Mix2 CPM in group B

Bdf = merged_ercc.loc[merged_ercc['ERCC group'] == 'B']
Bmix1df = Bdf.loc[Bdf['Mix']=='Mix 1']
Bmix1df['Mix1 CPM'] = Bmix1df[Bmix1df['Count per million'] > 0]['Count per million'].dropna()
Bmix1df = Bmix1df.groupby('ERCC ID')['Mix1 CPM'].agg(np.mean).rename('Avg Mix1 CPM')
Bmix1df = Bmix1df.to_frame()
Bmix2df = Bdf.loc[Bdf['Mix']=='Mix 2']
Bmix2df['Mix2 CPM'] = Bmix2df[Bmix2df['Count per million'] > 0]['Count per million'].dropna()
Bmix2df = Bmix2df.groupby('ERCC ID')['Mix2 CPM'].agg(np.mean).rename('Avg Mix2 CPM')
Bmix2df = Bmix2df.to_frame()

bdf = Bmix1df.merge(Bmix2df, on='ERCC ID')
bdf = bdf.reset_index()
bdf['Avg Mix1 CPM/ Avg Mix2 CPM'] = (bdf['Avg Mix1 CPM'] / bdf['Avg Mix2 CPM'])


# Filter Mix1 CPM and Mix2 CPM in group C

Cdf = merged_ercc.loc[merged_ercc['ERCC group'] == 'C']
Cmix1df = Cdf.loc[Cdf['Mix']=='Mix 1']
Cmix1df['Mix1 CPM'] = Cmix1df[Cmix1df['Count per million'] > 0]['Count per million'].dropna()
Cmix1df = Cmix1df.groupby('ERCC ID')['Mix1 CPM'].agg(np.mean).rename('Avg Mix1 CPM')
Cmix1df = Cmix1df.to_frame()
Cmix2df = Cdf.loc[Cdf['Mix']=='Mix 2']
Cmix2df['Mix2 CPM'] = Cmix2df[Cmix2df['Count per million'] > 0]['Count per million'].dropna()
Cmix2df = Cmix2df.groupby('ERCC ID')['Mix2 CPM'].agg(np.mean).rename('Avg Mix2 CPM')
Cmix2df = Cmix2df.to_frame()

cdf = Cmix1df.merge(Cmix2df, on='ERCC ID')
cdf = cdf.reset_index()
cdf['Avg Mix1 CPM/ Avg Mix2 CPM'] = (cdf['Avg Mix1 CPM'] / cdf['Avg Mix2 CPM'])


# Filter Mix1 CPM and Mix2 CPM in group D

Ddf = merged_ercc.loc[merged_ercc['ERCC group'] == 'D']
Dmix1df = Ddf.loc[Ddf['Mix']=='Mix 1']
Dmix1df['Mix1 CPM'] = Dmix1df[Dmix1df['Count per million'] > 0]['Count per million'].dropna()
Dmix1df = Dmix1df.groupby('ERCC ID')['Mix1 CPM'].agg(np.mean).rename('Avg Mix1 CPM')
Dmix1df = Dmix1df.to_frame()
Dmix2df = Ddf.loc[Ddf['Mix']=='Mix 2']
Dmix2df['Mix2 CPM'] = Dmix2df[Dmix2df['Count per million'] > 0]['Count per million'].dropna()
Dmix2df = Dmix2df.groupby('ERCC ID')['Mix2 CPM'].agg(np.mean).rename('Avg Mix2 CPM')
Dmix2df = Dmix2df.to_frame()

ddf = Dmix1df.merge(Dmix2df, on='ERCC ID')
ddf = ddf.reset_index()
ddf['Avg Mix1 CPM/ Avg Mix2 CPM'] = (ddf['Avg Mix1 CPM'] / ddf['Avg Mix2 CPM'])


##### Multi-sample ERCC analyses

### Create box and whisker plots of the log(2) CPM for each ERCC detected in group A in Mix 1 and Mix 2 spiked samples

a = sns.catplot(x="ERCC ID", y="log2 Count per million", order=groupA, hue="Mix",data=merged_ercc[merged_ercc['ERCC ID'].isin(groupA)], kind="box", col="ERCC group", height=5, aspect=1, palette=sns.color_palette(['blue', 'orange']))
a.set_xticklabels(rotation=90)
plt.text(23,2.5,"Mix1/ Mix2 = 4")

### Create bar plot of the average Mix1 CPM / average Mix 2 CPM for group A ERCC genes (for group A we expect Mix 1 CPM / Mix 2 CPM = 4)

a1 = sns.catplot(x="ERCC ID", y="Avg Mix1 CPM/ Avg Mix2 CPM", palette="rocket_r", data=adf, kind="bar", height=5, aspect=1, linewidth=0.5)
a1.set_xticklabels(rotation=90)
plt.title("ERCC Group A")
a1.set(ylim=(0, 6))
print('Number of ERCC detected in group A (out of 23) =', adf['Avg Mix1 CPM/ Avg Mix2 CPM'].count())


### Create box and whisker plots of the log(2) CPM for each ERCC detected in group B in Mix 1 and Mix 2 spiked samples

b = sns.catplot(x="ERCC ID", y="log2 Count per million", order=groupB, hue="Mix", data=merged_ercc[merged_ercc['ERCC ID'].isin(groupB)], kind="box", col="ERCC group", height=5, aspect=1, palette=sns.color_palette(['blue', 'orange']))
b.set_xticklabels(rotation=90)
plt.text(23,2.5,"Mix1/ Mix2 = 1")

### Create bar plot of the average Mix1 CPM / average Mix 2 CPM for group B ERCC genes (for group B we expect Mix 1 CPM / Mix 2 CPM = 1)

b = sns.catplot(x="ERCC ID", y="Avg Mix1 CPM/ Avg Mix2 CPM", palette="rocket_r", data=bdf, kind="bar", 
               height=5, aspect=1, linewidth=0.5)
b.set_xticklabels(rotation=90)
plt.title("ERCC Group B")
b.set(ylim=(0, 2))
print('Number of ERCC detected in group B (out of 23) =', bdf['Avg Mix1 CPM/ Avg Mix2 CPM'].count())


### Create box and whisker plots of the log(2) CPM for each ERCC detected in group C in Mix 1 and Mix 2 spiked samples

c = sns.catplot(x="ERCC ID", y="log2 Count per million", order=groupC, hue="Mix", data=merged_ercc[merged_ercc['ERCC ID'].isin(groupC)], kind="box", col="ERCC group", height=5, aspect=1, palette=sns.color_palette(['blue', 'orange']))
c.set_xticklabels(rotation=90)
plt.text(23,2.5,"Mix1/ Mix2 = 0.67")

### Create bar plot of the average Mix1 CPM / average Mix 2 CPM for group C ERCC genes (for group C we expect Mix 1 CPM / Mix 2 CPM = 0.67)

c = sns.catplot(x="ERCC ID", y="Avg Mix1 CPM/ Avg Mix2 CPM", palette="rocket_r", data=cdf, kind="bar", 
               height=5, aspect=1, linewidth=0.5)
c.set_xticklabels(rotation=90)
plt.title("ERCC Group C")
c.set(ylim=(0, 2))
print('Number of ERCC detected in group C (out of 23) =', cdf['Avg Mix1 CPM/ Avg Mix2 CPM'].count())


### Create box and whisker plots of the log(2) CPM for each ERCC detected in group D in Mix 1 and Mix 2 spiked samples

d = sns.catplot(x="ERCC ID", y="log2 Count per million", order=groupD, hue="Mix", data=merged_ercc[merged_ercc['ERCC ID'].isin(groupD)], col="ERCC group", kind="box", height=5, aspect=1, palette=sns.color_palette(['blue', 'orange']))
d.set_xticklabels(rotation=90)
plt.text(23,2.5,"Mix1/ Mix2 = 0.5")

### Create bar plot of the average Mix1 CPM / average Mix 2 CPM for group D ERCC genes (for group D we expect Mix 1 CPM / Mix 2 CPM = 0.5)

d = sns.catplot(x="ERCC ID", y="Avg Mix1 CPM/ Avg Mix2 CPM", palette="rocket_r", data=ddf, kind="bar", 
               height=5, aspect=1, linewidth=0.5)
d.set_xticklabels(rotation=90)
plt.title("ERCC Group D")
d.set(ylim=(0, 1))
print('Number of ERCC detected in group D (out of 23) =', ddf['Avg Mix1 CPM/ Avg Mix2 CPM'].count())



##### Individual sample ERCC analyses

#Calculate and plot ERCC metrics from individual samples, including limit of detection, dynamic range, and R^2 of counts vs. concentration.

#View the ERCC table
print(ercc_table.head(n=3))

# Check assay_table header to identify the 'Sample Name' column and the column title indicating the 'Spike-in Mix Nmber' if it's indicated in the metadata.

pd.set_option('display.max_columns', None)
print(assay_table.head(n=3))

# Get ERCC counts for all samples

ercc_counts = raw_counts_table[raw_counts_table.index.str.contains('^ERCC-')] 
ercc_counts = ercc_counts.sort_values(by=list(ercc_counts), ascending=False)
print(ercc_counts.head())


# Make a dictionary for ERCC concentrations for Mix 1
mix1_conc_dict = dict(zip(ercc_table['ERCC ID'], ercc_table['concentration in Mix 1 (attomoles/ul)']))

# Get samples spiked with Mix 1
mix1_samples = assay_table[assay_table['Parameter Value[Spike-in Mix Number]'] == 'Mix 1']['Sample Name']

# Get ERCC counts for Mix 1 spiked samples
ercc_counts_mix_1 = ercc_counts[mix1_samples]
ercc_counts_mix_1['ERCC conc (attomoles/ul)'] = ercc_counts_mix_1.index.map(mix1_conc_dict)
print(ercc_counts_mix_1.head(n=3))


# Make a dictionary for ERCC concentrations for Mix 2
mix2_conc_dict = dict(zip(ercc_table['ERCC ID'], ercc_table['concentration in Mix 2 (attomoles/ul)']))

# Get samples spiked with Mix 2
mix2_samples = assay_table[assay_table['Parameter Value[Spike-in Mix Number]'] == 'Mix 2']['Sample Name']

# Get ERCC counts for Mix 2 spiked samples

ercc_counts_mix_2 = ercc_counts[mix2_samples]
ercc_counts_mix_2['ERCC conc (attomoles/ul)'] = ercc_counts_mix_2.index.map(mix2_conc_dict)
print(ercc_counts_mix_2.head(n=3))


# Create a scatter plot of log(2) ERCC counts versus log(2) ERCC concentration for each sample

columns_mix_1 = ercc_counts_mix_1.columns.drop(['ERCC conc (attomoles/ul)'])
columns_mix_2 = ercc_counts_mix_2.columns.drop(['ERCC conc (attomoles/ul)'])
all_columns = columns_mix_1.to_list() + columns_mix_2.to_list()
total_columns = len(columns_mix_1) + len(columns_mix_2) 
side_size = np.int32(np.ceil(np.sqrt(total_columns)))# calculate grid side size. take sqrt of total plots and round up.
fig, axs = plt.subplots(side_size, side_size, figsize=(22,26), sharex='all', sharey='all'); #change figsize x,y labels if needed.
fig.tight_layout(pad=1, w_pad=2.5, h_pad=3.5)

counter = 0
for ax in axs.flat:
    
    if(counter < len(columns_mix_1)):
      ax.scatter(x=np.log2(ercc_counts_mix_1['ERCC conc (attomoles/ul)']), y=np.log2(ercc_counts_mix_1[all_columns[counter]]+1), s=7);
      ax.set_title(all_columns[counter][-45:], fontsize=9);
      ax.set_xlabel('log2 ERCC conc (attomoles/ ul)', fontsize=9);
      ax.set_ylabel('log2 Counts per million', fontsize=9);
      ax.tick_params(direction='in', axis='both', labelsize=9, labelleft=True, labelbottom=True);
      
    elif(counter >= len(columns_mix_1) and counter < total_columns):
      ax.scatter(x=np.log2(ercc_counts_mix_2['ERCC conc (attomoles/ul)']), y=np.log2(ercc_counts_mix_2[all_columns[counter]]+1), s=7);
      ax.set_title(all_columns[counter][-45:], fontsize=9);
      ax.set_xlabel('log2 ERCC conc (attomoles/ ul)', fontsize=9);
      ax.set_ylabel('log2 Counts per million', fontsize=9);
      ax.tick_params(direction='in', axis='both', labelsize=9, labelleft=True, labelbottom=True);
       
    else:
      pass

    counter = counter + 1


## Calculate and plot linear regression of log(2) ERCC counts versus log(2) ERCC concentration for each sample

# Filter counts > 0

nonzero_counts_list_1 = []
for i in range(0, len(ercc_counts_mix_1.columns)-1):
  counts = ercc_counts_mix_1[columns_mix_1[i]]
  counts.index.rename('Gene_ID', inplace=True)
  countsdf = pd.DataFrame(counts)
  nonzero_counts = countsdf[ercc_counts_mix_1[columns_mix_1[i]] > 0.0]
  nonzero_counts['Conc'] = nonzero_counts.index.map(mix1_conc_dict)
  nonzero_counts.columns = ['Counts','Conc']
  nonzero_counts_sorted = nonzero_counts.sort_values('Conc')
  nonzero_counts_list_1.append(nonzero_counts_sorted)

nonzero_counts_list_2 = []
for i in range(0, len(ercc_counts_mix_2.columns)-1):
  counts = ercc_counts_mix_2[columns_mix_2[i]]
  counts.index.rename('Gene_ID', inplace=True)
  countsdf = pd.DataFrame(counts)
  nonzero_counts = countsdf[ercc_counts_mix_2[columns_mix_2[i]] > 0.0]
  nonzero_counts['Conc'] = nonzero_counts.index.map(mix2_conc_dict)
  nonzero_counts.columns = ['Counts','Conc']
  nonzero_counts_sorted = nonzero_counts.sort_values('Conc')
  nonzero_counts_list_2.append(nonzero_counts_sorted)


# Plot each sample using linear regression of scatter plot with x = log2 Conc and y = log2 Counts.  Return min, max, R^2 and dynamic range (max / min) values.

samples = []
mins = []
maxs = []
dyranges = []
rs = []

fig, axs = plt.subplots(side_size, side_size, figsize=(22,26), sharex='all', sharey='all');
fig.tight_layout(pad=1, w_pad=2.5, h_pad=3.5)

counter = 0
list2counter = 0
for ax in axs.flat:
    
    if(counter < len(columns_mix_1)):

      nonzero_counts = nonzero_counts_list_1[counter]
      xvalues = nonzero_counts['Conc']
      yvalues = nonzero_counts['Counts']

      sns.regplot(x=np.log2(xvalues), y=np.log2(yvalues), ax=ax);
      ax.set_title(all_columns[counter][-47:], fontsize=9);
      ax.set_xlabel('log2 Conc (attomoles/ul)', fontsize=9);
      ax.set_ylabel('log2 Counts per million', fontsize=9);
      ax.tick_params(direction='in', axis='both', labelsize=9, labelleft=True, labelbottom=True)
      samples.append(all_columns[counter])

      if(len(xvalues) == 0):
        mins.append('NaN')
        maxs.append('NaN')
        dyranges.append('NaN')
        rs.append('NaN')

    
      else:
        min = xvalues[0];
        mins.append(min)
        minimum = f'Min:{min:.1f}';
        max = xvalues[-1];
        maxs.append(max)
        maximum = f'Max:{max:.1f}';
        dynamic_range = max / min;
        dyranges.append(dynamic_range)
        dyn_str = f'Dyn:{dynamic_range:.1f}';

        ax.text(0.02, 0.98, minimum,
        verticalalignment='top', horizontalalignment='left',
        transform=ax.transAxes,
        color='black', fontsize=10);
      
        ax.text(0.02, 0.88, maximum,
        verticalalignment='top', horizontalalignment='left',
        transform=ax.transAxes,
        color='black', fontsize=10);
      
        ax.text(0.02, 0.78, dyn_str,verticalalignment='top',
                horizontalalignment='left',transform=ax.transAxes,
                color='black', fontsize=10);
      
        if(len(xvalues) == 1):
          rs.append('NaN')

        else:
          slope, intercept, r, p, se = linregress(np.log2(xvalues), y=np.log2(yvalues))
          r_str = f'R:{r:.2f}'
          rs.append(r)

          ax.text(0.02, 0.68, r_str, verticalalignment='top',
                  horizontalalignment='left',transform=ax.transAxes,
                  color='black', fontsize=10);
    
    elif(counter >= len(columns_mix_1) and counter < total_columns):
      
      nonzero_counts = nonzero_counts_list_2[list2counter]
      xvalues = nonzero_counts['Conc']
      yvalues = nonzero_counts['Counts']

      sns.regplot(x=np.log2(xvalues), y=np.log2(yvalues), ax=ax);
      ax.set_title(all_columns[counter][-47:], fontsize=9);
      ax.set_xlabel('log2 Conc (attomoles/ul)', fontsize=9);
      ax.set_ylabel('log2 Counts per million', fontsize=9);
      ax.tick_params(direction='in', axis='both', labelsize=9, labelleft=True, labelbottom=True);
      samples.append(all_columns[counter])


      if(len(xvalues) == 0):
        mins.append('NaN')
        maxs.append('NaN')
        dyranges.append('NaN')
        rs.append('NaN')
    
      else:
        min = xvalues[0];
        mins.append(min)
        minimum = f'Min:{min:.1f}';
        max = xvalues[-1];
        maxs.append(max)
        maximum = f'Max:{max:.1f}';
        dynamic_range = max / min;
        dyranges.append(dynamic_range)
        dyn_str = f'Dyn:{dynamic_range:.1f}';

        ax.text(0.02, 0.98, minimum,
        verticalalignment='top', horizontalalignment='left',
        transform=ax.transAxes,
        color='black', fontsize=10);
      
        ax.text(0.02, 0.88, maximum,
        verticalalignment='top', horizontalalignment='left',
        transform=ax.transAxes,
        color='black', fontsize=10);
      
        ax.text(0.02, 0.78, dyn_str,verticalalignment='top',
                horizontalalignment='left',transform=ax.transAxes,
                color='black', fontsize=10);
      
        if(len(xvalues) == 1):
          rs.append('NaN')
          
        else:
          slope, intercept, r, p, se = linregress(np.log2(xvalues), y=np.log2(yvalues));
          r_str = f'R:{r:.2f}';
          rs.append(r)

          ax.text(0.02, 0.68, r_str, verticalalignment='top',
                  horizontalalignment='left',transform=ax.transAxes,
                  color='black', fontsize=10);

      list2counter = list2counter + 1
    
    else:
      pass

    counter = counter + 1


# Create directory for saved files

import os
os.makedirs(name="ERCC_analysis", exist_ok=True)

# Print tables containing the dynamic range and R^2 values for each sample.
# Remember to change file names to GLDS# analyzing

stats = pd.DataFrame(list(zip(samples, mins, maxs, dyranges, rs)))
stats.columns = ['Samples', 'Min', 'Max', 'Dynamic range', 'R']
stats.to_csv('ERCC_analysis/ERCC_stats_GLDS-NNN.csv', index = False) 
stats.filter(items = ['Samples', 'Dynamic range']).to_csv('ERCC_analysis/ERCC_dynrange_GLDS-NNN_mqc.csv', index = False) 
stats.filter(items = ['Samples', 'R']).to_csv('ERCC_analysis/ERCC_rsq_GLDS-NNN_mqc.csv', index = False) 



### Generate data and metadata files needed for ERCC DESeq2 analysis

# ERCC Mix 1 and Mix 2 are distributed so that half the samples receive Mix 1 spike-in and half receive Mix 2 spike-in. Transcripts in Mix 1 and Mix 2 are present at a known ratio, so we can determine how well these patterns are revealed in the dataset.

# Get sample table

combined = sample_table.merge(assay_table, on='Sample Name')
combined = combined.set_index(combined['Sample Name'])
pd.set_option('display.max_columns', None)
print(combined)

# Create metadata table containing samples and their respective ERCC spike-in Mix number
# Sometimes Number in [Spike-in Mix Number] is spelled 'number' and this could cause error in mismatch search 

ERCCmetadata = combined[['Parameter Value[Spike-in Mix Number]']]
ERCCmetadata.index = ERCCmetadata.index.str.replace('-','_')
ERCCmetadata.columns = ['Mix']
#ERCCmetadata = ERCCmetadata.rename(columns={'Parameter Value[Spike-in Mix Number]':'Mix'})
print(ERCCmetadata)

# Export ERCC sample metadata

ERCCmetadata.to_csv('ERCC_analysis/ERCCmetadata.csv') 

# Export ERCC count data

ercc_counts.columns = ercc_counts.columns.str.replace('-','_')
ERCCcounts = ercc_counts.loc[:,ERCCmetadata.index]
ERCCcounts.head()

ERCCcounts.to_csv('ERCC_analysis/ERCCcounts.csv') 
```


**Input Data:**

- *ISA.zip (compressed ISA directory containing Investigation, Study, and Assay (ISA) metadata files for the respective GLDS dataset, output from [Step 9a](#9a-create-sample-runsheet))
- RSEM_Unnormalized_Counts.csv (RSEM raw counts table, output from [Step 9](#ERCCspikeOut))

**Output Data:**

- ERCC_analysis/ERCC_stats_GLDS-*.csv (Samplewise counts statistics table containing 'Min', 'Max', 'Dynamic range', 'R')
- ERCC_analysis/ERCC_dynrange_GLDS-*.csv (Samplewise counts statistics subset table containing 'Dynamic range')
- ERCC_analysis/ERCC_rsq_GLDS-*.csv (Samplewise counts statistics subset table containing 'R')
- ERCC_analysis/ERCCmetadata.csv (Samplewise metadata table inlcuding ERCC mix number)
- ERCC_analysis/ERCCcounts.csv (Samplewise ERCC counts table)

<br>

### 10b. Perform DESeq2 Analysis of ERCC Counts in R

```R

## Install R packages if not already installed

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")

## Import DESeq2 library

library("DESeq2")

## Import and format ERCC count data and metadata

cts <- as.matrix(read.csv('ERCC_analysis/ERCCcounts.csv',sep=",",row.names="Gene_ID")) #INPUT
coldata <- read.csv('ERCC_analysis/ERCCmetadata.csv', row.names=1) #INPUT

coldata$Mix <- factor(coldata$Mix)
all(rownames(coldata) == colnames(cts))


## Make DESeqDataSet object

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ Mix)
dds


## Filter out ERCC genes with counts of less than 10 in all samples #####

keepGenes <- rowSums(counts(dds)) > 10
dds <- dds[keepGenes,]

dds


## Run DESeq2 analysis and calculate results

dds <- DESeq(dds)
res <- results(dds, contrast=c("Mix","Mix 1","Mix 2")) # remove space before mix number if needed
res

## Export DESeq2 results table and normalized ERCC counts table

write.csv(res, 'ERCC_analysis/ERCC_DESeq2.csv') #OUTPUT
normcounts = counts(dds, normalized=TRUE)
write.csv(normcounts, 'ERCC_analysis/ERCC_normcounts.csv') #OUTPUT
```

**Input Data:**

- ERCC_analysis/ERCCmetadata.csv (Samplewise metadata table inlcuding ERCC mix number, output from [Step 10a](#10a-evaluate-ercc-count-data-in-python))
- ERCC_analysis/ERCCcounts.csv (Samplewise ERCC counts table, output from [Step 10a](#10a-evaluate-ercc-count-data-in-python))

**Output Data:**

- ERCC_analysis/ERCC_DESeq2.csv (DESeq2 results table)
- ERCC_analysis/ERCC_normcounts.csv (Normalized ERCC Counts table)

<br>

### 10c. Analyze ERCC DESeq2 Results in Python

```python

# Import python packages

import pandas as pd
from urllib.request import urlopen, quote, urlretrieve
import seaborn as sns
import matplotlib.pyplot as plt


# Import ERCC DESeq2 results

deseq2out = pd.read_csv('ERCC_analysis/ERCC_DESeq2.csv', index_col=0) # INPUT
#deseq2out.index = deseq2out.index.str.replace('_','-')
deseq2out.rename(columns ={'baseMean' : 'meanNormCounts'}, inplace = True)
print(deseq2out.head())


# Get files containing ERCC gene concentrations and metadata

ercc_url = 'https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_095046.txt'
ercc_table = pd.read_csv(ercc_url, '\t', index_col='ERCC ID')
print(ercc_table.head(n=3))


# Combine ERCC DESeq2 results and ercc_table

combined = deseq2out.merge(ercc_table, left_index=True, right_index=True)
print(combined.head())


# Filter p-value and adj. p-value cutoff at 10^-3

combined['cleaned_padj'] = combined['padj']
combined.loc[(combined.cleaned_padj < 0.001),'cleaned_padj']=0.001

combined['cleaned_pvalue'] = combined['pvalue']
combined.loc[(combined.cleaned_pvalue < 0.001),'cleaned_pvalue']=0.001

print(combined.head())


# Export the filtered combined ERCC DESeq2 results and ercc_table
# Remember to change file name to GLDS# analyzing

combined.filter(items = ['ERCC ID', 'meanNormCounts', 'cleaned_pvalue','cleaned_padj']).to_csv('ERCC_analysis/ERCC_lodr_GLDS-NNN_mqc.csv') 


# Plot p-value vs. mean normalized ERCC counts

fig, ax = plt.subplots(figsize=(10, 7))

sns.scatterplot(data=combined, x="meanNormCounts", y="cleaned_pvalue",
            hue="expected fold-change ratio",
                palette=['red','green','black','blue'], ax=ax)

sns.lineplot(data=combined, x="meanNormCounts", y="cleaned_pvalue",
            hue="expected fold-change ratio",
                palette=['red','green','black','blue'], ax=ax)

#g.set_xscale("log", base=2)
ax.set_xscale("linear");
ax.set_yscale("log");


# Plot Adjp-value vs. mean normalized ERCC counts

fig, ax = plt.subplots(figsize=(10, 7))

sns.scatterplot(data=combined, x="meanNormCounts", y="cleaned_padj",
            hue="expected fold-change ratio",
                palette=['red','green','black','blue'], ax=ax)

sns.lineplot(data=combined, x="meanNormCounts", y="cleaned_padj",
            hue="expected fold-change ratio",
                palette=['red','green','black','blue'], ax=ax)

#g.set_xscale("log", base=2)
ax.set_xscale("linear");
ax.set_yscale("log");
```

**Input Data:**

- ERCC_analysis/ERCC_DESeq2.csv (ERCC DESeq2 results table, output from [Step 10b](#10b-perform-deseq2-analysis-of-ercc-counts-in-r))

**Output Data:**

- ERCC_analysis/ERCC_lodr_*.csv (ERCC Gene Table including mean counts, adjusted p-value and p-value, and filtered to genes with both adj. p-value and p-value < 0.001)

> All steps of the ERCC Spike-In Data Analysis are performed in a Jupyter Notebook (JN) and the completed JN is exported as an html file and published in the [GLDS repository](https://genelab-data.ndc.nasa.gov/genelab/projects) for the respective dataset.
