# GeneLab bioinformatics processing pipeline for 10X Chromium 3' single cell RNA-sequencing data

> **This page holds an overview and instructions for how GeneLab processes single cell RNA-sequencing (scRNAseq) datasets. Exact processing commands, GL-DPPD-7111 version used, and processed data output files for specific datasets are provided in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/).**  
>  
  > Note: Unlike [bulkRNAseq](../../../RNAseq), for scRNAseq trimming is performed during the [STARsolo alignment step](#3a-align-reads-to-reference-genome-with-starsolo) of the pipeline.

---

**Date:** August 18, 2022  
**Revision:** -  
**Document Number:** GL-DPPD-7111  

**Submitted by:**  
Lauren Sanders (GeneLab Data Processing Team)

**Approved by:**  
Amanda Saravia-Butler (GeneLab Data Processing Lead)  
Sylvain Costes (GeneLab Project Manager)  
Samrawit Gebre (GeneLab Deputy Project Manager and Interim GeneLab Configuration Manager)  
Jonathan Galazka (GeneLab Project Scientist)

---

# Table of contents  

- [**Software used**](#software-used)
- [**General processing overview with example commands**](#general-processing-overview-with-example-commands)
  - [**1. Raw Data QC**](#1-raw-data-qc)
    - [1a. Raw Data QC](#1a-raw-data-qc)
    - [1b. Compile Raw Data QC](#1b-compile-raw-data-qc)
  - [**2. Build STAR Reference**](#2-build-star-reference)
  - [**3. Align Reads to Reference Genome**](#3-align-reads-to-reference-genome)
    - [3a. Align Reads to Reference Genome with STARsolo](#3a-align-reads-to-reference-genome-with-starsolo)
    - [3b. Compile Alignment Logs](#3b-compile-alignment-logs)

---

# Software used  

|Program|Version|Relevant Links|
|:------|:------:|:-------------|
|FastQC|0.11.9|[https://www.bioinformatics.babraham.ac.uk/projects/fastqc/](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)|
|MultiQC|1.12|[https://multiqc.info/](https://multiqc.info/)|
|STAR|2.7.10a|[https://github.com/alexdobin/STAR](https://github.com/alexdobin/STAR)|

---

# General processing overview with example commands  

> Exact processing commands and output files listed in **bold** below are included with each scRNAseq processed dataset in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/).

---

## 1. Raw Data QC

<br>

### 1a. Raw Data QC  

```bash
fastqc -o /path/to/raw_fastqc/output/directory *.fastq.gz
```

**Parameter Definitions:**

- `-o` – the output directory to store results
- `*.fastq.gz` – the input reads are specified as a positional argument, and can be given all at once with wildcards like this, or as individual arguments with spaces in between them

**Input Data:**

- *fastq.gz (raw reads)

**Output Data:**

- *fastqc.html (FastQC report)
- *fastqc.zip (FastQC data)

<br>

### 1b. Compile Raw Data QC  

```bash
multiqc --interactive \ 
  -n raw_multiqc \ 
  -o /path/to/raw_multiqc/output/directory \   
  /path/to/directory/containing/raw_fastqc/files
```

**Parameter Definitions:**

- `--interactive` - force reports to use interactive plots
- `-n` - prefix name for output files
- `-o` – the output directory to store results
- `/path/to/directory/containing/raw_fastqc/files` – the directory holding the output data from the fastqc run, provided as a positional argument

**Input Data:**

- *fastqc.zip (FastQC data, output from [Step 1a](#1a-raw-data-qc))

**Output Data:**

- **raw_multiqc.html** (multiqc report)
- **/raw_multiqc_data** (directory containing multiqc data)

<br>

---


## 2. Build STAR Reference  

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
- `--genomeSAindexNbases` - length (in bases) of the SA pre-indexing string, usually between 10 and 15. Longer strings require more memory but allow for faster searches. This value should be scaled down for smaller genomes (like bacteria) to min(14, 
log2(GenomeLength)/2 - 1). For example, for a 1 megaBase genome this value would be 9.
- `--genomeDir` - specifies the path to the directory where the STAR reference will be stored. At least 100GB of available disk space is required for mammalian genomes.
- `--genomeFastaFiles` - specifies one or more fasta file(s) containing the genome reference sequences
- `--sjdbGTFfile` – specifies the file(s) containing annotated transcripts in the standard gtf format
- `--sjdbOverhang` - indicates the length of the genomic sequence around the annotated junction to be used in constructing the splice junctions database. The length should be one less than the maximum length of the reads.

**Input Data:**

- *.fasta (genome sequence, this scRCP version uses the Ensembl fasta file indicated in the `fasta` column of the [GL-DPPD-7110_annotations.csv](../../../GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110/GL-DPPD-7110_annotations.csv) GeneLab Annotations file)
- *.gtf (genome annotation, this scRCP version uses the Ensembl gtf file indicated in the `gtf` column of the [GL-DPPD-7110_annotations.csv](../../../GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110/GL-DPPD-7110_annotations.csv) GeneLab Annotations file)

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

## 3. Align Reads to Reference Genome

<br>

### 3a. Align Reads to Reference Genome with STARsolo

```bash
STAR --runThreadN <NumberOfThreads> \
  --genomeDir /path/to/STAR/genome/directory \
  --soloType CB_UMI_Simple \ # Used for 10X Chromium data 
  --clipAdapterType CellRanger4 \  
  --outFilterScoreMin 30 \  
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
  --soloUMIfiltering MultiGeneUMI_CR \
  --soloUMIdedup 1MM_CR \
  --soloUMIlen 12 \
  --soloCellFilter EmptyDrops_CR <ExpectedCells> 0.99 10 45000 90000 500 0.01 20000 0.01 10000 \
  --soloMultiMappers EM \
  --outSAMattributes NH HI nM AS CR UR GX GN sS sQ sM \
  --outSAMtype BAM Unsorted \
  --soloFeatures Gene GeneFull SJ Velocyto \
  --readFilesCommand zcat \
  --soloCBwhitelist CellBarcodeWhitelist \
  --outFileNamePrefix /path/to/STAR/output/directory/<sample_id> \
  --readFilesIn /path/to/cDNA_reads_file \
  /path/to/barcode_reads_file

```

**Parameter Definitions:**
> Note: Parameters selected to be consistent with the [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) pipeline for 10X Chromium data per the [STARsolo documentation](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md#running-starsolo-for-10x-chromium-scrna-seq-data) are indicated

- `--runThreadN` - indicates the number of threads to be used for STAR alignment and should be set to the number of available cores on the server node
- `--genomeDir` - specifies the path to the directory where the STAR reference is stored
- `--soloType CB_UMI_Simple` - Activates the STAR solo algorithm for 10X Chromium data 
- `--clipAdapterType` - specifies the type of trimming to perform, a value of `CellRanger4` is used to match CellRanger >= 4.0 such that the TSO adapter sequence is clipped from the 5' end of the cDNA read and the polyA-tail is trimmed from the 3' end
- `--outFilterScoreMin` - specifies the Q score to use for quality trimming, a value of Q`30` is used to match CellRanger >= 4.0
- `--soloCBmatchWLtype` - cell barcode and UMI collapsing parameter, a value of `1MM_multi_Nbase_pseudocounts` is used to get the best agreement between STARsolo and CellRanger >= 3.0
- `--soloUMIfiltering` - cell barcode and UMI collapsing parameter, a value of `MultiGeneUMI_CR` is used to get the best agreement between STARsolo and CellRanger >= 3.0 
- `--soloUMIdedup` - cell barcode and UMI collapsing parameter, a value of `1MM_CR` is used to get the best agreement between STARsolo and CellRanger >= 3.0
- `--soloUMIlen` - barcode length, a value of `12` is used to work for 10X Chromium V3 data
- `--soloCellFilter` - specifies the type of cell filtering to perform (ExpectedCells = number of expected cells), a value of `EmptyDrops_CR <ExpectedCells> 0.99 10 45000 90000 500 0.01 20000 0.01 10000` instructs STAR to use the CellRanger 3.0.0 advanced filtering based on the [EmptyDrop algorithm](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1662-y) 
- `--soloMultiMappers` - specifies the algorithm used to handle multi-mapped reads, a value of `EM` instructs STAR to use the Maximum Likelihood Estimation (MLE) to distribute multi-gene UMIs among their genes, taking into account other UMIs (both unique- and multi-gene) from the same cell (i.e. with the same CB).
- `--outSAMattributes` - list of desired sam attributes in the order desired for the output sam file; sam attribute descriptions can be found [here](https://samtools.github.io/hts-specs/SAMtags.pdf) 
- `--outSAMtype` - specifies desired output format, the `BAM Unsorted` options specify that the output file will be not be sorted and be in the bam format (required to output BAM tags in the BAM file) 
- `--soloFeatures` - specifies the genomic features collected for the UMI counts per Cell Barcode, the `Gene GeneFull SJ Velocyto` options indicate the following:
    - `Gene` - output gene counts 
    - `GeneFull` - output pre-mRNA counts (useful for single-nucleus RNA-seq), this option counts all reads that overlap gene loci, i.e. including both exonic and intronic reads
    - `SJ` - output counts for annotated and novel splice junctions
    - `Gene Velocyto` - output spliced, unspliced, and ambiguous counts per cell per gene, similar to the [velocyto.py](http://velocyto.org/) tool
- `--readFilesCommand` - specifies command needed to interpret input files, the `zcat` option indicates input files are compressed with gzip and zcat will be used to uncompress the gzipped input files
- `--soloCBwhitelist` - specifies the CellBarcode whitelist file, indicated here as `CellBarcodeWhitelist`, which contains the list of all known barcode sequences that have been included in the assay kit and are available during library preparation
- `--outFileNamePrefix` - specifies the path to and prefix for the output file names; for GeneLab the prefix is the sample id
- `--readFilesIn` - paths to the input reads files, the first file should contain the cDNA reads and second file should contain the respective barcode (cell+UMI) reads 


**Input Data:**

- STAR genome reference (output from [Step 2](#2-build-star-reference))
- *fastq.gz (raw reads)

**Output Data:**

- **\*Aligned.bam** (binary sequence alignment map with reads mapping to the genome)
- **\*Log.final.out** (log file containing alignment info/stats such as reads mapped, etc)
- *Log.out (main log file containing detailed info about the STAR run)
- *Log.progress.out (minute-by-minute report containing job progress statistics, such as the number of processed reads, %
of mapped reads etc.)
- **\*SJ.out.tab** (high confidence collapsed splice junctions in tab-delimited format)
- /*Solo.out (directory containing the following:)
  - **Barcodes.stats** (barcode statistics)
  - **/Gene** (sub-directory containing the following outputs using the `Gene` soloFeatures setting:)
    > Note: All files in the `/Gene` output directory are published in the [GLDS repository](https://genelab-data.ndc.nasa.gov/genelab/projects) as a `Gene.zip` file for each respective GLDS dataset.
    - Features.stats (quantitated features statistics)
    - Summary.csv (table containing summary statistics for filtered cells)
    - UMIperCellSorted.txt (list of the number of UMIs per cell sorted)
    - /raw (sub-directory containing output files for raw, unfiltered, expression data)
      - barcodes.tsv (table containing all barcodes processed)
      - features.tsv (table containing all gene IDs and symbols)
      - matrix.mtx (table containing the unique-gene UMI raw counts)
      - UniqueAndMult-EM.mtx (table containing the sum of unique+multi-gene UMI counts)
    - /filtered (sub-directory containing output files for filtered expression data)
      - barcodes.tsv (table containing remaining barcodes after filtering)
      - features.tsv (table containing remaining gene IDs and symbols after filtering)
      - matrix.mtx (table containing the unique-gene UMI filtered counts)
  - **/GeneFull** (sub-directory containing the following outputs using the `GeneFull` soloFeatures setting:)
     > Note: All files in the `/GeneFull` output directory are published in the [GLDS repository](https://genelab-data.ndc.nasa.gov/genelab/projects) as a `GeneFull.zip` file for each respective GLDS dataset.
    - Features.stats (quantitated features statistics)
    - Summary.csv (table containing summary statistics for filtered cells)
    - UMIperCellSorted.txt (list of the number of UMIs per cell sorted)
    - /raw (sub-directory containing output files for raw, unfiltered, expression data)
      - barcodes.tsv (table containing all barcodes processed)
      - features.tsv (table containing all gene IDs and symbols)
      - matrix.mtx (table containing the unique-gene UMI raw counts)
      - UniqueAndMult-EM.mtx (table containing the sum of unique+multi-gene UMI counts)
    - /filtered (sub-directory containing output files for filtered expression data)
      - barcodes.tsv (table containing remaining barcodes after filtering)
      - features.tsv (table containing remaining gene IDs and symbols after filtering)
      - matrix.mtx (table containing the unique-gene UMI filtered counts)
  - **/SJ** (sub-directory containing the following outputs using the `SJ` soloFeatures setting:)
     > Note: All files in the `/SJ` output directory are published in the [GLDS repository](https://genelab-data.ndc.nasa.gov/genelab/projects) as a `SJ.zip` file for each respective GLDS dataset.
    - Features.stats (quantitated features statistics)
    - Summary.csv (table containing summary statistics for filtered cells)
    - /raw (sub-directory containing output files for raw, unfiltered, expression data)
      - barcodes.tsv (table containing all barcodes processed)
      - features.tsv (table containing all gene IDs and symbols)
      - matrix.mtx (table containing the raw counts for annotated and novel splice junctions)
  - **/Velocyto** (sub-directory containing the following:)
    > Note: All files in the `/Velocyto` output directory are published in the [GLDS repository](https://genelab-data.ndc.nasa.gov/genelab/projects) as a `Velocyto.zip` file for each respective GLDS dataset.
    - Features.stats (quantitated features statistics)
    - Summary.csv (table containing summary statistics for filtered cells)
    - /raw (sub-directory containing output files for raw, unfiltered, expression data)
      - barcodes.tsv (table containing all barcodes processed)
      - features.tsv (table containing all gene IDs and symbols)
      - ambiguous.mtx (table containing the raw ambiguous counts per cell per gene)
      - spliced.mtx (table containing the raw spliced counts per cell per gene)
      - unspliced.mtx (table containing the raw unspliced counts per cell per gene)

<br>

### 3b. Compile Alignment Logs

```bash
multiqc --interactive \  
  -n align_multiqc \  
  -o /path/to/aligned_multiqc/output/directory \   
  /path/to/*Log.final.out/files
```

**Parameter Definitions:**

- `--interactive` - force reports to use interactive plots
- `-n` - prefix name for output files
- `-o` – the output directory to store results
- `/path/to/*Log.final.out/files` – the directory holding the *Log.final.out output files from the [STAR alignment step](#3a-align-reads-to-reference-genome-with-starsolo), provided as a positional argument

**Input Data:**

- *Log.final.out (log file conting alignment info/stats such as reads mapped, etc., output from [Step 3a](#3a-align-reads-to-reference-genome-with-starsolo))

**Output Data:**

- **align_multiqc.html** (multiqc report)
- **/align_multiqc_data** (directory containing multiqc data)

   
