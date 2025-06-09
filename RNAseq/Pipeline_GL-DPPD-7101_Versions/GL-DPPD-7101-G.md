# GeneLab bioinformatics processing pipeline for Illumina RNA-sequencing data from Eukaryotic organisms

> **This page holds an overview and instructions for how GeneLab processes RNAseq datasets derived from Eukaryotic organisms. Exact processing commands, GL-DPPD-7101 version used, and processed data output files for specific datasets are provided in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/).**  

---

**Date:** February 20, 2025  
**Revision:** G  
**Document Number:** GL-DPPD-7101-G  

**Submitted by:**  
Alexis Torres (GeneLab Data Processing Team)  
Crystal Han (GeneLab Data Processing Team)

**Approved by:**  
Barbara Novak (GeneLab Data Processing Lead)  
Amanda Saravia-Butler (GeneLab Science Lead)  
Samrawit Gebre (OSDR Project Manager)  
Danielle Lopez (OSDR Deputy Project Manager)  

---

## Updates from previous version  

Updated [Ensembl Reference Files](../../GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv) to the following releases:
- Animals: Ensembl release 112
- Plants: Ensembl plants release 59
- Bacteria: Ensembl bacteria release 59

Software Updates:

| Program | Previous Version | New Version    |
|:--------|:-----------------|:---------------|
| FastQC            | 0.11.9        | 0.12.1  |
| MultiQC           | 1.12          | 1.26    |
| Cutadapt          | 3.7           | 4.2     |
| TrimGalore!       | 0.6.7         | 0.6.10  |
| STAR              | 2.7.10a       | 2.7.11b |
| RSEM              | 1.3.1         | 1.3.3   |
| Samtools          | 1.15          | 1.21    |
| gtfToGenePred     | 377           | 469     |
| genePredToBed     | 377           | 469     |
| infer_experiment  | 4.0.0         | 5.0.4   |
| geneBody_coverage | 4.0.0         | 5.0.4   |
| inner_distance    | 4.0.0         | 5.0.4   |
| read_distribution | 4.0.0         | 5.0.4   |
| R                 | 4.1.3         | 4.4.2   |
| Bioconductor      | 3.14.0        | 3.20    |
| DESeq2            | 1.34          | 1.46.0  |
| tximport          | 1.27.1        | 1.34.0  |
| tidyverse         | 1.3.1         | 2.0.0   |
| stringr           | 1.4.1         | 1.5.1   |
| dp_tools          | 1.18*, 1.3.4* | 1.3.8   |
| pandas            | 1.5.0         | 2.2.3   |
| seaborn           | 0.12.0        | 0.13.2  |
| matplotlib        | 3.6.0         | 3.10.0  |
| numpy             | 1.23.3        | 2.2.1   |
| scipy             | 1.9.1         | 1.15.1  |

STAR Alignment
- Added unaligned reads FASTQ output file(s) via STAR `-outReadsUnmapped Fastx`

RSeQC Analysis
- Updated inner_distance.py invocation to use a lower minimum value to account for longer read lengths
  - Previously used fixed -150 minimum value 
  - Now uses -(max read length)

DESeq2 Analysis Workflow  
- Added variance-stabilizing transformation (VST) transformed counts output file, VST_Counts_GLbulkRNAseq.csv.
  
- Account for technical replicates
  - For datasets with uniform technical replicates (all samples have the same number of technical replicates):
    - Sum counts across technical replicates using DESeq2's collapseReplicates
  - For datasets with non-uniform technical replicates:
    - Keep only the first N technical replicates for each sample, where N is the smallest number of technical replicates among all samples
    - Sum counts across the kept technical replicates
      
- Removed DGE and PCA output tables previously used for GeneLab visualization (visualization_output_table_GLbulkRNAseq.csv and visualization_PCA_table_GLbulkRNAseq.csv)
  
- ERCC-normalized DGE analysis was removed. The following output files were removed:
  - ERCC_Normalized_Counts_GLbulkRNAseq.csv
  - ERCCnorm_differential_expression_GLbulkRNAseq.csv
  - ERCCnorm_contrasts_GLbulkRNAseq.csv
  - visualization_output_table_ERCCnorm_GLbulkRNAseq.csv
  - visualization_PCA_table_ERCCnorm_GLbulkRNAseq.csv

- Added parallel rRNA-removed DGE analysis:
  - Create filtered RSEM count files with rRNA features removed:
    - {sample}_rRNArm.genes.results
  - Normalize rRNA-removed counts
  - Perform DGE analysis using rRNA-removed counts
  - Output additional set of rRNA-removed counts and DGE results
  
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
    - [8d. Remove rRNA Genes from RSEM Genes Results](#8d-remove-rrna-genes-from-rsem-genes-results)
      - [8di. Extract rRNA Gene IDs from GTF](#8di-extract-rrna-gene-ids-from-gtf)
      - [8dii. Filter rRNA Genes from RSEM Genes Results](#8dii-filter-rrna-genes-from-rsem-genes-results)
  - [**9. Normalize Read Counts and Perform Differential Gene Expression Analysis**](#9-normalize-read-counts-and-perform-differential-gene-expression-analysis)
    - [9a. Create Sample RunSheet](#9a-create-sample-runsheet)
    - [9b. Environment Set Up](#9b-environment-set-up)
    - [9c. Configure Metadata, Sample Grouping, and Group Comparisons](#9c-configure-metadata-sample-grouping-and-group-comparisons)
    - [9d. Import RSEM GeneCounts](#9d-import-rsem-genecounts)
    - [9e. Perform DGE Analysis](#9e-perform-dge-analysis)
    - [9f. Add Statistics and Gene Annotations to DGE Results](#9f-add-statistics-and-gene-annotations-to-dge-results)
    - [9g. Export DGE Tables](#9g-export-dge-tables)

  - [**10. Evaluate ERCC Spike-In Data**](#10-evaluate-ercc-spike-in-data)
    - [10a. Evaluate ERCC Count Data in Python](#10a-evaluate-ercc-count-data-in-python)
    - [10b. Perform DESeq2 Analysis of ERCC Counts in R](#10b-perform-deseq2-analysis-of-ercc-counts-in-r)
    - [10c. Analyze ERCC DESeq2 Results in Python](#10c-analyze-ercc-deseq2-results-in-python)

---

# Software used  

|Program|Version|Relevant Links|
|:------|:------:|:-------------|
|FastQC|0.12.1|[https://www.bioinformatics.babraham.ac.uk/projects/fastqc/](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)|
|MultiQC|1.26|[https://multiqc.info/](https://multiqc.info/)|
|Cutadapt|4.2|[https://cutadapt.readthedocs.io/en/stable/](https://cutadapt.readthedocs.io/en/stable/)|
|TrimGalore!|0.6.10|[https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)|
|STAR|2.7.11b|[https://github.com/alexdobin/STAR](https://github.com/alexdobin/STAR)|
|RSEM|1.3.3|[https://github.com/deweylab/RSEM](https://github.com/deweylab/RSEM)|
|Samtools|1.21|[http://www.htslib.org/](http://www.htslib.org/)|
|gtfToGenePred|469|[http://hgdownload.cse.ucsc.edu/admin/exe/](http://hgdownload.cse.ucsc.edu/admin/exe/)|
|genePredToBed|469|[http://hgdownload.cse.ucsc.edu/admin/exe/](http://hgdownload.cse.ucsc.edu/admin/exe/)|
|infer_experiment|5.0.4|[http://rseqc.sourceforge.net/#infer-experiment-py](http://rseqc.sourceforge.net/#infer-experiment-py)|
|geneBody_coverage|5.0.4|[http://rseqc.sourceforge.net/#genebody-coverage-py](http://rseqc.sourceforge.net/#genebody-coverage-py)|
|inner_distance|5.0.4|[http://rseqc.sourceforge.net/#inner-distance-py](http://rseqc.sourceforge.net/#inner-distance-py)|
|read_distribution|5.0.4|[http://rseqc.sourceforge.net/#read-distribution-py](http://rseqc.sourceforge.net/#read-distribution-py)|
|R|4.4.2|[https://www.r-project.org/](https://www.r-project.org/)|
|Bioconductor|3.20|[https://bioconductor.org](https://bioconductor.org)|
|BiocParallel|1.40.0|[https://bioconductor.org/packages/release/bioc/html/BiocParallel.html](https://bioconductor.org/packages/release/bioc/html/BiocParallel.html)|
|DESeq2|1.46.0|[https://bioconductor.org/packages/release/bioc/html/DESeq2.html](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)|
|tximport|1.34.0|[https://github.com/mikelove/tximport](https://github.com/mikelove/tximport)|
|tidyverse|2.0.0|[https://www.tidyverse.org](https://www.tidyverse.org)|
|dplyr|1.1.4|[https://dplyr.tidyverse.org/](https://dplyr.tidyverse.org/)|
|knitr|1.49|[https://yihui.org/knitr/](https://yihui.org/knitr/)|
|stringr|1.5.1|[https://github.com/tidyverse/stringr](https://github.com/tidyverse/stringr)|
|yaml|2.3.10|[https://github.com/yaml/yaml](https://github.com/yaml/yaml)|
|dp_tools|1.3.8|[https://github.com/torres-alexis/dp_tools](https://github.com/torres-alexis/dp_tools)|
|pandas|2.2.3|[https://github.com/pandas-dev/pandas](https://github.com/pandas-dev/pandas)|
|seaborn|0.13.2|[https://seaborn.pydata.org/](https://seaborn.pydata.org/)|
|matplotlib|3.10.0|[https://matplotlib.org/stable](https://matplotlib.org/stable)|
|numpy|2.2.1|[https://numpy.org/](https://numpy.org/)|
|scipy|1.15.1|[https://scipy.org/](https://scipy.org/)|

---

# General processing overview with example commands  

> Exact processing commands and output files listed in **bold** below are included with each RNAseq processed dataset in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/). 

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

- *fastqc.html (FastQC output html summary)
- *fastqc.zip (FastQC output data)

<br>

### 1b. Compile Raw Data QC  

```bash
multiqc --interactive -n raw_multiqc_GLbulkRNAseq -o /path/to/raw_multiqc/output/directory /path/to/directory/containing/raw_fastqc/files

zip -r raw_multiqc_GLbulkRNAseq_data.zip raw_multiqc_GLbulkRNAseq_data
```

**Parameter Definitions:**

- `--interactive` – force reports to use interactive plots
- `-n` – prefix name for output files
- `-o` – the output directory to store results
- `/path/to/directory/containing/raw_fastqc/files` – the directory holding the output data from the FastQC run, provided as a positional argument

**Input Data:**

- *fastqc.zip (FastQC data, output from [Step 1a](#1a-raw-data-qc))

**Output Data:**

* **raw_multiqc_GLbulkRNAseq.html** (MultiQC output html summary)
* **raw_multiqc_GLbulkRNAseq_data.zip** (zipped directory containing MultiQC output data)

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
  --output_dir /path/to/TrimGalore/output/directory \
  --paired \ # only for PE studies, remove this parameter if raw data are SE
  sample1_R1_raw.fastq.gz sample1_R2_raw.fastq.gz sample2_R1_raw.fastq.gz sample2_R2_raw.fastq.gz
# if SE, replace the last line with only the forward reads (R1) of each sample

```

**Parameter Definitions:**

- `--gzip` – compress the output files with `gzip`
- `--path_to_cutadapt` – specify path to cutadapt software if it is not in your `$PATH`
- `--cores` – specify the number of threads available on the server node to perform trimming
- `--phred33` – instructs cutadapt to use ASCII+33 quality scores as Phred scores for quality trimming
- `--output_dir` – the output directory to store results
- `--paired` – indicates paired-end reads - both reads, forward (R1) and reverse (R2) must pass length threshold or else both reads are removed
- `sample1_R1_raw.fastq.gz sample1_R2_raw.fastq.gz sample2_R1_raw.fastq.gz sample2_R2_raw.fastq.gz` – the input reads are specified as a positional argument, paired-end read files are listed pairwise such that the forward reads (*R1_raw.fastq.gz) are immediately followed by the respective reverse reads (*R2_raw.fastq.gz) for each sample

**Input Data:**

- *fastq.gz (raw reads)

**Output Data:**

- **\*fastq.gz** (trimmed reads)
- **\*trimming_report.txt** (trimming report)

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

- *fastqc.html (FastQC output html summary)
- *fastqc.zip (FastQC output data)

<br>

### 2c. Compile Trimmed Data QC  

```bash
multiqc --interactive -n trimmed_multiqc_GLbulkRNAseq -o /path/to/trimmed_multiqc/output/directory /path/to/directory/containing/trimmed_fastqc/files /path/to/directory/containing/trimming_reports

zip -r trimmed_multiqc_GLbulkRNAseq_data.zip trimmed_multiqc_GLbulkRNAseq_data
```

**Parameter Definitions:**

- `--interactive` – force reports to use interactive plots
- `-n` – prefix name for output files
- `-o` – the output directory to store results
- `/path/to/directory/containing/trimmed_fastqc/files` – the directory holding the output data from the fastqc run, provided as a positional argument
- `/path/to/directory/containing/trimming_reports` – the directory containing the trimming reports from the trim/filter step, provided as a positional argument

**Input Data:**

- *fastqc.zip (FastQC data, output from [Step 2b](#2b-trimmed-data-qc))
- *trimming_report.txt (trimming report, output from [Step 2a](#2a-trimfilter-raw-data))

**Output Data:**

* **trimmed_multiqc_GLbulkRNAseq.html** (MultiQC output html summary)
* **trimmed_multiqc_GLbulkRNAseq_data.zip** (zipped directory containing MultiQC output data)

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
- `--runMode` – instructs STAR to run genome indices generation job
- `--limitGenomeGenerateRAM` – maximum RAM available (in bytes) to generate STAR reference, at least 35GB are needed for mouse and the example above shows 55GB
- `--genomeSAindexNbases` – length (in bases) of the SA pre-indexing string, usually between 10 and 15. Longer strings require more memory but allow for faster searches. This value should be scaled down for smaller genomes (like bacteria) to min(14, log2(GenomeLength)/2 - 1). For example, for a 1 megaBase genome this value would be 9.
- `--genomeDir` – specifies the path to the directory where the STAR reference will be stored. At least 100GB of available disk space is required for mammalian genomes.
- `--genomeFastaFiles` – specifies one or more fasta file(s) containing the genome reference sequences
- `--sjdbGTFfile` – specifies the file(s) containing annotated transcripts in the standard gtf format
- `--sjdbOverhang` – indicates the length of the genomic sequence around the annotated junction to be used in constructing the splice junctions database. The length should be one less than the maximum length of the reads.

**Input Data:**

- *.fasta (genome sequence, this pipeline version uses the Ensembl fasta file indicated in the `fasta` column of the [GL-DPPD-7110-A_annotations.csv](../../GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv) GeneLab Annotations file)
- *.gtf (genome annotation, this pipeline version uses the Ensembl gtf file indicated in the `gtf` column of the [GL-DPPD-7110-A_annotations.csv](../../GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv) GeneLab Annotations file)

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
 --outReadsUnmapped Fastx \
 --genomeLoad NoSharedMemory \
 --readFilesIn /path/to/trimmed_forward_reads \
 /path/to/trimmed_reverse_reads # only needed for PE studies

mv <sample_id>_Unmapped.out.mate1 <sample_id>_R1_unmapped.fastq  # Only needed for PE studies
mv <sample_id>_Unmapped.out.mate2 <sample_id>_R2_unmapped.fastq  # Only needed for PE studies
# mv <sample_id>_Unmapped.out.mate1 <sample_id>_unmapped.fastq    # Only needed for SE studies
gzip *_unmapped.fastq
```

**Parameter Definitions:**

- `--twopassMode` – specifies 2-pass mapping mode; the `Basic` option instructs STAR to perform the 1st pass mapping, then automatically extract junctions, insert them into the genome index, and re-map all reads in the 2nd mapping pass
- `--limitBAMsortRAM` – maximum RAM available (in bytes) to sort the bam files, the example above indicates 65GB
- `--genomeDir` – specifies the path to the directory where the STAR reference is stored
- `--outSAMunmapped` – specifies output of unmapped reads in the sam format; the `Within` option instructs STAR to output the unmapped reads within the main sam file
- `--outFilterType` – specifies the type of filtering; the `BySJout` option instructs STAR to keep only those reads that contain junctions that passed filtering in the SJ.out.tab output file
- `--outSAMattributes` – list of desired sam attributes in the order desired for the output sam file; sam attribute descriptions can be found [here](https://samtools.github.io/hts-specs/SAMtags.pdf)
- `--outFilterMultimapNmax` – specifies the maximum number of loci the read is allowed to map to; all alignments will be output only if the read maps to no more loci than this value
- `--outFilterMismatchNmax` – maximum number of mismatches allowed to be included in the alignment output
- `--outFilterMismatchNoverReadLmax` – ratio of mismatches to read length allowed to be included in the alignment output; the `0.04` value indicates that up to 4 mismatches are allowed per 100 bases
- `--alignIntronMin` – minimum intron size; a genomic gap is considered an intron if its length is equal to or greater than this value, otherwise it is considered a deletion
- `--alignIntronMax` – maximum intron size
- `--alignMatesGapMax` – maximum genomic distance (in bases) between two mates of paired-end reads; this option should be removed for single-end reads
- `--alignSJoverhangMin` – minimum overhang (i.e. block size) for unannotated spliced alignments
- `--alignSJDBoverhangMin` – minimum overhang (i.e. block size) for annotated spliced alignments
- `--sjdbScore` – additional alignment score for alignments that cross database junctions
- `--readFilesCommand` – specifies command needed to interpret input files; the `zcat` option indicates input files are compressed with gzip and zcat will be used to uncompress the gzipped input files
- `--runThreadN` – indicates the number of threads to be used for STAR alignment and should be set to the number of available cores on the server node
- `--outSAMtype` – specifies desired output format; the `BAM SortedByCoordinate` options specify that the output file will be sorted by coordinate and be in the bam format
- `--quantMode` – specifies the type(s) of quantification desired; the `TranscriptomeSAM` option instructs STAR to output a separate sam/bam file containing alignments to the transcriptome and the `GeneCounts` option instructs STAR to output a tab delimited file containing the number of reads per gene
- `--outSAMheaderHD` – indicates a header line for the sam/bam file
- `--outFileNamePrefix` – specifies the path to and prefix for the output file names; for GeneLab the prefix is the sample id
- `outReadsUnmapped` - specifies how to output unmapped and partially mapped reads (where only one mate of a paired-end read is mapped); the `Fastx` option outputs unmapped reads in separate fastq files 
- `--genomeLoad` – controls how the genome index is loaded into memory; `NoSharedMemory` specifies that each job will have its own private copy of the genome rather than using shared memory. This is the only option compatible with `--twopassMode Basic`.
- `--readFilesIn` – path to input read 1 (forward read) and read 2 (reverse read); for paired-end reads, read 1 and read 2 should be separated by a space; for single-end reads only read 1 should be indicated

**Input Data:**

- STAR genome reference (output from [Step 3](#3-build-star-reference))
- *fastq.gz (trimmed reads, output from [Step 2a](#2a-trimfilter-raw-data))

**Output Data:**

- *Aligned.sortedByCoord.out.bam (sorted mapping to genome)
- **\*Aligned.toTranscriptome.out.bam** (sorted mapping to transcriptome)
- **\*Log.final.out** (log file containing alignment info/stats such as reads mapped, etc)
- *ReadsPerGene.out.tab (tab delimitated file containing STAR read counts per gene with 4 columns that correspond to different strandedness options: column 1 = gene ID, column 2 = counts for unstranded RNAseq, column 3 = counts for 1st read strand aligned with RNA, column 4 = counts for 2nd read strand aligned with RNA)
- *Log.out (main log file containing detailed info about the STAR run)
- *Log.progress.out (minute-by-minute report containing job progress statistics, such as the number of processed reads, % of mapped reads etc.)
- **\*SJ.out.tab** (high confidence collapsed splice junctions in tab-delimited format)
- *_STARgenome (directory containing the following:)
  - sjdbInfo.txt
  - sjdbList.out.tab
- *_STARpass1 (directory containing the following:)
  - Log.final.out
  - SJ.out.tab
- *_STARtmp (directory containing the following:)
  - BAMsort (directory containing subdirectories that are empty – this was the location for temp files that were automatically removed after successful completion)
- **\*unmapped.fastq.gz** (unmapped and partially mapped reads)

<br>

### 4b. Compile Alignment Logs

```bash
multiqc --interactive -n align_multiqc_GLbulkRNAseq -o /path/to/align_multiqc/output/directory /path/to/*Log.final.out/files

zip -r align_multiqc_GLbulkRNAseq_data.zip align_multiqc_GLbulkRNAseq_data
```

**Parameter Definitions:**

- `--interactive` – force reports to use interactive plots
- `-n` – prefix name for output files
- `-o` – the output directory to store results
- `/path/to/*Log.final.out/files` – the directory holding the *Log.final.out output files from the [STAR alignment step](#4a-align-reads-to-reference-genome-with-star), provided as a positional argument

**Input Data:**

- *Log.final.out (log file containing alignment info/stats such as reads mapped, etc., output from [Step 4a](#4a-align-reads-to-reference-genome-with-star))

**Output Data:**

* **align_multiqc_GLbulkRNAseq.html** (MultiQC output html summary)
* **align_multiqc_GLbulkRNAseq_data.zip** (zipped directory containing MultiQC output data)

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
write.csv(counts,file='STAR_Unnormalized_Counts_GLbulkRNAseq.csv')


## print session info ##
print("Session Info below: ")
print("")
sessionInfo()
```

**Input Data:**

- samples.txt (a newline delimited list of sample IDs)
- *ReadsPerGene.out.tab (STAR counts per gene, output from [Step 4a](#4a-align-reads-to-reference-genome-with-star))

**Output Data:**

- **STAR_Unnormalized_Counts_GLbulkRNAseq.csv** (table containing raw STAR counts for each sample)

<br>

### 4d. Sort Aligned Reads

```bash
samtools sort -m 3G \
  --threads NumberOfThreads \
  -o /path/to/*Aligned.sortedByCoord_sorted.out.bam \
  /path/to/*Aligned.sortedByCoord.out.bam
```

**Parameter Definitions:**

- `-m` – memory available per thread, `3G` indicates 3 gigabytes, this can be changed based on user resources
- `--threads` – number of threads available on server node to sort genome alignment files
- `/path/to/*Aligned.sortedByCoord.out.bam` – path to the *Aligned.sortedByCoord.out.bam output files from the [STAR alignment step](#4a-align-reads-to-reference-genome-with-star), provided as a positional argument

**Input Data:**

- *Aligned.sortedByCoord.out.bam (sorted mapping to genome file, output from [Step 4a](#4a-align-reads-to-reference-genome-with-star))

**Output Data:**

- **\*Aligned.sortedByCoord_sorted.out.bam** (samtools sorted genome aligned bam file)

<br>

### 4e. Index Sorted Aligned Reads

```bash
samtools index -@ NumberOfThreads /path/to/*Aligned.sortedByCoord_sorted.out.bam
```

**Parameter Definitions:**

- `-@` – number of threads available on server node to index the sorted alignment files
- `/path/to/*Aligned.sortedByCoord_sorted.out.bam` – the path to the sorted *Aligned.sortedByCoord_sorted.out.bam output files from the [step 4d](#4d-sort-aligned-reads), provided as a positional argument

**Input Data:**

- *Aligned.sortedByCoord_sorted.out.bam (sorted mapping to genome file, output from [Step 4d](#4d-sort-aligned-reads))

**Output Data:**

- **\*Aligned.sortedByCoord_sorted.out.bam.bai** (index of sorted mapping to genome file)

<br>

---

## 5. Create Reference BED File

<br>

### 5a. Convert GTF to genePred File  

```bash
gtfToGenePred -geneNameAsName2 \
  /path/to/annotation/gtf/file \
  /path/to/output/genePred/file

```

**Parameter Definitions:**

- `--geneNameAsName2` – specifies that gene_name should be used as the name2 field in the genePred file
- `/path/to/annotation/gtf/file` – specifies the file(s) containing annotated reference transcripts in the standard gtf format, provided as a positional argument
- `/path/to/output/genePred/file` – specifies the location and name of the output genePred file(s), provided as a positional argument

**Input Data:**

- *.gtf (genome annotation, this pipeline version uses the Ensembl gtf file indicated in the `gtf` column of the [GL-DPPD-7110-A_annotations.csv](../../GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv) GeneLab Annotations file)

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
 -s 15000000 > /path/to/*.infer_expt.out
```

**Parameter Definitions:**

- `-r` – specifies the path to the reference annotation BED file
- `-i` – specifies the path to the input bam file(s)
- `-s` – specifies the number of reads to be sampled from the input bam file(s), 15M reads are sampled
- `>` – redirects standard output to specified file
- `/path/to/*.infer_expt.out` – specifies the location and name of the file containing the infer_experiment standard output

**Input Data:**

- *.bed (genome annotation in BED format, output from [Step 5b](#5b-convert-genepred-to-bed-file))
- *Aligned.sortedByCoord_sorted.out.bam (sorted mapping to genome file, output from [Step 4d](#4d-sort-aligned-reads))
- *Aligned.sortedByCoord_sorted.out.bam.bai (index of sorted mapping to genome file, output from [Step 4e](#4e-index-sorted-aligned-reads), although not indicated in the command, this file must be present in the same directory as the respective \*Aligned.sortedByCoord_sorted.out.bam file)

**Output Data:**

- *.infer_expt.out (file containing the infer_experiment standard output)

<br>

### 6b. Compile Strandedness Reports

```bash
multiqc --interactive -n infer_exp_multiqc_GLbulkRNAseq -o /path/to/infer_exp_multiqc/output/directory /path/to/*.infer_expt.out/files

zip -r infer_exp_multiqc_GLbulkRNAseq_data.zip infer_exp_multiqc_GLbulkRNAseq_data
```

**Parameter Definitions:**

- `--interactive` – force reports to use interactive plots
- `-n` – prefix name for output files
- `-o` – the output directory to store results
- `/path/to/*.infer_expt.out/files` – the directory holding the *.infer_expt.out output files from the [read strandedness step](#6a-determine-read-strandedness), provided as a positional argument

**Input Data:**

- *.infer_expt.out (file containing the infer_experiment standard output, output from [Step 6a](#6a-determine-read-strandedness))

**Output Data:**

* **infer_exp_multiqc_GLbulkRNAseq.html** (MultiQC output html summary)
* **infer_exp_multiqc_GLbulkRNAseq_data.zip** (zipped directory containing MultiQC output data)

<br>

### 6c. Evaluate GeneBody Coverage

```bash
geneBody_coverage.py -r /path/to/annotation/BED/file \
 -i /path/to/*Aligned.sortedByCoord_sorted.out.bam \
 -o /path/to/geneBody_coverage/output/directory/<sample_id>
```

**Parameter Definitions:**

- `-r` – specifies the path to the reference annotation BED file
- `-i` – specifies the path to the input bam file(s)
- `-o` – specifies the path to the output directory
- `/path/to/geneBody_coverage/output/directory/<sample_id>` – specifies the location and name of the directory containing the geneBody_coverage output files

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
multiqc --interactive -n genebody_cov_multiqc_GLbulkRNAseq -o /path/to/geneBody_cov_multiqc/output/directory /path/to/geneBody_coverage/output/files

zip -r genebody_cov_multiqc_GLbulkRNAseq_data.zip genebody_cov_multiqc_GLbulkRNAseq_data
```

**Parameter Definitions:**

- `--interactive` – force reports to use interactive plots
- `-n` – prefix name for output files
- `-o` – the output directory to store results
- `/path/to/geneBody_coverage/output/files` – the directory holding the geneBody_coverage output files from [step 6c](#6c-evaluate-genebody-coverage), provided as a positional argument

**Input Data:**

- *.geneBodyCoverage.txt (tab delimited file containing genebody coverage values, output from [Step 6c](#6c-evaluate-genebody-coverage))

**Output Data:**

* **genebody_cov_multiqc_GLbulkRNAseq.html** (MultiQC output html summary)
* **genebody_cov_multiqc_GLbulkRNAseq_data.zip** (zipped directory containing MultiQC output data)

<br>

### 6e. Determine Inner Distance (For Paired End Datasets ONLY)

```bash
inner_distance.py -r /path/to/annotation/BED/file \
 -i /path/to/*Aligned.sortedByCoord_sorted.out.bam \
 -k 15000000 \
 -l -(max read length) \
 -u 350 \ 
 -o  /path/to/inner_distance/output/directory
```

**Parameter Definitions:**

- `-r` – specifies the path to the reference annotation BED file
- `-i` – specifies the path to the input bam file(s)
- `-k` – specifies the number of reads to be sampled from the input bam file(s), 15M reads are sampled
- `-l` – specifies the lower bound of inner distance (bp), set to negative of the maximum read length
- `-u` – specifies the upper bound of inner distance (bp)
- `/path/to/inner_distance/output/directory` – specifies the location and name of the directory containing the inner_distance output files

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
multiqc --interactive -n inner_dist_multiqc_GLbulkRNAseq -o /path/to/inner_dist_multiqc/output/directory /path/to/inner_dist/output/files

zip -r inner_dist_multiqc_GLbulkRNAseq_data.zip inner_dist_multiqc_GLbulkRNAseq_data
```

**Parameter Definitions:**

- `--interactive` – force reports to use interactive plots
- `-n` – prefix name for output files
- `-o` – the output directory to store results
- `/path/to/inner_dist/output/files` – the directory holding the inner_distance output files from [Step 6e](#6e-determine-inner-distance-for-paired-end-datasets-only), provided as a positional argument

**Input Data:**

- *.inner_distance_freq.txt (tab delimited table of inner distances from [step 6e](#6e-determine-inner-distance-for-paired-end-datasets-only))

**Output Data:**

* **inner_dist_multiqc_GLbulkRNAseq.html** (MultiQC output html summary)
* **inner_dist_multiqc_GLbulkRNAseq_data.zip** (zipped directory containing MultiQC output data)

<br>

### 6g. Assess Read Distribution

```bash
read_distribution.py -r /path/to/annotation/BED/file \
 -i /path/to/*Aligned.sortedByCoord_sorted.out.bam > /path/to/*.read_dist.out
```

**Parameter Definitions:**

- `-r` – specifies the path to the reference annotation BED file
- `-i` – specifies the path to the input bam file(s)
- `>` – redirects standard output to specified file
- `/path/to/*.read_dist.out` – specifies the location and name of the file containing the read_distribution standard output

**Input Data:**

- *.bed (genome annotation in BED format, output from [Step 5b](#5b-convert-genepred-to-bed-file))
- *Aligned.sortedByCoord_sorted.out.bam (sorted mapping to genome file, output from [Step 4d](#4d-sort-aligned-reads))
- *Aligned.sortedByCoord_sorted.out.bam.bai (index of sorted mapping to genome file, output from [Step 4e](#4e-index-sorted-aligned-reads), although not indicated in the command, this file must be present in the same directory as the respective \*Aligned.sortedByCoord_sorted.out.bam file)

**Output Data:**

- *.read_dist.out (file containing the read distribution standard output)

<br>

### 6h. Compile Read Distribution Reports

```bash
multiqc --interactive -n read_dist_multiqc_GLbulkRNAseq -o /path/to/read_dist_multiqc/output/directory /path/to/*.read_dist.out/files

zip -r read_dist_multiqc_GLbulkRNAseq_data.zip read_dist_multiqc_GLbulkRNAseq_data
```

**Parameter Definitions:**

- `--interactive` - force reports to use interactive plots
- `-n` - prefix name for output files
- `-o` – the output directory to store results
- `/path/to/*.read_dist.out/files` – the directory holding the *.read_dist.out output files from [Step 6g](#6g-assess-read-distribution) provided as a positional argument

**Input Data:**

- *.read_dist.out (files containing the read_distribution standard output, output from [Step 6g](#6g-assess-read-distribution))

**Output Data:**

* **read_dist_multiqc_GLbulkRNAseq.html** (MultiQC output html summary)
* **read_dist_multiqc_GLbulkRNAseq_data.zip** (zipped directory containing MultiQC output data)

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
- `/path/to/RSEM/genome/directory/RSEM_ref_prefix` – specifies the path to the directory where the RSEM reference will be stored and the prefix desired for the RSEM reference files, provided as a positional argument

**Input Data:**

- *.fasta (genome sequence, this pipeline version uses the Ensembl fasta file indicated in the `fasta` column of the [GL-DPPD-7110-A_annotations.csv](../../GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv) GeneLab Annotations file)
- *.gtf (genome annotation, this pipeline version uses the Ensembl gtf file indicated in the `gtf` column of the [GL-DPPD-7110-A_annotations.csv](../../GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv) GeneLab Annotations file)

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
- `--alignments` – indicates that the input file contains alignments in sam, bam, or cram format
- `--bam` – specifies that the input alignments are in bam format
- `--paired-end` – indicates that the input reads are paired-end reads; this option should be removed if the input reads are single-end
- `--seed` – the seed for the random number generators used in calculating posterior mean estimates and credibility intervals; must be a non-negative 32-bit integer
- `--seed-length 20` – instructs RSEM to ignore any aligned read if it or its mates' (for paired-end reads) length is less than 20bp
- `--estimate-rspd` – instructs RSEM to estimate the read start position distribution (rspd) from the data
- `--no-bam-output` – instructs RSEM not to output any bam file
- `--strandedness` – defines the strandedness of the RNAseq reads; the `reverse` option is used if read strandedness (output from [step 6](#6a-determine-read-strandedness)) is antisense, `forward` is used with sense strandedness, and `none` is used if strandedness is half sense half antisense
- `/path/to/*Aligned.toTranscriptome.out.bam` – specifies path to input bam files, provided as a positional argument
- `/path/to/RSEM/genome/directory/RSEM_ref_prefix` – specifies the path to the directory where the RSEM reference is stored and its prefix, provided as a positional argument
- `/path/to/RSEM/counts/output/directory` – specifies the path to and prefix for the output file names; for GeneLab the prefix is the sample id

**Input Data:**

- RSEM genome reference (output from [Step 7](#7-build-rsem-reference))
- *Aligned.toTranscriptome.out.bam (sorted mapping to transcriptome, output from [Step 4a](#4a-align-reads-to-reference-genome-with-star))

**Output Data:**

- **\*genes.results** (counts per gene)
- **\*isoforms.results** (counts per isoform)
- *stat (directory containing the following stats files)
  - *cnt
  - *model
  - *theta

<br>

### 8b. Compile RSEM Count Logs

```bash
multiqc --interactive -n RSEM_count_multiqc_GLbulkRNAseq -o /path/to/RSEM_count_multiqc/output/directory /path/to/*stat/files

zip -r RSEM_count_multiqc_GLbulkRNAseq_data.zip RSEM_count_multiqc_GLbulkRNAseq_data
```

**Parameter Definitions:**

- `--interactive` – force reports to use interactive plots
- `-n` – prefix name for output files
- `-o` – the output directory to store results
- `/path/to/*stat/files` – the directories holding the *cnt file in the *stat output folder from the [RSEM Counts step](#8a-count-aligned-reads-with-rsem), provided as a positional argument

**Input Data:**

- *stat (directory containing the following stats files, output from [Step 8a](#8a-count-aligned-reads-with-rsem))
  - *cnt

**Output Data:**

* **RSEM_count_multiqc_GLbulkRNAseq.html** (MultiQC output html summary)
* **RSEM_count_multiqc_GLbulkRNAseq_data.zip** (zipped directory containing MultiQC output data)

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
write.csv(NumNonZeroGenes,file='NumNonZeroGenes_GLbulkRNAseq.csv')

## print session info ##
print("Session Info below: ")
print("")
sessionInfo()
```

**Input Data:**

- samples.txt (A newline delimited list of sample IDs)
- *genes.results (RSEM counts per gene, output from [Step 8a](#8a-count-aligned-reads-with-rsem))

**Output Data:**

- NumNonZeroGenes_GLbulkRNAseq.csv (A samplewise table of the number of genes expressed)

<br>

### 8d. Remove rRNA Genes from RSEM Genes Results  

<br>

#### 8di. Extract rRNA Gene IDs from GTF  

```bash
### Extract unique rRNA ENSEMBL gene IDs from GTF file ###
grep -E 'gene_biotype "rRNA"|gene_type "rRNA"|gbkey "rRNA"' /path/to/annotation/gtf/file \
    | grep -o 'gene_id "[^"]*"' \
    | sed 's/gene_id "\(.*\)"/\1/' \
    | sort -u > rrna_ids.txt
```

**Input Data:**
- *.gtf (genome annotation, this pipeline version uses the Ensembl gtf file indicated in the `gtf` column of the [GL-DPPD-7110-A_annotations.csv](../../GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv) GeneLab Annotations file)

**Output Data:**
- *rrna_ensembl_ids.txt (list of unique rRNA ENSEMBL gene IDs)

<br>

#### 8dii. Filter rRNA Genes from RSEM Genes Results 

```bash
### Filter out rRNA entries ###
awk 'NR==FNR {ids[$1]=1; next} !($1 in ids)' \
    *rrna_ensembl_ids.txt \
    *.genes.results > *_rRNArm.genes.results

### Count removed rRNA entries ###
rRNA_count=$(awk 'NR==FNR {ids[$1]=1; next} $1 in ids' \
    *rrna_ensembl_ids.txt \
    *.genes.results | wc -l)
echo "*: ${rRNA_count} rRNA entries removed." > *_rRNA_counts.txt
```

**Input Data:**
- *genes.results (RSEM counts per gene, output from [Step 8a](#8a-count-aligned-reads-with-rsem))
- *rrna_ensembl_ids.txt (file containing list of gene IDs with rRNA features, output from [Step 8di](#8di-extract-rrna-gene-ids-from-gtf))

**Output Data:**
- **\*rRNArm.genes.results** (RSEM gene counts with rRNA entries removed)
- *rRNA_counts.txt (Summary of number of rRNA entries removed)

<br>

---

## 9. Normalize Read Counts and Perform Differential Gene Expression Analysis

> Note: DGE Analysis is performed twice with different sets of input files:
> 1. Using RSEM genes.results files (*genes.results, output from [Step 8a](#8a-count-aligned-reads-with-rsem))
> 2. Using rRNA-removed RSEM genes.results files (*rRNArm.genes.results, output from [Step 8dii](#8dii-filter-rrna-genes-from-rsem-genes-results))

<br>

### 9a. Create Sample RunSheet

> Note: Rather than running the command below to create the runsheet needed for processing, the runsheet may also be created manually by following the [file specification](../Workflow_Documentation/NF_RCP/examples/runsheet/README.md).

```bash
### Download the *ISA.zip file from the Open Science Data Repository ###

dpt-get-isa-archive \
 --accession GLDS-###

### Parse the metadata from the *ISA.zip file to create a sample runsheet ###

dpt-isa-to-runsheet --accession GLDS-### \
 --config-type bulkRNASeq \
 --config-version Latest \
 --isa-archive *ISA.zip
```

**Parameter Definitions:**

- `--accession GLDS-###` – GLDS accession ID (replace ### with the GLDS number being processed), used to retrieve the urls for the ISA archive and raw reads hosted on the GeneLab Repository
  > *Note: you can also use the OSD identifier, e.g. `--accession OSD-###`* 
- `--config-type` – Instructs the script to extract the metadata required for `bulkRNAseq` processing from the ISA archive
- `--config-version` – Specifies the `dp-tools` configuration version to use, a value of `Latest` will specify the most recent version
- `--isa-archive` – Specifies the *ISA.zip file for the respective GLDS dataset, downloaded in the `dpt-get-isa-archive` command


**Input Data:**

- No input data required but the GLDS (or OSD) accession ID needs to be indicated, which is used to download the respective ISA archive 

**Output Data:**

- *ISA.zip (compressed ISA directory containing Investigation, Study, and Assay (ISA) metadata files for the respective GLDS dataset, used to define sample groups - the *ISA.zip file is located in the [OSDR repository](https://osdr.nasa.gov/bio/repo/) under 'Files' -> 'Study Metadata Files')

- **{GLDS-Accession-ID}_bulkRNASeq_v{version}_runsheet.csv** (table containing metadata required for processing, version denotes the dp_tools schema used to specify the metadata to extract from the ISA archive)

<br>

### 9b. Environment Set Up

```R
### Install and load required packages ###

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# List of required packages
cran_packages <- c("stringr", "knitr", "yaml", "dplyr")
bioc_packages <- c("tximport", "DESeq2", "BiocParallel")

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
library(stringr)
library(knitr)
library(yaml)
library(dplyr)
library(tximport)
library(DESeq2)
library(BiocParallel)

### Define a system-specific BiocParallelParam object to run DGE R script functions in parallel when applicable ###
BPPARAM <- SerialParam(RNGseed = 7)

### Set random seed for reproducibility ###
set.seed(7)

### Define which organism is used in the study - this should be consistent with the species name in the "species" column of the GL-DPPD-7110-A_annotations.csv file ###
organism <- "organism_that_samples_were_derived_from"

### Define the location of the input data and where the output data will be printed to ###
runsheet_path="/path/to/directory/containing/runsheet.csv/file" ## This is the runsheet created in Step 9a above
work_dir="/path/to/working/directory/where/script/is/executed/from" 
input_counts="/path/to/directory/containing/RSEM/counts/files"
norm_output="/path/to/normalized/counts/output/directory"
DGE_output="/path/to/DGE/output/directory"

### Pull in the GeneLab annotation table (GL-DPPD-7110-A_annotations.csv) file ###
org_table_link <- "https://raw.githubusercontent.com/nasa/GeneLab_Data_Processing/master/GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv"

org_table <- read.table(org_table_link, sep = ",", header = TRUE)

### Define the link to the GeneLab annotation table for the organism of interest ###
annotations_link <- org_table[org_table$species == organism, "genelab_annots_link"]

### Set your working directory to the directory where you will execute your DESeq2 script from ###
setwd(file.path(work_dir))
```

**Input Data:**

* {GLDS-Accession-ID}_bulkRNASeq_v{version}_runsheet.csv (runsheet, output from [Step 9a](#9a-create-sample-runsheet))
* `organism` (name of organism samples were derived from, found in the species column of [GL-DPPD-7110-A_annotations.csv](../../GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110-A/GL-DPPD-7110-A_annotations.csv) file)
- `BPPARAM` (system-specific BiocParallelParam object for parallel processing configuration, by default this is set to `SerialParam(RNGseed = 7)`)

**Output Data:**

* `runsheet_path` (variable containing path to runsheet created in [Step 9a](#9a-create-sample-runsheet)) 
* `annotations_link` (variable containing URL to GeneLab gene annotation table for the organism)
- `BPPARAM` (variable containing BiocParallelParam object for parallel processing configuration)

<br>

### 9c. Configure Metadata, Sample Grouping, and Group Comparisons

```R
### Pull all factors for each sample in the study from the runsheet created in Step 9a ###

compare_csv_from_runsheet <- function(runsheet_path) {
    df <- read.csv(runsheet_path)
    factors <- df %>%
        select(matches("Factor.Value", ignore.case = TRUE)) %>%
        rename_with(~ paste0("factor_", seq_along(.)))
    
    # Check if both Source.Name and Has.Tech.Reps columns exist
    if ("Source.Name" %in% colnames(df) && "Has.Tech.Reps" %in% colnames(df)) {
        result <- df %>%
            select(Sample.Name, Source.Name, Has.Tech.Reps) %>%
            bind_cols(factors)
    } else {
        result <- df %>%
            select(Sample.Name) %>%
            bind_cols(factors)
    }
    
    return(result)
}

### Load metadata from runsheet csv file ###
compare_csv <- compare_csv_from_runsheet(runsheet_path)

### Create data frame containing all samples and respective factors ###
study <- if ("Source.Name" %in% colnames(compare_csv) && "Has.Tech.Reps" %in% colnames(compare_csv)) {
    compare_csv[, -c(1, 2, 3), drop=FALSE]  # Exclude Sample.Name, Source.Name, and Has.Tech.Reps
} else {
    compare_csv[, -1, drop=FALSE]  # Exclude only Sample.Name
}
rownames(study) <- compare_csv$Sample.Name

### Format groups and indicate the group that each sample belongs to ###
group <- if (ncol(study) >= 2) {
    apply(study, 1, paste, collapse = " & ")
} else {
    study[[1]]
}
group_names <- paste0("(", group, ")") ## human readable group names
group <- sub("^BLOCKER_", "", make.names(paste0("BLOCKER_", group))) # group naming compatible with R models, this maintains the default behaviour of make.names with the exception that 'X' is never prepended to group names
names(group) <- group_names
rm(group_names)

### Format contrasts table, defining pairwise comparisons for all groups ###
contrast.names <- combn(levels(factor(names(group))),2)

### Generate matrix of pairwise group combinations for comparison ###
contrasts <- apply(contrast.names, MARGIN=2, function(col) sub("^BLOCKER_", "",  make.names(paste0("BLOCKER_", stringr::str_sub(col, 2, -2))))) # limited make.names call for each group (also removes leading parentheses)
contrast.names <- c(paste(contrast.names[1,],contrast.names[2,],sep = "v"),paste(contrast.names[2,],contrast.names[1,],sep = "v")) ## format combinations for output table files names
contrasts <- cbind(contrasts,contrasts[c(2,1),])
colnames(contrasts) <- contrast.names
rm(contrast.names) 
```

**Input Data:**

* `runsheet_path` (variable containing path to runsheet created in [Step 9a](#9a-create-sample-runsheet))

**Output Data:**

* `compare_csv` (data frame containing sample names, technical replicate information if provided, and factor levels from the runsheet)
* `study` (data frame specifying factor levels assigned to each sample)
* `group` (named vector specifying the group or set of factor levels for each sample)
* `contrasts` (matrix defining pairwise comparisons between groups)

<br>

### 9d. Import RSEM GeneCounts

```R
### Import RSEM gene count data ###
files <- list.files(
    path = input_counts, 
    pattern = ".genes.results", 
    full.names = TRUE
)

### Reorder files to match sample order ###
samples <- rownames(study)
reordering <- sapply(samples, function(x) {
    grep(paste0(x, ".genes.results$"), files, value = FALSE)
})

files <- files[reordering]
names(files) <- samples

### Import RSEM data ###
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)

### Verify sample count matches ###
if (dim(txi.rsem$counts)[2] != nrow(study)) {
    stop("Sample count mismatch between imported gene results and runsheet")
}

### Add 1 to genes with lengths of zero - needed to make DESeqDataSet object ###
txi.rsem$length[txi.rsem$length == 0] <- 1
```

**Input Data:**

* *genes.results (RSEM counts per gene, output from [Step 8a](#8a-count-aligned-reads-with-rsem) or from [Step 8dii](#8dii-filter-rrna-genes-from-rsem-genes-results) when using rRNA-removed count data)
* `study` (data frame containing sample condition values, output from [Step 9c](#9c-create-study-group-and-contrasts))

**Output Data:**

* `txi.rsem` (list containing imported RSEM data including counts, transcript lengths, and gene lengths)
* `samples` (vector of sample names from study data frame)

<br>

### 9e. Perform DGE Analysis

```R
### Create sample table ###
sampleTable <- data.frame(condition=factor(group))
rownames(sampleTable) <- colnames(txi.rsem$counts)

### Handle technical replicates in sample table - STEP 1/2: Filter which samples to retain ###
# Only if both `Source Name` and `Has Tech Reps` columns exist in runsheet
if ("Source.Name" %in% colnames(compare_csv) && "Has.Tech.Reps" %in% colnames(compare_csv)) {
    all_samples <- rownames(sampleTable)
    
    # Get source names and tech rep status for each sample
    sample_info <- data.frame(
        name = all_samples,
        source_name = compare_csv$Source.Name[match(all_samples, compare_csv$Sample.Name)],
        has_tech_reps = compare_csv$Has.Tech.Reps[match(all_samples, compare_csv$Sample.Name)],
        stringsAsFactors = FALSE
    )
    
    # Count samples per source name (all samples regardless of tech rep status)
    source_counts <- table(sample_info$source_name)
    min_samples_per_source <- min(source_counts)
    
    if (min_samples_per_source > 1) {
        # Use collapseReplicates approach
        # Create grouping variable: Source Name for tech reps, Sample Name for non-tech reps
        collapse_groups <- ifelse(
            toupper(sample_info$has_tech_reps) == "TRUE",
            sample_info$source_name,
            sample_info$name
        )
        
        # Get unique groups and count for balancing
        unique_groups <- unique(collapse_groups)
        group_counts <- sapply(unique_groups, function(grp) {
            sum(collapse_groups == grp)
        })
        
        min_group_size <- min(group_counts)
        
        # Keep first min_group_size samples from each group
        samples_to_keep <- character(0)
        for (grp in unique_groups) {
            indices <- which(collapse_groups == grp)
            samples_to_keep <- c(samples_to_keep, sample_info$name[indices[1:min_group_size]])
        }
        
        # Update sample table and counts
        sampleTable <- sampleTable[samples_to_keep, , drop=FALSE]
        
        if (params$microbes) {
            counts <- counts[, samples_to_keep]
        } else {
            txi.rsem$counts <- txi.rsem$counts[, samples_to_keep]
            txi.rsem$abundance <- txi.rsem$abundance[, samples_to_keep]
            txi.rsem$length <- txi.rsem$length[, samples_to_keep]
        }
        
    } else {
        # min_samples_per_source = 1, use manual filtering approach
        # Keep only first sample from each Source Name + Has Tech Reps = TRUE group
        # Keep all samples with Has Tech Reps = FALSE
        
        samples_to_keep <- character(0)
        
        # Group by source name and tech rep status
        for (src in unique(sample_info$source_name)) {
            src_samples <- sample_info[sample_info$source_name == src, ]
            
            # Separate tech reps from non-tech reps
            tech_rep_samples <- src_samples[toupper(src_samples$has_tech_reps) == "TRUE", ]
            non_tech_rep_samples <- src_samples[toupper(src_samples$has_tech_reps) == "FALSE", ]
            
            # Keep only first tech rep sample if any exist
            if (nrow(tech_rep_samples) > 0) {
                samples_to_keep <- c(samples_to_keep, tech_rep_samples$name[1])
            }
            
            # Keep all non-tech rep samples
            if (nrow(non_tech_rep_samples) > 0) {
                samples_to_keep <- c(samples_to_keep, non_tech_rep_samples$name)
            }
        }
        
        # Update sample table and counts
        sampleTable <- sampleTable[samples_to_keep, , drop=FALSE]
        
        if (params$microbes) {
            counts <- counts[, samples_to_keep]
        } else {
            txi.rsem$counts <- txi.rsem$counts[, samples_to_keep]
            txi.rsem$abundance <- txi.rsem$abundance[, samples_to_keep]
            txi.rsem$length <- txi.rsem$length[, samples_to_keep]
        }
    }
}

### Build dds object ###
dds <- DESeqDataSetFromTximport(
    txi = txi.rsem,
    colData = sampleTable,
    design = ~condition
)

### Handle technical replicates - STEP 2/2: Collapse retained tech reps in DESeq2 object ###
# Only if both `Source Name` and `Has Tech Reps` columns exist in runsheet
if ("Source.Name" %in% colnames(compare_csv) && "Has.Tech.Reps" %in% colnames(compare_csv)) {
    # Get info for remaining samples after filtering
    remaining_samples <- rownames(sampleTable)
    remaining_info <- data.frame(
        name = remaining_samples,
        source_name = compare_csv$Source.Name[match(remaining_samples, compare_csv$Sample.Name)],
        has_tech_reps = compare_csv$Has.Tech.Reps[match(remaining_samples, compare_csv$Sample.Name)],
        stringsAsFactors = FALSE
    )
    
    # Create collapse grouping: Source Name for tech reps, Sample Name for non-tech reps
    collapse_source_names <- ifelse(
        toupper(remaining_info$has_tech_reps) == "TRUE",
        remaining_info$source_name,
        remaining_info$name
    )
    
    # Only collapse if there are multiple samples with the same collapse group
    if (length(unique(collapse_source_names)) < length(collapse_source_names)) {
        # Collapse only if >1 replicate per group exists
        dds <- collapseReplicates(dds, groupby = collapse_source_names)
        
        collapsed_names <- unique(collapse_source_names)
        # Update sampleTable to match collapsed samples
        # For collapsed tech reps, use the source name; for non-tech reps, use sample name
        sampleTable <- sampleTable[match(collapsed_names, collapse_source_names), , drop = FALSE]
        rownames(sampleTable) <- collapsed_names
    }
}

### Filter low count genes ###
keep <- rowSums(counts(dds)) > 10
print(sprintf("Removed %d genes with dataset-wide count sum less than 10", sum(!keep)))
dds <- dds[keep,]

### Remove ERCC spike-in genes if present ###
if (length(grep("ERCC-", rownames(dds))) != 0) {
    dds <- dds[-c(grep("ERCC-", rownames(dds))), ]
}

### Perform DESeq analysis ###
dds <- DESeq(dds, parallel = TRUE, BPPARAM = BPPARAM)

### Generate normalized counts ###
normCounts <- as.data.frame(counts(dds, normalized = TRUE))
VSTCounts <- tryCatch({
  as.data.frame(assay(vst(dds)))
}, error = function(e) {
  # If vst() fails, use varianceStabilizingTransformation directly
  print(sprintf("DEBUG: %s: VST failed, falling back to direct varianceStabilizingTransformation call", Sys.time()))
  as.data.frame(assay(varianceStabilizingTransformation(dds)))
})

### Add 1 to normalized counts for log calculations ###
normCounts <- normCounts + 1

### Calculate LRT statistics ###
dds_lrt <- DESeq(dds, test = "LRT", reduced = ~1)
res_lrt <- results(dds_lrt)

```

**Input Data:**

* `group` (named vector specifying the group or set of factor levels for each sample, output from [Step 9c](#9c-configure-metadata-sample-grouping-and-group-comparisons))
* `compare_csv` (data frame containing sample names, technical replicates information if provided, and factor levels from the runsheet, output from [Step 9c](#9c-configure-metadata-sample-grouping-and-group-comparisons))
* `txi.rsem` (imported RSEM data containing counts matrix, output from [Step 9d](#9d-import-rsem-genecounts))
* `BPPARAM` (system-specific BiocParallelParam object for parallel processing configuration, output from [Step 9b](#9b-environment-set-up))

**Output Data:**

* `sampleTable` (data frame mapping samples to groups)
* `dds` (DESeq2 data object containing normalized counts, experimental design, and differential expression results)
* `normCounts` (data frame of normalized count values + 1)
* `VSTCounts` (data frame of variance stabilized transformed counts)
* `dds_lrt` (DESeq2 data object from likelihood ratio test)
* `res_lrt` (results object from likelihood ratio test)
* `output_table` (data frame containing normalized counts, DGE results)

<br>

### 9f. Add Statistics and Gene Annotations to DGE Results

```R
### Initialize output table with normalized counts ###
gene_id_type <- "ENSEMBL"
output_table <- tibble::rownames_to_column(normCounts, var = gene_id_type)

### Iterate through Wald Tests to generate pairwise comparisons of all groups ###
compute_contrast <- function(i) {
    res <- results(
        dds,
        contrast = c("condition", contrasts[1, i], contrasts[2, i]),
        parallel = FALSE  # Disable internal parallelization
    )
    res_df <- as.data.frame(res@listData)[, c(2, 4, 5, 6)]
    colnames(res_df) <- c(
        paste0("Log2fc_", colnames(contrasts)[i]),
        paste0("Stat_", colnames(contrasts)[i]),
        paste0("P.value_", colnames(contrasts)[i]),
        paste0("Adj.p.value_", colnames(contrasts)[i])
    )
    return(res_df)
}

### Use bplapply to compute results in parallel ###
res_list <- bplapply(1:dim(contrasts)[2], compute_contrast, BPPARAM = BPPARAM)

### Combine the list of data frames into a single data frame ###
res_df <- do.call(cbind, res_list)

### Combine with the existing output_table ###
output_table <- cbind(output_table, res_df)

### Add summary statistics ###
output_table$All.mean <- rowMeans(normCounts, na.rm = TRUE)
output_table$All.stdev <- rowSds(as.matrix(normCounts), na.rm = TRUE, useNames = FALSE)
output_table$LRT.p.value <- res_lrt@listData$padj

### Add group-wise statistics ###
tcounts <- as.data.frame(t(normCounts))
tcounts$group <- colData(dds)$condition[match(rownames(tcounts), rownames(colData(dds)))]

### Aggregate group means and standard deviations ###
agg_means <- aggregate(. ~ group, data = tcounts, FUN = mean, na.rm = TRUE)
agg_stdev <- aggregate(. ~ group, data = tcounts, FUN = sd, na.rm = TRUE)

### Save group names ###
group_ids <- agg_means$group

### Remove the 'group' column and transpose to match expected structure ###
group_means <- as.data.frame(t(agg_means[-1]))
group_stdev <- as.data.frame(t(agg_stdev[-1]))

### Get cleaned group names ###
orig_group_names <- group_names[match(group_ids, group)]

#### Interleave means and stdevs for each group ###
group_stats <- data.frame(matrix(ncol = 0, nrow = nrow(group_means)))
for (i in seq_along(group_ids)) {
    # Create column names with cleaned group names
    mean_colname <- paste0("Group.Mean_", orig_group_names[i])
    stdev_colname <- paste0("Group.Stdev_", orig_group_names[i])
    # Add columns
    group_stats[[mean_colname]] <- group_means[,i]
    group_stats[[stdev_colname]] <- group_stdev[,i]
}

### Add computed group means and standard deviations to output_table ###
output_table <- cbind(output_table, group_stats)

### Read in GeneLab annotation table for the organism of interest ###
annot <- read.table(annotations_link, 
    sep = "\t", 
    header = TRUE, 
    quote = "", 
    comment.char = ""
)

### Combine annotations table and the DGE table ###
# If gene ID column is missing from either table, just write the original DGE table
if (!(gene_id_type %in% colnames(annot)) || !(gene_id_type %in% colnames(output_table))) {
  warning(paste("Gene ID column", gene_id_type, "not found in both tables."))
} else {
  ### Combine annotations with data
  output_table <- annot %>%
    merge(output_table,
          by = gene_id_type,
          all.y = TRUE 
    ) %>%
    select(all_of(gene_id_type), everything())  # Make sure main gene ID is first column
}
```

**Input Data:**

- `gene_id_type` (Gene identifier type, e.g. ENSEMBL, used to merge the annotations with the DGE results)
* `normCounts` (data frame of normalized counts, output from [Step 9e](#9e-perform-dge-analysis))
* `res_lrt` (results object from likelihood ratio test, output from [Step 9e](#9e-perform-dge-analysis))
* `contrasts` (matrix defining pairwise comparisons, output from [Step 9c](#9c-configure-metadata-sample-grouping-and-group-comparisons))
* `dds` (DESeq2 data object containing normalized counts, experimental design, and differential expression results, output from [Step 9e](#9e-perform-dge-analysis))
* `annotations_link` (variable containing URL to GeneLab annotation table, output from [Step 9b](#9b-environment-set-up))
- `BPPARAM` (system-specific BiocParallelParam object for parallel processing configuration, output from [Step 9b](#9b-environment-set-up))

**Output Data:**

* `output_table` (data frame containing the following columns:
  - Gene identifier column (ENSEMBL or TAIR for plant studies)
  - Additional organism-specific gene annotations columns
  - Normalized counts for each sample
  - For each pairwise comparison:
    - Log2 fold change
    - Test statistic
    - P-value
    - Adjusted p-value
  - All.mean (mean across all samples)
  - All.stdev (standard deviation across all samples) 
  - LRT.p.value (likelihood ratio test adjusted p-value)
  - For each experimental group:
    - Group.Mean_(group) (mean within group)
    - Group.Stdev_(group) (standard deviation within group))

<br>

### 9g. Export DGE Tables

```R
### Export unnormalized and normalized counts tables ###
write.csv(txi.rsem$counts, 
    file.path(norm_output, "RSEM_Unnormalized_Counts_GLbulkRNAseq.csv"))

write.csv(normCounts,
    file.path(norm_output, "Normalized_Counts_GLbulkRNAseq.csv"))

write.csv(VSTCounts,
    file.path(norm_output, "VST_Counts_GLbulkRNAseq.csv"))

### Export sample grouping and contrasts tables ###
write.csv(sampleTable,
    file.path(DGE_output, "SampleTable_GLbulkRNAseq.csv"))

write.csv(contrasts,
    file.path(DGE_output, "contrasts_GLbulkRNAseq.csv"))

### Export DGE Results table ###
write.csv(output_table,
    file.path(DGE_output, "differential_expression_GLbulkRNAseq.csv"), 
    row.names = FALSE)

### print session info ###

print("Session Info below: ")
sessionInfo()

```


**Input Data:**

* `contrasts` (matrix defining pairwise comparisons between groups, output from [Step 9c](#9c-configure-metadata-sample-grouping-and-group-comparisons))
* `txi.rsem` (imported RSEM count data, output from [Step 9d](#9d-import-rsem-genecounts))
* `sampleTable` (data frame mapping samples to groups, output from [Step 9e](#9e-perform-dge-analysis))
* `normCounts` (normalized counts, output from [Step 9e](#9e-perform-dge-analysis))
* `VSTCounts` (variance stabilized transformed counts, output from [Step 9e](#9e-perform-dge-analysis)) 
* `output_table` (DGE output table, output from [Step 9f](#9f-add-statistics-and-gene-annotations-to-dge-results))

**Output Data:**

* **RSEM_Unnormalized_Counts_GLbulkRNAseq.csv** (raw RSEM gene counts for all samples, including technical replicates)
* **Normalized_Counts_GLbulkRNAseq.csv** (normalized gene counts)
* **VST_Counts_GLbulkRNAseq.csv** (variance stabilized transformed counts)
* **SampleTable_GLbulkRNAseq.csv** (table specifying the group or set of factor levels for each sample)
* **contrasts_GLbulkRNAseq.csv** (table listing all pairwise group comparisons)
* **differential_expression_GLbulkRNAseq.csv** (DGE results table containing the following columns:
  - Gene identifier column (ENSEMBL or TAIR for plant studies)
  - Additional organism-specific gene annotations columns
  - Normalized counts
  - For each pairwise group comparison:
    - Log2 fold change
    - Test statistic
    - P-value
    - Adjusted p-value
  - All.mean (mean across all samples)
  - All.stdev (standard deviation across all samples) 
  - LRT.p.value (likelihood ratio test adjusted p-value)
  - For each group:
    - Group.Mean_(group) (mean within group)
    - Group.Stdev_(group) (standard deviation within group))

> Note: Datasets with technical replicates are handled by collapsing them such that the minimum number of equal technical replicates is retained across all samples. Before normalization, the counts of technical replicates are summed to combine them into a single sample representing the biological replicate.

> Note: RNAseq processed data interactive tables and plots are found in the [GLDS visualization portal](https://visualization.genelab.nasa.gov/data/studies).

<br>

---

## 10. Evaluate ERCC Spike-In Data 

> Note: This is only applicable for datasets with ERCC spike-in

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
list_of_ISA_files = zip_file_object.namelist() # Print contents of zip file. Pick relevant one from list
UnnormalizedCountsPath = '/path/to/GLDS-NNN_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv'

GENE_ID_PREFIX = "ENSMU" # change according to a common prefix for all Gene IDs associated with the subject organism

# Print contents of ISA zip file to view file order
list_of_ISA_files

# There are datasets that have multiple assays (including microarray), so the RNAseq ISA files from the above output must be selected. 
# Txt files outputted above are indexed as 0, 1, 2, etc. Fill in the indexed number corresponding to the sample (s_*txt) and assay files for RNAseq (a_*_(RNA-Seq).txt) in the code block below.

# Extract metadata from the sample file (s_*txt)
sample_file = list_of_ISA_files[2] # replace [2] with index corresponding to the (s_*txt) file
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

### Get raw counts table

raw_counts_table = pd.read_csv(UnnormalizedCountsPath, index_col=0) 
raw_counts_table.index.rename('Gene_ID', inplace=True)
print(raw_counts_table.head(n=3))

raw_counts_transcripts = raw_counts_table[raw_counts_table.index.str.contains(f"^{GENE_ID_PREFIX}")]
assert len(raw_counts_transcripts) != 0, f"Looks like {GENE_ID_PREFIX} matched no genes, probably the wrong prefix"
raw_counts_transcripts = raw_counts_transcripts.sort_values(by=list(raw_counts_transcripts), ascending=False)
print(raw_counts_transcripts)

### Get ERCC counts

ercc_counts = raw_counts_table[raw_counts_table.index.str.contains('^ERCC-')] 
ercc_counts.reset_index(inplace=True)
ercc_counts = ercc_counts.rename(columns={'Gene_ID':'ERCC ID'})
ercc_counts = ercc_counts.sort_values(by=list(ercc_counts), ascending=False)
print(ercc_counts.head())

### Get files containing ERCC gene concentrations and metadata

ercc_url = 'https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_095046.txt'
ercc_table = pd.read_csv(ercc_url, sep = '\t')
print(ercc_table.head(n=3))

## Calculate the number of ERCC genes detected in each of the 4 (A, B, C and D) groups for each sample
### Extract ERCC counts and calculate the log(2)

meltERCC = ercc_counts.melt(id_vars=['ERCC ID'])
meltERCC['log2 Count'] = meltERCC['value']+1
meltERCC['log2 Count'] = np.log2(meltERCC['log2 Count'])
meltERCC = meltERCC.rename(columns={'variable':'Sample Name', 'value':'Count'})
print(meltERCC.head(n=3))

### Build Mix dictionary to link sample name to mix added and read depth using the assay table

mix_dict = assay_table.filter(['Sample Name','Parameter Value[Spike-in Mix Number]', 
                       'Parameter Value[Read Depth]'])
mix_dict = mix_dict.rename(columns={'Parameter Value[Spike-in Mix Number]':'Mix',
                                    'Parameter Value[Read Depth]':
                                    'Total Reads'})
print(mix_dict.head(n=3))

### Make combined ercc counts and assay table

merged_ercc = meltERCC.merge(mix_dict, on='Sample Name')
print(merged_ercc)

### Read ERCC info including concentrations from merged_ercc table

groupA = ercc_table.loc[ercc_table['subgroup'] == 'A']['ERCC ID']
groupB = ercc_table.loc[ercc_table['subgroup'] == 'B']['ERCC ID']
groupC = ercc_table.loc[ercc_table['subgroup'] == 'C']['ERCC ID']
groupD = ercc_table.loc[ercc_table['subgroup'] == 'D']['ERCC ID']

### Make a dictionary for ERCC groups

group_dict = dict(zip(ercc_table['ERCC ID'], ercc_table['subgroup']))

### Calculate ERCC counts per million and log(2) counts per million

merged_ercc['Count per million'] = merged_ercc['Count'] / (merged_ercc['Total Reads'] / 1000000.0)
merged_ercc['log2 Count per million'] = np.log2(merged_ercc['Count per million']+1)

### Add ERCC group column

merged_ercc['ERCC group'] = merged_ercc['ERCC ID'].map(group_dict)
merged_ercc = merged_ercc.sort_values(by=['Mix'], ascending=True)
print(merged_ercc)

## Filter and calculate mean counts per million of Mix1 and Mix2 spiked samples in each of the 4 groups
### Filter Mix1 CPM and Mix2 CPM in group A 

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

### Filter Mix1 CPM and Mix2 CPM in group B

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

### Filter Mix1 CPM and Mix2 CPM in group C

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

### Filter Mix1 CPM and Mix2 CPM in group D

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

## Multi-sample ERCC analyses
### Create box and whisker plots of the log(2) CPM for each ERCC detected in group A in Mix 1 and Mix 2 spiked samples

a = sns.catplot(x="ERCC ID", y="log2 Count per million", order=groupA[::-1], hue="Mix",data=merged_ercc[merged_ercc['ERCC ID'].isin(groupA)], kind="box", col="ERCC group", height=5, aspect=1, palette=sns.color_palette(['blue', 'orange']), showfliers=False)
a.set_xticklabels(rotation=90)
plt.text(23,2.5,"Mix1/ Mix2 = 4")

### Create bar plot of the average Mix1 CPM / average Mix 2 CPM for group A ERCC genes (for group A we expect Mix 1 CPM / Mix 2 CPM = 4)

a1 = sns.catplot(x="ERCC ID", y="Avg Mix1 CPM/ Avg Mix2 CPM", order=groupA[::-1], palette="rocket_r", data=adf, kind="bar", height=5, aspect=1, linewidth=0.5)
a1.set_xticklabels(rotation=90)
plt.title("ERCC Group A")
a1.set(ylim=(0, 6))
a1.set_axis_labels("ERCC genes ordered by concentration: low \u2192 high")
print('Number of ERCC detected in group A (out of 23) =', adf['Avg Mix1 CPM/ Avg Mix2 CPM'].count())

### Create box and whisker plots of the log(2) CPM for each ERCC detected in group B in Mix 1 and Mix 2 spiked samples

b = sns.catplot(x="ERCC ID", y="log2 Count per million", order=groupB[::-1], hue="Mix", data=merged_ercc[merged_ercc['ERCC ID'].isin(groupB)], kind="box", col="ERCC group", height=5, aspect=1, palette=sns.color_palette(['blue', 'orange']), showfliers=False)
b.set_xticklabels(rotation=90)
plt.text(23,2.5,"Mix1/ Mix2 = 1")

### Create bar plot of the average Mix1 CPM / average Mix 2 CPM for group B ERCC genes (for group B we expect Mix 1 CPM / Mix 2 CPM = 1)

b = sns.catplot(x="ERCC ID", y="Avg Mix1 CPM/ Avg Mix2 CPM", order=groupB[::-1], palette="rocket_r", data=bdf, kind="bar", 
               height=5, aspect=1, linewidth=0.5)
b.set_xticklabels(rotation=90)
plt.title("ERCC Group B")
b.set(ylim=(0, 2))
b.set_axis_labels("ERCC genes ordered by concentration: low \u2192 high")
print('Number of ERCC detected in group B (out of 23) =', bdf['Avg Mix1 CPM/ Avg Mix2 CPM'].count())

### Create box and whisker plots of the log(2) CPM for each ERCC detected in group C in Mix 1 and Mix 2 spiked samples

c = sns.catplot(x="ERCC ID", y="log2 Count per million", order=groupC[::-1], hue="Mix", data=merged_ercc[merged_ercc['ERCC ID'].isin(groupC)], kind="box", col="ERCC group", height=5, aspect=1, palette=sns.color_palette(['blue', 'orange']), showfliers=False)
c.set_xticklabels(rotation=90)
plt.text(23,2.5,"Mix1/ Mix2 = 0.67")

### Create bar plot of the average Mix1 CPM / average Mix 2 CPM for group C ERCC genes (for group C we expect Mix 1 CPM / Mix 2 CPM = 0.67)

c = sns.catplot(x="ERCC ID", y="Avg Mix1 CPM/ Avg Mix2 CPM", order=groupC[::-1], palette="rocket_r", data=cdf, kind="bar", 
               height=5, aspect=1, linewidth=0.5)
c.set_xticklabels(rotation=90)
plt.title("ERCC Group C")
c.set(ylim=(0, 2))
c.set_axis_labels("ERCC genes ordered by concentration: low \u2192 high")
print('Number of ERCC detected in group C (out of 23) =', cdf['Avg Mix1 CPM/ Avg Mix2 CPM'].count())

### Create box and whisker plots of the log(2) CPM for each ERCC detected in group D in Mix 1 and Mix 2 spiked samples

d = sns.catplot(x="ERCC ID", y="log2 Count per million", order=groupD[::-1], hue="Mix", data=merged_ercc[merged_ercc['ERCC ID'].isin(groupD)], col="ERCC group", kind="box", height=5, aspect=1, palette=sns.color_palette(['blue', 'orange']), showfliers=False)
d.set_xticklabels(rotation=90)
plt.text(23,2.5,"Mix1/ Mix2 = 0.5")

### Create bar plot of the average Mix1 CPM / average Mix 2 CPM for group D ERCC genes (for group D we expect Mix 1 CPM / Mix 2 CPM = 0.5)

d = sns.catplot(x="ERCC ID", y="Avg Mix1 CPM/ Avg Mix2 CPM", order=groupD[::-1], palette="rocket_r", data=ddf, kind="bar", 
               height=5, aspect=1, linewidth=0.5)
d.set_xticklabels(rotation=90)
plt.title("ERCC Group D")
d.set(ylim=(0, 1))
d.set_axis_labels("ERCC genes ordered by concentration: low \u2192 high")
print('Number of ERCC detected in group D (out of 23) =', ddf['Avg Mix1 CPM/ Avg Mix2 CPM'].count())

## Individual sample ERCC analyses
# Calculate and plot ERCC metrics from individual samples, including limit of detection, dynamic range, and R^2 of counts vs. concentration.

#Calculate and plot ERCC metrics from individual samples, including limit of detection, dynamic range, and R^2 of counts vs. concentration.
print(ercc_table.head(n=3))

# Make a dictionary for ERCC concentrations for each mix

mix1_conc_dict = dict(zip(ercc_table['ERCC ID'], ercc_table['concentration in Mix 1 (attomoles/ul)']))
mix2_conc_dict = dict(zip(ercc_table['ERCC ID'], ercc_table['concentration in Mix 2 (attomoles/ul)']))

# Check assay_table header to identify the 'Sample Name' column and the column title indicating the 'Spike-in Mix Nmber' if it's indicated in the metadata.

pd.set_option('display.max_columns', None)
print(assay_table.head(n=3))

# Get samples that use mix 1 and mix 2

mix1_samples = assay_table[assay_table['Parameter Value[Spike-in Mix Number]'] == 'Mix 1']['Sample Name']
mix2_samples = assay_table[assay_table['Parameter Value[Spike-in Mix Number]'] == 'Mix 2']['Sample Name']

# Get ERCC counts for all samples

ercc_counts = raw_counts_table[raw_counts_table.index.str.contains('^ERCC-')] 
ercc_counts = ercc_counts.sort_values(by=list(ercc_counts), ascending=False)
print(ercc_counts.head())

# Get ERCC counts for Mix 1 spiked samples

ercc_counts_mix_1 = ercc_counts[mix1_samples]
ercc_counts_mix_1['ERCC conc (attomoles/ul)'] = ercc_counts_mix_1.index.map(mix1_conc_dict)
print(ercc_counts_mix_1.head(n=3))

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

# Iterate over subplot positions, if subplot does not have data, hide the subplot with ax.set_visible(False)

for i in range(side_size):
    for j in range(side_size):
        ax = axs[i, j]
        index = i * side_size + j

        if index < total_columns:
            if index < len(columns_mix_1):
                ax.scatter(x=np.log2(ercc_counts_mix_1['ERCC conc (attomoles/ul)']), y=np.log2(ercc_counts_mix_1[all_columns[index]]+1), s=7)
                ax.set_title(all_columns[index][-45:], fontsize=9)
            else:
                ax.scatter(x=np.log2(ercc_counts_mix_2['ERCC conc (attomoles/ul)']), y=np.log2(ercc_counts_mix_2[all_columns[index]]+1), s=7)
                ax.set_title(all_columns[index][-45:], fontsize=9)

            ax.set_xlabel('log2 ERCC conc (attomoles/ ul)', fontsize=9)
            ax.set_ylabel('log2 Counts per million', fontsize=9)
            ax.tick_params(direction='in', axis='both', labelsize=9, labelleft=True, labelbottom=True)
        else:
            ax.set_visible(False)  # Hide the subplot if it's not needed

plt.show()

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


# Plot each sample using linear regression of scatter plot with x = log2 Conc and y = log2 Counts.
# Return min, max, R^2 and dynamic range (max / min) values.

samples = []
mins = []
maxs = []
dyranges = []
rs = []

fig, axs = plt.subplots(side_size, side_size, figsize=(22,26), sharex='all', sharey='all');
fig.tight_layout(pad=1, w_pad=2.5, h_pad=3.5)

counter = 0
list2counter = 0
for i in range(side_size):
    for j in range(side_size):
        ax = axs[i, j]
        index = i * side_size + j

        if index < len(columns_mix_1):
            nonzero_counts = nonzero_counts_list_1[index]
            xvalues = nonzero_counts['Conc']
            yvalues = nonzero_counts['Counts']

            if len(xvalues) > 0:
                sns.regplot(x=np.log2(xvalues), y=np.log2(yvalues), ax=ax)
                ax.set_title(all_columns[index][-47:], fontsize=9)
                ax.set_xlabel('log2 Conc (attomoles/ul)', fontsize=9)
                ax.set_ylabel('log2 Counts per million', fontsize=9)
                ax.tick_params(direction='in', axis='both', labelsize=9, labelleft=True, labelbottom=True)
                samples.append(all_columns[index])

                min_val = xvalues.min()
                max_val = xvalues.max()
                dynamic_range = max_val / min_val
                slope, intercept, r_value, p_value, std_err = linregress(np.log2(xvalues), np.log2(yvalues))

                ax.text(0.02, 0.98, f'Min:{min_val:.1f}', verticalalignment='top', horizontalalignment='left', transform=ax.transAxes, color='black', fontsize=10)
                ax.text(0.02, 0.88, f'Max:{max_val:.1f}', verticalalignment='top', horizontalalignment='left', transform=ax.transAxes, color='black', fontsize=10)
                ax.text(0.02, 0.78, f'Dyn:{dynamic_range:.1f}', verticalalignment='top', horizontalalignment='left', transform=ax.transAxes, color='black', fontsize=10)
                ax.text(0.02, 0.68, f'R:{r_value:.2f}', verticalalignment='top', horizontalalignment='left', transform=ax.transAxes, color='black', fontsize=10)

                mins.append(min_val)
                maxs.append(max_val)
                dyranges.append(dynamic_range)
                rs.append(r_value)
            else:
                ax.set_visible(False)

        elif index < total_columns:
            nonzero_counts = nonzero_counts_list_2[list2counter]
            xvalues = nonzero_counts['Conc']
            yvalues = nonzero_counts['Counts']

            if len(xvalues) > 0:
                sns.regplot(x=np.log2(xvalues), y=np.log2(yvalues), ax=ax)
                ax.set_title(all_columns[index][-47:], fontsize=9)
                ax.set_xlabel('log2 Conc (attomoles/ul)', fontsize=9)
                ax.set_ylabel('log2 Counts per million', fontsize=9)
                ax.tick_params(direction='in', axis='both', labelsize=9, labelleft=True, labelbottom=True)
                samples.append(all_columns[index])

                min_val = xvalues.min()
                max_val = xvalues.max()
                dynamic_range = max_val / min_val
                slope, intercept, r_value, p_value, std_err = linregress(np.log2(xvalues), np.log2(yvalues))

                ax.text(0.02, 0.98, f'Min:{min_val:.1f}', verticalalignment='top', horizontalalignment='left', transform=ax.transAxes, color='black', fontsize=10)
                ax.text(0.02, 0.88, f'Max:{max_val:.1f}', verticalalignment='top', horizontalalignment='left', transform=ax.transAxes, color='black', fontsize=10)
                ax.text(0.02, 0.78, f'Dyn:{dynamic_range:.1f}', verticalalignment='top', horizontalalignment='left', transform=ax.transAxes, color='black', fontsize=10)
                ax.text(0.02, 0.68, f'R:{r_value:.2f}', verticalalignment='top', horizontalalignment='left', transform=ax.transAxes, color='black', fontsize=10)

                mins.append(min_val)
                maxs.append(max_val)
                dyranges.append(dynamic_range)
                rs.append(r_value)
            else:
                ax.set_visible(False)

            list2counter += 1
        else:
            ax.set_visible(False)

plt.show()

# Create directory for saved files

import os
os.makedirs(name="ERCC_analysis", exist_ok=True)

# Print tables containing the dynamic range and R^2 values for each sample.
# Remember to change file names to GLDS# analyzing

stats = pd.DataFrame(list(zip(samples, mins, maxs, dyranges, rs)))
stats.columns = ['Samples', 'Min', 'Max', 'Dynamic range', 'R']
stats.to_csv(f'ERCC_analysis/ERCC_stats_{accession}_GLbulkRNAseq.csv', index = False)
stats.filter(items = ['Samples', 'Dynamic range']).to_csv(f'ERCC_analysis/ERCC_dynrange_{accession}_mqc_GLbulkRNAseq.csv', index = False)
stats.filter(items = ['Samples', 'R']).to_csv(f'ERCC_analysis/ERCC_rsq_{accession}_mqc_GLbulkRNAseq.csv', index = False)

## Generate data and metadata files needed for ERCC DESeq2 analysis

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

ERCCmetadata.to_csv('ERCC_analysis/ERCCmetadata_GLbulkRNAseq.csv') 

# Export ERCC count data

ercc_counts.columns = ercc_counts.columns.str.replace('-','_')
ERCCcounts = ercc_counts.loc[:,ERCCmetadata.index]
ERCCcounts.head()

ERCCcounts.to_csv('ERCC_analysis/ERCCcounts_GLbulkRNAseq.csv') 
```


**Input Data:**

- *ISA.zip (compressed ISA directory containing Investigation, Study, and Assay (ISA) metadata files for the respective GLDS dataset, output from [Step 9a](#9a-create-sample-runsheet))
- RSEM_Unnormalized_Counts.csv (RSEM raw counts table, output from [Step 9](#ERCCspikeOut))

**Output Data:**

- ERCC_analysis/ERCC_stats_GLDS-*_GLbulkRNAseq.csv (Samplewise counts statistics table containing 'Min', 'Max', 'Dynamic range', 'R')
- ERCC_analysis/ERCC_dynrange_GLDS-*_GLbulkRNAseq.csv (Samplewise counts statistics subset table containing 'Dynamic range')
- ERCC_analysis/ERCC_rsq_GLDS-*_GLbulkRNAseq.csv (Samplewise counts statistics subset table containing 'R')
- ERCC_analysis/ERCCmetadata_GLbulkRNAseq.csv (Samplewise metadata table inlcuding ERCC mix number)
- ERCC_analysis/ERCCcounts_GLbulkRNAseq.csv (Samplewise ERCC counts table)

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

cts <- as.matrix(read.csv('ERCC_analysis/ERCCcounts_GLbulkRNAseq.csv',sep=",",row.names="Gene_ID")) #INPUT
coldata <- read.csv('ERCC_analysis/ERCCmetadata_GLbulkRNAseq.csv', row.names=1) #INPUT

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

## Try first to use the default type="median", but if there is an error (usually due to zeros in genes), use type="poscounts"
## From DESeq2 manual: "The "poscounts" estimator deals with a gene with some zeros, by calculating a modified geometric mean by taking the n-th root of the product of the non-zero counts."
dds <- tryCatch(
      expr = { estimateSizeFactors(dds) },
      error = function(e) { estimateSizeFactors(dds, type="poscounts")}
)


dds <- DESeq(dds)
res <- results(dds, contrast=c("Mix","Mix 1","Mix 2")) # remove space before mix number if needed
res

## Export DESeq2 results table and normalized ERCC counts table

write.csv(res, 'ERCC_analysis/ERCC_DESeq2_GLbulkRNAseq.csv') #OUTPUT
normcounts = counts(dds, normalized=TRUE)
write.csv(normcounts, 'ERCC_analysis/ERCC_normcounts_GLbulkRNAseq.csv') #OUTPUT
```

**Input Data:**

- ERCC_analysis/ERCCmetadata_GLbulkRNAseq.csv (samplewise metadata table inlcuding ERCC mix number, output from [Step 10a](#10a-evaluate-ercc-count-data-in-python))
- ERCC_analysis/ERCCcounts_GLbulkRNAseq.csv (samplewise ERCC counts table, output from [Step 10a](#10a-evaluate-ercc-count-data-in-python))

**Output Data:**

- ERCC_analysis/ERCC_DESeq2_GLbulkRNAseq.csv (DESeq2 results table)
- ERCC_analysis/ERCC_normcounts_GLbulkRNAseq.csv (normalized ERCC Counts table)

<br>

### 10c. Analyze ERCC DESeq2 Results in Python

```python

# Import python packages

## Import python packages

import pandas as pd
from urllib.request import urlopen, quote, urlretrieve
import seaborn as sns
import matplotlib.pyplot as plt

accession = 'GLDS-NNN' # Replace Ns with GLDS number

## Import ERCC DESeq2 results

deseq2out = pd.read_csv('ERCC_analysis/ERCC_DESeq2_GLbulkRNAseq.csv', index_col=0) # INPUT
#deseq2out.index = deseq2out.index.str.replace('_','-')
deseq2out.rename(columns ={'baseMean' : 'meanNormCounts'}, inplace = True)
print(deseq2out.head())

## Get files containing ERCC gene concentrations and metadata

ercc_url = 'https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_095046.txt'
ercc_table = pd.read_csv(ercc_url, sep='\t', index_col='ERCC ID')
print(ercc_table.head(n=3))

## Combine ERCC DESeq2 results and ercc_table

combined = deseq2out.merge(ercc_table, left_index=True, right_index=True)
print(combined.head())

## Filter p-value and adj. p-value cutoff at 10^-3

combined['cleaned_padj'] = combined['padj']
combined.loc[(combined.cleaned_padj < 0.001),'cleaned_padj']=0.001

combined['cleaned_pvalue'] = combined['pvalue']
combined.loc[(combined.cleaned_pvalue < 0.001),'cleaned_pvalue']=0.001

print(combined.head())

## Export the filtered combined ERCC DESeq2 results and ercc_table
### Remember to change file name to GLDS# analyzing

combined.filter(items = ['ERCC ID', 'meanNormCounts', 'cleaned_pvalue','cleaned_padj']).to_csv(f'ERCC_analysis/ERCC_lodr_{accession}_mqc_GLbulkRNAseq.csv') 

## Plot p-value vs. mean normalized ERCC counts

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

## Plot Adjp-value vs. mean normalized ERCC counts

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

- ERCC_analysis/ERCC_DESeq2_GLbulkRNAseq.csv (ERCC DESeq2 results table, output from [Step 10b](#10b-perform-deseq2-analysis-of-ercc-counts-in-r))

**Output Data:**

- ERCC_analysis/ERCC_lodr_*_GLbulkRNAseq.csv (ERCC Gene Table including mean counts, adjusted p-value and p-value, and filtered to genes with both adj. p-value and p-value < 0.001)

> All steps of the ERCC Spike-In Data Analysis are performed in a Jupyter Notebook (JN) and the completed JN is exported as an html file (**ERCC_analysis_GLbulkRNAseq.html**) and published in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/) for the respective dataset.
