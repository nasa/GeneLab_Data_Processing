# GeneLab bioinformatics processing pipeline for miRNA-Seq sequencing data <!-- omit in toc -->

> **This page holds an overview and instructions for how GeneLab processes miRNA Sequencing (miRNA-Seq) datasets. Exact processing commands and GL-DPPD-XXXX version used for specific GeneLab datasets (GLDS) are provided with their processed data in the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo).**  
>

---

**Date:** Month DD, YYYY  
**Revision:** -  
**Document Number:** GL-DPPD-XXXX

**Submitted by:**  
Richard Barker (GeneLab Data Processing Team) and Barbara Novak (GeneLab Data Processing Lead)

**Approved by:**  

---

# Table of contents <!-- omit in toc -->

- [**Software used**](#software-used)
- [**General processing overview with example commands**](#general-processing-overview-with-example-commands)
  - [1. **Raw Data Quality Control (QC)**](#1-raw-data-quality-control-qc)
    - [1a. Raw Data QC](#1a-raw-data-qc)
    - [1b. Compile Raw Data QC](#1b-compile-raw-data-qc)
  - [2. **Trim/Filter Raw Data and Trimmed Data QC**](#2-trimfilter-raw-data-and-trimmed-data-qc)
    - [2a. Trim/Filter Raw Data](#2a-trimfilter-raw-data)
    - [2b. Trimmed Data QC](#2b-trimmed-data-qc)
    - [2c. Compile Trimmed Data QC](#2c-compile-trimmed-data-qc)
  - [3. **miRNA QC**](#3-mirna-qc)
  - [4. **Deduplicate reads**](#4-deduplicate-reads)
  - [5. **Prepare reference files**](#5-prepare-reference-files)
    - [5a. Remove Whitespace](#5a-remove-whitespace)
    - [5b. Build Bowtie Index for Reference Genome](#5b-build-bowtie-index-for-reference-genome)
    - [5c. Create miRBase species-specific references](#5c-create-mirbase-species-specific-references)
  - [6. **Predict Novel miRNAs**](#6-predict-novel-mirnas)
    - [6a. **Standard miRNA Prediction**](#6a-standard-mirna-prediction)
      - [6ai. Align Reads to Reference Genome](#6ai-align-reads-to-reference-genome)
      - [6aii. Convert Alignment File to mirdeep2 (ARF) format](#6aii-convert-alignment-file-to-mirdeep2-arf-format)
      - [6aiii. Predict Known and Novel miRNAs for Animals](#6aiii-predict-known-and-novel-mirnas)
      - [6aiv. Extract miRDeep2 counts](#6aiv-extract-mirdeep2-counts)
      - [6av. Prepare new species-specific reference](#6av-prepare-new-species-specific-reference)
    - [6b. **Plant miRNA Prediction**](#6b-plant-mirna-prediction)
      - [6bi. Align/process/predict known and novel miRNAs](#6bi-predict-known-and-novel-mirnas-for-plants)
      - [6bii. Convert miRNA predictions to FASTA format](#6bii-convert-mirna-predictions-to-fasta)
      - [6biii. Blast predicted to known mature sequences](#6biii-blast-predicted-to-known-mature-sequences)
      - [6biv. Prepare new species-specific reference](#6biv-prepare-new-species-specific-reference)
  - [7. Generate Read Counts for Known and Novel miRNAs](#7-generate-counts-for-known-and-predicted-mirna)
    - [7a. Align reads to species specific reference](#7a-align-reads-to-specifies-specific-reference)
    - [7b. Generate basic alignment statistics and read counts](#7b-generate-basic-alignment-statistics-and-read-counts)
    - [7c. Compile alignment QC](#7c-compile-alignment-qc)
    - [7d. Merge counts files](#7d-merge-read-counts)
  - [8. **Normalize Read Counts, Perform Differential Gene Expression Analysis, and Add Gene Annotations in R**](#8-differential-expression-analysis)
    - [8a Environment Setup](#8a-environment-setup)
    - [8b. Configure Metadata, Sample Grouping, and Group Comparisons](#8b-configure-metadata-sample-grouping-and-group-comparisons)
    - [8c. Import miRNA counts](#8c-import-mirna-counts)
    - [8d. Perform DGE](#8d-perform-dge)
    - [8e. Prepare GeneLab DGE Tables with Annotations](#8e-prepare-genelab-dge-tables-with-annotations)
    - [8f. Export GeneLab DGE Tables with Annotations](#8f-export-genelab-dge-tables-with-annotations)

---

# Software used

| Program            | Version | Relevant Links                                            |
|--------------------|---------|-----------------------------------------------------------|
| FastQC             | 0.12.1  | [https://www.bioinformatics.babraham.ac.uk/projects/fastqc/](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)  |
| Cutadapt           | 4.7     | [https://cutadapt.readthedocs.io/en/stable/](https://cutadapt.readthedocs.io/en/stable/) |
| TrimGalore!        | 0.6.10  | [https://github.com/FelixKrueger/TrimGalore](https://github.com/FelixKrueger/TrimGalore) |
| Bowtie             | 1.3.1   | [https://bowtie-bio.sourceforge.net/bowtie2/index.shtml](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) |
| MultiQC            | 1.21    | [https://multiqc.info/](https://multiqc.info/) |
| SAM tools          | 1.19.2  | [https://www.htslib.org/](https://www.htslib.org/) |
| Bioconductor       | 3.18    | [https://bioconductor.org/](https://bioconductor.org/)
| DESeq2             | 1.42.1  | [https://bioconductor.org/packages/release/bioc/html/DESeq2.html](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) |
| tidyverse          | 2.0.0   | [https://cran.r-project.org/web/packages/tidyverse/index.html](https://cran.r-project.org/web/packages/tidyverse/index.html)  |
| mirbase.db         | 1.2.0   | [https://bioconductor.org/packages/release/data/annotation/html/mirbase.db.html](https://bioconductor.org/packages/release/data/annotation/html/mirbase.db.html)  |
| MultiMiR           | 1.24.0  | [https://bioconductor.org/packages/release/bioc/html/multiMiR.html](https://bioconductor.org/packages/release/bioc/html/multiMiR.html)|
| R                  | 4.3.1   | [https://www.r-project.org/](https://www.r-project.org/) |
| mirDeep2           | 2.0.1.3-1 |[https://github.com/rajewsky-lab/mirdeep2](https://github.com/rajewsky-lab/mirdeep2)|
| miRDP2             | 1.1.4   | [https://sourceforge.net/projects/mirdp2/files/version%201.1.4/](https://sourceforge.net/projects/mirdp2/files/version%201.1.4/)|
| Perl               | 5.32.1  | [https://www.perl.org/](https://www.perl.org/)|
| Perl PDF::API2     | 2.045   | [https://metacpan.org/search?q=PDF%3A%3AAPI2](https://metacpan.org/search?q=PDF%3A%3AAPI2)|
| ViennaRNA          | 2.6.4   | [https://www.tbi.univie.ac.at/RNA/](https://www.tbi.univie.ac.at/RNA/)|
| Squid library      | 1.9g    | [http://eddylab.org/software.html](http://eddylab.org/software.html)|
| randfold           | 2.0.1-6 | [http://bioinformatics.psb.ugent.be/software/details/Randfold](http://bioinformatics.psb.ugent.be/software/details/Randfold)|
| seqkit             | v2.8.0  | [https://bioinf.shenwei.me/seqkit/](https://bioinf.shenwei.me/seqkit/)
| *BLAST+*             | 2.15.0  | [https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html#downloadblastdata](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html#downloadblastdata)|
| *Bio::Seq*           | 1.7.8   | [https://metacpan.org/pod/Bio::Seq](https://metacpan.org/pod/Bio::Seq) |
| *Bio::SeqIO*         | 1.7.8   | [https://metacpan.org/pod/Bio::SeqIO](https://metacpan.org/pod/Bio::SeqIO) |
| *Bio::SearchIO*      | 1.7.8   | [https://metacpan.org/pod/Bio::SearchIO](https://metacpan.org/pod/Bio::SearchIO) |
| *Switch*             | 2.17    | [https://metacpan.org/dist/Switch](https://metacpan.org/dist/Switch)
---


# Reference databases used

| Database | Version | Relevant Links |
|:---------|:-------:|---------------------------------------------------:|
| miRBase  | v22     | [https://www.mirbase.org](https://www.mirbase.org) |
| *ncRNA*  |         |                                                    |
---

# General processing overview with example commands  

The exact processing commands for a specific GLDS that has been released are provided with the processed data in the [OSDR](https://osdr.nasa.gov/bio/repo).

> All output files in **bold** are published with the miRNA-Seq processed data in the [OSDR](https://osdr.nasa.gov/bio/repo).

---
## 1. Raw Data Quality Control (QC)  

<br>

### 1a. Raw Data QC  

```bash
fastqc -t NumberOfThreads \
  -o /path/to/raw_fastqc_output \
  *.fastq.gz
```

**Parameter Definitions:**

- `-t` - the number of threads available to perform fastqc
- `-o` – the output directory to store results
- `*.fastq.gz` – the input reads are specified as a positional argument, and can be given all at once with wildcards like this, or as individual arguments with spaces in between them

**Input Data:**

- fastq, compressed or uncompressed (raw reads)

**Output Data:**

- fastqc.html (FastQC output html summary)
- fastqc.zip (FastQC output data)

<br>

### 1b. Compile Raw Data QC  

```bash
multiqc --interactive \
  -n raw_multiqc_GLmiRNAseq \
  -o /path/to/raw_multiqc/output/directory \
  /path/to/directory/containing/raw_fastqc_files
```

**Parameter Definitions:**

- `--interactive` - force reports to use interactive plots
- `-n` - prefix name for output files
- `-o` – the output directory to store results
- `/path/to/directory/containing/raw_fastqc_files` – the directory holding the output data from the fastqc run, provided as a positional argument

**Input Data:**

- *fastqc.zip (FastQC data, output from [Step 1a](#1a-raw-data-qc))

**Output Data:**

- **raw_multiqc_GLmiRNAseq.html** (multiqc output html summary)
- **raw_multiqc_GLmiRNAseq_data** (directory containing multiqc data)

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
  --length 17 \ 
  --quality 20 \
  --output_dir /path/to/TrimGalore/output/directory \
  /path/to/raw/fastq/files

for i in *_trimmed.fq.gz; do out=${i%fq.gz}_trimmed.fastq.gz; mv $i $out; done;
```

**Parameter Definitions:**

- `--gzip` – compress the output files with `gzip`
- `--path_to_cutadapt` - specify path to cutadapt software if it is not in your `$PATH`
- `--cores` - specify the number of threads available to perform trimming
- `--phred33` - instructs cutadapt to use ASCII+33 quality scores as Phred scores for quality trimming
- `--length` - specifies the minimum read length to keep, if reads become shorter than this value they are removed. 17 is used because it is the shortest length accepted by mirdeep2.
- `--quality` - Trim low-quality reads in addition to adapter removal. If data was generated on an Illumina 2-colour instrument (such as NextSeq or NovaSeq), use `--2colour` instead.
- `--output_dir` - the output directory to store results
- `/path/to/raw/fastq/files` - specify the input reads as a positional argument

**Input Data:**

- *fastq.gz (raw reads)

**Output Data:**

- **\*trimmed.fastq.gz** (trimmed reads)
- **\*trimming_report.txt** (trimming report)

<br>

### 2b. Trimmed Data QC  

```bash
fastqc -t NumberOfThreads \
  -o /path/to/trimmed_fastqc/output/directory \
  *.fastq.gz
```

**Parameter Definitions:**

- `-t` - specify the number of threads available to perform fastqc
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
multiqc --interactive \
  --zip-data-dir \
  -n trimmed_multiqc_GLmiRNAseq \
  -o /path/to/trimmed_multiqc/output/trimmed_multiqc_GLmiRNAseq_report \
  /path/to/directory/containing/trimmed_fastqc/files
```

**Parameter Definitions:**

- `--interactive` - force reports to use interactive plots
- `--zip-data-dir` - compress the output data folder
- `-n` - prefix name for output files
- `-o` – the output directory to store results
- `/path/to/directory/containing/trimmed_fastqc/files` – the directory holding the output data from the fastqc run, provided as a positional argument

**Input Data:**

- *fastqc.zip (FastQC data, output from [Step 2b](#2b-trimmed-data-qc))

**Output Data:**

- **trimmed_multiqc_GLmiRNAseq.html** (multiqc report)
- **/trimmed_multiqc_GLmiRNAseq_data** (directory containing multiqc data)

<br>

---
## 3. miRNA QC  

```bash
mirtrace qc --species ${species} \
  --protocol ${protocol} \
  --num-threads ${threads} \
  --output-dir /path/to/mirtrace/output/directory \
  /path/to/directory/containing/trimmed_fastq/files
```

**Parameter Definitions:**

- `qc` - positional argument specifying "qc" mode. Runs all qc modules and generates full report.
- `--species ` - miRBase species code
- `--protocol` - sample prep protocol defining the read structure 
  * `illumina` - (miRNA--3'-adapter--index) [DEFAULT]
  * `qiaseq` - (miRNA--3'-adapter--UMI--3'-adapter--index), NOTE: Only the first (leftmost) 3' adapter should be specified.
  * `cats` - (NNN--miRNA--poly-A--3'-adapter--index), NOTE: It's not possible to specify an adapter for -p cats.
  * `nextflex` - (NNNN--miRNA--NNNN--3'-adapter--index)
- `--num-threads - specify the number of threads available to perform qc
- `--output-dir` – the output directory to store results
- `/path/to/directory/containing/trimmed_fastq/files` – the trimmed fastq.gz files generated in [Step 2a](#2a-trimfilter-raw-data)

**Input Data:**

- *trimmed.fastq.gz (trimmed fastq files, output from [Step 2a](#2a-trimfilter-raw-data))

**Output Data:**

- **mirtrace-report.html**				
- mirtrace-stats-length.tsv			
- mirtrace-stats-rnatype.tsv
- mirtrace-results.json				
- mirtrace-stats-mirna-complexity.tsv		
- mirtrace-stats-contamination_basic.tsv		
- mirtrace-stats-phred.tsv			
- mirtrace-stats-contamination_detailed.tsv	
- mirtrace-stats-qcstatus.tsv

<br>

---
## 4. Deduplicate Reads


```bash
  zcat ${in_dir}/${sample}_trimmed.fastq.gz | perl -s -E '\$c=1; while(<>){if(/^@/){\$seq=<>;chomp(\$seq); <>;<>; \$hash{\$seq}++}} foreach my \$s(sort keys %hash){printf \">\".\$name.\"_%08d_x\$hash{\$s}\n\$s\n\",\$c; \$c++;}' -- -name=${prefix} > ${out_dir}/${sample}_collapsed_reads.fasta
```

**Parameter Definitions:**

- `${in_dir}` - the input directory, specified as an environment variable
- `${sample}` - the sample name specified as an environment variable
- `| perl` - directs each sample to be run through the perl command indicated to remove duplicates
- `-name=${prefix}` - specifies the three letter string to use as the base for forming the output readnames, specified as an environment variable
- `${out_dir}` - output directory, specified as an environment variable

**Input Data:**

- *trimmed.fastq.gz (trimmed reads, output from [Step 2a](#2a-trimfilter-raw-data))

**Output Data:**

- *collapsed_reads.fa (deduplicated (collapsed) trimmed reads)

<br>

--- 
## 5. Prepare reference files
### 5a. Remove whitespace

```bash
seqkit seq -i -o /path/to/genome_nows.fasta /path/to/genome.fasta
```

**Parameter Definitions:**

- `-i` - only print ID portion of sequence name (everything up to the first space)
- `-o` - specifies output file name
- `/path/to/genome.fasta` - specify the input genome file, provided as a positional argument

**Input Data:**

- genome.fasta - genome sequence, this miRNA-Seq pipeline version uses the Ensembl fasta file indicated in the fasta column of the [GL-DPPD-7110_annotations.csv](https://github.com/nasa/GeneLab_Data_Processing/blob/master/GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110/GL-DPPD-7110_annotations.csv) GeneLab Annotations file)

**Output Data:**

- `/path/to/genome_nows.fasta` - genome sequence, with fullnames truncated at first whitespace character

<br>

### 5b. Build Bowtie Index for Reference Genome 

```bash
bowtie-build --threads NumberOfThreads \
  /path/to/genome_nows.fasta \
  {index_prefix}
```

**Parameter Definitions:**

- `--threads` - specify the number of threads available to build the index
- `/path/to/genome_nows.fasta` - specifies the genome reference fasta file, provided as a positional argument
- `{index_prefix}` - specifies the prefix to use for the Bowtie index files, provided as a positional argument

**Input Data:**

- genome_nows.fasta - genome sequence, with fullnames truncated at first whitespace character from [Step 5a.](#5a-remove-whitespace)

**Output Data:**

Bowtie index, which consists of the following files:

- {index_prefix}.1.ebwt 
- {index_prefix}.2.ebwt 
- {index_prefix}.3.ebwt 
- {index_prefix}.4.ebwt 
- {index_prefix}.rev.1.ebwt 
- {index_prefix}.rev.2.ebwt 

<br>

### 5c. Create miRBase species-specific references

Using mirdeep2's `extract_miRNAs.pl` script, convert miRBase reference from RNA to DNA and create 2 sets of files:
- species-specific mature and hairpin reference fastas
- related species-specific mature miRNA reference fasta

```bash
extract_miRNAs.pl /path/to/mirbase/mature.fa {species_code} > mature_ref.fa
extract_miRNAs.pl /path/to/mirbase/hairpin.fa {species_code} > hairpin_ref.fa
extract_miRNAs.pl /path/to/mirbase/mature.fa {near_species_codes} > mature_other.fa
```

**Parameter Definitions:**

- `/path/to/mirbase/mature.fa` - specifies the miRBase mature fasta file, provided as a positional argument
- `/path/to/mirbase/hairpin.fa` - specifies the miRBase hairpin fasta file, provided as a positional argument
- `{species_code}` - 3 letter miRBase species code for the species of interest, provided as a positional argument
- `{near_species_codes}` - comma separated list of 3 letter miRBase species codes for related species, provided as a positional argument

**Input Data:**

- mature.fa - miRBase mature sequences
- hairpin.fa - miRBase hairpin sequences


**Output Data:**
Species specific reference fasta files in fasta format
- mature_ref.fa
- hairpin_ref.fa 
- mature_other.fa

<br>

---
## 6. Predict novel miRNAs
### 6a. Standard miRNA prediction
### 6ai. Align Reads to Reference Genome

```bash
bowtie --threads NumberOfThreads \
  --seed 12345 \ 
  --seedlen 16 \
  --seedmms 0 \
  --maqerr 80 \
  --all \
  -m 5 \
  --best \
  --strata \
  -x /path/to/bowtie/reference/{index_prefix} \
  -f /path/to/collapsed_reads \
  > /path/to/bowtie/output/directory/{sample_id}_mapped.bwt
```

**Parameter Definitions:**

- `--threads` - specify the number of threads available to perform alignment
- `--seed` - specify the seed to use for pseudo-random number generator to ensure reproducibility 
- `--seedlen` - specifies the seed length, 16 is used due to the shorter length of miRNA reads
- `--seedmms` - specifies the number of mismatches allowed in the seed, a value of 0 means no mismatches are allowed in the seed
- `--maqerr` - specifies the maximum sum of quality values at each mismatch position. When input files are in fasta format, bowtie assigns a default value of 40 for the quality. With a quality value of 40 for each base (the default when input is a fasta file), "-e 80" allows up to 2 mismatches to occur after the seed region of a read alignment.
- `--all` - Report all valid alignments per read
- `-m` - suppresses all alignments for a read if more than the indicated number of reportable alignments exist for it
- `--best` - guarantees that reported singleton alignments are the "best" in terms of stratum, i.e. fewest number of mis-matches and highest quality score
- `--strata` - report only those alignemnts that fall into the best stratum if multiple alignments exist, are reportable, and fall into multiple alignment "stratum"
- `-x` - specifies the path to the Bowtie index and the index_prefix used
- `-f` - indicates that the input read file(s) are in fasta format
- `> {sample_id}__mapped.bwt` - re-directs Bowtie standard output to a file; for GeneLab the file prefix is the sample id 

**Input Data:**

- /path/to/collapsed_reads - trimmed and collapsed reads in fasta format, output from [Step 4.](#4-deduplicate-reads)
- Bowtie reference genome index - Bowtie index, output from [Step 5b.](#5b-build-bowtie-index-for-reference-genome)

**Output Data:**

- {sample_id}_mapped.bwt - alignment file in Bowtie standard output format


<br>

---
### 6aii. Convert Alignment File to mirdeep2 (ARF) format

Converts bowtie output to miRDeep2 ARF format.

```bash
convert_bowtie_output.pl {sample_id}_mapped.bwt > {sample_id}_mapped.arf
parse_mappings.pl {sample_id}_mapped.arf -j > {sample_id}_mapped_parsed.arf
```

**Parameter Definitions:**

*convert_bowtie_output.pl*
- `{sample_id}_mapped.bwt` - input alignment file, provided as a positional argument
- `> {sample_id}_mapped.arf` - redirects standard output to an *.arf file

*parse_mappings.pl*
- `{sample_id}_mapped.arf` - alignment file in ARF format, provided as a positional argument
- `-j` - remove any unmatched nucleotides at the 3' end
- `> {sample_id}_mapped_parsed.arf` - redirects standard output to a *parsed.arf file


**Input Data:**

- {sample_id}_mapped.bwt - Alignment file in Bowtie standard output format, output from [Step 6ai.](#6ai-align-reads-to-reference-genome)

**Output Data:**

- {sample_id}_mapped_parsed.arf - parsed alignment file in ARF format for mirdeep2

<br>

---
### 6aiii. Predict Known and Novel miRNAs

```bash
miRDeep2.pl ${sample_id}_collapsed_reads.fa \
${genome}.fa \
${sample_id}_mapped_parsed.arf \
/path/to/mature_ref_miRNA \
/path/to/hairpin_ref_miRNA \
-t ${mirbase_species_tag} \
-d \
-P \
-z _${sample_id} \
2> ${sample_id}_report.log

```

**Parameter Definitions:**

*miRDeep2.pl*
- `${sample_id}_collapsed_reads.fa` -- input trimmed and collapsed reads in fasta format, provided as a positional argument
- `${genome}.fa` -- specifies the genome refrence fasta file, provided as a positional argument
- `${sample_id}_mapped_parsed.arf` -- parsed alignment file in ARF format, provided as a positional argument
- `/path/to/mature_ref_miRNA` -- mature miRNA reference, provided as a positional argument
- `/path/to/hairpin_ref_miRNA` -- hairpin miRNA precursor reference, provided as a positional argument
- `-t [Animal Species]` -- specifies the 3 letter code from miRBase for the species being analyzed
- `-P` -- indicates that mature_ref_miRNA identifiers follow miRBase format >= v18 (include 5p and 3p designation).
- `-d` -- disable pdf generation
- `-z _${sample_id}` -- append sample_id to end of timestamp in output file names
- `2> ${sample_id}_report.log` -- redirects all progress output to a log file


**Input Data:**

- ${sample_id}_collapsed_reads.fa -- trimmed and collapsed reads in fasta format, output from [Step 4](#4-deduplicate-reads)
- ${genome}.fa -- genome sequence, this miRNA-Seq pipeline version uses the Ensembl fasta file indicated in the fasta column of the [GL-DPPD-7110_annotations.csv](https://github.com/nasa/GeneLab_Data_Processing/blob/master/GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110/GL-DPPD-7110_annotations.csv) GeneLab Annotations file
- ${sample_id}_mapped_parsed.arf -- parsed alignment file in ARF format, output from [Step 6aii](#6aii-convert-alignment-file-to-mirdeep2-arf-format)
- mature_ref_miRNA -- miRBase mature miRNA reference for the respective organism [Step 5c](#5c-create-mirbase-species-specific-references)
- hairpin_ref_miRNA -- miRBase hairpin miRNA precursor reference for the respective organism [Step 5c](#5c-create-mirbase-species-specific-references)

**Output Data:**
- dir_prepare_signature_${timestring}	-- Prepared signature files used by miRDeep2			
- error_${timestamp}_${sample_id}.log -- miRDeep2 error log
- expression_${timestamp}_${sample_id}.html -- an overview of all miRNAs in the input data
- expression_analyses	-- folder containing intermediate files from miRDeep2 quantification step 
- **miRNAs_expressed_all_samples_${timestamp}_${sample_id}.csv** -- tab separated file with miRNA identifiers and their read count
- mirdeep_runs/run_${timestamp}_${sample_id}/output.mrd -- miRNA read signatures in text format		
- mirdeep_runs/run_${timestamp}_${sample_id}/survey.csv -- miRDeep2 performance survey in text format
- mirdeep_runs/run_${timestamp}_${sample_id}/run_{timestamp}_{sample_id}_parameters -- miRDeep2 run paramters
- mirna_results_${timestamp}_${sample_id} -- fasta and BED files representing all detected miRNAs, split into separate files by miRDeep score and miRNA type
- **result_${timestamp}_${sample_id}.html** -- an html overview of all detected miRNAs in the deep sequencing input data.
- **result_${timestamp}_${sample_id}.csv** -- a copy of the html overview in tab-separated text format
- **result_${timestamp}_${sample_id}.bed** -- all detected miRNAs in BED format
- ${sample_id}_report.log -- miRDeep2 progress output
-- 
<br>

---
### 6aiv. Extract miRDeep2 counts

Extract read counts from miRDeep2 results csv file

```bash
mkdir results counts
awk -v NAME=${sample_id} -v RS="" '{
  if (NR==1){
    TYPE="mirdeep2_stats"
  } else if(NR==2) {
    TYPE="novel"
  } else if(NR==3) {
    TYPE="known"
  } else {
    TYPE="mirbase_undetected"
  }
  FILE=sprintf("results/%s_%s.csv", NAME, TYPE);
  print > FILE
}' result_*_${sample_id}.csv

# Remove old header, extract counts, and re-header file for novel miRNAs
tail -n +3 results/result_${sample_id}_novel.csv | awk -F\t -v OFS='\t' -v SAMPLE=${sample_id} \
'{
  print SAMPLE"_"$1,$5,SAMPLE"_"$1,$5;
}' > counts/${sample_id}_novel_counts.csv

awk -v OFS="\t" '{print $1,$2,$3}' mirdeep2_${sample_id}_*/miRNAs_expressed_all_samples_*_${sample_id}.csv | cat < counts/${sample_id}_novel_counts.csv  > counts/${sample_id}_all_counts.csv
```

**Parameter Definitions:**

- `${sample_id}` -- the sample ID to include in the output file

**Input Data:**

- `result_${timestamp}_${sample_id}.csv -- result CSV file generated by miRDeep2 [Step 6aiii.](#6aiii-predict-known-and-novel-mirnas)

**Output Data:**

results/
- ${sample_id}_stats.csv -- Prediction result statistics extracted from miRDeep2 results file.
- ${sample_id}_novel.csv -- Prediction results for novel miRNAs extracted from miRDeep2 results file.
- ${sample_id}_known.csv -- Prediction results matching known miRNAs extracted from miRDeep2 results file.
- ${sample_id}_mirbase_undetected.csv -- Prediction results for known miRNAs undetected by miRDeep2 extracted from miRDeep2 results file.

counts/
- ${sample_id}_novel_counts.csv -- Count data for novel miRDeep predicted sequences.
- **${sample_id}_all_counts.csv** -- Count data for known mature miRNAs from miRBase combined with counts for novel miRDeep2 
                                     predicted sequences and known mature miRNAs from miRBase. 3 columns: miRNA ID, total count, precursor ID
                                     miRNA ID for predicted miRNAs consists of source sampleName concatenated with miRDeep2 ID.
-- 

### 6av. Prepare new species-specific reference

```bash
for RESULTDIR in mirdeep2_*/mirna_results_*; do
  sed 's/\([0-9]\)$/\1_s/' ${RESULTDIR}/novel_star_*.fa | cat < ${RESULTDIR}/novel_mature_*.fa  >> novel_predictions.fa;
done;
cat seqkit rmdup -s -P novel_predictions.fa > unique_novel_predictions.fa
cat mature_ref.fa unique_novel_predictions.fa > mature_ref_plus_novel.fa

bowtie-build mature_ref_plus_novel.fa mature_ref_plus_novel
```

**Parameters:**

- `RESULTDIR` -- the directory holding all prediction results, output from [Step 6aiii.](#6aiii-predict-known-and-novel-mirnas)
*seqkit*
- `rmdup` -- select duplicate sequence removal mode
- `-s` -- deduplicate by sequence
- `-P` -- only consider positive strand when comparing sequences
*bowtie-build*
- `mature_ref_plus_novel.fa` --  reference FASTA
- `mature_ref_plus_novel` -- name of bowtie index

**Input Data:**

- novel_star_*.fa -- novel predicted miRNA star sequences
- novel_star_*.fa -- novel predicted miRNA mature sequences

**Output Data:**

- `novel_predictions.fa` -- FASTA file containing all predicted mature and star miRNA sequences not already found in reference database.
- **mature_ref_plus_novel.fa** - FASTA file containing all known mature and novel predicted mature and star miRNA sequences
*bowtie index files for new reference*
- mature_ref_plus_novel.1.ebwt 
- mature_ref_plus_novel.2.ebwt 
- mature_ref_plus_novel.3.ebwt 
- mature_ref_plus_novel.4.ebwt 
- mature_ref_plus_novel.rev.1.ebwt 
- mature_ref_plus_novel.rev.2.ebwt 
---

<br>

---
### 6b. Plant miRNA prediction
### 6bi. Predict known and novel miRNAs for plants

```bash
miRDP2-v1.1.4_pipeline.bash \
  --genome {genome}.fa \
  --index {index_prefix} \
  --fasta \
  --input {sample_id}_collapsed_reads.fa \
  --output {output_directory} \
  --thread 1
```

**Parameter Definitions:**

- `--thread` -- specify the number of threads available to perform alignment
- `--genome {genome}.fa` -- the genome reference fasta file
- `--index {index_prefix}` -- Bowtie index prefix for the 
- `--input {sample_id}_collapsed_reads.fa` -- input read file
- `--fasta` -- flag to indicate that the sample file is in fasta format

**Input Data:**

- `{sample_id}_collapsed_reads.fa` -- (trimmed and collapsed reads in fasta format, output from [Step 4](#4-deduplicate-reads))
- `{index_prefix}` -- (Bowtie index, output from [Step 4](#5b-build-bowtie-index-for-reference-genome)
- `{genome}.fa` -- (genome sequence, this miRNA-Seq pipeline version uses the Ensembl fasta file indicated in the fasta column of the [GL-DPPD-7110_annotations.csv](https://github.com/nasa/GeneLab_Data_Processing/blob/master/GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110/GL-DPPD-7110_annotations.csv) GeneLab Annotations file)

**Output Data:**

- **{output_directory}/{sample_id}/{sample_id}_filter_P_prediction** - A tab delimited file containing information on final predicted miRNAs. Columns are: *chromosome id, strand direction, representative reads id, precursor id, mature miRNA location, precursor location, mature sequence, and precursor sequence*.
- **{output_directory}/{sample_id}/{sample_id}.bed** - a BED file of the miRNA predictions


### 6bii. Convert miRNA predictions to fasta

Convert predictions to fasta
```bash
awk '{print ">"$4"_"$3"\n"$7}' *_filter_P_prediction | seqkit rmdup -s -P > all_samples_filter_P_prediction.fa
```

**Parameter Definitions:**
- `*_filter_P_prediction` -- list of miRDeep-P2 prediction outputs, provided as a positional argument
- `all_samples_filter_P_prediction.fa` -- output filename, provided as a positional argument

**Input Data:**
- `{*_filter_P_prediction}` -- a list of file paths to the *filter_P_prediction files produced by [Step 6bi.](#6bi-predict-known-and-novel-mirnas-for-plants)


**Output Data:**
- **all_samples_filter_P_prediction.fa** - fasta file containing unique predicted miRNA mature sequences. The sequence names are derived from the miRDeep-P2 sequence identifer and the collapsed read name from the input fasta files.


### 6biii. BLAST predicted to known mature sequences

```bash
makeblastdb -in ${mature_ref.fa} -dbtype nucl

blastn -task blastn-short -db ${mature_ref.fa} -query all_samples_filter_P_prediction.fa\
  -out predictions_v_mature_blastn.tsv -strand plus\
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs"
```

**Parameter Definitions:**

*makeblastdb*
- `in` -- reference fasta
- `-dbtype` -- sequence type (nucl: DNA/RNA nucleotide sequences)

*blastn*
- `-db` -- BLAST database
- `-query` -- query sequences in fasta format
- `-out` -- output filename
- `-word_size` -- Length of initial exact match. Set to shorter size for short sequences.
- `-num_alignments` -- Maximum number of alignments reported.

**Input**
- `mature_ref.fa` -- species specific mature miRNA from miRBase from [Step 5c.](#5c-create-mirbase-species-specific-references)
- `all_samples_filter_P_prediction.fa` -- mirdeep-P2 mature miRNA predictions from [Step 6bii](#6bii-convert-mirna-predictions-to-fasta)

**Output**
- `predictions_v_mature_blastn.tsv` -- BLAST output

### 6biv. Prepare new species-specific reference

```bash
awk '{if ($NF>=90) {print $1} }' predictions_v_mature_blastn.tsv | sort --unique > known_ids.txt
seqkit grep --pattern-file known_ids.txt --invert-match all_samples_filter_P_prediction.fa > novel_predictions.fa

cat mature_ref.fa novel_prediction.fasta > mature_ref_plus_novel.fa

bowtie-build mature_ref_plus_novel.fa mature_ref_plus_novel
```

**Parameter Definitions:**

*sort*
- `--unique` -- return only unique entries
*seqkit*
- `grep` -- seqkit tool specification
- `--pattern-file` -- file containing IDs to search for (1 ID per line)
- `--invert-match` -- select the sequences whose IDs are not in the provided pattern file.

*bowtie-build*
- `mature_ref_plus_novel.fa` --  reference FASTA
- `mature_ref_plus_novel` -- name of bowtie index

**Input**
- `predictions_mature_blastn.tsv` -- BLAST output file
- `all_samples_filter_P_prediction.fa` -- fasta file of all mirDeep-P2 predicted miRNA from [Step 6bi.](#6bi-predict-known-and-novel-mirnas-for-plants)
- `mature_ref.fa` -- known mature sequences from miRBase [Step 5c](#5c-create-mirbase-species-specific-references)

**Output**
- `novel_predictions.fa` -- FASTA file containing all predicted mature miRNA sequences not already found in reference database.
- **mature_ref_plus_novel.fa** - FASTA file containing all known mature and novel predicted mature and star miRNA sequences
*bowtie index files for new reference*
- mature_ref_plus_novel.1.ebwt 
- mature_ref_plus_novel.2.ebwt 
- mature_ref_plus_novel.3.ebwt 
- mature_ref_plus_novel.4.ebwt 
- mature_ref_plus_novel.rev.1.ebwt 
- mature_ref_plus_novel.rev.2.ebwt 


## 7. Generate counts for known and predicted miRNA

### 7a. Align reads to specifies-specific reference

```bash


bowtie --threads ${num_threads} \
  --seed 12345 \
  --seedlen 16 \
  --seedmms 0 \
  --maqerr 80 \
  --all \
  -k 100 \
  --best \
  --strata \
  --sam \
  --un ${sample_id}_unmapped.fq \
  -x ${mature_ref_plus_novel} \
  -q <(zcat < ${sample}.fq.gz) \
  --sam > ${sample_id}_mapped.sam
  
samtools sort --write-index -o ${sample_id}_mapped.bam##idx##${sample_id}_mapped.bam.bai ${sample_id}_mapped.sam
```

**Parameter Definitions:**

*bowtie*
- `--threads` -- specify the number of threads available to perform alignment
- `--seed` -- specify the seed to use for pseudo-random number generator to ensure reproducibility ## test with and without setting the seed to make sure there is no bias ##
- `--seedlen` -- specifies the seed length, 16 is used due to the shorter length of miRNA reads
- `--seedmms` -- specifies the number of mismatches allowed in the seed, a value of 0 means no mismatches are allowed in the seed
- `--maqerr` -- specifies the maximum sum of quality values at each mismatch position. When input files are in fasta format, bowtie assigns a default value of 40 for the quality. With a quality value of 40 for each base (the default when input is a fasta file), "-e 80" allows up to 2 mismatches to occur after the seed region of a read alignment.
- `--all` -- Report all valid alignments per read
- `-m` -- suppresses all alignments for a read if more than the indicated number of reportable alignments exist for it
- `--best` -- guarantees that reported singleton alignments are the "best" in terms of stratum, i.e. fewest number of mis-matches and highest quality score
- `--strata` -- report only those alignemnts that fall into the best stratum if multiple alignments exist, are reportable, and fall into multiple alignment "stratum"
- `--sam` -- output alignments in SAM format
- `--un` -- output unmapped reads to a fastq file
- `-x` -- specifies the path to the Bowtie index and the index_prefix used
- `-q` -- indicates that the input read file(s) are in fastq format
- `> {sample_id}_mapped.bwt` -- re-directs Bowtie standard output to a file; for GeneLab the file prefix is the sample id 

*samtools sort*
- `--write-index` -- automatically index the output file
- `-o` -- Specify the output file
- `${sample_id}_mapped.bam##idx##${sample_id}_mapped.bam.bai` -- Indicates the BAM filename followed by the name of the BAM index, separated by "##idx##"

**Input Data:**
- `${sample}.fq.gz` -- Trimmed fastq file from [Step 2a](#2a-trimfilter-raw-data)
- `${mature_ref_plus_novel}` -- Reference miRNA mature sequences from [Step 6av.](#6av-prepare-new-species-specific-reference) or [Step 6b](#6biv-prepare-new-species-specific-reference)

**Output Data:**
- `${sample_id}_mapped.bam` -- Sorted BAM output file
- `${sample_id}_mapped.bam.bai` -- BAM index
  

### 7b. Generate basic alignment statistics and read counts

idxstats provides counts, normalized counts, and total counts for each contig in a BAM file.


```bash
samtools idxstat ${sample_id}_mapped.bam > ${sample_id}_idxstat.txt
echo -e "mirname\tlength\tmapped\tunmapped" > header.txt
cat header.txt ${sample_id}_idxstat.txt | grep -v "^\*" > ${sample_id}_counts.txt
rm header.txt

samtools stats ${sample_id}_mapped.bam > ${sample_id}_stats.txt

samtools flagstat ${sample_id}_mapped.bam > ${sample_id}_flagstats.txt

```

**Parameter Definitions:**
*samtools*
  - 'stats' -- produces comprehensive statistics from alignment file, provided as a positional argument
  - 'idxstats' -- reports alignment summary statistics, provided as a positional argument
  - 'flagstat' -- counts the number of alignments for each FLAG type, provided as a positional argument
*echo*
  - `-e` -- enable interpretation of backslash escapes (to correctly interpret '\t' as tab character)

**Input Data:**
  - `${sample_id}_mapped.bam` files from [Step 7a](#7a-align-reads-to-specifies-specific-reference)

**Output Data:**
  - **${sample_id}_flagstat.txt** - SAM flag alignment statistics
  - **${sample_id}_stats.txt** - Comprehensive alignment statistics produced by `samtools stats`. For additional information, see the [samtools documentation](http://www.htslib.org/doc/samtools-stats.html)
  - **${sample_id}_idxstat.txt** - A tab delimited file containing reference sequence alignment statistics consisting of the following columns: sequence name, sequence length, number of mapped reads, number of unmapped reads. Note that this will count reads multiple times if they are mapped more than once. 
- **${sample_id}_counts.txt** - Adds a header to the idxstat file and removes the last row of the file containing counts of unmapped reads.


--- 

### 7c. Compile alignment QC
**Code:**

```bash
multiqc --interactive --filename alignment_stats_multiqc_report.html \
--outdir ${output_directory} --zip-data-dir \
/path/to/alignment_stats/
```
**Parameter Definitions:**

- `--interactive` -- force reports to use interactive plots
- `--filename` -- report filename
- `--outdir` -– the output directory to store results
- `--zip-data-dir` -- compress the data directory
- `/path/to/alignment_stats/` -– the directory holding the output data from the samtools statistics files, provided as a positional argument

**Input Data:**

- *.txt -- Samtools alignment statistics, output from [Step 7b.](#7b-generate-basic-alignment-statistics-and-read-counts)

**Output Data:**
- **alignment_stats_multiqc_report.html** -- multiqc report
- **alignment_stats_multiqc_report_data.zip** -- zip file of directory containing multiqc data


<br>

---
### 7d. Merge read counts

**Code:**

```bash 
paste /path/to/alignment_stats/*_counts.txt | awk -v OFS=',' '{OUTL=$1;for(i=0;i<NF;i+=4){c=i+3; OUTL=OUTL","$c}print OUTL}' > all_counts.txt
```

***Parameter Definition:**

- `/path/to/alignment_stats/` - the directory holding the output data from the samtools statistics files, provided as a positional argument

**Input Data:**

- *_counts.txt - Counts file for each sample, output from [Step 7b.](#7b-generate-basic-alignment-statistics-and-read-counts)

**Output Data:**
- **all_counts.txt** -- Comma-separated file containing read counts for each sample.

<br>

## 8. Differential expression analysis

Perform differential expression analysis using DESeq2. Adapted from GeneLab bulk RNASeq DESeq2 R-script.

## 8a. Environment setup

```R
### Install R packages if not already installed ###

install.packages("tidyverse")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.14")
BiocManager::install("DESeq2")
BiocManager::install("mirbase.db")

### Import libraries (tximport, DESeq2, tidyverse, stringr) ###

library(DESeq2)
library(tidyverse)
library(mirbase.db)

analysis_dir <- "/path/to/analysis/location"
metadata_file <- file.path(analysis_dir, "metadata.csv")
all_cts_file <- "/path/to/alignment/stats/output/all_counts.txt"
work_dir <- analysis_dir
counts_dir <- "/path/to/alignment/stats/output"
norm_output <- file.path(analysis_dir, "NormCounts")
dge_output <-file.path(analysis_dir, "DGEOutput")
interaction_output <- file.path(analysis_dir, "TargetInteractions")
setwd(file.path(work_dir))
```

### 8b. Configure Metadata, Sample Grouping, and Group Comparisons

```R
### Pull all factors for each sample in the study from the metadata in the ISA files ###
study <- read.csv(metadata_file, header = TRUE, row.names = 1, stringsAsFactors = TRUE)

##### Format groups and indicate the group that each sample belongs to #####
if (ncol(study) >= 2) {
  group <- apply(study, 1, paste, collapse = " & ") # concatenate multiple factors into one condition per sample
} else {
  group <- study[, 1]
}
group_names <- paste0("(", group, ")", sep = "") # human readable group names
group <- make.names(group) # group naming compatible with R models
names(group) <- group_names
rm(group_names)

##### Format contrasts table, defining pairwise comparisons for all groups #####

# generate matrix of pairwise group combinations for comparison
contrast_names <- combn(levels(factor(names(group))), 2)
# limited make.names call for each group (also removes leading parentheses)
contrasts <- apply(contrast_names, MARGIN = 2,
  function(col) sub("^BLOCKER_", "", make.names(paste0("BLOCKER_", stringr::str_sub(col, 2, -2))))
)
# format combinations for output table files names
contrast_names <- c(paste(contrast_names[1, ], contrast_names[2, ], sep = "v"),
                    paste(contrast_names[2, ], contrast_names[1, ], sep = "v"))
contrasts <- cbind(contrasts, contrasts[c(2, 1), ])
colnames(contrasts) <- contrast_names
rm(contrast_names)

```

### 8c. Import miRNA counts

```R

# Import counts data as matrix
read_counts <- as.matrix(read.csv(all_cts_file, sep = ",", row.names = "mirname"))

## Create data frame defining which group each sample belongs to
sample_table <- data.frame(condition = factor(group))
rownames(sample_table) <- rownames(study)
# reorder counts columns to match the ordering of samples in the metadata
read_counts <- read_counts[, rownames(sample_table)]

```

### 8d. Perform DGE
```R
# make DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = read_counts,
                              colData = sample_table,
                              design = ~condition)

# filter out miRNAs with counts of less than 10 in all conditions
keep <- rowSums(counts(dds)) > 10
dds <- dds[keep, ]
summary(dds)
dim(dds)

#### Perform DESeq analysis
dds <- DESeq(dds)
```

### 8e. Prepare GeneLab DGE Tables with Annotations

```R
##### Generate F statistic p-value (similar to ANOVA p-value) using DESeq2 likelihood ratio test (LRT) design #####
dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ 1)
res_lrt <- results(dds_lrt)

### Create a data frames containing normalized counts ###

norm_counts <- as.data.frame(counts(dds, normalized = TRUE))

## Add 1 to all counts to avoid issues with log transformation
norm_counts <- norm_counts + 1


## output table 1 will be used to generate computer-readable DGE table,
## which is used to create GeneLab visualization plots
output_table_1 <- norm_counts

## reduced output table 1 will be used to generate human-readable DGE table
reduced_output_table_1 <- norm_counts

## Iterate through Wald Tests to generate pairwise comparisons of all groups
for (i in seq_len(dim(contrasts)[2])) {
  res_1 <- results(dds, contrast = c("condition", contrasts[1, i], contrasts[2, i]))
  res_1 <- as.data.frame(res_1@listData)[, c(2, 4, 5, 6)]
  colnames(res_1) <- c(paste0("Log2fc_", colnames(contrasts)[i]),
                       paste0("Stat_", colnames(contrasts)[i]),
                       paste0("P.value_", colnames(contrasts)[i]),
                       paste0("Adj.p.value_", colnames(contrasts)[i]))
  output_table_1 <- cbind(output_table_1, res_1)
  reduced_output_table_1 <- cbind(reduced_output_table_1, res_1)
  rm(res_1)
}

## Generate and add all sample mean column to the normalized counts table

output_table_1$All.mean <- rowMeans(norm_counts, na.rm = TRUE, dims = 1)
reduced_output_table_1$All.mean <- rowMeans(norm_counts, na.rm = TRUE, dims = 1)

## Generate and add all sample stdev column to the normalized counts table
output_table_1$All.stdev <- rowSds(as.matrix(norm_counts), na.rm = TRUE, dims = 1)
reduced_output_table_1$All.stdev <- rowSds(as.matrix(norm_counts), na.rm = TRUE, dims = 1)

## Add F statistic p-value (similar to ANOVA p-value) column to the normalized counts table
output_table_1$LRT.p.value <- res_lrt@listData$padj
reduced_output_table_1$LRT.p.value <- res_lrt@listData$padj

## Generate and add group mean and stdev columns to the normalized counts table
tcounts <- as.data.frame(t(norm_counts))
## Use final table group name formatting (e.g. '( Space Flight & Blue Light )' )
tcounts$group <- names(group)

# Compute group name group-wise means
group_means <- as.data.frame(t(aggregate(. ~ group, data = tcounts, mean)))
# assign group name as column names
colnames(group_means) <- paste0("Group.Mean_", group_means["group", ])

# Compute group name group-wise standard deviation
group_stdev <- as.data.frame(t(aggregate(. ~ group, data = tcounts, sd)))
# assign group name as column names
colnames(group_stdev) <- paste0("Group.Stdev_", group_stdev["group", ])

# Drop group name row from data rows (now present as column names)
group_means <- group_means[-c(1), ]
group_stdev <- group_stdev[-c(1), ]

# Column bind the group-wise data
output_table_1 <- cbind(output_table_1, group_means, group_stdev)
reduced_output_table_1 <- cbind(reduced_output_table_1, group_means, group_stdev)

rm(group_stdev, group_means, tcounts)

### Add columns needed to generate GeneLab visulaization plots to the normalized counts table

## Add column to indicate the sign (positive/negative) of log2fc for each pairwise comparison
updown_table <- sign(output_table_1[, grep("Log2fc_", colnames(output_table_1))])
colnames(updown_table) <- gsub("Log2fc", "Updown", grep("Log2fc_", colnames(output_table_1), value = TRUE))
output_table_1 <- cbind(output_table_1, updown_table)
rm(updown_table)

## Add column to indicate contrast significance with p <= 0.1
sig.1_table <- output_table_1[, grep("P.value_", colnames(output_table_1))] <= .1
colnames(sig.1_table) <- gsub("P.value", "Sig.1", grep("P.value_", colnames(output_table_1), value = TRUE))
output_table_1 <- cbind(output_table_1, sig.1_table)
rm(sig.1_table)

## Add column to indicate contrast significance with p <= 0.05
sig.05_table <- output_table_1[, grep("P.value_",colnames(output_table_1))] <= .05
colnames(sig.05_table) <- gsub("P.value", "Sig.05",grep("P.value_", colnames(output_table_1), value = TRUE))
output_table_1 <- cbind(output_table_1, sig.05_table)
rm(sig.05_table)

## Add columns for the volcano plot with p-value and adjusted p-value
log_pval_table <- log2(output_table_1[, grep("P.value_", colnames(output_table_1))])
colnames(log_pval_table) <- paste0("Log2_", colnames(log_pval_table))
output_table_1 <- cbind(output_table_1, log_pval_table)
rm(log_pval_table)
log_adj_pval_table <- log2(output_table_1[, grep("Adj.p.value_", colnames(output_table_1))])
colnames(log_adj_pval_table) <- paste0("Log2_", colnames(log_adj_pval_table))
output_table_1 <- cbind(output_table_1, log_adj_pval_table)
rm(log_adj_pval_table)

## Prepare PCA table for GeneLab visualization plots ##

exp_raw <- log2(norm_counts)
PCA_raw <- prcomp(t(exp_raw), scale = FALSE)

```

### 8f. Export GeneLab DGE Tables with Annotations

```r
# add mirbase.db annotation columns for microRNA family and description
mir_desc <- as.data.frame(toTable(mirbaseDESCRIPTION))
mir_fam <- as.data.frame(toTable(mirbaseFAMILY))
mir_mature <- as.data.frame(toTable(mirbaseMATURE))
mir_annot <- merge(mir_mature, mir_fam, by = "mirna_id", all.x = TRUE)
mir_annot <- merge(mir_annot, mir_desc, by = "mirna_id", all.x = TRUE)


### Combine annotations table and the DGE tables
output_table_1 <- merge(mir_annot, output_table_1, by.x = "mature_name", by.y = "row.names", all.y = TRUE)
reduced_output_table_1 <- merge(mir_annot, reduced_output_table_1, by.x = "mature_name", by.y = "row.names", all.y = TRUE)

### Export counts tables ###

norm_counts_exp <- as.data.frame(counts(dds, normalized = TRUE))

write.csv(read_counts, file.path(norm_output, "Unnormalized_Counts.csv"))
write.csv(norm_counts_exp, file.path(norm_output, "Normalized_Counts.csv"))


### Export sample grouping and contrasts tables for normalized data ###

write.csv(sample_table, file.path(dge_output, "SampleTable.csv"))
write.csv(contrasts, file.path(dge_output, "contrasts.csv"))


### Export human-readable normalized DGE tables ###

write.csv(reduced_output_table_1, file.path(dge_output, "differential_expression.csv"), row.names = FALSE)

### Export computer-readable DGE and PCA tables used for GeneLab visualization ###

write.csv(output_table_1, file.path(dge_output, "visualization_output_table.csv"), row.names = FALSE)
write.csv(PCA_raw$x, file.path(dge_output, "visualization_PCA_table.csv"), row.names = TRUE)

### print session info ###

print("Session Info below: ")
sessionInfo()
```
	
**Input Data:**

 - all_counts.txt -- Counts from [Step 7d.](#7d-merge-read-counts)
 - runsheet.txt -- table containing metadata required for analysis

**Output Data:**

 - **Unnormalized_counts.csv** -- table containing raw miRNA counts for each sample
 - **Normalized_counts.csv** -- table containing normalized gene counts for each sample
 - **SampleTable.csv** -- table containing samples and their respective groups
 - **visualization_output_table.csv** -- file used to generate GeneLab DGE visualizations
 - **visualization_PCA_table.csv** -- file used to generate GeneLab PCA plots
 - **differential_expression.csv** -- table containing normalized counts for each sample, group statistics, DESeq2 DGE results for each pairwise comparison, and miRNA annotations
 - **contrasts.csv** -- table containing all pairwise comparisons
<br>
