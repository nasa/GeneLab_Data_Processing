# QC Metrics Summary Specification

## Description

The qc_metrics summary file is a comma-separated file that lists a summary of metadata and qc 
metrics for a dataset created after processing is complete. [Metadata](#metadata-fields) for the 
sample and assay are pulled from the ISA.zip file. QC metrics include 
[raw and trimmed read metrics](#read-qc-metrics), [alignment metrics](#alignment-metrics), 
[gene count metrics](#gene-count-metrics), and [RSeQC metrics](#rseqc-metrics) are pulled from the 
MultiQC reports generated during processing.

> **NOTE:** Eukaryotic and prokaryotic data use different tools for alignment and gene counting, so 
the metrics reported will differ for those data types. Similarly, paired-end data will include 
metrics for the reverse read, while single-end data will not. All data columns will be present in 
the output files regardless of data type. Any fields that are not relevant for a particular data 
type will be left empty.
 
<br>
<br>

-----------
### Metadata fields
*Source: ISA.zip*   
List selected metadata fields for each sample

| Column Name           | Data Source         | Description                                                                                               | Example(s)                               |
|:----------------------|:--------------------|:----------------------------------------------------------------------------------------------------------|:-----------------------------------------|
| osd_num               | investigation table | Study identifier                                                                                          | "OSD-515"                                |
| sample                | sample table        | Sample Identifier                                                                                         | "RR23_TMS_FLT_F1"                        |
| organism              | sample table        | Genus and Species of primary organism under study                                                         | "Mus musculus"                           |
| tissue                | sample table        | 'Material Type' listed in sample table (e.g. organ, tissue, cell culture, whole organism)                 | "Cells, cultured", "thymus"              |
| sequencing_instrument | assay table         | Sequencer used to sequence the samples                                                                    | "Illumina NovaSeq 6000"                  |
| library_selection     | assay table         | Biomolecule selection (rRNA depletion, mRNA enrichment, target amplification)                             | "Ribo-depletion", "polyA enrichment"     |
| library_layout        | assay table         | RNA library construction parameter indicating whether the library was sequenced single- or paired-end.    | "SINGLE", "PAIRED"                       |
| strandedness          | assay table         | RNAseq library construction parameter indicating whether the library preserves the transcript orientation | "STRANDED", "UNSTRANDED"                 |
| read_depth            | assay table         | Total number of sequenced fragments                                                                       | 82186286                                 |
| read_length           | assay table         | Raw data read length in bases                                                                             | 151                                      |
| rrna_contamination    | assay table         | Percent rRNA found in the data                                                                            | 2.74                                     |
| rin                   | assay table         | RNA integrity number                                                                                      | 7.4                                      |
| organism_part         | sample table        | For plant samples, the part of the plant from which the sample was derived.                               | "Plant Roots"                            |
| cell_line             | sample table        | Cell line identifier                                                                                      | "IMR90 hiPSC"                            |
| cell_type             | sample table        | The morphological or functional form of the cells used.                                                   | "Myocytes, Cardiac"                      |
| secondary_organism    | sample table        | Additional organism(s) which may confound assay measurements (e.g. food source, host, symbiote)           | "Vibrio fischeri ES114"                  |
| strain                | sample table        | Strain, breed, ecotype, etc.                                                                              | "C57BL/6J"                               |
| animal_source         | sample table        | The source of the animal(s) used in the study                                                             | "Jackson Laboratory"                     |
| seed_source           | sample table        | The source of the seed(s) used in the study                                                               | "Arabidopsis Biological Resource Center" |
| source_accession      | sample table        | The source accession number for the animal(s)/seed(s)/cell(s) used in the study                           | "SALK_027956C"                           |
| mix                   | assay table         | The ERCC spike-in mix number used                                                                         | "Mix 1"                                  |

<br>
<br>

-----------
### Read QC metrics
*Source: Raw and Trimmed MultiQC*   
QC metrics describing the read quality.
> *NOTE:* the same fields are extracted for both raw and trimmed reads, the prefixes "raw_" and 
"trimmed_" indicate the source of the metric. Similarly, the same fields are also extracted for both 
forward and reverse reads with the suffixes "_f" and "_r" denoting the read type. The table below 
lists each metric only once.

| Column Name            | Type  | Description                                                                                                       | Example(s)  |
|:-----------------------|:------|:------------------------------------------------------------------------------------------------------------------|:------------|
| total_sequences        | int   | Total number of sequences for this read type                                                                      | 82186286    |
| avg_sequence_length    | float | Average sequence length for this read type                                                                        | 138.7187431 |
| median_sequence_length | float | Median sequence length for this read type                                                                         | 151         |
| quality_score_mean     | float | Average quality score                                                                                             | 35.48372691 |
| quality_score_median   | float | Median quality score                                                                                              | 35.61967608 |
| percent_duplicates     | float | Percentage estimated sequence duplication, measured based on first 50bp of the first 100,000 reads in the dataset | 67.01577646 |
| percent_gc             | float | Overall %GC of all bases in all sequences for this read type                                                      | 48          |
| gc_min_1pct            | int   | Minimum %GC value reached by at least 1% of the total reads of this read type                                     | 31          |
| gc_max_1pct            | int   | Maximum %GC value reached by at least 1% of the total reads of this read type                                     | 65          |
| gc_auc_25pct           | int   | %GC value at the 25 quartile of total reads of this read                                                          | 42          |
| gc_auc_50pct           | int   | %GC value at the 50th quartile of total reads of this read type                                                   | 49          |
| gc_auc_75pct           | int   | %GC value at the 75th quartile of total reads of this read type                                                   | 57          |
| n_content_sum          | float | %N base calls summed across all base positions in the read for this read type                                     | 1.510487284 |

<br>
<br>

-----------
### Alignment metrics

#### STAR alignment metrics 
*Source: STAR MultiQC*  
Alignment metrics generated by STAR.
> *Present only for data processed with the eukaryotic pipeline.*

| Column Name                 | Type  | Description                                                                         | Example |
|:----------------------------|:------|:------------------------------------------------------------------------------------|:--------|
| uniquely_mapped_percent     | float | Percent of reads mapped uniquely in the genome                                      | 79.95   |
| multimapped_percent         | float | Percent of reads mapped to multiple loci in the genome                              | 13.14   |
| multimapped_toomany_percent | float | Percent of reads mapped to 20 or more loci in the genome                            | 0.13	  |
| unmapped_tooshort_percent   | float | Percent of reads where the best alignment is shorter than the allowed mapped length | 6.53	  |  
| unmapped_other_percent      | float | Percent of reads that couldn't be mapped at all                                     | 0.25    |  

#### Bowtie2 alignment metrics 
*Source: bowtie2 MultiQC*   
Alignment metrics generated by bowtie2.
> *Present only for data processed with the prokaryotic pipeline.*

| Column Name            | Type  | Description                                                                                                      | Example  |
|:-----------------------|:------|:-----------------------------------------------------------------------------------------------------------------|:---------|
| total_reads            | int   | total input reads (or read-pairs for paired-end data)                                                            | 15066949 |
| overall_alignment_rate | float | percentage of input reads that mapped to the genome (includes discordant and partial alignments).                | 98.03    |
| aligned_none           | int   | number of reads (or read-pairs) that aligned 0 times to the reference genome (concordantly if paired-end)        | 516325   |
| aligned_one            | int   | number of reads (or read-pairs) that aligned exactly 1 time to the reference genome (concordantly if paired-end) | 11294617 |
| aligned_multi          | int   | number of reads (or read-pairs) that aligned > 1 times to the reference genome (concordantly if paired-end)      | 3256007  |

<br>
<br>

-----------
### Gene count metrics

#### Common gene count metrics
*Source: RSEM_Unnormalized_Counts_GLbulkRNAseq.csv or FeatureCounts_Unnormalized_Counts_GLbulkRNAseq.csv*  
Summarizes the data in the Unnormalized_Counts_GLbulkRNASeq.csv files.
> *Present for data processed with either the eukaryotic or prokaryotic pipelines.*

| Column Name            | Type  | Description                                       | Example    |
|:-----------------------|:------|:--------------------------------------------------|:-----------|
| gene_total             | int   | total number of genes detected                    | 57186      |
| gene_detected_gt10     | int   | number of genes detected with > 10 read depth     | 19207      |
| gene_detected_gt10_pct | float | percentage of genes detected with > 10 read depth | 33.5868919 |


#### RSEM gene count metrics
*Source: RSEM MultiQC*  
Summarizes the alignment rates generated by RSEM during gene expression estimation.
> *Present only for data processed with the eukaryotic pipeline.*

| Column Name            | Type  | Description                                             | Example     |
|:-----------------------|:------|:--------------------------------------------------------|:------------|
| num_uniquely_aligned   | int   | number of reads aligned uniquely to a gene              | 41152050    |
| pct_uniquely_aligned   | float | percentage of reads aligned unique to a gene            | 74.48868329 |
| pct_multi_aligned      | float | percentage of reads aligned to multiple genes           | 15.27982918 |
| pct_filtered           | float | percentage of reads filtered due to too many alignments | 0           |
| pct_unalignable        | float | percentage of reads unalignable to any gene             | 10.23148753 |

#### FeatureCounts gene count metrics
*Source: FeatureCounts MultiQC*  
Summarizes the feature counts assigned to genome features.
> *Present only for data processed with the prokaryotic pipeline.*

| Column Name               | Type  | Description                                              | Example     |
|:--------------------------|:------|:---------------------------------------------------------|:------------|
| total_count               | int   | total number of input alignments                         | 496         |
| num_assigned              | int   | number of alignments assigned to a feature               | 431         |
| pct_assigned              | float | percentage of alignments assigned to a feature           | 86.89516129 |
| num_unassigned_nofeatures | int   | number of alignments that overlap no features            | 19          |
| num_unassigned_ambiguity  | int   | number of alignments that overlap 2 or more features     | 43          |
| pct_unassigned_nofeatures | float | percentage of alignments that overlap no features        | 3.830645161 |	
| pct_unassigned_ambiguity  | float | percentage of alignments that overlap 2 or more features | 8.669354839 |

<br>
<br>

-----------
### RSeQC metrics

#### Genebody coverage metrics
*Source: RSeQC Gene Body Coverage MultiQC*  
Summarizes the read coverage over gene bodies. Used to check if the coverage is uniform and if 
there is any 5' or 3' bias.

| Column Name               | Type  | Description                                                                          | Example     |
|:--------------------------|:------|:-------------------------------------------------------------------------------------|:------------|
| mean_genebody_cov_5_20    | float | average read coverage between the 5th and 20th gene body percentile from the 5' end         | 70.60419705 |
| mean_genebody_cov_40_60   | float | average read coverage between the 40th and 60th gene body percentile from the 5' end        | 99.38086546 |
| mean_genebody_cov_80_95   | float | average read coverage between the 80th and 95th gene body percentile from the 5' end        | 87.1712315  |
| ratio_genebody_cov_3_to_5 | float | ratio of the 3' (mean_genebody_cov_80_95) to 5' (mean_genebody_cov_5_20) gene body coverage | 1.234646595 |

#### Infer experiment metrics
*Source: RSeQC Infer Experiment MultiQC*    
Summarizes the percentage of reads and read pairs that match the strandedness of the overlapping 
transcripts. Can be used to infer whether the RNAseq library prep was stranded or unstranded.

| Column Name               | Type  | Description                                                  | Example |
|:--------------------------|:------|:-------------------------------------------------------------|:--------|
| pct_sense                 | float | percentage of reads aligned to the sense strand              | 91.16   |
| pct_antisense             | float | percentage of reads aligned to the antisense strand          | 1.76    |
| pct_undetermined          | float | percentage of reads where the strand could not be determined | 7.08    |

#### Inner distance metrics
*Source: RSeQC Inner Distance MultiQC*  
Summarizes the inner distance (or insert size) between two paired RNA reads. Note that this can be 
negative if the two reads overlap.
> *Present only for paired-end data.*

| Column Name               | Type  | Description                                         | Example     |
|:--------------------------|:------|:----------------------------------------------------|:------------|
| peak_inner_dist           | float | The inner distance at the peak of the distribution  | -123.5      |
| peak_inner_dist_pct_reads | float | percentage of reads at the peak of the distribution | 5.112384053 |

#### Read distribution metrics
*Source: RSeQC Read Distribution MultiQC*   
Summarizes the distribution of reads over genome features.

| Column Name               | Type  | Description                                                                | Example.    |
|:--------------------------|:------|:---------------------------------------------------------------------------|:------------|
| cds_exons_pct             | float | percentage of reads aligned to CDS exons                                   | 43.5637117  |
| 5_utr_exons_pct           | float | percentage of reads aligned to 5' UTRs                                     | 8.886329966 |
| 3_utr_exons_pct           | float | percentage of reads aligned to 3' UTRs                                     | 21.41150152 |
| introns_pct               | float | percentage of reads aligned to introns                                     | 20.55254242 |
| tss_up_1kb_pct            | float | percentage of reads aligned 1 kb upstream of a transcription start site    | 0.086612962 |
| tss_up_1kb_5kb_pct        | float | percentage of reads aligned 1-5 kb upstream of a transcription start site  | 0.389163859 |
| tss_up_5kb_10kb_pct       | float | percentage of reads aligned 5-10 kb upstream of a transcription start site | 0.116168783 |
| tes_down_1kb_pct          | float | percentage of reads aligned 1 kb downstream of a transcription end site    | 0.276327172 |
| tss_down_1kb_5kb_pct      | float | percentage of reads aligned 1-5 kb downstream of a transcription end site  | 0.419487364 |
| tss_down_5kb_10kb_pct     | float | percentage of reads aligned 5-10 kb downstream of a transcription end site | 0.137565531 |
| other_intergenic_pct      | float | percentage of reads aligned to other intergenic regions                    | 4.160588723 |

