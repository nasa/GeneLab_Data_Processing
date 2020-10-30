####################################################################################
## Command-line processing file for 16S and ITS data of NASA GLDS-126             ##
## https://genelab-data.ndc.nasa.gov/genelab/accession/GLDS-126/                  ##
##                                                                                ##
## This code as written expects to be run within the "processing_info/" directory ##
## with the raw starting fastq files present in a "Raw_Data/"" directory          ##
## Processed by Michael D. Lee (Mike.Lee@nasa.gov)                                ##
####################################################################################

fastqc -v # v0.11.8
multiqc -v # v1.7
cutadapt --version # 2.3


  # running fastqc on raw data
fastqc ../Raw_Data/*.gz -t 15
  # summarizing fastqc results with multiqc
multiqc -z -o ../FastQC_Outputs/ ../Raw_Data/
  # renaming the outputs
mv ../FastQC_Outputs/multiqc_data.zip ../FastQC_Outputs/raw_multiqc_data.zip
mv ../FastQC_Outputs/multiqc_report.html ../FastQC_Outputs/raw_multiqc_report.html
  # removing the individual fastqc files now that summarized in multiqc output
rm ../Raw_Data/*fastqc*

  # running cutadapt to remove primers:
    # 16S-universal
for sample in $(cat 16S-universal-unique-sample-IDs.txt); do cutadapt -a ^GTGCCAGCMGCCGCGGTAA...ATTAGAWACCCBDGTAGTCC -A ^GGACTACHVGGGTWTCTAAT...TTACCGCGGCKGCTGGCAC -o ../Trimmed_Sequence_Data/${sample}-R1-primer-trimmed.fastq.gz -p ../Trimmed_Sequence_Data/${sample}-R2-primer-trimmed.fastq.gz ../Raw_Data/${sample}_R1_raw.fastq.gz ../Raw_Data/${sample}_R2_raw.fastq.gz --discard-untrimmed -m 100; done > ../Trimmed_Reads/16S-universal-cutadapt.log 2>&1
      # keeping track of reads maintained (for some reason there are some very short sequences mixed in with the "16S_dust_FLT_uni_noINCd_VAC_2008_RISS2" sample. these are being dropped here. not sure how they are here as these should all be the same length if not yet manipulated)
cat <( printf "sample\traw_reads\tcutadapt_trimmed\n" ) <( paste 16S-universal-unique-sample-IDs.txt <( grep "read pairs processed" ../Trimmed_Sequence_Data/16S-universal-cutadapt.log | tr -s " " "\t" | cut -f 5 | tr -d "," ) <( grep "Pairs written" ../Trimmed_Sequence_Data/16S-universal-cutadapt.log | tr -s " " "\t" | cut -f 5 | tr -d "," ) ) > ../Trimmed_Sequence_Data/16S-universal-trimmed-read-counts.tsv
    # 16S-archaeal
for sample in $(cat 16S-archaeal-unique-sample-IDs.txt); do cutadapt -a ^GYGCASCAGKCGMGAAW...CAGCMGCCGCGGTAA -A ^TTACCGCGGCKGCTG...WTTCKCGMCTGSTGCRC -o ../Trimmed_Sequence_Data/${sample}-R1-primer-trimmed.fastq.gz -p ../Trimmed_Sequence_Data/${sample}-R2-primer-trimmed.fastq.gz ../Raw_Data/${sample}_R1_raw.fastq.gz ../Raw_Data/${sample}_R2_raw.fastq.gz --discard-untrimmed -m 100; done > ../Trimmed_Sequence_Data/16S-archaeal-cutadapt.log 2>&1
      # keeping track of reads maintained
cat <( printf "sample\traw_reads\tcutadapt_trimmed\n" ) <( paste 16S-archaeal-unique-sample-IDs.txt <( grep "read pairs processed" ../Trimmed_Sequence_Data/16S-archaeal-cutadapt.log | tr -s " " "\t" | cut -f 5 | tr -d "," ) <( grep "Pairs written" ../Trimmed_Sequence_Data/16S-archaeal-cutadapt.log | tr -s " " "\t" | cut -f 5 | tr -d "," ) ) > ../Trimmed_Sequence_Data/16S-archaeal-trimmed-read-counts.tsv

  # processing in R (quality filtering, ASV processing, taxonomic classification, and generating primary output files)
    # 16S-universal
Rscript 16S-universal-full-R-processing.R
    # 16S-archaeal
Rscript 16S-archaeal-full-R-processing.R

  # running fastqc on filtered reads
fastqc ../Filtered_Sequence_Data/*.gz -t 15
  # summarizing fastqc results with multiqc
multiqc -z -o ../FastQC_Outputs/ ../Filtered_Sequence_Data/
  # renaming the outputs
mv ../FastQC_Outputs/multiqc_data.zip ../FastQC_Outputs/filtered_multiqc_data.zip
mv ../FastQC_Outputs/multiqc_report.html ../FastQC_Outputs/filtered_multiqc_report.html
  # removing the individual fastqc files now that summarized in multiqc output
rm ../Filtered_Sequence_Data/*fastqc*
