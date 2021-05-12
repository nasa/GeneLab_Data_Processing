####################################################################################
## Command-line processing file for 16S and ITS data of NASA GLDS-200             ##
## https://genelab-data.ndc.nasa.gov/genelab/accession/GLDS-200/                  ##
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
    # 16S
for sample in $(cat 16S-unique-sample-IDs.txt); do cutadapt -a ^GTGCCAGCMGCCGCGGTAA...ATTAGATACCCSBGTAGTCC -A ^GGACTACVSGGGTATCTAAT...TTACCGCGGCKGCTGGCAC -o ../Trimmed_Sequence_Data/${sample}-R1-primer-trimmed.fastq.gz -p ../Trimmed_Sequence_Data/${sample}-R2-primer-trimmed.fastq.gz ../Raw_Data/${sample}_R1_raw.fastq.gz ../Raw_Data/${sample}_R2_raw.fastq.gz --discard-untrimmed; done > ../Trimmed_Sequence_Data/16S-cutadapt.log 2>&1
      # keeping track of reads maintained
cat <( printf "sample\traw_reads\tcutadapt_trimmed\n" ) <( paste 16S-unique-sample-IDs.txt <( grep "read pairs processed" ../Trimmed_Sequence_Data/16S-cutadapt.log | tr -s " " "\t" | cut -f 5 | tr -d "," ) <( grep "Pairs written" ../Trimmed_Sequence_Data/16S-cutadapt.log | tr -s " " "\t" | cut -f 5 | tr -d "," ) ) > ../Trimmed_Sequence_Data/16S-trimmed-read-counts.tsv
    # ITS
for sample in $(cat ITS-unique-sample-IDs.txt); do cutadapt -a ^CTTGGTCATTTAGAGGAAGTAA...GCATCGATGAAGAACGCAGC -A ^GCTGCGTTCTTCATCGATGC...TTACTTCCTCTAAATGACCAAG -o ../Trimmed_Sequence_Data/${sample}-R1-primer-trimmed.fastq.gz -p ../Trimmed_Sequence_Data/${sample}-R2-primer-trimmed.fastq.gz ../Raw_Data/${sample}_R1_raw.fastq.gz ../Raw_Data/${sample}_R2_raw.fastq.gz --discard-untrimmed; done > ../Trimmed_Sequence_Data/ITS-cutadapt.log 2>&1
      # keeping track of reads maintained
cat <( printf "sample\traw_reads\tcutadapt_trimmed\n" ) <( paste ITS-unique-sample-IDs.txt <( grep "read pairs processed" ../Trimmed_Sequence_Data/ITS-cutadapt.log | tr -s " " "\t" | cut -f 5 | tr -d "," ) <( grep "Pairs written" ../Trimmed_Sequence_Data/ITS-cutadapt.log | tr -s " " "\t" | cut -f 5 | tr -d "," ) ) > ../Trimmed_Sequence_Data/ITS-trimmed-read-counts.tsv

  # processing in R (quality filtering, ASV processing, taxonomic classification, and generating primary output files)
    # 16S
Rscript 16S-full-R-processing.R
    # ITS
Rscript ITS-full-R-processing.R

  # running fastqc on filtered reads
fastqc ../Filtered_Sequence_Data/*.gz -t 15
  # summarizing fastqc results with multiqc
multiqc -z -o ../FastQC_Outputs/ ../Filtered_Sequence_Data/
  # renaming the outputs
mv ../FastQC_Outputs/multiqc_data.zip ../FastQC_Outputs/filtered_multiqc_data.zip
mv ../FastQC_Outputs/multiqc_report.html ../FastQC_Outputs/filtered_multiqc_report.html
  # removing the individual fastqc files now that summarized in multiqc output
rm ../Filtered_Sequence_Data/*fastqc*
