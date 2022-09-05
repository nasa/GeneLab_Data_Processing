############################################################################################
## Command-line processing file for 16S data of NASA GLDS-191                             ##
## https://genelab-data.ndc.nasa.gov/genelab/accession/GLDS-191/                          ##
##                                                                                        ##
## This code as written expects to be run within the "processing_info/" directory         ##
## with the raw starting fastq files present in a "Raw_Data/"" directory                  ##
## Processed by Michael D. Lee (Mike.Lee@nasa.gov)                                        ##
############################################################################################

fastqc -v # v0.11.8
multiqc -v # v1.7

  # running fastqc on raw data
fastqc ../Raw_Data/*.gz -t 20
  # summarizing fastqc results with multiqc
multiqc -z -o ../FastQC_Outputs/ ../Raw_Data/ --interactive
  # renaming the outputs
mv ../FastQC_Outputs/multiqc_data.zip ../FastQC_Outputs/raw_multiqc_data.zip
mv ../FastQC_Outputs/multiqc_report.html ../FastQC_Outputs/raw_multiqc_report.html
  # removing the individual fastqc files now that summarized in multiqc output
rm ../Raw_Data/*fastqc*


 # running cutadapt to remove primers:
for sample in $(cat unique-sample-IDs.txt); do cutadapt -g CCTACGGGNGGCWGCAG -G GACTACHVGGGTATCTAATCC -o ../Trimmed_Sequence_Data/${sample}-R1-primer-trimmed.fastq.gz -p ../Trimmed_Sequence_Data/${sample}-R2-primer-trimmed.fastq.gz --discard-untrimmed -m 250 ../Raw_Data/${sample}_R1_raw.fastq.gz ../Raw_Data/${sample}_R2_raw.fastq.gz; done > ../Trimmed_Sequence_Data/cutadapt.log 2>&1

cat <( printf "sample\traw_reads\tcutadapt_trimmed\n" ) <( paste unique-sample-IDs.txt <( grep "read pairs processed" ../Trimmed_Sequence_Data/cutadapt.log | tr -s " " "\t" | cut -f 5 | tr -d "," ) <( grep "Pairs written" ../Trimmed_Sequence_Data/cutadapt.log | tr -s " " "\t" | cut -f 5 | tr -d "," ) ) > ../Trimmed_Sequence_Data/trimmed-read-counts.tsv

  # processing in R (quality filtering, ASV processing, taxonomic classification, and generating primary output files)
Rscript full-R-processing.R

  # running fastqc on filtered reads
fastqc ../Filtered_Sequence_Data/*.gz -t 20
  # summarizing fastqc results with multiqc
multiqc -z -o ../FastQC_Outputs/ ../Filtered_Sequence_Data/ --interactive
  # renaming the outputs
mv ../FastQC_Outputs/multiqc_data.zip ../FastQC_Outputs/filtered_multiqc_data.zip
mv ../FastQC_Outputs/multiqc_report.html ../FastQC_Outputs/filtered_multiqc_report.html
  # removing the individual fastqc files now that summarized in multiqc output
rm ../Filtered_Sequence_Data/*fastqc*
