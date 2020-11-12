############################################################################################
## Command-line processing file for 16S data of NASA GLDS-146                             ##
## https://genelab-data.ndc.nasa.gov/genelab/accession/GLDS-146/                          ##
##                                                                                        ##
## This code as written expects to be run within the "processing_info/" directory         ##
## with the raw starting fastq files present in a "Raw_Data/"" directory                  ##
## Processed by Michael D. Lee (Mike.Lee@nasa.gov)                                        ##
## The submitted data already had primers removed and gone through some quality filtering ##
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
