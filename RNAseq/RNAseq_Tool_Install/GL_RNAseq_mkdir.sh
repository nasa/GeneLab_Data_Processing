mkdir 00-RawData  01-TG_Preproc  02-STAR_Alignment  03-RSEM_Counts  04-DESeq2_NormCounts  05-DESeq2_DGE  Metadata processing_scripts
mkdir 00-RawData/Fastq
mkdir 00-RawData/FastQC_Reports
mkdir 00-RawData/FastQC_Reports/raw_multiqc_report

mkdir 01-TG_Preproc/Fastq
mkdir 01-TG_Preproc/FastQC_Reports
mkdir 01-TG_Preproc/FastQC_Reports/trimmed_multiqc_report
mkdir 01-TG_Preproc/Trimming_Reports

mkdir 05-DESeq2_DGE/ERCC_NormDGE

mkdir processing_scripts/00-RawData
mkdir processing_scripts/01-TG_Preproc
mkdir processing_scripts/02-STAR_Alignment
mkdir processing_scripts/03-RSEM_Counts
mkdir processing_scripts/04-05-DESeq2_NormCounts_DGE

mkdir processing_scripts/00-RawData/raw_fastqc_out_logs
mkdir processing_scripts/01-TG_Preproc/TG_out_logs
mkdir processing_scripts/01-TG_Preproc/trimmed_fastqc_out_logs
mkdir processing_scripts/02-STAR_Alignment/STAR_align_out_logs
mkdir processing_scripts/03-RSEM_Counts/rsem_count_out_logs
