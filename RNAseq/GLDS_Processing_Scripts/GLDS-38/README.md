# GeneLab processing commands for GLDS-38
This dataset was processed with [GL-DPPD-7101-A](../../Pipeline_GL-DPPD-7101_Versions/GL-DPPD-7101-A.md).

The raw and processed data are available from: [https://genelab-data.ndc.nasa.gov/genelab/accession/GLDS-38/](https://genelab-data.ndc.nasa.gov/genelab/accession/GLDS-38/)

These samples were processed in batch using the SLURM job scheduler. The following sub-directories contain the scripts used to generate each level of processed data.
  - [00-RawData](00-RawData): Slurm scripts used to generate raw fastQC and multiQC reports
  - [01-TG_PreProc](01-TG_Preproc): Slurm scripts used to generate trimmed reads, trimmed fastQC and multiQC reports
  - [02-STAR_Alignment](02-STAR_Alignment): Slurm scripts used to generate STAR reference and STAR alignment files
  - [03-RSEM_Counts](03-RSEM_Counts): Slurm scripts used to generate RSEM reference and RSEM count files
  - [04-05-DESeq2_NormCounts_DGE](04-05-DESeq2_NormCounts_DGE): Slurm and R scripts used to generate raw and normlaized counts tables, differential gene expression (DGE), and DGE tables with annotations

## Software used  
|Program|Version|Relevant Links|
|:------|:------:|:-------------|
|FastQC|0.11.8|[https://www.bioinformatics.babraham.ac.uk/projects/fastqc/](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)|
|MultiQC|1.7|[https://multiqc.info/](https://multiqc.info/)|
|Cutadapt|2.3|[https://cutadapt.readthedocs.io/en/stable/](https://cutadapt.readthedocs.io/en/stable/)|
|TrimGalore!|0.6.2|[https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)|
|STAR|2.7.1a|[https://github.com/alexdobin/STAR](https://github.com/alexdobin/STAR)|
|RSEM|1.3.1|[https://github.com/deweylab/RSEM](https://github.com/deweylab/RSEM)|
|Bioconductor|3.9|[https://bioconductor.org](https://bioconductor.org)|
|DESeq2|1.26.0|[https://bioconductor.org/packages/release/bioc/html/DESeq2.html](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)|
|tximport|1.14.0|[https://bioconductor.org/packages/release/bioc/html/tximport.html](https://bioconductor.org/packages/release/bioc/html/tximport.html)|
|tidyverse|1.2.1|[https://www.tidyverse.org](https://www.tidyverse.org)|
|STRINGdb|1.24.0|[https://www.bioconductor.org/packages/release/bioc/html/STRINGdb.html](https://www.bioconductor.org/packages/release/bioc/html/STRINGdb.html)|
|PANTHER.db|1.0.4|[https://bioconductor.org/packages/release/data/annotation/html/PANTHER.db.html](https://bioconductor.org/packages/release/data/annotation/html/PANTHER.db.html)|
|org.At.tair.db|3.8.2|[https://bioconductor.org/packages/release/data/annotation/html/org.At.tair.db.html](https://bioconductor.org/packages/release/data/annotation/html/org.At.tair.db.html)|
