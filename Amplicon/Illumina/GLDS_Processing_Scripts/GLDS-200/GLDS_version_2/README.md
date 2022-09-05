# GeneLab processing commands for GLDS-200

This dataset was processed with [GL-DPPD-7104-A](../../../Pipeline_GL-DPPD-7104_Versions/GL-DPPD-7104-A.md).

The raw and processed data are available from: [https://genelab-data.ndc.nasa.gov/genelab/accession/GLDS-200/](https://genelab-data.ndc.nasa.gov/genelab/accession/GLDS-200/)

The [full-command-line-log.sh](full-command-line-log.sh) file holds all commands as the processing was performed, which calls the two R scripts [16S-full-R-processing.R](16S-full-R-processing.R) and [ITS-full-R-processing.R](ITS-full-R-processing.R). 


## Software used  

|Program|Version|Relevant Links|
|:------|:-----:|:-------------|
|FastQC|0.11.8|[https://www.bioinformatics.babraham.ac.uk/projects/fastqc/](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)|
|MultiQC|1.7|[https://multiqc.info/](https://multiqc.info/)|
|Cutadapt|2.3|[https://cutadapt.readthedocs.io/en/stable/](https://cutadapt.readthedocs.io/en/stable/)|
|DADA2|1.12|[https://www.bioconductor.org/packages/release/bioc/html/dada2.html](https://www.bioconductor.org/packages/release/bioc/html/dada2.html)|
|DECIPHER|2.14|[https://bioconductor.org/packages/release/bioc/html/DECIPHER.html](https://bioconductor.org/packages/release/bioc/html/DECIPHER.html)|
|biomformat|1.12.0|[https://github.com/joey711/biomformat](https://github.com/joey711/biomformat)|

