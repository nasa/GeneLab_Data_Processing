# GeneLab processing commands for GLDS-120
The incorrect RSEM `--strandedness` setting was used to generate the original processed data, therefore this dataset was re-processed 
starting with the RSEM counts step using pipeline version [GL-DPPD-7101-D](../../Pipeline_GL-DPPD-7101_Versions/GL-DPPD-7101-D.md).
> For the re-processed data, the RSEM `--strandedness` was set to none.

The raw and processed data are available from: [https://osdr.nasa.gov/bio/repo/data/studies/OSD-120](https://osdr.nasa.gov/bio/repo/data/studies/OSD-120)

These samples were processed in batch using the SLURM job scheduler. The following sub-directories contain the scripts used to generate each level of re-processed data.
  - [RSeQC_analyses](RSeQC_analyses): Scripts used to generate the RSeQC analyses: infer experiment
  - [03-RSEM_Counts](03-RSEM_Counts): Slurm scripts used to generate RSEM reference and RSEM count files
  - [Metadata](Metadata): Metadata table used for DGE analysis
  - [04-05-DESeq2_NormCounts_DGE](04-05-DESeq2_NormCounts_DGE): Slurm and R scripts used to generate raw and normlaized counts tables, differential gene expression (DGE), and DGE tables with annotations

## Software used  
|Program|Version|Relevant Links|
|:------|:------:|:-------------|
|MultiQC|1.9|[https://multiqc.info/](https://multiqc.info/)|
|infer_experiment.py|4.0.0|[https://sourceforge.net/projects/rseqc](https://sourceforge.net/projects/rseqc)|
|RSEM|1.3.1|[https://github.com/deweylab/RSEM](https://github.com/deweylab/RSEM)|
|R|4.0.3|[https://www.r-project.org](https://www.r-project.org)|
|Bioconductor|3.12|[https://bioconductor.org](https://bioconductor.org)|
|DESeq2|1.30.0|[https://bioconductor.org/packages/release/bioc/html/DESeq2.html](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)|
|tximport|1.18.0|[https://bioconductor.org/packages/release/bioc/html/tximport.html](https://bioconductor.org/packages/release/bioc/html/tximport.html)|
|tidyverse|1.3.0|[https://www.tidyverse.org](https://www.tidyverse.org)|
|STRINGdb|2.2.0|[https://www.bioconductor.org/packages/release/bioc/html/STRINGdb.html](https://www.bioconductor.org/packages/release/bioc/html/STRINGdb.html)|
|PANTHER.db|1.0.10|[https://bioconductor.org/packages/release/data/annotation/html/PANTHER.db.html](https://bioconductor.org/packages/release/data/annotation/html/PANTHER.db.html)|
|org.At.tair.db|3.12.0|[https://bioconductor.org/packages/release/data/annotation/html/org.At.tair.db.html](https://bioconductor.org/packages/release/data/annotation/html/org.At.tair.db.html)|
