# GeneLab processing commands for GLDS-120
The metadata for this dataset was updated to be consistent with the updated GeneLab plant metadata standards, therefore the RSEM raw 
counts data was used as input to re-run differential gene expression (DGE) analysis using pipeline version 
[GL-DPPD-7101-F](../../Pipeline_GL-DPPD-7101_Versions/GL-DPPD-7101-F.md).

The raw and processed data are available from: [https://osdr.nasa.gov/bio/repo/data/studies/OSD-120](https://osdr.nasa.gov/bio/repo/data/studies/OSD-120)

These samples were processed using the SLURM job scheduler. The following sub-directories contain the scripts used to generate each level of re-processed data.
  - [Metadata](Metadata): Scripts used to generate the runsheet used for DGE analysis
  - [04-05-DESeq2_NormCounts_DGE](04-05-DESeq2_NormCounts_DGE): Slurm and R scripts used to generate raw and normlaized counts tables, differential gene expression (DGE), and DGE tables with annotations

## Software used  
| Program              | Version          | Relevant Links                                                    |
|:---------------------|:-----------------|:------------------------------------------------------------------|
| FastQC               | 0.11.8           | https://www.bioinformatics.babraham.ac.uk/projects/fastqc/        |
| MultiQC              | 1.7              | https://multiqc.info/                                             |
| Cutadapt             | 2.3              | https://cutadapt.readthedocs.io/en/stable/                        |
| TrimGalore!          | 0.6.2            | https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/   |
| STAR                 | 2.7.1a           | https://github.com/alexdobin/STAR                                 |
| RSEM                 | 1.3.1            | https://github.com/deweylab/RSEM                                  |
| infer_experiment.py  | 4.0.0            | https://sourceforge.net/projects/rseqc                            |
| R                    | 4.1.2            | https://www.r-project.org                                         |
| Bioconductor         | 3.14             | https://bioconductor.org                                          |
| DESeq2               | 1.34.0           | https://bioconductor.org/packages/release/bioc/html/DESeq2.html   |
| tximport             | 1.22.0           | https://bioconductor.org/packages/release/bioc/html/tximport.html |
| dp_tools             | 1.1.6            | https://github.com/J-81/dp_tools                                  |
