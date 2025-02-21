#!/usr/bin/env python
import click
from pathlib import Path
import json
import os

# Output structure config
STRUCTURE = {
    "rnaseq": {
        "microbes": {
            "components": {
                "raw_reads": {
                    "outputs": {
                        "raw_fastq": "00-RawData/Fastq",
                        "raw_fastqc": "00-RawData/FastQC_Reports",
                        "raw_multiqc": "00-RawData/FastQC_Reports"
                    }
                },
                "trimmed_reads": {
                    "outputs": {
                        "fastq": "01-TG_Preproc/Fastq",
                        "fastqc": "01-TG_Preproc/FastQC_Reports",
                        "trimming_reports": "01-TG_Preproc/Trimming_Reports"
                    }
                },
                "alignments": {
                    "02-Bowtie2_Alignment": {
                        "{sample_name}": {}  # Sample-specific subdirectories
                    }
                },
                "counts": {
                    "03-FeatureCounts": {
                    }
                },
                "dge": {
                    "04-DESeq2_NormCounts": {},
                    "04-rRNArm_DESeq2_NormCounts": {},
                    "05-rRNArm_DESeq2_DGE": {},
                    "05-DESeq2_DGE": {}
                }
            }
        }
    }
}

# Future tissue-specific structure could look like:
TISSUE_STRUCTURE = {
    "dge": {
        "04-{tissue}-DESeq2_NormCounts": {},
        "04-{tissue}-rRNArm_DESeq2_NormCounts": {},
        "05-{tissue}-rRNArm_DESeq2_DGE": {},
        "05-{tissue}-DESeq2_DGE": {}
    }
}

@click.command()
@click.option('--assay-type', type=click.Choice(['rnaseq', 'scrna']), default='rnaseq')
@click.option('--assay-suffix', type=click.STRING, default="_GLbulkRNAseq")
@click.option('--runsheet-path', type=click.Path(exists=True), help="Path to runsheet")
@click.option('--outdir', type=click.Path(), default=Path.cwd(), help="Output directory")
@click.option('--paired-end', type=click.STRING, help="Paired end setting")
@click.option('--mode', type=click.Choice(['microbes', 'default']), default='default')
@click.option('--run-components', type=click.STRING, help="Component to validate (e.g. raw_reads)")
@click.option('--raw-fastq', type=click.Path(exists=True), help="Path to raw fastq directory")
@click.option('--raw-fastqc', type=click.Path(exists=True), help="Path to raw fastqc directory")
@click.option('--raw-multiqc', type=click.Path(exists=True), help="Path to raw multiqc directory")
def vv(assay_type, assay_suffix, runsheet_path, outdir, paired_end, mode, run_components, raw_fastq, raw_fastqc, raw_multiqc):
    """Organize pipeline outputs and optionally validate"""
    outdir = Path(outdir)
    
    # Stage files if inputs provided
    if any([raw_fastq, raw_fastqc, raw_multiqc]):
        file_paths = {
            'raw_fastq': raw_fastq,
            'raw_fastqc': raw_fastqc,
            'raw_multiqc': raw_multiqc
        }
        stage_files(assay_type, 'raw_reads', **file_paths)
    
    # Run validation if component specified
    if run_components:
        with open("VV_log.tsv", "w") as f:
            f.write(f"Stub validation log for {run_components}\n")

def stage_files(assay_type, section, **file_paths):
    """
    Stage files either by component or direct paths
    
    Args:
        assay_type (str): e.g. 'rnaseq'
        section (str): e.g. 'raw_reads'
        **file_paths: Keyword args for direct file paths (raw_fastq, raw_fastqc, etc)
    """
    structure = STRUCTURE[assay_type]['microbes']['components'][section]['outputs']
    
    # Direct path staging
    for file_type, path in file_paths.items():
        if path:  # Only process if path was provided
            target_dir = structure[file_type]
            stage_to_location(path, target_dir)

def stage_to_location(source_path, target_dir):
    """Helper to stage files to their target location"""
    os.makedirs(target_dir, exist_ok=True)
    
    # Get the ultimate source by following all symlinks
    ultimate_source = os.path.realpath(source_path)
    
    if os.path.isdir(source_path):
        # For directories, link their contents directly into target_dir
        for item in os.listdir(source_path):
            src = os.path.realpath(os.path.join(source_path, item))  # Get ultimate source for each file
            dst = os.path.join(target_dir, item)
            os.symlink(src, dst)
    else:
        # For single files
        dst = os.path.join(target_dir, os.path.basename(source_path))
        os.symlink(ultimate_source, dst)

def get_target_dir(structure, file_type):
    """Traverse structure to find target directory for file type"""
    # Implementation depends on your exact structure format
    pass

if __name__ == '__main__':
    vv()

# Component based:
stage_files('rnaseq', 'raw_reads', 
           components=['raw_reads'],
           raw_fastq='path/to/fastq',
           raw_fastqc='path/to/fastqc')

# Direct path based:
stage_files('rnaseq', 'raw_reads',
           raw_fastq='path/to/fastq',
           raw_fastqc='path/to/fastqc')
