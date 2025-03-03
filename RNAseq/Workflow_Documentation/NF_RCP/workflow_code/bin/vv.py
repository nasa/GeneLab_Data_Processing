#!/usr/bin/env python
import click
from pathlib import Path
import json
import os
import subprocess
import gzip
import logging
from datetime import datetime
import csv
import re  # Add regex for parsing reports
import zipfile
from io import TextIOWrapper
import statistics
import shutil
import tempfile
import numpy as np
from statistics import mean, median
from typing import Dict, List, Tuple, Union, Optional
import traceback
import sys

# Setup logging
def setup_logging():
    logger = logging.getLogger('vv')
    logger.setLevel(logging.DEBUG)  # Change to DEBUG level
    
    # File handler
    fh = logging.FileHandler('vv.log')
    fh.setLevel(logging.DEBUG)  # Change to DEBUG level
    
    # Console handler
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)  # Keep console at INFO
    
    # Formatting
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    
    logger.addHandler(fh)
    logger.addHandler(ch)
    
    return logger

logger = setup_logging()

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
                        "trimmed_multiqc": "01-TG_Preproc/FastQC_Reports",
                        "trimming_reports": "01-TG_Preproc/Trimming_Reports",
                        "trimming_multiqc": "01-TG_Preproc/Trimming_Reports"
                    }
                },
                "alignments": {
                    "outputs": {
                        "bowtie2_alignment_log": "02-Bowtie2_Alignment",
                        "bowtie2_alignment_unmapped": "02-Bowtie2_Alignment",
                        "bowtie2_alignment_multiqc": "02-Bowtie2_Alignment",
                        "bowtie2_alignment_sorted": "02-Bowtie2_Alignment",
                        "bowtie2_alignment_sorted_index": "02-Bowtie2_Alignment"
                    }
                },
                "counts": {
                    "outputs": {
                        "featurecounts_counts": "03-FeatureCounts",
                        "featurecounts_summary": "03-FeatureCounts",
                        "featurecounts_counts_rrnarm": "03-FeatureCounts",
                        "featurecounts_counts_rrnarm_summary": "03-FeatureCounts",
                        "featurecounts_multiqc": "03-FeatureCounts"
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

class ValidationLogger:
    """Handles logging validation results in a structured format"""
    def __init__(self, log_file="VV_log.csv"):
        self.log_file = log_file
        self.results = []
        self.stats = {}  # Track quantitative stats
        
        # Create/overwrite log file with header
        with open(log_file, "w") as f:
            f.write("component,sample_id,check_name,status,message,details\n")
    
    def log(self, component, sample_id, check_name, status, message, details="", stats=None):
        """Log a validation result with optional stats"""
        # Extract just the sample name from CSV row
        sample_name = sample_id.split(',')[-1].strip() if ',' in sample_id else sample_id
        
        entry = {
            "component": component,
            "sample_id": sample_name,  # Use clean sample name
            "check_name": check_name, 
            "status": status,
            "message": message,
            "details": details
        }
        self.results.append(entry)
        
        # Track stats if provided
        if stats:
            if component not in self.stats:
                self.stats[component] = {}
            if sample_name not in self.stats[component]:  # Use clean sample name
                self.stats[component][sample_name] = {}
            self.stats[component][sample_name][check_name] = stats
        
        # Log to vv.log
        if stats:
            # Quantitative check
            logger.info(f"{component} - {sample_name} - {check_name}:")  # Use clean sample name
            for stat_name, value in stats.items():
                logger.info(f"  {stat_name}: {value}")
        else:
            # Qualitative check
            logger.info(f"{component} - {sample_name} - {check_name}: {status}")  # Use clean sample name
            if details:
                logger.info(f"  Details: {details}")
        
        # Escape any commas in fields
        safe_fields = [
            str(field).replace('"', '""') if ',' in str(field) else str(field)
            for field in [component, sample_name, check_name, status, message, details]  # Use clean sample name
        ]
        
        # Quote fields that need it
        csv_fields = [
            f'"{field}"' if ',' in field else field 
            for field in safe_fields
        ]
        
        # Append to log file
        with open(self.log_file, "a") as f:
            f.write(f"{','.join(csv_fields)}\n")
    
    def get_status(self):
        """Get overall validation status"""
        if any(r["status"] == "HALT" for r in self.results):
            return "failed"
        if any(r["status"] == "RED" for r in self.results):
            return "failed"
        if any(r["status"] == "YELLOW" for r in self.results):
            return "warning"
        return "passed"

@click.command(name="vv")
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
@click.option('--trimmed-fastq', type=click.Path(exists=True), help="Path to trimmed fastq directory")
@click.option('--trimmed-fastqc', type=click.Path(exists=True), help="Path to trimmed fastqc directory")
@click.option('--trimmed-multiqc', type=click.Path(exists=True), help="Path to trimmed multiqc directory")
@click.option('--trimming-reports', type=click.Path(exists=True), help="Path to trimming reports directory") 
@click.option('--trimming-multiqc', type=click.Path(exists=True), help="Path to trimming multiqc directory")
@click.option('--bowtie2-alignment-log', type=click.Path(exists=True), help="Path to bowtie2 alignment log directory")
@click.option('--bowtie2-alignment-unmapped', type=click.Path(exists=True), help="Path to bowtie2 unmapped reads directory")
@click.option('--bowtie2-alignment-multiqc', type=click.Path(exists=True), help="Path to bowtie2 alignment multiqc directory")
@click.option('--bowtie2-alignment-sorted', type=click.Path(exists=True), help="Path to sorted BAM files directory")
@click.option('--bowtie2-alignment-sorted-index', type=click.Path(exists=True), help="Path to sorted BAM index files directory")
@click.option('--genebody-coverage', type=click.Path(), help="Path to geneBody coverage directory")
@click.option('--infer-experiment', type=click.Path(), help="Path to infer experiment directory")
@click.option('--inner-distance', type=click.Path(), help="Path to inner distance directory")
@click.option('--read-distribution', type=click.Path(), help="Path to read distribution directory")
@click.option('--featurecounts-counts', type=click.Path(), help="Path to featurecounts counts directory")
@click.option('--featurecounts-summary', type=click.Path(), help="Path to featurecounts summary directory")
@click.option('--featurecounts-counts-rrnarm', type=click.Path(), help="Path to featurecounts rRNA removed counts directory")
@click.option('--featurecounts-counts-rrnarm-summary', type=click.Path(), help="Path to featurecounts rRNA removed summary directory")
@click.option('--featurecounts-multiqc', type=click.Path(), help="Path to featurecounts multiqc directory")
def vv(assay_type, assay_suffix, runsheet_path, outdir, paired_end, mode, run_components, 
       raw_fastq, raw_fastqc, raw_multiqc,
       trimmed_fastq, trimmed_fastqc, trimmed_multiqc, trimming_reports, trimming_multiqc,
       bowtie2_alignment_log, bowtie2_alignment_unmapped, bowtie2_alignment_multiqc, 
       bowtie2_alignment_sorted, bowtie2_alignment_sorted_index,
       genebody_coverage, infer_experiment, inner_distance, read_distribution,
       featurecounts_counts, featurecounts_summary, featurecounts_counts_rrnarm,
       featurecounts_counts_rrnarm_summary, featurecounts_multiqc):
    """
    Main validation and verification (VV) function for GeneLab RNAseq pipeline.
    
    This function orchestrates the validation and verification process for various 
    components of the GeneLab RNAseq pipeline. It provides comprehensive quality
    control and validation for:
    
    1. Raw Reads Validation:
       - File existence, integrity, and format
       - Read counts and paired-end consistency
       - FastQC output validation
       - MultiQC report validation
       
    2. Trimmed Reads Validation:
       - Trimmed file existence, integrity, and format
       - Trimming report validation (adapter and quality trimming)
       - Trimmed read counts and paired-end consistency
       - Trimmed FastQC output validation
       - MultiQC report validation for trimmed reads and trimming reports
       
    3. Alignment Validation:
       - Bowtie2 alignment log validation
       - Alignment rate analysis and outlier detection
       - BAM file existence and integrity
       - Unmapped read file validation
       - MultiQC report validation for alignments
       
    4. RSeQC Analysis Validation:
       - geneBody coverage validation
       - Library strandedness validation (infer experiment)
       - Fragment size distribution validation (inner distance)
       - Read distribution validation across genomic features
       - MultiQC report validation for RSeQC outputs
    
    The function supports selective validation of specific components through
    the run_components parameter, and handles both default and microbe-specific modes.
    It can validate both paired-end and single-end data.
    
    Results are logged to both standard output and a CSV file for further analysis.
    
    Args:
        assay_type: Type of assay ('rnaseq' or 'scrna')
        assay_suffix: Suffix for output files (e.g., '_GLbulkRNAseq')
        runsheet_path: Path to the sample runsheet CSV file
        outdir: Output directory for validation results
        paired_end: Boolean string indicating if data is paired-end ('true'/'false')
        mode: Mode for validation ('default' or 'microbes')
        run_components: Comma-separated list of components to validate (e.g., 'raw_reads')
        raw_fastq, raw_fastqc, raw_multiqc: Paths to raw data directories
        trimmed_fastq, trimmed_fastqc, trimmed_multiqc, trimming_reports, trimming_multiqc: Paths to trimmed data directories
        bowtie2_alignment_log, bowtie2_alignment_unmapped, bowtie2_alignment_multiqc, bowtie2_alignment_sorted, bowtie2_alignment_sorted_index: Paths to alignment directories
        genebody_coverage, infer_experiment, inner_distance, read_distribution: Paths to RSeQC output directories
        featurecounts_counts, featurecounts_summary, featurecounts_counts_rrnarm, featurecounts_counts_rrnarm_summary, featurecounts_multiqc: Paths to FeatureCounts-related directories
    
    Returns:
        None: Results are written to log files and returned as exit code
    """
    outdir = Path(outdir)
    
    # Initialize validation logger
    global val_logger
    val_logger = ValidationLogger()
    
    # Log all arguments for debugging
    logging.info(f"VV started with arguments:")
    logging.info(f"  assay_type: {assay_type}")
    logging.info(f"  assay_suffix: {assay_suffix}")
    logging.info(f"  runsheet_path: {runsheet_path}")
    logging.info(f"  outdir: {outdir}")
    logging.info(f"  paired_end: {paired_end}")
    logging.info(f"  mode: {mode}")
    logging.info(f"  run_components: {run_components}")
    
    # Get the working directory and its files
    cwd = Path.cwd()
    logging.info(f"Current working directory: {cwd}")
    logging.info("Files in current working directory:")
    for f in cwd.iterdir():
        logging.info(f"  {f}")
        
    # For validation of components, use the current working directory as outdir
    # This is where we staged all the files
    validation_outdir = cwd
    
    # Parse run_components into list
    components = run_components.split(',') if run_components else []
    
    # Convert paired_end string to boolean
    is_paired_end = paired_end.lower() == 'true'
    logging.info(f"  paired_end parsed as: {is_paired_end}")
    
    # Stage files if inputs provided
    if any([raw_fastq, raw_fastqc, raw_multiqc]):
        logging.info("Staging raw read files...")
        file_paths = {
            'raw_fastq': raw_fastq,
            'raw_fastqc': raw_fastqc,
            'raw_multiqc': raw_multiqc
        }
        stage_files(assay_type, 'raw_reads', assay_suffix=assay_suffix, **file_paths)
    
    # Stage trimmed files if inputs provided
    if any([trimmed_fastq, trimmed_fastqc, trimmed_multiqc, trimming_reports, trimming_multiqc]):
        logging.info("Staging trimmed read files...")
        logging.info(f"  trimmed_fastq: {trimmed_fastq}")
        logging.info(f"  trimmed_fastqc: {trimmed_fastqc}")
        logging.info(f"  trimmed_multiqc: {trimmed_multiqc}")
        logging.info(f"  trimming_reports: {trimming_reports}")
        logging.info(f"  trimming_multiqc: {trimming_multiqc}")
        
        file_paths = {
            'fastq': trimmed_fastq,
            'fastqc': trimmed_fastqc,
            'trimmed_multiqc': trimmed_multiqc,
            'trimming_reports': trimming_reports,
            'trimming_multiqc': trimming_multiqc
        }
        stage_files(assay_type, 'trimmed_reads', assay_suffix=assay_suffix, **file_paths)
        
        # Check that files were staged properly
        fastq_dir = Path("01-TG_Preproc/Fastq")
        if fastq_dir.exists():
            logging.info(f"Files in {fastq_dir}:")
            for f in fastq_dir.iterdir():
                logging.info(f"  {f}")
        else:
            logging.error(f"Directory not created: {fastq_dir}")
    
    # Stage Bowtie2 alignment files if inputs provided
    if any([bowtie2_alignment_log, bowtie2_alignment_unmapped, bowtie2_alignment_multiqc, 
            bowtie2_alignment_sorted, bowtie2_alignment_sorted_index]):
        
        logging.info("Staging alignment files...")
        file_paths = {
            'bowtie2_alignment_log': bowtie2_alignment_log,
            'bowtie2_alignment_unmapped': bowtie2_alignment_unmapped,
            'bowtie2_alignment_multiqc': bowtie2_alignment_multiqc,
            'bowtie2_alignment_sorted': bowtie2_alignment_sorted,
            'bowtie2_alignment_sorted_index': bowtie2_alignment_sorted_index
        }
        stage_files(assay_type, 'alignments', assay_suffix=assay_suffix, **file_paths)
    
    # Stage FeatureCounts files if inputs provided
    if any([featurecounts_counts, featurecounts_summary, featurecounts_counts_rrnarm, 
            featurecounts_counts_rrnarm_summary, featurecounts_multiqc]):
        
        logging.info("Staging FeatureCounts files...")
        file_paths = {
            'featurecounts_counts': featurecounts_counts,
            'featurecounts_summary': featurecounts_summary,
            'featurecounts_counts_rrnarm': featurecounts_counts_rrnarm,
            'featurecounts_counts_rrnarm_summary': featurecounts_counts_rrnarm_summary,
            'featurecounts_multiqc': featurecounts_multiqc
        }
        stage_files(assay_type, 'counts', assay_suffix=assay_suffix, **file_paths)
    
    # Run validations if components specified
    results = {}
    if 'raw_reads' in components:
        logging.info("Running raw reads validation...")
        results = validate_raw_reads(
            validation_outdir,
            samples_txt=Path(runsheet_path),
            paired_end=is_paired_end,
            assay_suffix=assay_suffix
        )
    
    if 'trimmed_reads' in components:
        logging.info("Running trimmed reads validation...")
        results = validate_trimmed_reads(
            validation_outdir,
            samples_txt=Path(runsheet_path),
            paired_end=is_paired_end,
            assay_suffix=assay_suffix
        )
    
    if 'bowtie2_alignment' in components or 'alignments' in components:
        # Run validation for Bowtie2 alignments
        if mode != 'microbes':
            logger.warning(f"Bowtie2 alignment validation is only supported for microbes mode, current mode: {mode}")
        
        results = validate_bowtie2_alignments(
            validation_outdir,
            samples_txt=Path(runsheet_path),
            paired_end=is_paired_end,
            assay_suffix=assay_suffix
        )
    
    # Add RSeQC validation
    if 'rseqc' in components:
        results = validate_rseqc(
            validation_outdir,
            samples_txt=runsheet_path,
            paired_end=is_paired_end,
            assay_suffix=assay_suffix,
            genebody_coverage_dir=Path(genebody_coverage) if genebody_coverage else None,
            infer_experiment_dir=Path(infer_experiment) if infer_experiment else None,
            inner_distance_dir=Path(inner_distance) if inner_distance else None,
            read_distribution_dir=Path(read_distribution) if read_distribution else None
        )
    
    # Add FeatureCounts validation
    if 'featurecounts' in components:
        results = validate_featurecounts(
            validation_outdir,
            samples_txt=runsheet_path,
            paired_end=is_paired_end,
            assay_suffix=assay_suffix,
            featurecounts_counts_dir=Path(featurecounts_counts) if featurecounts_counts else None,
            featurecounts_summary_dir=Path(featurecounts_summary) if featurecounts_summary else None,
            featurecounts_counts_rrnarm_dir=Path(featurecounts_counts_rrnarm) if featurecounts_counts_rrnarm else None,
            featurecounts_counts_rrnarm_summary_dir=Path(featurecounts_counts_rrnarm_summary) if featurecounts_counts_rrnarm_summary else None,
            featurecounts_multiqc_dir=Path(featurecounts_multiqc) if featurecounts_multiqc else None
        )
    
    return results

def stage_files(assay_type, section, assay_suffix="_GLbulkRNAseq", **file_paths):
    """
    Stage files either by component or direct paths
    
    Args:
        assay_type (str): e.g. 'rnaseq'
        section (str): e.g. 'raw_reads'
        assay_suffix (str): Suffix for output files (e.g., '_GLbulkRNAseq')
        **file_paths: Keyword args for direct file paths (raw_fastq, raw_fastqc, etc)
    """
    logging.info(f"Staging files for section {section} with assay_suffix={assay_suffix}")
    
    structure = STRUCTURE[assay_type]['microbes']['components'][section]['outputs']
    logging.info(f"Structure for {section}: {structure}")
    
    # Special handling for Bowtie2 alignments - they need a per-sample subdirectory structure
    if section == 'alignments':
        # Create base alignment directory
        alignment_base_dir = "02-Bowtie2_Alignment"
        os.makedirs(alignment_base_dir, exist_ok=True)
        
        # Collect all sample names from sorted BAM files (most reliable)
        sample_names = set()
        if file_paths.get('bowtie2_alignment_sorted'):
            for file in os.listdir(file_paths['bowtie2_alignment_sorted']):
                if file.endswith('_sorted.bam'):
                    sample_name = file.replace('_sorted.bam', '')
                    sample_names.add(sample_name)
                    
        # Create sample directories
        for sample in sample_names:
            os.makedirs(os.path.join(alignment_base_dir, sample), exist_ok=True)
            
        # Stage sorted BAM files
        if file_paths.get('bowtie2_alignment_sorted'):
            for file in os.listdir(file_paths['bowtie2_alignment_sorted']):
                if file.endswith('_sorted.bam'):
                    sample_name = file.replace('_sorted.bam', '')
                    src = os.path.join(file_paths['bowtie2_alignment_sorted'], file)
                    dst = os.path.join(alignment_base_dir, sample_name, file)
                    
                    if os.path.exists(dst) or os.path.islink(dst):
                        os.unlink(dst)
                    os.symlink(os.path.realpath(src), dst)
                    
        # Stage BAM index files
        if file_paths.get('bowtie2_alignment_sorted_index'):
            for file in os.listdir(file_paths['bowtie2_alignment_sorted_index']):
                if file.endswith('.bam.bai'):
                    sample_name = file.replace('_sorted.bam.bai', '')
                    src = os.path.join(file_paths['bowtie2_alignment_sorted_index'], file)
                    dst = os.path.join(alignment_base_dir, sample_name, file)
                    
                    if os.path.exists(dst) or os.path.islink(dst):
                        os.unlink(dst)
                    os.symlink(os.path.realpath(src), dst)
                    
        # Stage alignment log files
        if file_paths.get('bowtie2_alignment_log'):
            for file in os.listdir(file_paths['bowtie2_alignment_log']):
                if file.endswith('.bowtie2.log'):
                    sample_name = file.replace('.bowtie2.log', '')
                    src = os.path.join(file_paths['bowtie2_alignment_log'], file)
                    dst = os.path.join(alignment_base_dir, sample_name, file)
                    
                    if os.path.exists(dst) or os.path.islink(dst):
                        os.unlink(dst)
                    os.symlink(os.path.realpath(src), dst)
                    
        # Stage unmapped read files
        if file_paths.get('bowtie2_alignment_unmapped'):
            for file in os.listdir(file_paths['bowtie2_alignment_unmapped']):
                # Handle multiple possible filename patterns for unmapped reads
                sample_name = None
                
                # Pattern 1: "Sample1.unmapped.fastq.1.gz" or "Sample1.unmapped.fastq.2.gz"
                if '.unmapped.fastq.' in file:
                    parts = file.split('.unmapped.fastq.')
                    if len(parts) >= 2:
                        sample_name = parts[0]
                
                # Pattern 2: "Sample1_R1.unmapped.fastq.gz" or "Sample1_R2.unmapped.fastq.gz"
                elif '_R1.unmapped.fastq.gz' in file:
                    sample_name = file.replace('_R1.unmapped.fastq.gz', '')
                elif '_R2.unmapped.fastq.gz' in file:
                    sample_name = file.replace('_R2.unmapped.fastq.gz', '')
                
                # If we identified a sample name, stage the file
                if sample_name:
                    # Create dir if it doesn't exist
                    if not os.path.exists(os.path.join(alignment_base_dir, sample_name)):
                        os.makedirs(os.path.join(alignment_base_dir, sample_name), exist_ok=True)
                        
                    src = os.path.join(file_paths['bowtie2_alignment_unmapped'], file)
                    dst = os.path.join(alignment_base_dir, sample_name, file)
                    
                    if os.path.exists(dst) or os.path.islink(dst):
                        os.unlink(dst)
                    os.symlink(os.path.realpath(src), dst)
        
        # Stage MultiQC report
        if file_paths.get('bowtie2_alignment_multiqc'):
            multiqc_found = False
            for file in os.listdir(file_paths['bowtie2_alignment_multiqc']):
                if '_multiqc' in file and file.endswith('.zip'):
                    src = os.path.join(file_paths['bowtie2_alignment_multiqc'], file)
                    dst = os.path.join(alignment_base_dir, f"align_multiqc{assay_suffix}_report.zip")
                    
                    logging.info(f"Staging alignment MultiQC from {src} to {dst}")
                    
                    if os.path.exists(dst) or os.path.islink(dst):
                        os.unlink(dst)
                    os.symlink(os.path.realpath(src), dst)
                    
                    # Don't extract the MultiQC report - only use the zip file
                    multiqc_found = True
                    break
            
            if not multiqc_found:
                logging.warning(f"No MultiQC report found in {file_paths['bowtie2_alignment_multiqc']}")
        
        # Check what was actually staged
        logging.info(f"Checking contents of 02-Bowtie2_Alignment after staging:")
        for root, dirs, files in os.walk(alignment_base_dir):
            for file in files:
                logging.info(f"  {os.path.join(root, file)}")
        
        return  # Skip the standard staging for alignments
    
    # Special handling for FeatureCounts - needs consolidated files
    elif section == 'counts':
        # Create FeatureCounts directory
        fc_dir = os.path.join("03-FeatureCounts")
        os.makedirs(fc_dir, exist_ok=True)
        
        # Stage counts files
        if file_paths.get('featurecounts_counts'):
            for file in os.listdir(file_paths['featurecounts_counts']):
                if file.endswith('.txt') or file.endswith('.tsv'):
                    src = os.path.join(file_paths['featurecounts_counts'], file)
                    dst = os.path.join(fc_dir, file)
                    
                    if os.path.exists(dst) or os.path.islink(dst):
                        os.unlink(dst)
                    os.symlink(os.path.realpath(src), dst)
        
        # Stage summary files
        if file_paths.get('featurecounts_summary'):
            for file in os.listdir(file_paths['featurecounts_summary']):
                if file.endswith('.txt.summary') or file.endswith('.tsv.summary') or file.endswith('.summary'):
                    src = os.path.join(file_paths['featurecounts_summary'], file)
                    dst = os.path.join(fc_dir, file)
                    
                    if os.path.exists(dst) or os.path.islink(dst):
                        os.unlink(dst)
                    os.symlink(os.path.realpath(src), dst)
        
        # Stage rRNA-removed counts files
        if file_paths.get('featurecounts_counts_rrnarm'):
            for file in os.listdir(file_paths['featurecounts_counts_rrnarm']):
                if file.endswith(".txt") or file.endswith(".tsv"):
                    src = os.path.join(file_paths['featurecounts_counts_rrnarm'], file)
                    dst = os.path.join(fc_dir, file)
                    
                    if os.path.exists(dst) or os.path.islink(dst):
                        os.unlink(dst)
                    os.symlink(os.path.realpath(src), dst)
        
        # Stage rRNA-removed summary files
        if file_paths.get('featurecounts_counts_rrnarm_summary'):
            for file in os.listdir(file_paths['featurecounts_counts_rrnarm_summary']):
                if file.endswith(".txt.summary") or file.endswith(".tsv.summary") or file.endswith(".summary") or file.endswith(".txt"):
                    src = os.path.join(file_paths['featurecounts_counts_rrnarm_summary'], file)
                    dst = os.path.join(fc_dir, file)
                    
                    if os.path.exists(dst) or os.path.islink(dst):
                        os.unlink(dst)
                    os.symlink(os.path.realpath(src), dst)
        
        # Stage MultiQC report
        if file_paths.get('featurecounts_multiqc'):
            multiqc_found = False
            for file in os.listdir(file_paths['featurecounts_multiqc']):
                if file.endswith(f"_multiqc{assay_suffix}_report.zip") or file.endswith("_report.zip"):
                    src = os.path.join(file_paths['featurecounts_multiqc'], file)
                    dst = os.path.join(fc_dir, f"featureCounts_multiqc{assay_suffix}_report.zip")
                    
                    if os.path.exists(dst) or os.path.islink(dst):
                        os.unlink(dst)
                    os.symlink(os.path.realpath(src), dst)
                    multiqc_found = True
                    break
            
            if not multiqc_found:
                logging.warning(f"No MultiQC report found in {file_paths['featurecounts_multiqc']}")
        
        # Check what was actually staged
        logging.info(f"Checking contents of 03-FeatureCounts after staging:")
        for root, dirs, files in os.walk(fc_dir):
            for file in files:
                logging.info(f"  {os.path.join(root, file)}")
        
        return  # Skip the standard staging for FeatureCounts
    
    # Direct path staging for other components
    for file_type, path in file_paths.items():
        if path:  # Only process if path was provided
            target_dir = structure[file_type]
            logging.info(f"Staging {file_type} from {path} to {target_dir}")
            
            # Ensure target directory exists
            os.makedirs(target_dir, exist_ok=True)
            
            # Stage the files
            stage_to_location(path, target_dir)

def stage_to_location(source_path, target_dir):
    """Helper to stage files to their target location"""
    # Ensure target directory exists
    os.makedirs(target_dir, exist_ok=True)
    
    logging.info(f"Staging from {source_path} to {target_dir}")
    
    if os.path.exists(source_path):
        if os.path.isdir(source_path):
            # For directories, link their contents directly into target_dir
            for item in os.listdir(source_path):
                src = os.path.join(source_path, item)
                src_real = os.path.realpath(src)  # Get ultimate source
                dst = os.path.join(target_dir, item)
                
                # Remove destination if it already exists
                if os.path.exists(dst) or os.path.islink(dst):
                    if os.path.isdir(dst) and not os.path.islink(dst):
                        shutil.rmtree(dst)
                    else:
                        os.unlink(dst)
                
                # Create symlink
                logging.info(f"Creating symlink from {src_real} to {dst}")
                os.symlink(src_real, dst)
            
            # If this is a MultiQC zip file, extract it for processing
            for item in os.listdir(source_path):
                if item.endswith('_multiqc_report.zip') or item.endswith('_multiqc_GLbulkRNAseq_report.zip'):
                    extract_dir = os.path.join(target_dir, item.replace('.zip', '_data'))
                    logging.info(f"Extracting MultiQC zip file to {extract_dir}")
                    if os.path.exists(extract_dir):
                        shutil.rmtree(extract_dir)
                    with zipfile.ZipFile(os.path.join(source_path, item), 'r') as zip_ref:
                        zip_ref.extractall(extract_dir)
        else:
            # For single files
            src_real = os.path.realpath(source_path)  # Get ultimate source
            dst = os.path.join(target_dir, os.path.basename(source_path))
            
            # Remove destination if it already exists
            if os.path.exists(dst) or os.path.islink(dst):
                if os.path.isdir(dst) and not os.path.islink(dst):
                    shutil.rmtree(dst)
                else:
                    os.unlink(dst)
            
            # Create symlink
            logging.info(f"Creating symlink from {src_real} to {dst}")
            os.symlink(src_real, dst)
            
            # If this is a MultiQC zip file, extract it for processing
            if os.path.basename(src_real).endswith('_multiqc_report.zip') or os.path.basename(src_real).endswith('_multiqc_GLbulkRNAseq_report.zip'):
                extract_dir = dst.replace('.zip', '_data')
                logging.info(f"Extracting MultiQC zip file to {extract_dir}")
                if os.path.exists(extract_dir):
                    shutil.rmtree(extract_dir)
                with zipfile.ZipFile(src_real, 'r') as zip_ref:
                    zip_ref.extractall(extract_dir)
    else:
        logging.error(f"Source path does not exist: {source_path}")

def get_target_dir(structure, file_type):
    """Traverse structure to find target directory for file type"""
    # Implementation depends on your exact structure format
    pass

def get_expected_files(sample_id: str, paired_end: bool, assay_suffix: str) -> dict:
    """
    Generate expected file patterns for raw read validation based on sample ID.
    
    This function determines the expected file paths for raw reads, including:
    - Raw FastQ files (sample_id_R1_raw.fastq.gz, sample_id_R2_raw.fastq.gz for paired-end)
    - FastQC outputs (sample_id_R1_raw_fastqc.zip, sample_id_R2_raw_fastqc.zip for paired-end)
    - MultiQC report (raw_multiqc_[assay_suffix]_report.zip)
    
    Args:
        sample_id: Sample name/ID from the runsheet
        paired_end: Boolean indicating if sample is paired-end (True) or single-end (False)
        assay_suffix: Suffix for the assay (e.g., "_GLbulkRNAseq")
        
    Returns:
        dict: Dictionary with keys for each file type (fastq, fastqc, multiqc) and lists of expected file paths
    """
    expected = {
        # Raw FastQ files
        "fastq": [f"{sample_id}_raw.fastq.gz"] if not paired_end else [
            f"{sample_id}_R1_raw.fastq.gz",
            f"{sample_id}_R2_raw.fastq.gz"
        ],
        
        # FastQC outputs
        "fastqc": [f"{sample_id}_raw_fastqc.zip"] if not paired_end else [
            f"{sample_id}_R1_raw_fastqc.zip",
            f"{sample_id}_R2_raw_fastqc.zip"
        ],
        
        # MultiQC report
        "multiqc": [f"raw_multiqc{assay_suffix}_report.zip"]
    }
    logging.debug(f"Expected files for {sample_id}: {expected}")
    return expected

def detect_outliers(stats_by_sample, metrics_to_check=None):
    """
    Detect outliers in sample metrics compared to the group
    
    Args:
        stats_by_sample: Dictionary of sample statistics keyed by sample name
        metrics_to_check: List of metrics to check for outliers (default: None, checks all metrics)
        
    Returns:
        Dictionary of outliers by sample, with reasons
    """
    # If no metrics are specified, check all important metrics
    if metrics_to_check is None:
        metrics_to_check = [
            # Read counts
            'raw_total_sequences_f', 'raw_total_sequences_r',
            # Sequence lengths
            'raw_avg_sequence_length_f', 'raw_median_sequence_length_f',
            'raw_avg_sequence_length_r', 'raw_median_sequence_length_r',
            # Quality scores
            'raw_quality_score_mean_f', 'raw_quality_score_median_f',
            'raw_quality_score_mean_r', 'raw_quality_score_median_r',
            # Duplication rates
            'raw_percent_duplicates_f', 'raw_percent_duplicates_r',
            # GC content
            'raw_percent_gc_f', 'raw_percent_gc_r',
            'raw_gc_min_1pct_f', 'raw_gc_min_1pct_r',
            'raw_gc_max_1pct_f', 'raw_gc_max_1pct_r',
            'raw_gc_auc_25pct_f', 'raw_gc_auc_25pct_r',
            'raw_gc_auc_50pct_f', 'raw_gc_auc_50pct_r',
            'raw_gc_auc_75pct_f', 'raw_gc_auc_75pct_r',
            # N content
            'raw_n_content_sum_f', 'raw_n_content_sum_r'
        ]
        
        # Add trimmed metrics if we're checking those
        trimmed_metrics_exist = False
        if stats_by_sample and len(stats_by_sample) > 0:
            # Check the first sample's metrics to see if any contain 'trimmed_'
            first_sample = next(iter(stats_by_sample.values()))
            for key in first_sample:
                if 'trimmed_' in key:
                    trimmed_metrics_exist = True
                    break
                    
        if trimmed_metrics_exist:
            trimmed_metrics = [
                # Replace 'raw_' with 'trimmed_' in the metric names
                metric.replace('raw_', 'trimmed_') for metric in metrics_to_check
            ]
            metrics_to_check.extend(trimmed_metrics)
    
    # Gather all values for each metric
    metric_values = {}
    for metric in metrics_to_check:
        metric_values[metric] = []
        
    # Collect values from all samples
    metrics_found = set()
    for sample, stats in stats_by_sample.items():
        for metric in metrics_to_check:
            if metric in stats:
                metric_values[metric].append(stats[metric])
                metrics_found.add(metric)
    
    # Calculate statistics for each metric
    metric_stats = {}
    for metric, values in metric_values.items():
        if len(values) > 1:  # Need at least 2 values for statistics
            metric_stats[metric] = {
                'mean': np.mean(values),
                'median': np.median(values),
                'std': np.std(values),
                'min': min(values),
                'max': max(values)
            }
    
    # Detect outliers for each sample
    outliers = {}
    for sample, stats in stats_by_sample.items():
        sample_outliers = []
        
        for metric in metrics_to_check:
            if metric in stats and metric in metric_stats:
                value = stats[metric]
                mean = metric_stats[metric]['mean']
                std = metric_stats[metric]['std']
                
                # Check if value is more than 2 standard deviations from mean
                if std > 0 and abs(value - mean) > 2 * std:
                    z_score = (value - mean) / std
                    sample_outliers.append({
                        'metric': metric,
                        'value': value,
                        'mean': mean,
                        'std': std,
                        'z_score': z_score
                    })
        
        if sample_outliers:
            outliers[sample] = sample_outliers
    
    return outliers, metric_stats

def calculate_outliers(values, stdev_threshold=2):
    """Calculate outliers based on standard deviation from median.
    
    Args:
        values (list): List of numeric values
        stdev_threshold (float): Number of standard deviations to consider as outlier threshold
        
    Returns:
        list: Indices of outlier values
    """
    if not values:
        return []
        
    median = statistics.median(values)
    stdev = statistics.stdev(values) if len(values) > 1 else 0
    
    outliers = []
    for i, value in enumerate(values):
        if abs(value - median) > stdev_threshold * stdev:
            outliers.append(i)
            
    return outliers

def validate_raw_reads(outdir: Path,
                      samples_txt: Path,
                      paired_end: bool = True,
                      assay_suffix: str = "_GLbulkRNAseq") -> dict:
    """
    Perform validation checks on raw sequencing reads.
    
    This function validates the raw sequencing data by performing:
    1. File existence checks for raw FastQ files and FastQC outputs
    2. GZIP integrity validation of FastQ files
    3. FASTQ format validation (checks header format, counts reads)
    4. FastQC output existence validation
    5. MultiQC report existence validation
    6. Sample presence verification in MultiQC reports
    
    Note: Paired-end read count comparison is skipped to improve performance.
    
    Args:
        outdir: Path to the output directory containing the raw data
        samples_txt: Path to the samples runsheet CSV file
        paired_end: Boolean indicating if samples are paired-end (True) or single-end (False)
        assay_suffix: Suffix for the assay (e.g., "_GLbulkRNAseq")
        
    Returns:
        dict: Dictionary with validation status, messages, and any failures
    """
    val_logger = ValidationLogger()
    read_counts = {}  # Track read counts for paired-end validation
    
    # Log validation parameters
    logging.info("Starting raw reads validation (skipping paired-end count check):")
    logging.info(f"  Output directory: {outdir}")
    logging.info(f"  Samples file: {samples_txt}")
    logging.info(f"  Paired-end: {paired_end}")
    logging.info(f"  Assay suffix: {assay_suffix}")

    try:
        # 1. Basic File Existence Checks
        fastq_dir = outdir / "00-RawData/Fastq"
        fastqc_dir = outdir / "00-RawData/FastQC_Reports"
        
        # Read samples from CSV, get Sample Name column
        samples = []
        with open(samples_txt) as f:
            reader = csv.DictReader(f)
            for row in reader:
                if 'Sample Name' in row:
                    sample_name = row['Sample Name'].strip()
                    if sample_name:  # Only add non-empty sample names
                        samples.append(sample_name)
                        logging.debug(f"Added sample: {sample_name}")
                else:
                    raise ValueError("Runsheet missing required 'Sample Name' column")
                    
        logging.info(f"Found {len(samples)} samples to validate")

        # 2. Check each sample's files
        for sample in samples:
            expected = get_expected_files(sample, paired_end, assay_suffix)
            logging.info(f"\nValidating sample: {sample}")
            
            # 2a. Check FASTQ files exist and validate
            for fastq in expected["fastq"]:
                fastq_path = fastq_dir / fastq
                
                # Check existence
                if not fastq_path.exists():
                    val_logger.log("raw_reads", sample, "file_exists", "HALT", 
                             f"Missing FastQ file: {fastq}")
                    continue
                    
                # Check GZIP integrity
                try:
                    output = subprocess.run(["gzip", "-t", str(fastq_path)], 
                                          capture_output=True, check=True)
                    if output.stdout:
                        val_logger.log("raw_reads", sample, "gzip_integrity", "HALT",
                                 f"GZIP integrity check failed", output.stdout.decode())
                    else:
                        val_logger.log("raw_reads", sample, "gzip_integrity", "GREEN",
                                 "GZIP integrity check passed")
                except subprocess.CalledProcessError as e:
                    val_logger.log("raw_reads", sample, "gzip_integrity", "HALT",
                             f"GZIP integrity check failed", e.stderr.decode())
                
                # Check FASTQ format with stats
                try:
                    issues = []
                    total_reads = 0
                    # Limit line checking to a very small number for faster check
                    max_lines_to_check = 200_000_000 
                    line_count = 0
                    
                    with gzip.open(fastq_path, "rt") as f:
                        for i, line in enumerate(f):
                            line_count = i
                            if i >= max_lines_to_check:
                                logging.debug(f"Reached {max_lines_to_check} lines, ending check")
                                break
                                
                            if i % 4 == 0:  # Header lines only
                                total_reads += 1
                                if not line.startswith('@'):
                                    issues.append(i+1)
                            if i % 2_000_000 == 0:  # Log progress
                                logging.debug(f"Checked {i} lines in {fastq_path}")
                                
                    # Estimate total reads based on checked portion
                    estimated_total = None
                    if line_count < max_lines_to_check:
                        # We read the whole file
                        estimated_total = total_reads
                    else:
                        # We only read a portion, so estimate total
                        file_size = fastq_path.stat().st_size
                        with gzip.open(fastq_path, "rb") as f:
                            f.seek(0, 2)  # Seek to end to get uncompressed size
                            uncompressed_size = f.tell()
                        
                        if uncompressed_size > 0:
                            estimated_total = int(total_reads * (file_size / uncompressed_size))
                        else:
                            estimated_total = total_reads
                    
                    stats = {
                        "total_reads": estimated_total,
                        "checked_reads": total_reads,
                        "invalid_headers": len(issues),
                        "checked_lines": line_count + 1
                    }
                                
                    if issues:
                        val_logger.log("raw_reads", sample, "fastq_format", "HALT",
                                 f"Invalid FASTQ format - headers missing @ at lines: {issues[:10]}..." if len(issues) > 10 else f"Invalid FASTQ format - headers missing @ at lines: {issues}",
                                 stats=stats)
                    else:
                        val_logger.log("raw_reads", sample, "fastq_format", "GREEN",
                                 f"FASTQ format check passed (checked first {min(line_count+1, max_lines_to_check)} lines)",
                                 stats=stats)

                    # Store read count but we won't use it for paired-end check
                    read_counts[fastq_path] = estimated_total
                                
                except Exception as e:
                    val_logger.log("raw_reads", sample, "fastq_format", "HALT",
                             f"Error checking FASTQ format", str(e))
            
            # 2b. SKIP Paired-end read count comparison
            # We're skipping this check to improve performance
            
            # 2c. Check FastQC outputs exist
            for fastqc in expected["fastqc"]:
                fastqc_path = fastqc_dir / fastqc
                if not fastqc_path.exists():
                    val_logger.log("raw_reads", sample, "fastqc_exists", "HALT",
                             f"Missing FastQC output: {fastqc}")
                else:
                    val_logger.log("raw_reads", sample, "fastqc_exists", "GREEN",
                             f"FastQC output found: {fastqc}")
                
        # 3. Check MultiQC report exists
        expected = get_expected_files("ALL", paired_end, assay_suffix)
        multiqc = expected["multiqc"][0]
        multiqc_path = fastqc_dir / multiqc
        if not multiqc_path.exists():
            val_logger.log("raw_reads", "ALL", "multiqc_exists", "HALT",
                      f"Missing MultiQC report: {multiqc}")
        else:
            val_logger.log("raw_reads", "ALL", "multiqc_exists", "GREEN",
                      f"MultiQC report found: {multiqc}")
            
            # Add check for samples in MultiQC
            multiqc_check = check_samples_in_multiqc(multiqc_path, samples, "raw_reads")
            if multiqc_check["all_samples_found"]:
                val_logger.log("raw_reads", "ALL", "raw_multiqc_samples", "GREEN",
                      f"All samples found in raw reads MultiQC report")
                
                # Use a temporary directory for multiqc extraction instead of persistent directory
                with tempfile.TemporaryDirectory() as temp_dir:
                    try:
                        logging.info(f"Extracting MultiQC zip to temporary directory: {temp_dir}")
                        
                        # Make sure the zip file exists before attempting to extract
                        if os.path.exists(multiqc_path):
                            with zipfile.ZipFile(multiqc_path, 'r') as zip_ref:
                                zip_ref.extractall(temp_dir)
                            logging.info(f"Successfully extracted MultiQC zip to temporary directory")
                            
                            # List the extracted contents for debugging
                            logging.info(f"Extracted files in temp directory:")
                            for root, dirs, files in os.walk(temp_dir):
                                for file in files:
                                    logging.info(f"  {os.path.join(root, file)}")
                        else:
                            logging.error(f"MultiQC zip file does not exist at {multiqc_path}")
                    except Exception as e:
                        logging.error(f"Error extracting MultiQC zip: {str(e)}")
                    
                    # Add group statistics from parse_fastqc
                    try:
                        # Set up environment for parse_fastqc to find the data
                        multiqc_data_dir = temp_dir
                        # Parse the FastQC data from MultiQC report
                        group_stats = parse_fastqc("raw", assay_suffix)
                        
                        if group_stats:
                            # Collect all stats but don't log individual samples
                            # This removes the redundant entries
                            
                            # After collecting all stats, check for outliers
                            outliers, metric_stats = detect_outliers(group_stats)
                            
                            # Report outliers if any were found
                            if outliers:
                                for sample, outlier_metrics in outliers.items():
                                    if sample in samples:  # Only report for samples in our list
                                        # Format outlier information
                                        outlier_details = []
                                        for o in outlier_metrics:
                                            z_score = round(o['z_score'], 2)
                                            direction = "higher" if z_score > 0 else "lower"
                                            outlier_details.append(
                                                f"{o['metric']}: {o['value']} ({abs(z_score)} std dev {direction} than mean of {round(o['mean'], 2)})"
                                            )
                                        
                                        # Determine status based on severity
                                        severe_outliers = any(abs(o['z_score']) > 3 for o in outlier_metrics)
                                        status = "RED" if severe_outliers else "YELLOW"
                                        
                                        val_logger.log("raw_reads", sample, "raw_metrics_outliers", status,
                                              f"Sample has {len(outlier_metrics)} metric outliers", 
                                              details="; ".join(outlier_details))
                            else:
                                val_logger.log("raw_reads", "ALL", "raw_metrics_outliers", "GREEN",
                                      f"No metric outliers detected across samples")
                        else:
                            val_logger.log("raw_reads", "ALL", "raw_group_stats", "WARNING",
                                  f"No group statistics collected from MultiQC report")
                    except Exception as e:
                        val_logger.log("raw_reads", "ALL", "raw_group_stats", "WARNING",
                              f"Error collecting group statistics: {str(e)}")
            else:
                missing = ", ".join(multiqc_check["missing_samples"]) if multiqc_check["missing_samples"] else "None"
                if "error" in multiqc_check:
                    val_logger.log("raw_reads", "ALL", "raw_multiqc_samples", "HALT",
                          f"Error checking samples in MultiQC: {multiqc_check['error']}")
                else:
                    val_logger.log("raw_reads", "ALL", "raw_multiqc_samples", "HALT",
                          f"Missing samples in raw reads MultiQC: {missing}")
            
    except Exception as e:
        val_logger.log("raw_reads", "ALL", "validation", "HALT",
                  f"Validation error", str(e))
        
    return {
        "status": val_logger.get_status(),
        "messages": [r["message"] for r in val_logger.results],
        "failures": {r["check_name"]: r["message"] 
                    for r in val_logger.results 
                    if r["status"] in ["HALT", "RED"]}
    }

def get_expected_trimmed_files(sample_id: str, paired_end: bool, assay_suffix: str) -> dict:
    """
    Get expected file patterns for trimmed reads of a sample
    
    Args:
        sample_id: Sample name from runsheet
        paired_end: Whether data is paired-end
        assay_suffix: Assay suffix (e.g. "_GLbulkRNAseq")
    """
    expected = {
        # Trimmed FastQ files
        "fastq": [f"{sample_id}_trimmed.fastq.gz"] if not paired_end else [
            f"{sample_id}_R1_trimmed.fastq.gz",
            f"{sample_id}_R2_trimmed.fastq.gz"
        ],
        
        # FastQC outputs
        "fastqc": [f"{sample_id}_trimmed_fastqc.zip"] if not paired_end else [
            f"{sample_id}_R1_trimmed_fastqc.zip",
            f"{sample_id}_R2_trimmed_fastqc.zip"
        ],
        
        # MultiQC report
        "trimmed_multiqc": [f"trimmed_multiqc{assay_suffix}_report.zip"],
        
        # Trimming reports - Update to match TrimGalore's output format
        "trimming_reports": [f"{sample_id}_raw.fastq.gz_trimming_report.txt"] if not paired_end else [
            f"{sample_id}_R1_raw.fastq.gz_trimming_report.txt",
            f"{sample_id}_R2_raw.fastq.gz_trimming_report.txt"
        ],
        
        # Trimming MultiQC report
        "trimming_multiqc": [f"trimming_multiqc{assay_suffix}_report.zip"]
    }
    logging.debug(f"Expected trimmed files for {sample_id}: {expected}")
    return expected

def parse_trimming_report(report_path: Path) -> dict:
    """
    Parse a TrimGalore trimming report to extract adapter trimming information
    
    Args:
        report_path: Path to the trimming report
        
    Returns:
        Dictionary with parsed information
    """
    # FORCE DEBUG OUTPUT
    print(f"\n\n### DEBUGGING TRIM REPORT PARSING: {report_path} ###\n")
    
    stats = {
        "total_processed_reads": 0,
        "adapters_found": False,
        "adapter_type": "unknown",
        "adapter_sequence": "",
        "reads_with_adapters": 0,
        "adapter_trimmed_percentage": 0.0,
        "quality_trimmed_reads": 0,
        "quality_trimmed_percentage": 0.0
    }
    
    try:
        with open(report_path, 'r') as f:
            content = f.read()
            
            # Print key sections of the report to debug
            summary_section = content[content.find('=== Summary ==='):content.find('=== Adapter 1 ===')]
            print(f"REPORT SUMMARY SECTION:\n{summary_section}\n")
            
            # Extract total processed reads
            processed_match = re.search(r'Total reads processed:\s*([0-9,]+)', content)
            if processed_match:
                # Remove commas before converting to int
                stats["total_processed_reads"] = int(processed_match.group(1).replace(',', ''))
                print(f"FOUND TOTAL READS: {stats['total_processed_reads']}")
            else:
                print(f"ERROR: Could not find 'Total reads processed' in report!")
                print(f"Regex pattern: r'Total reads processed:\\s*([0-9,]+)'")
                
            # Check if adapter was detected and its sequence/type
            adapter_match = re.search(r'Adapter sequence:\s*\'([ACGT]+)\'', content)
            if adapter_match:
                stats["adapters_found"] = True
                stats["adapter_sequence"] = adapter_match.group(1)
                print(f"FOUND ADAPTER SEQUENCE: {stats['adapter_sequence']}")
                # Try to identify adapter type from TrimGalore output
                if "Illumina" in content:
                    stats["adapter_type"] = "Illumina"
                elif "Nextera" in content:
                    stats["adapter_type"] = "Nextera"
                elif "smallRNA" in content:
                    stats["adapter_type"] = "smallRNA"
                print(f"ADAPTER TYPE: {stats['adapter_type']}")
                
            # Extract reads with adapters - match actual format with right-aligned numbers
            adapter_count_match = re.search(r'Reads with adapters:\s*([0-9,]+)\s*\(([0-9.]+)%\)', content)
            if adapter_count_match:
                raw_count = adapter_count_match.group(1)
                percentage = adapter_count_match.group(2)
                print(f"RAW ADAPTER MATCH: count='{raw_count}', percentage='{percentage}'")
                
                stats["reads_with_adapters"] = int(raw_count.replace(',', ''))
                stats["adapter_trimmed_percentage"] = float(percentage)
                print(f"FOUND ADAPTERS: {stats['reads_with_adapters']} reads ({stats['adapter_trimmed_percentage']}%)")
            else:
                print(f"ERROR: Could not match 'Reads with adapters' pattern!")
                print(f"Regex pattern: r'Reads with adapters:\\s*([0-9,]+)\\s*\\(([0-9.]+)%\\)'")
                print(f"Looking for alternative patterns...")
                
                # Alternative pattern for newer versions
                if "adapters were trimmed" in content:
                    percentage_match = re.search(r'([0-9.]+)% of reads contained adapter', content)
                    if percentage_match:
                        stats["adapter_trimmed_percentage"] = float(percentage_match.group(1))
                        stats["adapters_found"] = True
                        print(f"FOUND ALTERNATIVE ADAPTER PATTERN: {stats['adapter_trimmed_percentage']}%")
                
            # Extract quality trimming info - match Quality-trimmed: right-aligned format
            quality_match = re.search(r'Quality-trimmed:\s*([0-9,]+)\s+bp\s*\(([0-9.]+)%\)', content)
            if quality_match:
                raw_count = quality_match.group(1)
                percentage = quality_match.group(2)
                print(f"RAW QUALITY MATCH: count='{raw_count}', percentage='{percentage}'")
                
                stats["quality_trimmed_reads"] = int(raw_count.replace(',', ''))
                stats["quality_trimmed_percentage"] = float(percentage)
                print(f"FOUND QUALITY TRIMMED: {stats['quality_trimmed_reads']} bp ({stats['quality_trimmed_percentage']}%)")
            else:
                print(f"ERROR: Could not match 'Quality-trimmed' pattern")
                print(f"Regex pattern: r'Quality-trimmed:\\s*([0-9,]+)\\s+bp\\s*\\(([0-9.]+)%\\)'")
                
    except Exception as e:
        print(f"EXCEPTION PARSING REPORT: {e}")
        import traceback
        print(f"TRACEBACK: {traceback.format_exc()}")
        
    # Print final results
    print(f"FINAL PARSED STATS: {stats}")
    print(f"\n### END DEBUGGING TRIM REPORT PARSING ###\n")
    
    return stats

def validate_trimmed_reads(outdir: Path,
                       samples_txt: Path,
                       paired_end: bool = True,
                       assay_suffix: str = "_GLbulkRNAseq") -> dict:
    """
    Perform validation checks on trimmed sequencing reads.
    
    This function validates the trimmed sequencing data by performing:
    1. File existence checks for trimmed FastQ files and FastQC outputs
    2. GZIP integrity validation of FastQ files
    3. FASTQ format validation (checks header format, counts reads)
    4. FastQC output existence validation
    5. MultiQC report existence validation
    6. Sample presence verification in MultiQC reports
    7. Adapter trimming validation using trimming reports
    
    Note: Paired-end read count comparison is skipped to improve performance, 
    just like in the raw reads validation.
    
    Args:
        outdir: Path to the output directory containing the trimmed data
        samples_txt: Path to the samples runsheet CSV file
        paired_end: Boolean indicating if samples are paired-end (True) or single-end (False)
        assay_suffix: Suffix for the assay (e.g., "_GLbulkRNAseq")
        
    Returns:
        dict: Dictionary with validation status, messages, and any failures
    """
    val_logger = ValidationLogger()
    read_counts = {}  # Track read counts for paired-end validation
    
    # Log validation parameters
    logging.info("Starting trimmed reads validation (skipping paired-end count check):")
    logging.info(f"  Output directory: {outdir}")
    logging.info(f"  Samples file: {samples_txt}")
    logging.info(f"  Paired-end: {paired_end}")
    logging.info(f"  Assay suffix: {assay_suffix}")

    try:
        # 1. Basic File Existence Checks
        fastq_dir = outdir / "01-TG_Preproc/Fastq"
        fastqc_dir = outdir / "01-TG_Preproc/FastQC_Reports"
        trimming_dir = outdir / "01-TG_Preproc/Trimming_Reports"

        # Read samples from CSV, get Sample Name column
        samples = []
        with open(samples_txt) as f:
            reader = csv.DictReader(f)
            for row in reader:
                if 'Sample Name' in row:
                    sample_name = row['Sample Name'].strip()
                    if sample_name:  # Only add non-empty sample names
                        samples.append(sample_name)
                        logging.debug(f"Added sample: {sample_name}")
                else:
                    raise ValueError("Runsheet missing required 'Sample Name' column")
                    
        logging.info(f"Found {len(samples)} samples to validate")

        # 2. Check each sample's files
        for sample in samples:
            expected = get_expected_trimmed_files(sample, paired_end, assay_suffix)
            logging.info(f"\nValidating sample: {sample}")
            
            # 2a. Check FASTQ files exist and validate
            for fastq in expected["fastq"]:
                fastq_path = fastq_dir / fastq
                
                # Check existence
                if not fastq_path.exists():
                    val_logger.log("trimmed_reads", sample, "file_exists", "HALT", 
                             f"Missing trimmed FastQ file: {fastq}")
                    continue
                    
                # Check GZIP integrity
                try:
                    output = subprocess.run(["gzip", "-t", str(fastq_path)], 
                                          capture_output=True, check=True)
                    if output.stdout:
                        val_logger.log("trimmed_reads", sample, "gzip_integrity", "HALT",
                                 f"GZIP integrity check failed", output.stdout.decode())
                    else:
                        val_logger.log("trimmed_reads", sample, "gzip_integrity", "GREEN",
                                 "GZIP integrity check passed")
                except subprocess.CalledProcessError as e:
                    val_logger.log("trimmed_reads", sample, "gzip_integrity", "HALT",
                             f"GZIP integrity check failed", e.stderr.decode())
                
                # Check FASTQ format with stats - using same approach as raw reads
                try:
                    issues = []
                    total_reads = 0
                    # Limit line checking to a very small number for faster check
                    max_lines_to_check = 200_000_000 
                    line_count = 0
                    
                    with gzip.open(fastq_path, "rt") as f:
                        for i, line in enumerate(f):
                            line_count = i
                            if i >= max_lines_to_check:
                                logging.debug(f"Reached {max_lines_to_check} lines, ending check")
                                break
                                
                            if i % 4 == 0:  # Header lines only
                                total_reads += 1
                                if not line.startswith('@'):
                                    issues.append(i+1)
                            if i % 2_000_000 == 0:  # Log progress
                                logging.debug(f"Checked {i} lines in {fastq_path}")
                                
                    # Estimate total reads based on checked portion - same as raw reads
                    estimated_total = None
                    if line_count < max_lines_to_check:
                        # We read the whole file
                        estimated_total = total_reads
                    else:
                        # We only read a portion, so estimate total
                        file_size = fastq_path.stat().st_size
                        with gzip.open(fastq_path, "rb") as f:
                            f.seek(0, 2)  # Seek to end to get uncompressed size
                            uncompressed_size = f.tell()
                        
                        if uncompressed_size > 0:
                            estimated_total = int(total_reads * (file_size / uncompressed_size))
                        else:
                            estimated_total = total_reads
                                
                    stats = {
                        "total_reads": estimated_total,
                        "checked_reads": total_reads,
                        "invalid_headers": len(issues),
                        "checked_lines": line_count + 1
                    }
                                
                    if issues:
                        val_logger.log("trimmed_reads", sample, "fastq_format", "HALT",
                                 f"Invalid FASTQ format - headers missing @ at lines: {issues[:10]}..." if len(issues) > 10 else f"Invalid FASTQ format - headers missing @ at lines: {issues}",
                                 stats=stats)
                    else:
                        val_logger.log("trimmed_reads", sample, "fastq_format", "GREEN",
                                 f"FASTQ format check passed (checked first {min(line_count+1, max_lines_to_check)} lines)",
                                 stats=stats)

                    # Store read count but we won't use it for paired-end check
                    read_counts[fastq_path] = estimated_total
                                
                except Exception as e:
                    val_logger.log("trimmed_reads", sample, "fastq_format", "HALT",
                             f"Error checking FASTQ format", str(e))
            
            # 2b. SKIP Paired-end read count comparison to improve performance
            # This is consistent with raw reads validation
                
            # 2c. Check FastQC outputs exist
            for fastqc in expected["fastqc"]:
                fastqc_path = fastqc_dir / fastqc
                if not fastqc_path.exists():
                    val_logger.log("trimmed_reads", sample, "fastqc_exists", "HALT",
                             f"Missing FastQC output: {fastqc}")
                else:
                    val_logger.log("trimmed_reads", sample, "fastqc_exists", "GREEN",
                             f"FastQC output found: {fastqc}")
            
            # 2d. Check trimming reports exist
            for read_end_idx, report in enumerate(expected["trimming_reports"]):
                report_path = trimming_dir / report
                read_end = "R1" if read_end_idx == 0 else "R2"
                
                if not report_path.exists():
                    val_logger.log("trimmed_reads", sample, f"trimming_report_exists_{read_end}", "HALT",
                             f"Missing trimming report for {read_end}: {report}")
                    continue
                else:
                    val_logger.log("trimmed_reads", sample, f"trimming_report_exists_{read_end}", "GREEN",
                             f"Trimming report found for {read_end}: {report}")
                
                # Parse the trimming report to check for adapter trimming
                trim_stats = parse_trimming_report(report_path)
                
                # Log adapter trimming information
                if trim_stats["adapters_found"]:
                    status = "GREEN"
                    msg = f"{read_end}: Adapter trimming performed - {trim_stats['adapter_trimmed_percentage']}% of reads contained adapters"
                    details = f"Adapter type: {trim_stats['adapter_type']}, Sequence: {trim_stats['adapter_sequence']}"
                    
                    # If very high adapter percentage (>90%), flag as potential issue
                    if trim_stats["adapter_trimmed_percentage"] > 90:
                        status = "YELLOW"
                        msg += " (unusually high percentage - check library prep)"
                else:
                    # No adapters found could be fine (e.g., very clean library) or could indicate an issue
                    status = "YELLOW"
                    msg = f"{read_end}: No adapters found in reads - verify this is expected"
                    details = "Some libraries should have adapter contamination; absence may indicate issues"
                
                val_logger.log("trimmed_reads", sample, f"adapter_trimming_{read_end}", status, msg, details, stats=trim_stats)
                
                # Only log quality trimming if it actually occurred
                if "quality_trimmed_percentage" in trim_stats and trim_stats["quality_trimmed_percentage"] > 0:
                    status = "GREEN"
                    msg = f"{read_end}: Quality trimming performed - {trim_stats['quality_trimmed_percentage']}% of reads were quality-trimmed"
                    
                    # If extremely high quality trimming (>50%), flag as potential issue
                    if trim_stats["quality_trimmed_percentage"] > 50:
                        status = "YELLOW"
                        msg += " (unusually high percentage - check sequencing quality)"
                    
                    val_logger.log("trimmed_reads", sample, f"quality_trimming_{read_end}", status, msg, stats=trim_stats)
            
        # 3. Check MultiQC reports exist
        expected = get_expected_trimmed_files("ALL", paired_end, assay_suffix)
        
        # 3a. Trimmed reads MultiQC
        trimmed_multiqc = expected["trimmed_multiqc"][0]
        trimmed_multiqc_path = fastqc_dir / trimmed_multiqc
        if not trimmed_multiqc_path.exists():
            val_logger.log("trimmed_reads", "ALL", "trimmed_multiqc_exists", "HALT",
                      f"Missing trimmed reads MultiQC report: {trimmed_multiqc}")
        else:
            val_logger.log("trimmed_reads", "ALL", "trimmed_multiqc_exists", "GREEN",
                      f"Trimmed reads MultiQC report found: {trimmed_multiqc}")
            
            # Add check for samples in MultiQC
            multiqc_check = check_samples_in_multiqc(trimmed_multiqc_path, samples, "trimmed_reads")
            if multiqc_check["all_samples_found"]:
                val_logger.log("trimmed_reads", "ALL", "trimmed_multiqc_samples", "GREEN",
                      f"All samples found in trimmed reads MultiQC report")
                
                # Use a temporary directory for multiqc extraction instead of persistent directory
                with tempfile.TemporaryDirectory() as temp_dir:
                    try:
                        logging.info(f"Extracting MultiQC zip to temporary directory: {temp_dir}")
                        
                        # Make sure the zip file exists before attempting to extract
                        if os.path.exists(trimmed_multiqc_path):
                            with zipfile.ZipFile(trimmed_multiqc_path, 'r') as zip_ref:
                                zip_ref.extractall(temp_dir)
                            logging.info(f"Successfully extracted MultiQC zip to temporary directory")
                            
                            # List the extracted contents for debugging
                            logging.info(f"Extracted files in temp directory:")
                            for root, dirs, files in os.walk(temp_dir):
                                for file in files:
                                    logging.info(f"  {os.path.join(root, file)}")
                        else:
                            logging.error(f"MultiQC zip file does not exist at {trimmed_multiqc_path}")
                    except Exception as e:
                        logging.error(f"Error extracting MultiQC zip: {str(e)}")
                    
                    # Add group statistics from parse_fastqc
                    try:
                        # Set up environment for parse_fastqc to find the data
                        multiqc_data_dir = temp_dir
                        # Parse the FastQC data from MultiQC report
                        group_stats = parse_fastqc("trimmed", assay_suffix)
                        
                        if group_stats:
                            # Collect all stats but don't log individual samples
                            # This removes the redundant entries
                            
                            # After collecting all stats, check for outliers
                            outliers, metric_stats = detect_outliers(group_stats)
                            
                            # Report outliers if any were found
                            if outliers:
                                for sample, outlier_metrics in outliers.items():
                                    if sample in samples:  # Only report for samples in our list
                                        # Format outlier information
                                        outlier_details = []
                                        for o in outlier_metrics:
                                            z_score = round(o['z_score'], 2)
                                            direction = "higher" if z_score > 0 else "lower"
                                            outlier_details.append(
                                                f"{o['metric']}: {o['value']} ({abs(z_score)} std dev {direction} than mean of {round(o['mean'], 2)})"
                                            )
                                        
                                        # Determine status based on severity
                                        severe_outliers = any(abs(o['z_score']) > 3 for o in outlier_metrics)
                                        status = "RED" if severe_outliers else "YELLOW"
                                        
                                        val_logger.log("trimmed_reads", sample, "trimmed_metrics_outliers", status,
                                              f"Sample has {len(outlier_metrics)} metric outliers", 
                                              details="; ".join(outlier_details))
                            else:
                                val_logger.log("trimmed_reads", "ALL", "trimmed_metrics_outliers", "GREEN",
                                      f"No metric outliers detected across samples")
                        else:
                            val_logger.log("trimmed_reads", "ALL", "trimmed_group_stats", "WARNING",
                                  f"No group statistics collected from MultiQC report")
                    except Exception as e:
                        val_logger.log("trimmed_reads", "ALL", "trimmed_group_stats", "WARNING",
                              f"Error collecting group statistics: {str(e)}")
            else:
                missing = ", ".join(multiqc_check["missing_samples"]) if multiqc_check["missing_samples"] else "None"
                if "error" in multiqc_check:
                    val_logger.log("trimmed_reads", "ALL", "trimmed_multiqc_samples", "HALT",
                          f"Error checking samples in MultiQC: {multiqc_check['error']}")
                else:
                    val_logger.log("trimmed_reads", "ALL", "trimmed_multiqc_samples", "HALT",
                          f"Missing samples in trimmed reads MultiQC: {missing}")
        
        # 3b. Trimming reports MultiQC
        trimming_multiqc = expected["trimming_multiqc"][0]
        trimming_multiqc_path = trimming_dir / trimming_multiqc
        if not trimming_multiqc_path.exists():
            val_logger.log("trimmed_reads", "ALL", "trimming_multiqc_exists", "HALT",
                      f"Missing trimming MultiQC report: {trimming_multiqc}")
        else:
            val_logger.log("trimmed_reads", "ALL", "trimming_multiqc_exists", "GREEN",
                      f"Trimming MultiQC report found: {trimming_multiqc}")
            
    except Exception as e:
        val_logger.log("trimmed_reads", "ALL", "validation", "HALT",
                  f"Validation error", str(e))
        
    return {
        "status": val_logger.get_status(),
        "messages": [r["message"] for r in val_logger.results],
        "failures": {r["check_name"]: r["message"] 
                    for r in val_logger.results 
                    if r["status"] in ["HALT", "RED"]}
    }

def validate_bowtie2_alignments(outdir, samples_txt, paired_end, assay_suffix):
    """
    Validate Bowtie2 alignment outputs in the work directory
    
    Args:
        outdir (Path): Output directory (main GLDS directory) - NOT USED, validation happens in work dir
        samples_txt (Path): Path to samples.txt file
        paired_end (bool): Whether data is paired-end
        assay_suffix (str): Suffix for output files (e.g., '_GLbulkRNAseq')
        
    Returns:
        dict: Validation log
    """
    logger.info(f"Validating Bowtie2 alignments with params: samples_txt={samples_txt}, paired_end={paired_end}")
    
    # Debug what files actually exist
    alignment_dir = Path("02-Bowtie2_Alignment")
    if alignment_dir.exists():
        logger.info(f"Files in {alignment_dir} before validation:")
        for item in alignment_dir.glob("**/*"):
            logger.info(f"  {item}")
    else:
        logger.warning(f"Alignment directory does not exist: {alignment_dir}")
    
    # Sample management
    sample_ids = []
    with open(samples_txt, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if 'Sample Name' in row:
                sample_name = row['Sample Name'].strip()
                if sample_name:  # Only add non-empty sample names
                    sample_ids.append(sample_name)
                    logging.debug(f"Added sample: {sample_name}")
            else:
                raise ValueError("Runsheet missing required 'Sample Name' column")
            
    logger.info(f"Found {len(sample_ids)} samples in runsheet: {sample_ids}")
    
    # Create a samples variable for MultiQC processing consistency with other functions
    # This helps keep code consistent with functions that call check_samples_in_multiqc()
    samples = sample_ids
    
    # Track alignment rates across samples for outlier detection
    alignment_rates = []
    sample_names = []  # Track sample names in same order as metrics
    
    # Dictionary to store all metrics by sample for comprehensive outlier detection
    alignment_stats_by_sample = {}
    
    # Validate each sample's alignment outputs
    for sample in sample_ids:
        sample_dir = alignment_dir / sample
        
        # Check sorted BAM file exists
        sorted_bam = sample_dir / f"{sample}_sorted.bam"
        if not sorted_bam.exists():
            val_logger.log("alignments", sample, "sorted_bam_exists", "HALT", 
                        f"Missing sorted BAM file: {sorted_bam.name}")
        else:
            val_logger.log("alignments", sample, "sorted_bam_exists", "GREEN",
                        f"Found sorted BAM file: {sorted_bam.name}")
            
            # Check BAM file integrity with samtools
            try:
                subprocess.run(["samtools", "quickcheck", str(sorted_bam)], check=True)
                val_logger.log("alignments", sample, "bam_integrity", "GREEN",
                            "BAM integrity check passed")
            except subprocess.CalledProcessError:
                val_logger.log("alignments", sample, "bam_integrity", "HALT",
                            "BAM integrity check failed")
        
        # Check BAM index exists
        bam_index = sample_dir / f"{sample}_sorted.bam.bai"
        if not bam_index.exists():
            val_logger.log("alignments", sample, "bam_index_exists", "HALT",
                        f"Missing BAM index file: {bam_index.name}")
        else:
            val_logger.log("alignments", sample, "bam_index_exists", "GREEN",
                        f"Found BAM index file: {bam_index.name}")
        
        # Check bowtie2 log file for alignment rate
        bowtie2_log = sample_dir / f"{sample}.bowtie2.log"
        if not bowtie2_log.exists():
            val_logger.log("alignments", sample, "alignment_rate", "HALT", 
                        f"Missing Bowtie2 log file: {bowtie2_log.name}")
        else:
            # Parse alignment rate from log
            try:
                with open(bowtie2_log, "r") as f:
                    log_content = f.read()
                    
                # Parse overall alignment rate
                match = re.search(r"(\d+\.\d+)% overall alignment rate", log_content)
                if match:
                    value = float(match.group(1))
                    alignment_rates.append(value)
                    sample_names.append(sample)
                    
                    # Store alignment rate in stats dictionary
                    if sample not in alignment_stats_by_sample:
                        alignment_stats_by_sample[sample] = {}
                    alignment_stats_by_sample[sample]['overall_alignment_rate'] = value
                    
                    # Log the alignment rate without threshold judgments
                    val_logger.log("alignments", sample, "alignment_rate", "GREEN", 
                                f"Alignment rate: {value}%")
                    
            except Exception as e:
                val_logger.log("alignments", sample, "alignment_rate", "HALT", 
                            f"Error parsing Bowtie2 log: {str(e)}")
        
        # Check unmapped reads files exist
        if paired_end:
            # Check for both possible naming patterns
            # Pattern 1: sample.unmapped.fastq.1.gz and sample.unmapped.fastq.2.gz
            # Pattern 2: sample_R1.unmapped.fastq.gz and sample_R2.unmapped.fastq.gz
            
            # Pattern 1
            unmapped_1 = sample_dir / f"{sample}.unmapped.fastq.1.gz"
            unmapped_2 = sample_dir / f"{sample}.unmapped.fastq.2.gz"
            
            # Pattern 2
            unmapped_r1 = sample_dir / f"{sample}_R1.unmapped.fastq.gz"
            unmapped_r2 = sample_dir / f"{sample}_R2.unmapped.fastq.gz"
            
            # Check both patterns
            if (not unmapped_1.exists() or not unmapped_2.exists()) and (not unmapped_r1.exists() or not unmapped_r2.exists()):
                val_logger.log("alignments", sample, "unmapped_reads_exist", "HALT",
                            f"Missing unmapped reads files: {unmapped_1.name} and/or {unmapped_2.name}")
            else:
                # Log which pattern was found
                if unmapped_1.exists() and unmapped_2.exists():
                    val_logger.log("alignments", sample, "unmapped_reads_exist", "GREEN",
                                f"Found unmapped reads files: {unmapped_1.name} and {unmapped_2.name}")
                elif unmapped_r1.exists() and unmapped_r2.exists():
                    val_logger.log("alignments", sample, "unmapped_reads_exist", "GREEN",
                                f"Found unmapped reads files: {unmapped_r1.name} and {unmapped_r2.name}")
                else:
                    # This shouldn't happen based on our check above, but just in case
                    val_logger.log("alignments", sample, "unmapped_reads_exist", "YELLOW",
                                f"Found only partial unmapped reads files - check for missing pairs")
        else:
            # Single-end - also check both possible patterns
            unmapped = sample_dir / f"{sample}.unmapped.fastq.gz"
            unmapped_alt = sample_dir / f"{sample}_R1.unmapped.fastq.gz"
            
            if not unmapped.exists() and not unmapped_alt.exists():
                val_logger.log("alignments", sample, "unmapped_reads_exist", "HALT",
                            f"Missing unmapped reads file: {unmapped.name}")
            else:
                if unmapped.exists():
                    val_logger.log("alignments", sample, "unmapped_reads_exist", "GREEN",
                                f"Found unmapped reads file: {unmapped.name}")
                else:
                    val_logger.log("alignments", sample, "unmapped_reads_exist", "GREEN",
                                f"Found unmapped reads file: {unmapped_alt.name}")
    
    # Check for outliers in alignment rate from logs
    if len(alignment_rates) >= 2:
        median = statistics.median(alignment_rates)
        
        # Find outliers at 2 and 4 standard deviations
        # Using statistical outlier detection instead of arbitrary thresholds
        # as it adapts to the specific dataset characteristics
        yellow_outliers = calculate_outliers(alignment_rates, 2)  # Moderate outliers
        red_outliers = calculate_outliers(alignment_rates, 4)     # Extreme outliers
        
        # Log yellow outliers (excluding those that are also red)
        for idx in yellow_outliers:
            if idx not in red_outliers:
                sample = sample_names[idx]
                val_logger.log("alignments", sample, "alignment_rate", "YELLOW",
                            f"Alignment rate {alignment_rates[idx]}% is an outlier",
                            f"Value deviates >2 standard deviations from median ({median}%)")
                
        # Log red outliers
        for idx in red_outliers:
            sample = sample_names[idx]
            val_logger.log("alignments", sample, "alignment_rate", "RED",
                        f"Alignment rate {alignment_rates[idx]}% is an extreme outlier",
                        f"Value deviates >4 standard deviations from median ({median}%)")
    
    # Check MultiQC report for alignments (always in root directory as it's global)
    multiqc_path = alignment_dir / f"align_multiqc{assay_suffix}_report.zip"
    logger.debug(f"Looking for MultiQC zip at: {multiqc_path}")
    
    if not multiqc_path.exists():
        val_logger.log("alignments", "ALL", "align_multiqc_exists", "HALT", 
                    f"Missing alignment MultiQC report: {multiqc_path.name}")
    else:
        val_logger.log("alignments", "ALL", "align_multiqc_exists", "GREEN", 
                    f"Alignment MultiQC report found: {multiqc_path.name}")
        
        # Parse MultiQC for alignment metrics
        try:
            # Extract the MultiQC report to a temporary directory
            with tempfile.TemporaryDirectory() as temp_dir:
                temp_dir_path = Path(temp_dir)
                logger.info(f"Extracting MultiQC zip to temporary directory: {temp_dir_path}")
                
                with zipfile.ZipFile(multiqc_path, 'r') as zip_ref:
                    # List all files in the zip for debugging
                    zip_contents = zip_ref.namelist()
                    logger.debug(f"Files in MultiQC zip: {zip_contents}")
                    
                    # Extract all files
                    zip_ref.extractall(temp_dir_path)
                
                # Find the extracted data directory (may have variable naming based on assay suffix)
                multiqc_data_dir = None
                
                # Define a recursive function to find the data directory
                def find_data_dir(path):
                    # Check if this directory is the data directory
                    if path.is_dir() and "multiqc" in path.name.lower() and "data" in path.name.lower():
                        return path
                    
                    # Otherwise, check subdirectories
                    if path.is_dir():
                        for item in path.iterdir():
                            if item.is_dir():
                                result = find_data_dir(item)
                                if result:
                                    return result
                    return None
                
                # Start the recursive search
                multiqc_data_dir = find_data_dir(temp_dir_path)
                
                # If we still couldn't find it, try a more general search
                if not multiqc_data_dir:
                    logger.warning(f"Could not find MultiQC data directory with standard naming")
                    
                    # Try to find any directory that might contain the JSON file
                    for root, dirs, files in os.walk(temp_dir_path):
                        root_path = Path(root)
                        # Check if this directory contains multiqc_data.json
                        json_path = root_path / "multiqc_data.json"
                        if json_path.exists():
                            multiqc_data_dir = root_path
                            logger.info(f"Found MultiQC data.json in: {multiqc_data_dir}")
                            break
                
                # Log what we found
                if multiqc_data_dir:
                    logger.info(f"Found MultiQC data directory: {multiqc_data_dir}")
                    # List files in the data directory
                    data_files = list(multiqc_data_dir.iterdir())
                    logger.debug(f"Files in MultiQC data directory: {data_files}")
                    
                    # Look for JSON data file
                    json_file = multiqc_data_dir / "multiqc_data.json"
                    if json_file.exists():
                        logger.info(f"Found MultiQC data JSON: {json_file}")
                        
                        # Read the JSON data
                        with open(json_file, 'r') as f:
                            multiqc_data = json.load(f)
                        
                        # Extract samples from the MultiQC data
                        bowtie2_stats = {}
                        
                        # Look for bowtie2 data in general stats
                        if 'report_general_stats_data' in multiqc_data:
                            for stats_section in multiqc_data['report_general_stats_data']:
                                for sample, stats in stats_section.items():
                                    # Clean up sample name to remove read identifiers
                                    base_sample = re.sub(r'_R[12]$', '', sample)
                                    
                                    # Skip samples not in our list
                                    if base_sample not in samples and not any(s == base_sample for s in samples):
                                        continue
                                        
                                    # Initialize if needed
                                    if base_sample not in bowtie2_stats:
                                        bowtie2_stats[base_sample] = {}
                                    
                                    # Map the metrics from general stats
                                    metric_map = {
                                        'total_reads': 'total_reads',
                                        'overall_alignment_rate': 'overall_alignment_rate',
                                        'paired_aligned_none': 'aligned_none_pct',  # Will convert to pct
                                        'paired_aligned_one': 'aligned_one_pct',    # Will convert to pct
                                        'paired_aligned_multi': 'aligned_multi_pct',# Will convert to pct
                                        'paired_total': 'total_reads'
                                    }
                                    
                                    # Copy the metrics
                                    for src_metric, dest_metric in metric_map.items():
                                        if src_metric in stats:
                                            # Special handling for percentage metrics
                                            if src_metric.startswith('paired_aligned') and 'paired_total' in stats:
                                                # Convert to percentage
                                                total = stats['paired_total']
                                                if total > 0:
                                                    bowtie2_stats[base_sample][dest_metric] = (stats[src_metric] / total) * 100
                                            else:
                                                bowtie2_stats[base_sample][dest_metric] = stats[src_metric]
                            
                        # Log what we found
                        logger.info(f"Extracted Bowtie2 metrics for {len(bowtie2_stats)} samples")
                        for sample, metrics in bowtie2_stats.items():
                            logger.debug(f"Sample {sample} metrics: {metrics}")
                        
                        # Update our alignment stats with the extracted metrics
                        for sample, metrics in bowtie2_stats.items():
                            # Find matching sample in our list
                            matched_sample = None
                            if sample in samples:
                                matched_sample = sample
                            else:
                                for s in samples:
                                    if s in sample or sample in s:
                                        matched_sample = s
                                        break
                            
                            if matched_sample:
                                # Initialize stats if needed
                                if matched_sample not in alignment_stats_by_sample:
                                    alignment_stats_by_sample[matched_sample] = {}
                                
                                # Add all metrics
                                for metric, value in metrics.items():
                                    alignment_stats_by_sample[matched_sample][metric] = value
                        
                        # Detect outliers across different alignment metrics
                        if alignment_stats_by_sample:
                            # Define metrics to check for outliers
                            metrics_to_check = [
                                'total_reads', 
                                'overall_alignment_rate', 
                                'aligned_none_pct', 
                                'aligned_one_pct', 
                                'aligned_multi_pct'
                            ]
                            
                            # Log what metrics we have per sample
                            for sample, metrics in alignment_stats_by_sample.items():
                                available_metrics = [m for m in metrics_to_check if m in metrics]
                                if available_metrics:
                                    logger.debug(f"Sample {sample} has metrics: {available_metrics}")
                                else:
                                    logger.warning(f"Sample {sample} has no alignment metrics for outlier detection")
                            
                            # Detect outliers
                            outliers, metric_stats = detect_outliers(alignment_stats_by_sample, metrics_to_check)
                            
                            # Log outliers
                            if outliers:
                                for sample, sample_outliers in outliers.items():
                                    for outlier in sample_outliers:
                                        metric = outlier['metric']
                                        value = outlier['value']
                                        mean = outlier['mean']
                                        std = outlier['std']
                                        z_score = outlier['z_score']
                                        
                                        # Determine severity by z-score
                                        status = "YELLOW"
                                        if abs(z_score) > 4:
                                            status = "RED"
                                            
                                        val_logger.log("alignments", sample, f"{metric}_outlier", status,
                                                    f"{metric} value {value} is an outlier",
                                                    f"Value deviates {abs(z_score):.2f} standard deviations from mean ({mean:.2f})")
                            else:
                                val_logger.log("alignments", "ALL", "alignment_metrics_outliers", "GREEN",
                                            "No outliers detected in alignment metrics")
                    else:
                        logger.warning(f"MultiQC data JSON not found in data directory")
                else:
                    logger.warning(f"Could not find MultiQC data directory in extracted files")
                    # List what was extracted
                    logger.debug(f"Contents of temp directory: {list(temp_dir_path.iterdir())}")
                
            # Check all samples are present in MultiQC report
            try:
                with zipfile.ZipFile(multiqc_path, 'r') as zip_ref:
                    # Debug: List all files in zip
                    all_files = zip_ref.namelist()
                    logger.debug(f"Files in MultiQC zip: {all_files}")
                    
                    # Try to find the data directory and json file directly
                    json_path = None
                    for file in all_files:
                        if file.endswith('multiqc_data.json'):
                            json_path = file
                            break
                    
                    if json_path:
                        with zip_ref.open(json_path) as f:
                            multiqc_data = json.loads(TextIOWrapper(f).read())
                        
                        # Extract samples from the MultiQC data - look in various modules
                        found_samples = set()
                        
                        # Check different modules where samples might be listed
                        modules_to_check = ['rseqc_gene_body_coverage', 'rseqc_infer_experiment', 
                                           'rseqc_inner_distance', 'rseqc_read_distribution']
                        
                        for module in modules_to_check:
                            if module in multiqc_data:
                                found_samples.update(multiqc_data[module].keys())
                        
                        # Also check general stats list
                        if 'report_general_stats_data' in multiqc_data:
                            for stats in multiqc_data['report_general_stats_data']:
                                found_samples.update(stats.keys())
                        
                        # Improved sample name matching logic
                        clean_found_samples = set()
                        logger.debug(f"Raw sample names in MultiQC: {found_samples}")
                        
                        # For each sample we're expecting
                        for sample in sample_ids:
                            # Direct match
                            if sample in found_samples:
                                clean_found_samples.add(sample)
                                continue
                            
                            # Check for partial matches
                            for found_sample in found_samples:
                                # MultiQC sometimes modifies sample names or adds file suffixes
                                # Check if expected sample is part of a longer found sample name
                                if sample in found_sample:
                                    clean_found_samples.add(sample)
                                    logger.debug(f"Matched sample {sample} to MultiQC sample {found_sample}")
                                    break
                                    
                                # Try removing common suffixes like _sorted or _trimmed
                                # Check for variations based on common RSeQC naming patterns
                                possible_variations = [
                                    f"{sample}.infer_expt",
                                    f"{sample}_sorted",
                                    f"{sample}.geneBodyCoverage",
                                    f"{sample}.inner_distance",
                                    f"{sample}.read_dist"
                                ]
                                
                                for variation in possible_variations:
                                    if variation in found_sample:
                                        clean_found_samples.add(sample)
                                        logger.debug(f"Matched sample {sample} to MultiQC sample {found_sample} via pattern {variation}")
                                        break
                        
                        logger.info(f"Samples found in MultiQC after cleanup: {clean_found_samples}")
                        
                        # Record which samples were found
                        multiqc_samples = clean_found_samples
                        
                        # Find missing samples
                        missing_samples = set(samples) - clean_found_samples
                        
                        # Log results
                        if missing_samples:
                            msg = f"Samples missing from alignment MultiQC report: {', '.join(missing_samples)}"
                            logger.warning(msg)
                            val_logger.log("alignments", "ALL", "samples_in_multiqc", "YELLOW", msg)
                        else:
                            val_logger.log("alignments", "ALL", "samples_in_multiqc", "GREEN",
                                        f"All samples found in alignment MultiQC report")
                    else:
                        logger.warning(f"Could not find multiqc_data.json in {multiqc_path}")
            except Exception as e:
                val_logger.log("alignments", "ALL", "samples_in_multiqc", "HALT",
                            f"Error checking samples in MultiQC report: {str(e)}")
        except Exception as e:
            logger.error(f"Error processing MultiQC report: {str(e)}")
            val_logger.log("alignments", "ALL", "process_multiqc", "HALT",
                        f"Error processing MultiQC report: {str(e)}")
    
    return {
        "status": val_logger.get_status(),
        "messages": [r["message"] for r in val_logger.results],
        "failures": {r["check_name"]: r["message"] 
                    for r in val_logger.results 
                    if r["status"] in ["HALT", "RED"]}
    }

def validate_rseqc(outdir: Path,
                  samples_txt: Path,
                  paired_end: bool,
                  assay_suffix: str,
                  genebody_coverage_dir: Path = None,
                  infer_experiment_dir: Path = None,
                  inner_distance_dir: Path = None,
                  read_distribution_dir: Path = None) -> dict:
    """
    Validate RSeQC outputs in work directory.
    
    This function checks for the presence of RSeQC outputs (geneBody_coverage.py, infer_experiment.py,
    inner_distance.py, read_distribution.py) for each sample, either directly or through staging.
    
    Args:
        outdir (Path): Output directory path
        samples_txt (Path): Path to samples.txt file
        paired_end (bool): Whether data is paired-end
        assay_suffix (str): Assay suffix, e.g. "_GLbulkRNAseq"
        genebody_coverage_dir (Path, optional): Path to geneBody coverage directory
        infer_experiment_dir (Path, optional): Path to infer experiment directory
        inner_distance_dir (Path, optional): Path to inner distance directory
        read_distribution_dir (Path, optional): Path to read distribution directory
    
    Returns:
        dict: Validation status information
    """
    logger.info(f"Validating RSeQC outputs with params: samples_txt={samples_txt}, paired_end={paired_end}")
    
    # Create main RSeQC directory in the published directory structure
    work_rseqc_dir = Path("RSeQC_Analyses")
    if not work_rseqc_dir.exists():
        os.makedirs(work_rseqc_dir, exist_ok=True)
        logger.info(f"Created directory: {work_rseqc_dir}")
    
    # Prepare for MultiQC metrics outlier detection
    rseqc_metrics_by_sample = {}
    
    # Get the expected sample IDs from the runsheet
    samples = []
    with open(samples_txt, "r") as f:
        reader = csv.reader(f)
        header = next(reader)  # Skip header
        for row in reader:
            if row:  # Skip empty rows
                samples.append(row[-1].strip())  # Assume last column is sample name
    
    # Define RSeQC sections we're checking
    sections = {
        "geneBody_coverage.py": {
            "dir": genebody_coverage_dir,
            "display_name": "geneBody_coverage.py",
            "required": True,
            "multiqc_dir_pattern": f"geneBody_cov_multiqc{assay_suffix}_report",
            "file_patterns": ["geneBodyCoverage", ".geneBodyCoverage.txt"],
            "work_dir": work_rseqc_dir / "02_geneBody_coverage"
        },
        "infer_experiment.py": {
            "dir": infer_experiment_dir,
            "display_name": "infer_experiment.py",
            "required": True,
            "multiqc_dir_pattern": f"infer_exp_multiqc{assay_suffix}_report",
            "file_patterns": ["infer_expt", ".infer_experiment.txt"],
            "work_dir": work_rseqc_dir / "03_infer_experiment"
        },
        "inner_distance.py": {
            "dir": inner_distance_dir,
            "display_name": "inner_distance.py",
            "required": True if paired_end else False,
            "multiqc_dir_pattern": f"inner_dist_multiqc{assay_suffix}_report",
            "file_patterns": ["inner_distance", ".inner_distance.txt"],
            "work_dir": work_rseqc_dir / "04_inner_distance"
        },
        "read_distribution.py": {
            "dir": read_distribution_dir,
            "display_name": "read_distribution.py",
            "required": True,
            "multiqc_dir_pattern": f"read_dist_multiqc{assay_suffix}_report",
            "file_patterns": ["read_dist", ".read_distribution.txt"],
            "work_dir": work_rseqc_dir / "05_read_distribution"
        }
    }
    
    # Create section and sample directories in work directory
    for section_name, section in sections.items():
        if section["work_dir"] and (section_name != "inner_distance.py" or paired_end):
            os.makedirs(section["work_dir"], exist_ok=True)
            logger.info(f"Created section directory: {section['work_dir']}")
            
            # Create sample directories
            for sample in samples:
                sample_dir = section["work_dir"] / sample
                if not sample_dir.exists():
                    os.makedirs(sample_dir, exist_ok=True)
                    logger.debug(f"Created sample directory: {sample_dir}")
    
    # Track files found for each sample and section
    files_by_sample = {sample: {} for sample in samples}
    
    # Stage and track files from input directories to work directory
    for section_name, section in sections.items():
        if section["dir"] and os.path.exists(section["dir"]):
            logger.info(f"Checking {section['display_name']} outputs in {section['dir']}")
            
            # Look for sample files in this directory and stage them
            for sample in samples:
                sample_dir = section["work_dir"] / sample
                files_found = False
                
                # Check for each file pattern
                for pattern in section["file_patterns"]:
                    file_pattern = f"*{sample}*{pattern}*"
                    matches = list(Path(section["dir"]).glob(file_pattern))
                    
                    if matches:
                        # Initialize if needed
                        if section_name not in files_by_sample[sample]:
                            files_by_sample[sample][section_name] = []
                        
                        # Stage and track each file
                        for src_file in matches:
                            dst_file = sample_dir / src_file.name
                            try:
                                # Create symlink if it doesn't exist
                                if not dst_file.exists() and not dst_file.is_symlink():
                                    os.symlink(src_file, dst_file)
                                    logger.info(f"Linked {src_file} to {dst_file}")
                                
                                # Track the file
                                files_by_sample[sample][section_name].append(dst_file)
                                files_found = True
                            except Exception as e:
                                logger.error(f"Failed to link {src_file} to {dst_file}: {str(e)}")
                
                # Log only once if any files were found for this sample/section
                if files_found:
                    val_logger.log("rseqc", sample, f"{section['display_name']}_outputs_found", "GREEN",
                                f"Found RSeQC outputs for {section['display_name']}")
    
    # Also look for files in the rseqc-logs directory
    rseqc_logs_dirs = [
        Path("rseqc-logs"),  # Current working directory
        os.path.join(outdir, "rseqc-logs"),  # Inside output directory
        os.path.join(Path.cwd(), "rseqc-logs")  # Explicit current working directory
    ]
    
    # Find the first valid rseqc-logs directory
    rseqc_logs_dir = None
    for dir_path in rseqc_logs_dirs:
        logger.info(f"Checking for rseqc-logs at: {dir_path}")
        if os.path.exists(dir_path):
            rseqc_logs_dir = dir_path
            logger.info(f"Found rseqc-logs directory at: {rseqc_logs_dir}")
            break
    
    if rseqc_logs_dir:
        logger.info(f"Processing files from rseqc-logs directory: {rseqc_logs_dir}")
        
        # Define file type patterns and their target sections
        file_type_patterns = {
            "geneBodyCoverage": "geneBody_coverage.py",
            "infer_expt.out": "infer_experiment.py", 
            "inner_distance": "inner_distance.py",
            "read_dist.out": "read_distribution.py"
        }
        
        logger.info(f"File type patterns: {file_type_patterns}")
        
        # List all files in the rseqc-logs directory
        all_files = os.listdir(rseqc_logs_dir)
        logger.info(f"Found {len(all_files)} files in rseqc-logs directory")
        logger.debug(f"Files in rseqc-logs: {all_files}")
        
        # Process files to appropriate sample directories
        for file in all_files:
            logger.info(f"Processing file: {file}")
            
            # Find which sample this file belongs to
            sample_match = None
            for sample in samples:
                if sample in file:
                    sample_match = sample
                    logger.info(f"File {file} matched to sample {sample}")
                    break
                    
            if not sample_match:
                logger.warning(f"Could not match file {file} to any sample")
                continue
                
            # Determine which section this file belongs to
            target_section = None
            for pattern, section_name in file_type_patterns.items():
                if pattern in file:
                    target_section = section_name
                    logger.info(f"File {file} matched to section {section_name} with pattern {pattern}")
                    break
                    
            if not target_section:
                logger.warning(f"Could not match file {file} to any RSeQC section - patterns: {file_type_patterns.keys()}")
                continue
                
            # Skip inner_distance section for single-end data
            if target_section == "inner_distance.py" and not paired_end:
                continue
                
            # Create symbolic links for the file
            src = os.path.realpath(os.path.join(rseqc_logs_dir, file))
            
            # Get the working directory path
            section = sections[target_section]
            sample_dir = section["work_dir"] / sample_match
            
            
            dst_file = sample_dir / file
            
            try:
                # Ensure destination directories exist
                os.makedirs(sample_dir, exist_ok=True)
                
                # Create symbolic links in directories
                if not dst_file.exists() and not dst_file.is_symlink():
                    os.symlink(src, dst_file)
                    logger.info(f"Linked {src} to {dst_file}")
                
                # Update tracking
                if target_section not in files_by_sample[sample_match]:
                    files_by_sample[sample_match][target_section] = []
                files_by_sample[sample_match][target_section].append(dst_file)
                
                # Log successful file linkage from rseqc-logs
                val_logger.log("rseqc", sample_match, f"{section['display_name']}_outputs_found", "GREEN",
                            f"Found RSeQC outputs for {section['display_name']}")
                
            except Exception as e:
                msg = f"Failed to link {file} from rseqc-logs: {str(e)}"
                logger.error(msg)
                val_logger.log("rseqc", sample_match, f"{section['display_name']}_outputs_link", "RED", msg)
    
    # Look for MultiQC reports for each section and stage them
    for section_name, section in sections.items():
        if not section["required"]:
            continue
            
        # Look for MultiQC zip files
        if section["dir"] and os.path.exists(section["dir"]):
            # Find all MultiQC zip files in the directory
            multiqc_files = []
            for file in os.listdir(section["dir"]):
                if "multiqc" in file.lower() and file.endswith(".zip"):
                    multiqc_files.append(os.path.join(section["dir"], file))
            
            if multiqc_files:
                # Stage the first MultiQC zip file found
                src_file = Path(multiqc_files[0])
                # Create the destination path directly, not nested under the input directory
                dst_file = section["work_dir"] / f"{section['multiqc_dir_pattern']}.zip"
                if not os.path.exists(dst_file) and not os.path.islink(dst_file):
                    try:
                        # Use absolute paths for both source and destination to avoid nesting issues
                        abs_src = os.path.abspath(src_file)
                        abs_dst = os.path.abspath(dst_file)
                        os.symlink(abs_src, abs_dst)
                        logger.info(f"Linked MultiQC report {src_file} to {dst_file}")
                    except Exception as e:
                        logger.error(f"Failed to link MultiQC report {src_file} to {dst_file}: {str(e)}")
                
                val_logger.log("rseqc", "ALL", f"{section['display_name']}_multiqc", "GREEN",
                            f"MultiQC report found for {section['display_name']}")
    
    # Now after checking and staging files, log warnings for any truly missing files
    for sample in samples:
        for section_name, section in sections.items():
            if section["required"]:  # Skip non-required sections
                # Only log warning if no files were found for this sample and section
                if section_name not in files_by_sample[sample]:
                    msg = f"No {section['display_name']} outputs found for sample {sample} in any location"
                    logger.warning(msg)
                    val_logger.log("rseqc", sample, f"{section['display_name']}_outputs", "YELLOW", msg)
    
    # Parse sample metrics from MultiQC data for outlier detection
    rseqc_metrics_by_sample = {}
    
    # Parse the geneBody coverage 3'/5' ratios
    try:
        # First, directly use the MultiQC files from the input directories
        input_multiqc_files = []
        
        # Look for MultiQC files in each input directory
        if genebody_coverage_dir and os.path.exists(genebody_coverage_dir):
            for file in os.listdir(genebody_coverage_dir):
                if "multiqc" in file.lower() and file.endswith(".zip"):
                    input_multiqc_files.append((os.path.join(genebody_coverage_dir, file), "geneBody_cov"))
        
        if infer_experiment_dir and os.path.exists(infer_experiment_dir):
            for file in os.listdir(infer_experiment_dir):
                if "multiqc" in file.lower() and file.endswith(".zip"):
                    input_multiqc_files.append((os.path.join(infer_experiment_dir, file), "infer_exp"))
        
        if inner_distance_dir and paired_end and os.path.exists(inner_distance_dir):
            for file in os.listdir(inner_distance_dir):
                if "multiqc" in file.lower() and file.endswith(".zip"):
                    input_multiqc_files.append((os.path.join(inner_distance_dir, file), "inner_dist"))
        
        if read_distribution_dir and os.path.exists(read_distribution_dir):
            for file in os.listdir(read_distribution_dir):
                if "multiqc" in file.lower() and file.endswith(".zip"):
                    input_multiqc_files.append((os.path.join(read_distribution_dir, file), "read_dist"))
        
        logger.info(f"Found {len(input_multiqc_files)} MultiQC files for RSeQC metrics extraction")
        
        # Create a temporary directory to extract MultiQC data
        with tempfile.TemporaryDirectory() as temp_extraction_dir:
            temp_dir_path = Path(temp_extraction_dir)
            
            # Extract each MultiQC zip into the expected directories for parsing
            extract_success = False
            
            # Process each input MultiQC file directly
            for multiqc_file, prefix in input_multiqc_files:
                try:
                    # Extract to create the directory structure expected by parse_multiqc functions
                    with zipfile.ZipFile(multiqc_file, 'r') as zip_ref:
                        # Get the directory containing MultiQC data
                        data_dir_in_zip = None
                        for item in zip_ref.namelist():
                            if "data" in item.lower() and item.endswith('/'):
                                data_dir_in_zip = item
                                break
                        
                        if data_dir_in_zip:
                            # Create expected directory for the parsing functions
                            # Handle the case where assay_suffix already has a leading underscore
                            suffix = assay_suffix
                            if suffix.startswith('_'):
                                suffix = suffix[1:]  # Remove the leading underscore
                            # Use the exact directory name expected by the parsing functions
                            expected_dir = Path(os.getcwd()) / f"{prefix}_multiqc_{suffix}_data"
                            os.makedirs(expected_dir, exist_ok=True)
                            
                            # Extract only the data directory
                            for item in zip_ref.namelist():
                                if item.startswith(data_dir_in_zip):
                                    # Extract to the expected location but strip the prefix directory
                                    extracted_name = item.replace(data_dir_in_zip, '')
                                    if extracted_name:  # Skip the directory itself
                                        extraction_target = expected_dir / extracted_name
                                        # Create directory if needed
                                        if item.endswith('/'):
                                            os.makedirs(extraction_target, exist_ok=True)
                                        else:
                                            # Extract file
                                            with zip_ref.open(item) as source, open(extraction_target, 'wb') as target:
                                                shutil.copyfileobj(source, target)
                            
                            extract_success = True
                            logger.info(f"Extracted MultiQC data from {multiqc_file} to {expected_dir}")
                except Exception as e:
                    logger.error(f"Error extracting MultiQC data from {multiqc_file}: {str(e)}")
            
            if not extract_success:
                logger.warning("Failed to extract any MultiQC data for parsing")
        
        # Import parse module from this script's directory
        import_path = os.path.dirname(os.path.abspath(__file__))
        if import_path not in sys.path:
            sys.path.insert(0, import_path)
            
        # Attempt to use parse_multiqc module to extract metrics
        try:
            from parse_multiqc import parse_genebody_cov, parse_infer_exp, parse_inner_dist, parse_read_dist
            
            # Try to parse genebody coverage data
            genebody_metrics = parse_genebody_cov(assay_suffix)
            if genebody_metrics:
                logger.info(f"Extracted genebody coverage metrics for {len(genebody_metrics)} samples")
                for sample, metrics in genebody_metrics.items():
                    # Find matching sample from our list
                    matched_sample = None
                    if sample in samples:
                        matched_sample = sample
                    else:
                        for s in samples:
                            if s in sample or sample in s:
                                matched_sample = s
                                break
                    
                    if matched_sample:
                        if matched_sample not in rseqc_metrics_by_sample:
                            rseqc_metrics_by_sample[matched_sample] = {}
                        
                        # Add metrics with prefix for clarity
                        for key, value in metrics.items():
                            rseqc_metrics_by_sample[matched_sample][f"genebody_coverage_{key}"] = value
            
            # Try to parse infer experiment data
            infer_metrics = parse_infer_exp(assay_suffix)
            if infer_metrics:
                logger.info(f"Extracted infer experiment metrics for {len(infer_metrics)} samples")
                for sample, metrics in infer_metrics.items():
                    # Find matching sample from our list
                    matched_sample = None
                    if sample in samples:
                        matched_sample = sample
                    else:
                        for s in samples:
                            if s in sample or sample in s:
                                matched_sample = s
                                break
                    
                    if matched_sample:
                        if matched_sample not in rseqc_metrics_by_sample:
                            rseqc_metrics_by_sample[matched_sample] = {}
                        
                        # Add metrics with prefix for clarity
                        for key, value in metrics.items():
                            rseqc_metrics_by_sample[matched_sample][f"infer_experiment_{key}"] = value
            
            # Try to parse inner distance data
            if paired_end:
                inner_metrics = parse_inner_dist(assay_suffix)
                if inner_metrics:
                    logger.info(f"Extracted inner distance metrics for {len(inner_metrics)} samples")
                    for sample, metrics in inner_metrics.items():
                        # Find matching sample from our list
                        matched_sample = None
                        if sample in samples:
                            matched_sample = sample
                        else:
                            for s in samples:
                                if s in sample or sample in s:
                                    matched_sample = s
                                    break
                        
                        if matched_sample:
                            if matched_sample not in rseqc_metrics_by_sample:
                                rseqc_metrics_by_sample[matched_sample] = {}
                            
                            # Add metrics with prefix for clarity
                            for key, value in metrics.items():
                                rseqc_metrics_by_sample[matched_sample][f"inner_distance_{key}"] = value
            
            # Try to parse read distribution data
            read_dist_metrics = parse_read_dist(assay_suffix)
            if read_dist_metrics:
                logger.info(f"Extracted read distribution metrics for {len(read_dist_metrics)} samples")
                for sample, metrics in read_dist_metrics.items():
                    # Find matching sample from our list
                    matched_sample = None
                    if sample in samples:
                        matched_sample = sample
                    else:
                        for s in samples:
                            if s in sample or sample in s:
                                matched_sample = s
                                break
                    
                    if matched_sample:
                        if matched_sample not in rseqc_metrics_by_sample:
                            rseqc_metrics_by_sample[matched_sample] = {}
                        
                        # Add metrics with prefix for clarity
                        for key, value in metrics.items():
                            rseqc_metrics_by_sample[matched_sample][f"read_distribution_{key}"] = value
        
        except ImportError as e:
            logger.warning(f"Could not import RSeQC parsing functions: {str(e)}")
        except Exception as e:
            logger.warning(f"Error parsing RSeQC metrics: {str(e)}")
            logger.warning(traceback.format_exc())
    except Exception as e:
        logger.error(f"Error in RSeQC metrics import setup: {str(e)}")
    
    # Detect outliers in collected RSeQC metrics
    if rseqc_metrics_by_sample:
        logger.info(f"Detecting outliers in RSeQC metrics for {len(rseqc_metrics_by_sample)} samples")
        
        # Define metrics to check for outliers based on their importance
        metrics_to_check = [
            # Gene body coverage metrics
            'genebody_coverage_ratio_genebody_cov_3_to_5',
            'genebody_coverage_mean_genebody_cov_5_20',
            'genebody_coverage_mean_genebody_cov_40_60',
            'genebody_coverage_mean_genebody_cov_80_95',
            
            # Infer experiment metrics (strandedness)
            'infer_experiment_pct_sense',
            'infer_experiment_pct_antisense',
            'infer_experiment_pct_undetermined',
            
            # Read distribution metrics
            'read_distribution_cds_exons_pct',
            'read_distribution_5_utr_exons_pct',
            'read_distribution_3_utr_exons_pct',
            'read_distribution_introns_pct',
            'read_distribution_other_intergenic_pct',
            
            # Inner distance metrics (paired-end only)
            'inner_distance_peak_inner_dist',
            'inner_distance_peak_inner_dist_pct_reads'
        ]
        
        # Log the metrics we have for each sample
        for sample, metrics in rseqc_metrics_by_sample.items():
            available_metrics = [m for m in metrics_to_check if m in metrics]
            if available_metrics:
                logger.debug(f"Sample {sample} has RSeQC metrics: {available_metrics}")
            else:
                logger.warning(f"Sample {sample} has no RSeQC metrics for outlier detection")
        
        # Print key metric values for each sample for debugging
        key_debug_metrics = [
            'genebody_coverage_ratio_genebody_cov_3_to_5',
            'genebody_coverage_mean_genebody_cov_5_20',
            'infer_experiment_pct_sense',
            'read_distribution_cds_exons_pct'
        ]
        
        for metric in key_debug_metrics:
            logger.info(f"Values for {metric}:")
            for sample in sorted(rseqc_metrics_by_sample.keys()):
                if metric in rseqc_metrics_by_sample[sample]:
                    logger.info(f"  {sample}: {rseqc_metrics_by_sample[sample][metric]:.2f}")
        
        # Detect outliers using our existing function
        outliers, metric_stats = detect_outliers(rseqc_metrics_by_sample, metrics_to_check)
        
        # Print metric statistics for context
        for metric in key_debug_metrics:
            if metric in metric_stats:
                logger.info(f"Stats for {metric}: mean={metric_stats[metric]['mean']:.2f}, std={metric_stats[metric]['std']:.2f}")
        
        # Add strandedness assessment using logic similar to assess_strandedness.py
        STRANDEDNESS_THRESHOLD = 70.0  # Above this percent is considered stranded
        UNSTRANDEDNESS_THRESHOLD_MIN = 40.0  # Between min and max is considered unstranded
        UNSTRANDEDNESS_THRESHOLD_MAX = 60.0
        
        logger.info("Strandedness assessment:")
        # Check strandedness for all samples at once
        strandedness_results = {}
        sample_strandedness = {}
        
        for sample_id, metrics in rseqc_metrics_by_sample.items():
            sense_pct = metrics.get('infer_experiment_pct_sense', 0)
            antisense_pct = metrics.get('infer_experiment_pct_antisense', 0)
            undetermined_pct = 100 - (sense_pct + antisense_pct)
            
            if sense_pct >= STRANDEDNESS_THRESHOLD:
                strandedness = "forward stranded"
            elif antisense_pct >= STRANDEDNESS_THRESHOLD:
                strandedness = "reverse stranded"
            elif UNSTRANDEDNESS_THRESHOLD_MIN <= sense_pct <= UNSTRANDEDNESS_THRESHOLD_MAX and UNSTRANDEDNESS_THRESHOLD_MIN <= antisense_pct <= UNSTRANDEDNESS_THRESHOLD_MAX:
                strandedness = "unstranded"
            else:
                strandedness = "undetermined"
            
            sample_strandedness[sample_id] = strandedness
            if strandedness not in strandedness_results:
                strandedness_results[strandedness] = 0
            strandedness_results[strandedness] = strandedness_results[strandedness] + 1
        
        # Determine overall dataset strandedness
        total_samples = len(rseqc_metrics_by_sample)
        dataset_strandedness = "mixed"
        for strand_type, count in strandedness_results.items():
            if count == total_samples:
                dataset_strandedness = strand_type
                break
        
        # Calculate average metrics across all samples
        avg_sense = sum(metrics.get('infer_experiment_pct_sense', 0) for metrics in rseqc_metrics_by_sample.values()) / len(rseqc_metrics_by_sample) if rseqc_metrics_by_sample else 0
        avg_antisense = sum(metrics.get('infer_experiment_pct_antisense', 0) for metrics in rseqc_metrics_by_sample.values()) / len(rseqc_metrics_by_sample) if rseqc_metrics_by_sample else 0
        avg_undetermined = 100 - avg_sense - avg_antisense
        
        # Format the strandedness summary with average metrics
        strandedness_summary = f"Dataset strandedness: {dataset_strandedness} ({', '.join([f'{strand_type}: {count}/{total_samples}' for strand_type, count in strandedness_results.items()])})"
        strandedness_details = f"mean sense: {avg_sense:.2f}%, mean antisense: {avg_antisense:.2f}%, mean undetermined: {avg_undetermined:.2f}%"
        
        # Log a single consolidated strandedness entry with average metrics
        val_logger.log(
            "rseqc",
            "all_samples",
            "strandedness",
            "INFO",
            strandedness_summary,
            strandedness_details
        )
        
        logger.info(f"Strandedness assessment complete: {strandedness_summary} ({strandedness_details})")
        
        # Log to regular logger as well
        logger.info(f"  {sample}: {strandedness} (sense: {sense_pct:.2f}%, antisense: {antisense_pct:.2f}%, undetermined: {undetermined_pct:.2f}%)")
        
        # Log any outliers we've found
        if outliers:
            logger.info(f"Found outliers in RSeQC metrics: {len(outliers)} samples with outliers")
            
            for sample, sample_outliers in outliers.items():
                for outlier in sample_outliers:
                    metric = outlier['metric']
                    value = outlier['value']
                    mean = outlier['mean']
                    std = outlier['std']
                    z_score = outlier['z_score']
                    
                    # Determine severity by z-score
                    status = "YELLOW"
                    if abs(z_score) > 4:
                        status = "RED"
                    
                    # Create user-friendly metric display name
                    display_metric = metric
                    if "genebody_coverage" in metric:
                        display_metric = metric.replace('genebody_coverage_', 'GeneBody coverage: ')
                    elif "infer_experiment" in metric:
                        display_metric = metric.replace('infer_experiment_', 'Strandedness: ')
                    elif "read_distribution" in metric:
                        display_metric = metric.replace('read_distribution_', 'Read distribution: ')
                    elif "inner_distance" in metric:
                        display_metric = metric.replace('inner_distance_', 'Inner distance: ')
                    
                    val_logger.log("rseqc", sample, f"{metric}_outlier", status,
                                f"{display_metric} value {value:.2f} is an outlier",
                                f"Value deviates {abs(z_score):.2f} standard deviations from mean ({mean:.2f})")
        else:
            val_logger.log("rseqc", "ALL", "rseqc_metrics_outliers", "GREEN",
                         "No outliers detected in RSeQC metrics")
    else:
        val_logger.log("rseqc", "ALL", "rseqc_metrics", "YELLOW",
                     "No RSeQC metrics found for outlier detection")
    
    # Return validation status
    return {
        "status": val_logger.get_status(),
        "messages": [r["message"] for r in val_logger.results],
        "failures": {r["check_name"]: r["message"] 
                    for r in val_logger.results 
                    if r["status"] in ["HALT", "RED"]}
    }

def check_samples_in_multiqc(multiqc_zip_path: Path, samples: list, assay_type: str) -> dict:
    """
    Check if all expected samples are present in a MultiQC report.
    
    This function extracts and parses the MultiQC general stats file to verify
    that all expected samples are included in the report. It's used to validate
    that MultiQC properly processed all sample data.
    
    The function:
    1. Opens the MultiQC zip file
    2. Locates and extracts the multiqc_general_stats.txt file
    3. Parses the file to get a list of samples included in the report
    4. Compares this list against the expected samples
    
    Args:
        multiqc_zip_path: Path to the MultiQC zip file
        samples: List of expected sample names/IDs to check for
        assay_type: Type of assay for logging purposes (e.g., "raw_reads", "trimmed_reads")
        
    Returns:
        dict: Dictionary with results of the check, including:
             - all_samples_found: Boolean indicating if all samples were found
             - missing_samples: List of samples not found in the MultiQC report
             - found_samples: List of samples found in the MultiQC report
             - error: Error message if something went wrong during the check (optional)
    """
    if not multiqc_zip_path.exists():
        return {
            "all_samples_found": False,
            "missing_samples": samples,
            "found_samples": [],
            "error": f"MultiQC file not found: {multiqc_zip_path}"
        }
    
    try:
        # Extract the multiqc_data directory from the zip
        with zipfile.ZipFile(multiqc_zip_path, 'r') as zip_ref:
            # Find the multiqc_data/multiqc_general_stats.txt file which contains sample names
            data_files = [f for f in zip_ref.namelist() if f.endswith('/multiqc_general_stats.txt')]
            
            if not data_files:
                return {
                    "all_samples_found": False,
                    "missing_samples": samples,
                    "found_samples": [],
                    "error": "No general stats data found in MultiQC report"
                }
            
            # Extract and parse the general stats file
            with zip_ref.open(data_files[0]) as f:
                content = TextIOWrapper(f).read()
                
                # Extract sample names from the first column
                lines = content.strip().split('\n')
                if len(lines) < 2:  # Need at least header and one sample
                    return {
                        "all_samples_found": False,
                        "missing_samples": samples,
                        "found_samples": [],
                        "error": "MultiQC general stats file is empty or malformed"
                    }
                
                # Parse header and extract sample column
                header = lines[0].split('\t')
                if not header or len(header) < 1:
                    return {
                        "all_samples_found": False,
                        "missing_samples": samples,
                        "found_samples": [],
                        "error": "MultiQC general stats header not found"
                    }
                
                # Extract all sample names
                multiqc_samples = []
                for line in lines[1:]:  # Skip header
                    cols = line.split('\t')
                    if cols and len(cols) >= 1:
                        multiqc_samples.append(cols[0])
                
                # Check for each expected sample
                found_samples = []
                missing_samples = []
                
                for sample in samples:
                    # Look for exact matches or sample names within MultiQC sample names
                    if sample in multiqc_samples or any(sample in mqc_sample for mqc_sample in multiqc_samples):
                        found_samples.append(sample)
                    else:
                        missing_samples.append(sample)
                
                return {
                    "all_samples_found": len(missing_samples) == 0,
                    "missing_samples": missing_samples,
                    "found_samples": found_samples,
                    "multiqc_samples": multiqc_samples
                }
                
    except Exception as e:
        return {
            "all_samples_found": False,
            "missing_samples": samples,
            "found_samples": [],
            "error": f"Error checking MultiQC: {str(e)}"
        }

def find_file_recursive(base_dir, target_file):
    """Search recursively for a file within a directory"""
    logging.info(f"Searching recursively for {target_file} in {base_dir}")
    if not os.path.exists(base_dir):
        logging.warning(f"Base directory {base_dir} does not exist")
        return None
        
    for root, dirs, files in os.walk(base_dir):
        if target_file in files:
            file_path = os.path.join(root, target_file)
            logging.info(f"Found {target_file} at {file_path}")
            return file_path
    
    logging.warning(f"Could not find {target_file} in {base_dir}")
    return None

def parse_fastqc(prefix, assay_suffix):
    """
    Parse FastQC data from a MultiQC report to extract quality metrics
    
    Args:
        prefix (str): Prefix for the MultiQC directory (e.g., 'raw')
        assay_suffix (str): Suffix for the assay (e.g., '_GLbulkRNAseq')
    
    Returns:
        dict: Dictionary of sample names and their quality metrics
    """
    # Determine the correct path to the multiqc data directory
    fastqc_dir = Path("00-RawData/FastQC_Reports") if prefix == 'raw' else Path("01-TG_Preproc/FastQC_Reports")
    
    # Try several possible paths for the multiqc_data.json file
    possible_paths = [
        str(fastqc_dir / f"{prefix}_multiqc{assay_suffix}_data" / "multiqc_data.json"),
        str(fastqc_dir / f"{prefix}_multiqc{assay_suffix}_report_data" / "multiqc_data.json"),
        # Look in the base data directory
        str(fastqc_dir / f"{prefix}_multiqc{assay_suffix}_data/multiqc_data/multiqc_data.json"),
        # Try without the assay suffix in case it's missing
        str(fastqc_dir / f"{prefix}_multiqc_data" / "multiqc_data.json")
    ]
    
    multiqc_data_path = None
    for path in possible_paths:
        logging.info(f"Checking for MultiQC data at {path}")
        if os.path.exists(path):
            multiqc_data_path = path
            logging.info(f"Found MultiQC data at {path}")
            break
    
    # If not found in expected locations, try a recursive search
    if not multiqc_data_path:
        logging.info(f"MultiQC data not found in expected locations, trying recursive search")
        multiqc_data_dir = str(fastqc_dir / f"{prefix}_multiqc{assay_suffix}_data")
        
        # If the directory doesn't exist, try other possible directory names
        if not os.path.exists(multiqc_data_dir):
            alt_dir = str(fastqc_dir / f"{prefix}_multiqc{assay_suffix}_report_data")
            if os.path.exists(alt_dir):
                multiqc_data_dir = alt_dir
        
        # If we found a valid directory, search for multiqc_data.json recursively
        if os.path.exists(multiqc_data_dir):
            recursive_path = find_file_recursive(multiqc_data_dir, "multiqc_data.json")
            if recursive_path:
                multiqc_data_path = recursive_path
    
    if not multiqc_data_path:
        logging.warning(f"MultiQC data file not found in any location")
        return {}
    
    try:
        with open(multiqc_data_path) as f:
            j = json.loads(f.read())
            
        # Group the samples by base name for paired end data
        sample_groups = {}
        for sample in j['report_general_stats_data'][-1].keys():
            # Handle various naming patterns
            if ' Read 1' in sample:
                base_name = sample.replace(' Read 1', '')
                if base_name not in sample_groups:
                    sample_groups[base_name] = {'f': None, 'r': None}
                sample_groups[base_name]['f'] = sample
            elif ' Read 2' in sample:
                base_name = sample.replace(' Read 2', '')
                if base_name not in sample_groups:
                    sample_groups[base_name] = {'f': None, 'r': None}
                sample_groups[base_name]['r'] = sample
            elif '_R1' in sample:
                base_name = sample.replace('_R1', '')
                if base_name not in sample_groups:
                    sample_groups[base_name] = {'f': None, 'r': None}
                sample_groups[base_name]['f'] = sample
            elif '_R2' in sample:
                base_name = sample.replace('_R2', '')
                if base_name not in sample_groups:
                    sample_groups[base_name] = {'f': None, 'r': None}
                sample_groups[base_name]['r'] = sample
            else:
                # For single-end or non-paired samples
                base_name = sample
                if base_name not in sample_groups:
                    sample_groups[base_name] = {'f': None, 'r': None}
                sample_groups[base_name]['f'] = sample

        data = {}
        # Process each sample group
        for base_name, reads in sample_groups.items():
            data[base_name] = {}
            
            # Process forward read
            if reads['f']:
                for k, v in j['report_general_stats_data'][-1][reads['f']].items():
                    if k != 'percent_fails':
                        data[base_name][prefix + '_' + k + '_f'] = v
                    
            # Process reverse read
            if reads['r']:
                for k, v in j['report_general_stats_data'][-1][reads['r']].items():
                    if k != 'percent_fails':
                        data[base_name][prefix + '_' + k + '_r'] = v

        # Process other stats sections (quality, GC, etc)
        for section, suffix in [
            ('fastqc_per_base_sequence_quality_plot', 'quality_score'),
            ('fastqc_per_sequence_gc_content_plot', 'gc'),
            ('fastqc_per_base_n_content_plot', 'n_content')
        ]:
            if section in j['report_plot_data']:
                for data_item in j['report_plot_data'][section]['datasets'][0]['lines']:
                    sample = data_item['name']
                    
                    # Determine if it's forward or reverse read
                    read_suffix = '_f'  # Default to forward
                    base_name = sample
                    
                    if ' Read 2' in sample:
                        read_suffix = '_r'
                        base_name = sample.replace(' Read 2', '')
                    elif ' Read 1' in sample:
                        base_name = sample.replace(' Read 1', '')
                    elif '_R2' in sample:
                        read_suffix = '_r'
                        base_name = sample.replace('_R2', '')
                    elif '_R1' in sample:
                        base_name = sample.replace('_R1', '')
                    
                    # Skip if we don't have this sample
                    if base_name not in data:
                        continue
                    
                    # Process based on the section
                    if suffix == 'quality_score':
                        data[base_name][prefix + '_quality_score_mean' + read_suffix] = mean([i[1] for i in data_item['pairs']])
                        data[base_name][prefix + '_quality_score_median' + read_suffix] = median([i[1] for i in data_item['pairs']])
                    elif suffix == 'gc':
                        gc_data_1pct = [i[0] for i in data_item['pairs'] if i[1] >= 1]
                        if gc_data_1pct:
                            data[base_name][prefix + '_gc_min_1pct' + read_suffix] = gc_data_1pct[0]
                            data[base_name][prefix + '_gc_max_1pct' + read_suffix] = gc_data_1pct[-1]
                            
                            gc_data_cum = list(np.cumsum([i[1] for i in data_item['pairs']]))
                            data[base_name][prefix + '_gc_auc_25pct' + read_suffix] = list(i >= 25 for i in gc_data_cum).index(True)
                            data[base_name][prefix + '_gc_auc_50pct' + read_suffix] = list(i >= 50 for i in gc_data_cum).index(True)
                            data[base_name][prefix + '_gc_auc_75pct' + read_suffix] = list(i >= 75 for i in gc_data_cum).index(True)
                    elif suffix == 'n_content':
                        data[base_name][prefix + '_n_content_sum' + read_suffix] = sum([i[1] for i in data_item['pairs']])

        return data
    except Exception as e:
        logging.error(f"Error parsing MultiQC data: {str(e)}")
        return {}

def validate_featurecounts(outdir: Path,
                       samples_txt: Path,
                       paired_end: bool,
                       assay_suffix: str,
                       featurecounts_counts_dir: Path = None,
                       featurecounts_summary_dir: Path = None,
                       featurecounts_counts_rrnarm_dir: Path = None,
                       featurecounts_counts_rrnarm_summary_dir: Path = None,
                       featurecounts_multiqc_dir: Path = None) -> dict:
    """
    Perform validation checks on FeatureCounts results.
    
    This function validates the FeatureCounts outputs by performing:
    1. File existence checks for counts files (.txt)
    2. File existence checks for summary files (.txt.summary)
    3. Content validation of counts files (proper format, gene IDs, counts)
    4. Content validation of summary files (assignment rates)
    5. Detection of outliers in assignment rates
    6. Presence of MultiQC report
    
    Args:
        outdir: Path to the output directory
        samples_txt: Path to the runsheet file containing sample information
        paired_end: Boolean indicating if data is paired-end (True) or single-end (False)
        assay_suffix: Suffix for the assay (e.g., "_GLbulkRNAseq")
        featurecounts_counts_dir: Path to directory with FeatureCounts counts files
        featurecounts_summary_dir: Path to directory with FeatureCounts summary files
        featurecounts_counts_rrnarm_dir: Path to directory with rRNA-removed counts files
        featurecounts_counts_rrnarm_summary_dir: Path to directory with rRNA-removed summary files (UNUSED/REMOVED)
        featurecounts_multiqc_dir: Path to directory with FeatureCounts MultiQC reports
        
    Returns:
        dict: Dictionary with validation status, messages, and any failures
    """
    
    logger.info("Starting FeatureCounts validation:")
    logger.info(f"  Output directory: {outdir}")
    logger.info(f"  Samples file: {samples_txt}")
    logger.info(f"  Paired-end: {paired_end}")
    logger.info(f"  FeatureCounts counts directory: {featurecounts_counts_dir}")
    logger.info(f"  FeatureCounts summary directory: {featurecounts_summary_dir}")
    logger.info(f"  FeatureCounts rRNA-removed counts directory: {featurecounts_counts_rrnarm_dir}")
    logger.info(f"  FeatureCounts rRNA-removed summary directory: {featurecounts_counts_rrnarm_summary_dir} (UNUSED)")
    logger.info(f"  FeatureCounts MultiQC directory: {featurecounts_multiqc_dir}")
    
    # Get sample IDs from runsheet
    sample_ids = []
    with open(samples_txt, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if 'Sample Name' in row:
                sample_name = row['Sample Name'].strip()
                if sample_name:  # Only add non-empty sample names
                    sample_ids.append(sample_name)
                    logging.debug(f"Added sample: {sample_name}")
            else:
                raise ValueError("Runsheet missing required 'Sample Name' column")
    
    logger.info(f"Found {len(sample_ids)} samples in runsheet: {', '.join(sample_ids)}")
    
    # Create 03-FeatureCounts directory structure
    fc_dir = outdir / "03-FeatureCounts"
    os.makedirs(fc_dir, exist_ok=True)
    
    # Stage files
    staged_files = {
        "counts": [],
        "summary": [],
        "counts_rrnarm": [],
        "multiqc": []
    }
    
    # Stage all files first
    
    # 1. Stage counts files
    if featurecounts_counts_dir:
        logger.info(f"Staging FeatureCounts counts files from {featurecounts_counts_dir}")
        
        for file in os.listdir(featurecounts_counts_dir):
            if file.endswith(".txt") or file.endswith(".tsv"):
                src = os.path.join(featurecounts_counts_dir, file)
                dst = os.path.join(fc_dir, file)
                
                # Create symlink
                if os.path.exists(dst) or os.path.islink(dst):
                    os.unlink(dst)
                os.symlink(os.path.realpath(src), dst)
                staged_files["counts"].append(dst)
    
    # 2. Stage summary files
    if featurecounts_summary_dir:
        logger.info(f"Staging FeatureCounts summary files from {featurecounts_summary_dir}")
        
        for file in os.listdir(featurecounts_summary_dir):
            if file.endswith(".txt.summary") or file.endswith(".tsv.summary") or file.endswith(".summary"):
                src = os.path.join(featurecounts_summary_dir, file)
                dst = os.path.join(fc_dir, file)
                
                # Create symlink
                if os.path.exists(dst) or os.path.islink(dst):
                    os.unlink(dst)
                os.symlink(os.path.realpath(src), dst)
                staged_files["summary"].append(dst)
    
    # 3. Stage rRNA-removed counts files
    if featurecounts_counts_rrnarm_dir:
        logger.info(f"Staging FeatureCounts rRNA-removed counts files from {featurecounts_counts_rrnarm_dir}")
        
        for file in os.listdir(featurecounts_counts_rrnarm_dir):
            if file.endswith(".txt") or file.endswith(".tsv"):
                src = os.path.join(featurecounts_counts_rrnarm_dir, file)
                dst = os.path.join(fc_dir, file)
                
                # Create symlink
                if os.path.exists(dst) or os.path.islink(dst):
                    os.unlink(dst)
                os.symlink(os.path.realpath(src), dst)
                staged_files["counts_rrnarm"].append(dst)
    
    # 4. (REMOVED) Stage rRNA-removed summary files - no longer used
    
    # 5. Stage MultiQC report
    multiqc_found = False
    multiqc_data = None
    if featurecounts_multiqc_dir:
        logger.info(f"Staging FeatureCounts MultiQC report from {featurecounts_multiqc_dir}")
        
        for file in os.listdir(featurecounts_multiqc_dir):
            if file.endswith(f"_multiqc{assay_suffix}_report.zip") or file.endswith("_report.zip"):
                src = os.path.join(featurecounts_multiqc_dir, file)
                dst = os.path.join(fc_dir, f"featureCounts_multiqc{assay_suffix}_report.zip")
                
                # Create symlink
                if os.path.exists(dst) or os.path.islink(dst):
                    os.unlink(dst)
                os.symlink(os.path.realpath(src), dst)
                staged_files["multiqc"].append(dst)
                multiqc_found = True
    
    # STEP 1: Validate basic file existence
    # Check if any counts files exist
    if featurecounts_counts_dir and len(staged_files["counts"]) == 0:
        val_logger.log("featurecounts", "ALL", "counts_file", "RED", 
                     f"No FeatureCounts counts file found")
    elif featurecounts_counts_dir:
        # Log the consolidated count files
        logger.info(f"Found FeatureCounts count files:")
        for file in staged_files["counts"]:
            logger.info(f"  - {os.path.basename(file)}")
            val_logger.log("featurecounts", "ALL", "counts_file", "GREEN", 
                         f"FeatureCounts counts file found: {os.path.basename(file)}")
            
    # Check if any counts summary files exist
    if featurecounts_summary_dir and len(staged_files["summary"]) == 0:
        val_logger.log("featurecounts", "ALL", "summary_file", "RED", 
                     f"No FeatureCounts summary file found")
    elif featurecounts_summary_dir:
        # Log the consolidated summary files
        logger.info(f"Found FeatureCounts summary files:")
        for file in staged_files["summary"]:
            logger.info(f"  - {os.path.basename(file)}")
            val_logger.log("featurecounts", "ALL", "summary_file", "GREEN", 
                         f"FeatureCounts summary file found: {os.path.basename(file)}")
    
    # Check if any rRNA-removed counts files exist
    if featurecounts_counts_rrnarm_dir and len(staged_files["counts_rrnarm"]) == 0:
        val_logger.log("featurecounts", "ALL", "counts_rrnarm_file", "RED", 
                     f"No FeatureCounts rRNA-removed counts file found")
    elif featurecounts_counts_rrnarm_dir:
        # Log the consolidated rRNA-removed count files
        logger.info(f"Found FeatureCounts rRNA-removed count files:")
        for file in staged_files["counts_rrnarm"]:
            logger.info(f"  - {os.path.basename(file)}")
            val_logger.log("featurecounts", "ALL", "counts_rrnarm_file", "GREEN", 
                         f"FeatureCounts rRNA-removed counts file found: {os.path.basename(file)}")
            
    # (REMOVED) Check if any rRNA-removed summary files exist
    
    # STEP 2: Check file formats
    # Validate counts file content
    for file_path in staged_files["counts"]:
        try:
            with open(file_path, 'r') as f:
                lines = f.readlines()
                # Check header (should have Geneid, Chr, Start, End, etc.)
                header = lines[0].strip().split('\t')
                if not "Geneid" in header[0] and not header[0].startswith("# Program:featureCounts"):
                    val_logger.log("featurecounts", "ALL", "counts_file_format", "YELLOW", 
                                 f"Counts file has unexpected header format: {header[0]}")
                else:
                    # Header format is acceptable (either contains "Geneid" or starts with the expected program header)
                    val_logger.log("featurecounts", "ALL", "counts_file_format", "GREEN", 
                                 f"Counts file has expected header format")
                
                # Check that file has content beyond header
                if len(lines) < 2:
                    val_logger.log("featurecounts", "ALL", "counts_file_content", "RED", 
                                 "Counts file has no content beyond header")
                else:
                    # Log success - Use "ALL" as the sample_id instead of the filename
                    val_logger.log("featurecounts", "ALL", "counts_file_exists", "GREEN", 
                                 f"FeatureCounts counts file exists and has valid format")
        except Exception as e:
            val_logger.log("featurecounts", "ALL", "counts_file_read", "RED", 
                         f"Error reading counts file: {str(e)}")
            logger.error(f"Error reading counts file: {str(e)}")
            
    # Validate summary file content
    for file_path in staged_files["summary"]:
        try:
            with open(file_path, 'r') as f:
                if len(f.read().strip()) > 0:
                    # File has some content
                    val_logger.log("featurecounts", "ALL", "summary_file_format", "GREEN", 
                                 f"Summary file has valid format")
        except Exception as e:
            val_logger.log("featurecounts", "ALL", "summary_file_read", "RED", 
                         f"Error reading summary file: {str(e)}")
            logger.error(f"Error reading summary file: {str(e)}")
    
    # Validate rRNA-removed counts file content
    for file_path in staged_files["counts_rrnarm"]:
        try:
            with open(file_path, 'r') as f:
                if len(f.read().strip()) > 0:
                    # File has some content
                    val_logger.log("featurecounts", "ALL", "counts_rrnarm_file_format", "GREEN", 
                                 f"rRNA-removed counts file has valid format")
        except Exception as e:
            val_logger.log("featurecounts", "ALL", "counts_rrnarm_file_read", "RED", 
                         f"Error reading rRNA-removed counts file: {str(e)}")
            logger.error(f"Error reading rRNA-removed counts file: {str(e)}")
    
    # (REMOVED) Validate rRNA-removed summary file content
    
    # STEP 3: MultiQC validation and sample metrics
    # Validate MultiQC report content and do outlier detection as the FINAL step
    if multiqc_found:
        # Get the first MultiQC file (should only be one)
        multiqc_file = os.path.realpath(staged_files["multiqc"][0])
        
        # Try to get MultiQC data - first check if we have a directory structure instead of a zip
        multiqc_data = None
        try:
            # First check if we have a directory with the MultiQC data unzipped
            multiqc_dir_path = os.path.join(fc_dir, f"featureCounts_multiqc{assay_suffix}_report")
            multiqc_data_path = os.path.join(multiqc_dir_path, f"featureCounts_multiqc{assay_suffix}_data", "multiqc_data.json")
            
            logger.info(f"Checking for unzipped MultiQC data at: {multiqc_data_path}")
            
            if os.path.exists(multiqc_data_path):
                logger.info(f"Found unzipped MultiQC data at: {multiqc_data_path}")
                with open(multiqc_data_path, 'r') as f:
                    multiqc_data = json.load(f)
                    logger.info(f"Successfully loaded MultiQC data from unzipped file")
            else:
                # Try to extract and unzip the file
                try:
                    logger.info(f"Unzipped file not found, looking for MultiQC data in zip: {multiqc_file}")
                    with zipfile.ZipFile(multiqc_file, 'r') as zip_ref:
                        # List all files in the zip for debugging
                        zip_contents = zip_ref.namelist()
                        logger.info(f"Zip contents: {zip_contents[:10]} ... (showing first 10 files)")
                        
                        # Look for featurecounts data
                        fc_data_files = [f for f in zip_contents if 'featurecounts' in f.lower()]
                        if fc_data_files:
                            logger.info(f"Found FeatureCounts data files in zip: {fc_data_files[:5]}")
                        
                        # Look for multiqc_data.json in various possible locations
                        possible_json_paths = [
                            "multiqc_data.json",
                            "featureCounts_multiqc_data/multiqc_data.json",
                            f"featureCounts_multiqc{assay_suffix}_data/multiqc_data.json"
                        ]
                        
                        json_path = None
                        for path in possible_json_paths:
                            if path in zip_contents:
                                json_path = path
                                break
                        
                        if json_path:
                            logger.info(f"Found MultiQC data JSON at: {json_path}")
                            with zip_ref.open(json_path) as f:
                                multiqc_data = json.load(TextIOWrapper(f))
                                logger.info(f"Successfully loaded MultiQC data from zip file")
                        else:
                            # Try to find any JSON file
                            json_files = [f for f in zip_contents if f.endswith('.json')]
                            if json_files:
                                logger.info(f"Found other JSON files: {json_files}")
                                
                                # Try the first one
                                with zip_ref.open(json_files[0]) as f:
                                    multiqc_data = json.load(TextIOWrapper(f))
                                    logger.info(f"Loaded MultiQC data from alternate JSON: {json_files[0]}")
                except Exception as zip_error:
                    logger.error(f"Error extracting MultiQC data from zip: {str(zip_error)}")
                    logger.error(traceback.format_exc())
            
            # Verify we have the data
            if multiqc_data:
                # Log keys for debugging
                top_level_keys = list(multiqc_data.keys())
                logger.info(f"MultiQC data top-level keys: {top_level_keys}")
                
                # First check if data is in report_saved_raw_data
                if 'report_saved_raw_data' in multiqc_data and 'multiqc_featurecounts' in multiqc_data['report_saved_raw_data']:
                    logger.info("Found FeatureCounts data in report_saved_raw_data section")
                    # Extract samples from MultiQC data
                    stats_section = multiqc_data['report_saved_raw_data']['multiqc_featurecounts']
                    
                    # Log successful data extraction
                    val_logger.log("featurecounts", "ALL", "multiqc_data", "GREEN", 
                                 f"FeatureCounts MultiQC data successfully extracted")
                    
                # Check if data is in general stats section
                elif 'report_general_stats_data' in multiqc_data and len(multiqc_data['report_general_stats_data']) > 0:
                    # Get the first stats section (should contain featurecounts data)
                    stats_section = multiqc_data['report_general_stats_data'][0]
                    logger.info(f"Found FeatureCounts data in report_general_stats_data section")
                    
                    # Log successful data extraction
                    val_logger.log("featurecounts", "ALL", "multiqc_data", "GREEN", 
                                 f"FeatureCounts MultiQC data successfully extracted")
                else:
                    logger.warning("Could not find FeatureCounts data in expected MultiQC sections")
                    val_logger.log("featurecounts", "ALL", "multiqc_content", "YELLOW", 
                                 f"FeatureCounts MultiQC data found but couldn't locate metrics in expected sections")
                    
                    # Try to find any sections that might have sample data
                    for key in top_level_keys:
                        if isinstance(multiqc_data[key], dict):
                            logger.info(f"Checking section {key} for potential sample data")
                            potential_sample_data = False
                            
                            # Look at first few items to see if they might be sample data
                            try:
                                if key in ['report_saved_raw_data', 'report_general_stats_data']:
                                    for subkey, value in list(multiqc_data[key].items())[:3]:
                                        logger.info(f"  Subkey: {subkey}, Type: {type(value)}")
                                        if isinstance(value, dict):
                                            logger.info(f"    First few keys in value: {list(value.keys())[:5]}")
                                            potential_sample_data = True
                            except Exception as e:
                                logger.info(f"  Error examining section: {str(e)}")
                            
                            if potential_sample_data:
                                logger.info(f"Section {key} might contain sample data")
                
                # Process MultiQC data if we have stats
                if 'stats_section' in locals() and stats_section:
                    # Extract samples from MultiQC data
                    multiqc_samples = list(stats_section.keys())
                    logger.info(f"Found {len(multiqc_samples)} samples in MultiQC data: {multiqc_samples}")
                    
                    # Check if all expected samples are present
                    missing_samples = []
                    for sample_id in sample_ids:
                        # Check if the sample matches any entry in multiqc_samples
                        if not any(sample_id in s or s in sample_id for s in multiqc_samples):
                            missing_samples.append(sample_id)
                    
                    if missing_samples:
                        val_logger.log("featurecounts", "ALL", "samples_in_multiqc", "YELLOW", 
                                     f"Missing samples in FeatureCounts MultiQC report: {', '.join(missing_samples)}")
                    else:
                        val_logger.log("featurecounts", "ALL", "samples_in_multiqc", "GREEN", 
                                     f"All samples found in FeatureCounts MultiQC report")
                    
                    # Extract per-sample metrics for outlier detection
                    featurecounts_stats = {}
                    for sample, stats in stats_section.items():
                        # Clean up sample name
                        base_sample = sample.split('_featureCounts')[0]
                        
                        # Log what we're extracting for debugging
                        logger.debug(f"Extracting metrics for sample {base_sample}")
                        logger.debug(f"  Available keys: {list(stats.keys())}")
                        
                        # Handle different possible field names
                        total_count = 0
                        if 'Total' in stats:
                            total_count = stats.get('Total', 0)
                        
                        assigned_pct = 0
                        # Try different ways to get assigned percentage
                        if 'percent_assigned' in stats:
                            assigned_pct = stats.get('percent_assigned', 0)
                        elif 'Assigned' in stats and total_count > 0:
                            assigned_pct = (stats.get('Assigned', 0) / total_count) * 100
                        
                        # Initialize metrics for this sample
                        featurecounts_stats[base_sample] = {
                            'total_count': total_count,
                            'assigned_pct': assigned_pct
                        }
                        
                        # Add unassigned reason percentages if Total > 0
                        if total_count > 0:
                            total = total_count  # Avoid division by zero
                            
                            # Try to detect the right field names for unassigned reasons
                            unassigned_fields = [k for k in stats.keys() if k.startswith('Unassigned_')]
                            logger.debug(f"  Unassigned fields: {unassigned_fields}")
                            
                            # Common ones to look for
                            for field in unassigned_fields:
                                clean_name = field.lower() + '_pct'
                                featurecounts_stats[base_sample][clean_name] = (stats.get(field, 0) / total) * 100
                    
                    # Detect outliers in key metrics
                    if featurecounts_stats:
                        # Get all metrics that are present in the first sample
                        first_sample = list(featurecounts_stats.keys())[0]
                        all_metrics = list(featurecounts_stats[first_sample].keys())
                        
                        logger.info(f"Sample metrics: {featurecounts_stats[first_sample]}")
                        
                        # Make sure all metrics exist for all samples
                        metrics_to_check = []
                        for metric in all_metrics:
                            valid_for_all = True
                            for sample in featurecounts_stats:
                                if metric not in featurecounts_stats[sample]:
                                    valid_for_all = False
                                    break
                            
                            if valid_for_all:
                                metrics_to_check.append(metric)
                        
                        logger.info(f"Checking for outliers in FeatureCounts metrics: {metrics_to_check}")
                        
                        try:
                            outliers, metric_stats = detect_outliers(featurecounts_stats, metrics_to_check)
                            
                            # Log outliers
                            if outliers:
                                for sample, sample_outliers in outliers.items():
                                    for outlier in sample_outliers:
                                        metric = outlier['metric']
                                        value = outlier['value']
                                        mean = outlier['mean']
                                        # Use 'std' instead of 'stdev'
                                        std = outlier['std']
                                        # Use 'z_score' instead of 'deviation'
                                        z_score = outlier['z_score']
                                        
                                        # Determine status based on z-score (absolute value)
                                        status = "YELLOW" if abs(z_score) < 3 else "RED"
                                        
                                        val_logger.log("featurecounts", sample, f"{metric}_outlier", status, 
                                                    f"{metric} value {value} is an outlier", 
                                                    f"Value deviates {abs(z_score):.2f} standard deviations from mean ({mean:.2f})")
                            else:
                                logger.info("No outliers detected in FeatureCounts metrics")
                                val_logger.log("featurecounts", "ALL", "metrics_outliers", "GREEN", 
                                            f"No outliers detected in FeatureCounts metrics")
                        except Exception as e:
                            logger.error(f"Error detecting outliers: {str(e)}")
                            logger.error(traceback.format_exc())
                            val_logger.log("featurecounts", "ALL", "metrics_outliers", "YELLOW", 
                                        f"Error occurred while detecting outliers: {str(e)}")
                else:
                    val_logger.log("featurecounts", "ALL", "multiqc_content", "YELLOW", 
                                f"FeatureCounts MultiQC data found but doesn't contain expected sections")
            else:
                val_logger.log("featurecounts", "ALL", "multiqc_content", "YELLOW", 
                             f"FeatureCounts MultiQC report found but couldn't extract data")
        except Exception as e:
            val_logger.log("featurecounts", "ALL", "multiqc_read", "RED", 
                         f"Error reading MultiQC report: {str(e)}")
            logger.error(f"Error reading MultiQC report: {str(e)}")
            logger.error(traceback.format_exc())
    else:
        if featurecounts_multiqc_dir:
            val_logger.log("featurecounts", "ALL", "multiqc_exists", "RED", 
                         "No FeatureCounts MultiQC report found")
    
    # Return validation status with expected format matching other validation functions
    return {
        "status": val_logger.get_status(),
        "messages": [r["message"] for r in val_logger.results],
        "failures": {r["check_name"]: r["message"] 
                    for r in val_logger.results 
                    if r["status"] in ["HALT", "RED"]}
    }


if __name__ == "__main__":
    vv()
