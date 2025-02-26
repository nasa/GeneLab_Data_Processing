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
def vv(assay_type, assay_suffix, runsheet_path, outdir, paired_end, mode, run_components, 
       raw_fastq, raw_fastqc, raw_multiqc,
       trimmed_fastq, trimmed_fastqc, trimmed_multiqc, trimming_reports, trimming_multiqc,
       bowtie2_alignment_log, bowtie2_alignment_unmapped, bowtie2_alignment_multiqc, 
       bowtie2_alignment_sorted, bowtie2_alignment_sorted_index,
       genebody_coverage, infer_experiment, inner_distance, read_distribution):
    """Organize pipeline outputs and optionally validate"""
    outdir = Path(outdir)
    
    # Initialize validation logger
    global val_logger
    val_logger = ValidationLogger()
    
    # Stage files if inputs provided
    if any([raw_fastq, raw_fastqc, raw_multiqc]):
        file_paths = {
            'raw_fastq': raw_fastq,
            'raw_fastqc': raw_fastqc,
            'raw_multiqc': raw_multiqc
        }
        stage_files(assay_type, 'raw_reads', **file_paths)
    
    # Stage trimmed files if inputs provided
    if any([trimmed_fastq, trimmed_fastqc, trimmed_multiqc, trimming_reports, trimming_multiqc]):
        file_paths = {
            'fastq': trimmed_fastq,
            'fastqc': trimmed_fastqc,
            'trimmed_multiqc': trimmed_multiqc,
            'trimming_reports': trimming_reports,
            'trimming_multiqc': trimming_multiqc
        }
        stage_files(assay_type, 'trimmed_reads', **file_paths)
    
    # Stage Bowtie2 alignment files if inputs provided
    if any([bowtie2_alignment_log, bowtie2_alignment_unmapped, bowtie2_alignment_multiqc, 
            bowtie2_alignment_sorted, bowtie2_alignment_sorted_index]):
        
        # Create 02-Bowtie2_Alignment directory in CURRENT WORKING DIRECTORY for Nextflow
        work_alignment_dir = Path("02-Bowtie2_Alignment")
        os.makedirs(work_alignment_dir, exist_ok=True)
        
        # Keep track of all sample names to create directories
        sample_names = set()
        
        # Get sample names from all files
        # From sorted BAM files (these have the most reliable pattern)
        if bowtie2_alignment_sorted:
            for file in os.listdir(bowtie2_alignment_sorted):
                if file.endswith('_sorted.bam'):
                    # Extract sample name from file name (e.g., "Sample1_sorted.bam" -> "Sample1")
                    sample_name = file.replace('_sorted.bam', '')
                    sample_names.add(sample_name)
        
        # Create sample directories in working directory
        for sample in sample_names:
            # Create in working directory for Nextflow
            work_sample_dir = work_alignment_dir / sample
            os.makedirs(work_sample_dir, exist_ok=True)
        
        # Stage alignment logs
        if bowtie2_alignment_log:
            for file in os.listdir(bowtie2_alignment_log):
                src = os.path.realpath(os.path.join(bowtie2_alignment_log, file))
                
                # Special handling for sample-specific files
                if any(sample in file for sample in sample_names):
                    # Find which sample this file belongs to
                    for sample in sample_names:
                        if sample in file:
                            # Link to work directory for Nextflow
                            work_sample_dir = work_alignment_dir / sample
                            work_dst = os.path.join(work_sample_dir, file)
                            if not os.path.exists(work_dst) and not os.path.islink(work_dst):
                                os.symlink(src, work_dst)
                            break
                else:
                    # Handle non-sample specific files (keep at root level)
                    # Link to work directory for Nextflow
                    work_dst = os.path.join(work_alignment_dir, file)
                    if not os.path.exists(work_dst) and not os.path.islink(work_dst):
                        os.symlink(src, work_dst)
        
        # Stage unmapped reads
        if bowtie2_alignment_unmapped:
            for file in os.listdir(bowtie2_alignment_unmapped):
                src = os.path.realpath(os.path.join(bowtie2_alignment_unmapped, file))
                
                # Find which sample this file belongs to
                for sample in sample_names:
                    if sample in file:
                        # Link to work directory for Nextflow
                        work_sample_dir = work_alignment_dir / sample
                        work_dst = os.path.join(work_sample_dir, file)
                        if not os.path.exists(work_dst) and not os.path.islink(work_dst):
                            os.symlink(src, work_dst)
                        break
                else:
                    # If no sample match found, place at root level
                    # Link to work directory for Nextflow
                    work_dst = os.path.join(work_alignment_dir, file)
                    if not os.path.exists(work_dst) and not os.path.islink(work_dst):
                        os.symlink(src, work_dst)
        
        # Stage alignment MultiQC report (at root level since it's not sample-specific)
        if bowtie2_alignment_multiqc:
            for file in os.listdir(bowtie2_alignment_multiqc):
                src = os.path.realpath(os.path.join(bowtie2_alignment_multiqc, file))
                
                # Link to work directory for Nextflow
                work_dst = os.path.join(work_alignment_dir, file)
                if not os.path.exists(work_dst) and not os.path.islink(work_dst):
                    os.symlink(src, work_dst)
        
        # Stage sorted BAM files
        if bowtie2_alignment_sorted:
            for file in os.listdir(bowtie2_alignment_sorted):
                src = os.path.realpath(os.path.join(bowtie2_alignment_sorted, file))
                
                # Find which sample this file belongs to
                for sample in sample_names:
                    if sample in file:
                        # Link to work directory for Nextflow
                        work_sample_dir = work_alignment_dir / sample
                        work_dst = os.path.join(work_sample_dir, file)
                        if not os.path.exists(work_dst) and not os.path.islink(work_dst):
                            os.symlink(src, work_dst)
                        break
                else:
                    # If no sample match found, place at root level
                    # Link to work directory for Nextflow
                    work_dst = os.path.join(work_alignment_dir, file)
                    if not os.path.exists(work_dst) and not os.path.islink(work_dst):
                        os.symlink(src, work_dst)
        
        # Stage sorted BAM index files
        if bowtie2_alignment_sorted_index:
            for file in os.listdir(bowtie2_alignment_sorted_index):
                src = os.path.realpath(os.path.join(bowtie2_alignment_sorted_index, file))
                
                # Find which sample this file belongs to
                for sample in sample_names:
                    if sample in file:
                        # Link to work directory for Nextflow
                        work_sample_dir = work_alignment_dir / sample
                        work_dst = os.path.join(work_sample_dir, file)
                        if not os.path.exists(work_dst) and not os.path.islink(work_dst):
                            os.symlink(src, work_dst)
                        break
                else:
                    # If no sample match found, place at root level
                    # Link to work directory for Nextflow
                    work_dst = os.path.join(work_alignment_dir, file)
                    if not os.path.exists(work_dst) and not os.path.islink(work_dst):
                        os.symlink(src, work_dst)
    
    # Run validation if component specified
    if run_components:
        # Convert paired_end string to bool
        is_paired = paired_end.lower() == 'true' if isinstance(paired_end, str) else paired_end
        
        if 'raw_reads' in run_components:
            # Run validation
            results = validate_raw_reads(
                outdir=outdir,
                samples_txt=Path(runsheet_path),  # Using runsheet as samples file for now
                paired_end=is_paired,
                assay_suffix=assay_suffix
            )
        
        elif 'trimmed_reads' in run_components:
            # Run validation for trimmed reads
            results = validate_trimmed_reads(
                outdir=outdir,
                samples_txt=Path(runsheet_path),
                paired_end=is_paired,
                assay_suffix=assay_suffix
            )
        
        elif 'bowtie2_alignment' in run_components or 'alignments' in run_components:
            # Run validation for Bowtie2 alignments
            if mode != 'microbes':
                logger.warning(f"Bowtie2 alignment validation is only supported for microbes mode, current mode: {mode}")
            
            results = validate_bowtie2_alignments(
                outdir=outdir,
                samples_txt=Path(runsheet_path),
                paired_end=is_paired,
                assay_suffix=assay_suffix
            )
        
        # Add RSeQC validation
        if 'rseqc' in run_components:
            results = validate_rseqc(
                outdir=outdir,
                samples_txt=runsheet_path,
                paired_end=is_paired,
                assay_suffix=assay_suffix,
                genebody_coverage_dir=Path(genebody_coverage) if genebody_coverage else None,
                infer_experiment_dir=Path(infer_experiment) if infer_experiment else None,
                inner_distance_dir=Path(inner_distance) if inner_distance else None,
                read_distribution_dir=Path(read_distribution) if read_distribution else None
            )
        
        return results

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

def get_expected_files(sample_id: str, paired_end: bool, assay_suffix: str) -> dict:
    """
    Get expected file patterns for a sample
    
    Args:
        sample_id: Sample name from runsheet
        paired_end: Whether data is paired-end
        assay_suffix: Assay suffix (e.g. "_GLbulkRNAseq")
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

def validate_raw_reads(outdir: Path,
                      samples_txt: Path,
                      paired_end: bool = True,
                      assay_suffix: str = "_GLbulkRNAseq") -> dict:
    """
    Original dp_tools raw reads validation checks
    """
    val_logger = ValidationLogger()
    read_counts = {}  # Track read counts for paired-end validation
    
    # Log validation parameters
    logging.info("Starting raw reads validation:")
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
                    with gzip.open(fastq_path, "rt") as f:
                        for i, line in enumerate(f):
                            if i % 4 == 0:  # Header lines only
                                total_reads += 1
                                if not line.startswith('@'):
                                    issues.append(i+1)
                            if i % 2_000_000 == 0:  # Log progress
                                logging.debug(f"Checked {i} lines in {fastq_path}")
                                
                    stats = {
                        "total_reads": total_reads,
                        "invalid_headers": len(issues),
                        "checked_lines": i + 1
                    }
                                
                    if issues:
                        val_logger.log("raw_reads", sample, "fastq_format", "HALT",
                                 f"Invalid FASTQ format - headers missing @ at lines: {issues[:10]}..." if len(issues) > 10 else f"Invalid FASTQ format - headers missing @ at lines: {issues}",
                                 stats=stats)
                    else:
                        val_logger.log("raw_reads", sample, "fastq_format", "GREEN",
                                 "FASTQ format check passed",
                                 stats=stats)

                    # Store read count for paired-end check later
                    read_counts[fastq_path] = total_reads
                                
                except Exception as e:
                    val_logger.log("raw_reads", sample, "fastq_format", "HALT",
                             f"Error checking FASTQ format", str(e))
            
            # 2b. If paired-end, check read counts match with stats
            if paired_end:
                try:
                    r1_path = fastq_dir / expected["fastq"][0]
                    r2_path = fastq_dir / expected["fastq"][1]
                    
                    # Use counts from format checking
                    r1_count = read_counts.get(r1_path, 0)
                    r2_count = read_counts.get(r2_path, 0)
                    
                    stats = {
                        "r1_reads": r1_count,
                        "r2_reads": r2_count,
                        "difference": abs(r1_count - r2_count)
                    }
                    
                    if r1_count != r2_count:
                        val_logger.log("raw_reads", sample, "paired_counts", "HALT",
                                 f"Paired read counts don't match",
                                 f"R1={r1_count}, R2={r2_count}",
                                 stats=stats)
                    else:
                        val_logger.log("raw_reads", sample, "paired_counts", "GREEN",
                                 f"Paired read counts match ({r1_count} reads)",
                                 stats=stats)
                except Exception as e:
                    val_logger.log("raw_reads", sample, "paired_counts", "HALT",
                             f"Error checking read counts", str(e))
                
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
        multiqc = expected["multiqc"][0]
        multiqc_path = fastqc_dir / multiqc
        if not multiqc_path.exists():
            val_logger.log("raw_reads", "ALL", "multiqc_exists", "HALT",
                      f"Missing MultiQC report: {multiqc}")
        else:
            val_logger.log("raw_reads", "ALL", "multiqc_exists", "GREEN",
                      f"MultiQC report found: {multiqc}")
            
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
            
            # Extract total processed reads
            processed_match = re.search(r'Processed reads:\s+(\d+)', content)
            if processed_match:
                stats["total_processed_reads"] = int(processed_match.group(1))
                
            # Check if adapter was detected and its sequence/type
            adapter_match = re.search(r'Adapter sequence:\s+\'([ACGT]+)\'', content)
            if adapter_match:
                stats["adapters_found"] = True
                stats["adapter_sequence"] = adapter_match.group(1)
                # Try to identify adapter type from TrimGalore output
                if "Illumina" in content:
                    stats["adapter_type"] = "Illumina"
                elif "Nextera" in content:
                    stats["adapter_type"] = "Nextera"
                elif "smallRNA" in content:
                    stats["adapter_type"] = "smallRNA"
                
            # Extract reads with adapters (pattern may vary based on TrimGalore version)
            adapter_count_match = re.search(r'Reads with adapters:\s+(\d+)\s+\(([0-9.]+)%\)', content)
            if adapter_count_match:
                stats["reads_with_adapters"] = int(adapter_count_match.group(1))
                stats["adapter_trimmed_percentage"] = float(adapter_count_match.group(2))
            # Alternative pattern
            elif "adapters were trimmed" in content:
                percentage_match = re.search(r'([0-9.]+)% of reads contained adapter', content)
                if percentage_match:
                    stats["adapter_trimmed_percentage"] = float(percentage_match.group(1))
                    stats["adapters_found"] = True
                
            # Extract quality trimming info
            quality_match = re.search(r'Quality-trimmed:\s+(\d+)\s+\(([0-9.]+)%\)', content)
            if quality_match:
                stats["quality_trimmed_reads"] = int(quality_match.group(1))
                stats["quality_trimmed_percentage"] = float(quality_match.group(2))
                
    except Exception as e:
        logging.error(f"Error parsing trimming report {report_path}: {e}")
        
    return stats

def validate_trimmed_reads(outdir: Path,
                       samples_txt: Path,
                       paired_end: bool = True,
                       assay_suffix: str = "_GLbulkRNAseq") -> dict:
    """
    Validation checks for trimmed reads and QC outputs
    """
    val_logger = ValidationLogger()
    read_counts = {}  # Track read counts for paired-end validation
    
    # Log validation parameters
    logging.info("Starting trimmed reads validation:")
    logging.info(f"  Output directory: {outdir}")
    logging.info(f"  Samples file: {samples_txt}")
    logging.info(f"  Paired-end: {paired_end}")
    logging.info(f"  Assay suffix: {assay_suffix}")

    try:
        # 1. Basic File Existence Checks
        fastq_dir = outdir / "01-TG_Preproc/Fastq"
        fastqc_dir = outdir / "01-TG_Preproc/FastQC_Reports"
        trimming_dir = outdir / "01-TG_Preproc/Trimming_Reports"
        raw_fastqc_dir = outdir / "00-RawData/FastQC_Reports"
        
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
            logging.info(f"\nValidating trimmed sample: {sample}")
            
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
                
                # Check FASTQ format with stats
                try:
                    issues = []
                    total_reads = 0
                    with gzip.open(fastq_path, "rt") as f:
                        for i, line in enumerate(f):
                            if i % 4 == 0:  # Header lines only
                                total_reads += 1
                                if not line.startswith('@'):
                                    issues.append(i+1)
                            if i % 2_000_000 == 0:  # Log progress
                                logging.debug(f"Checked {i} lines in {fastq_path}")
                                
                    stats = {
                        "total_reads": total_reads,
                        "invalid_headers": len(issues),
                        "checked_lines": i + 1
                    }
                                
                    if issues:
                        val_logger.log("trimmed_reads", sample, "fastq_format", "HALT",
                                 f"Invalid FASTQ format - headers missing @ at lines: {issues[:10]}..." if len(issues) > 10 else f"Invalid FASTQ format - headers missing @ at lines: {issues}",
                                 stats=stats)
                    else:
                        val_logger.log("trimmed_reads", sample, "fastq_format", "GREEN",
                                 "FASTQ format check passed",
                                 stats=stats)

                    # Store read count for paired-end check later
                    read_counts[fastq_path] = total_reads
                                
                except Exception as e:
                    val_logger.log("trimmed_reads", sample, "fastq_format", "HALT",
                             f"Error checking FASTQ format", str(e))
            
            # 2b. If paired-end, check read counts match with stats
            if paired_end:
                try:
                    r1_path = fastq_dir / expected["fastq"][0]
                    r2_path = fastq_dir / expected["fastq"][1]
                    
                    # Use counts from format checking
                    r1_count = read_counts.get(r1_path, 0)
                    r2_count = read_counts.get(r2_path, 0)
                    
                    stats = {
                        "r1_reads": r1_count,
                        "r2_reads": r2_count,
                        "difference": abs(r1_count - r2_count)
                    }
                    
                    if r1_count != r2_count:
                        val_logger.log("trimmed_reads", sample, "paired_counts", "HALT",
                                 f"Paired read counts don't match",
                                 f"R1={r1_count}, R2={r2_count}",
                                 stats=stats)
                    else:
                        val_logger.log("trimmed_reads", sample, "paired_counts", "GREEN",
                                 f"Paired read counts match ({r1_count} reads)",
                                 stats=stats)
                except Exception as e:
                    val_logger.log("trimmed_reads", sample, "paired_counts", "HALT",
                             f"Error checking read counts", str(e))
                
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
        # 3a. Trimmed reads MultiQC
        trimmed_multiqc = expected["trimmed_multiqc"][0]
        trimmed_multiqc_path = fastqc_dir / trimmed_multiqc
        if not trimmed_multiqc_path.exists():
            val_logger.log("trimmed_reads", "ALL", "trimmed_multiqc_exists", "HALT",
                      f"Missing trimmed reads MultiQC report: {trimmed_multiqc}")
        else:
            val_logger.log("trimmed_reads", "ALL", "trimmed_multiqc_exists", "GREEN",
                      f"Trimmed reads MultiQC report found: {trimmed_multiqc}")
        
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
    
    # Calculate directory paths - use work directory
    alignment_dir = Path("02-Bowtie2_Alignment")
    if not alignment_dir.exists():
        logger.error(f"Work directory {alignment_dir} does not exist")
        return {
            "status": "failed",
            "messages": [f"Work directory {alignment_dir} does not exist"],
            "failures": {"directory_exists": f"Work directory {alignment_dir} does not exist"}
        }
    
    # Get sample names from runsheet (assumes CSV format for now)
    samples = []
    with open(samples_txt, "r") as f:
        reader = csv.reader(f)
        header = next(reader)  # Skip header
        for row in reader:
            if row:  # Skip empty rows
                samples.append(row[-1].strip())  # Assume last column is sample name
    
    # Track alignment rates across samples for outlier detection
    alignment_rates = []
    sample_names = []  # Track sample names in same order as metrics
    
    # Validate each sample's alignment outputs
    for sample in samples:
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
                    
                    # Log the alignment rate without threshold judgments
                    val_logger.log("alignments", sample, "alignment_rate", "GREEN", 
                                f"Alignment rate: {value}%")
                    
            except Exception as e:
                val_logger.log("alignments", sample, "alignment_rate", "HALT", 
                            f"Error parsing Bowtie2 log: {str(e)}")
        
        # Check unmapped reads files exist
        if paired_end:
            unmapped_1 = sample_dir / f"{sample}.unmapped.fastq.1.gz"
            unmapped_2 = sample_dir / f"{sample}.unmapped.fastq.2.gz"
            if not unmapped_1.exists() or not unmapped_2.exists():
                val_logger.log("alignments", sample, "unmapped_reads_exist", "HALT",
                            f"Missing unmapped reads files: {unmapped_1.name} and/or {unmapped_2.name}")
            else:
                val_logger.log("alignments", sample, "unmapped_reads_exist", "GREEN",
                            f"Found unmapped reads files: {unmapped_1.name} and {unmapped_2.name}")
        else:
            unmapped = sample_dir / f"{sample}.unmapped.fastq.gz"
            if not unmapped.exists():
                val_logger.log("alignments", sample, "unmapped_reads_exist", "HALT",
                            f"Missing unmapped reads file: {unmapped.name}")
            else:
                val_logger.log("alignments", sample, "unmapped_reads_exist", "GREEN",
                            f"Found unmapped reads file: {unmapped.name}")
    
    # Check for outliers in alignment rate
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
                    
                    # Also check general samples list
                    if 'report_general_stats_data' in multiqc_data:
                        for stats in multiqc_data['report_general_stats_data']:
                            found_samples.update(stats.keys())
                    
                    # Improved sample name matching logic
                    clean_found_samples = set()
                    logger.debug(f"Raw sample names in MultiQC: {found_samples}")
                    
                    # For each sample we're expecting
                    for sample in samples:
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
                        msg = f"Samples missing from MultiQC report for {section_name}: {', '.join(missing_samples)}"
                        logger.warning(msg)
                        val_logger.log("rseqc", "ALL", f"{section_name}_samples_in_multiqc", "YELLOW", msg)
                    else:
                        val_logger.log("rseqc", "ALL", f"{section_name}_samples_in_multiqc", "GREEN",
                                    f"All samples found in MultiQC report for {section_name}")
                else:
                    logger.warning(f"Could not find multiqc_data.json in {multiqc_path}")
        except Exception as e:
            val_logger.log("alignments", "ALL", "samples_in_multiqc", "HALT",
                        f"Error checking samples in MultiQC report: {str(e)}")
    
    return {
        "status": val_logger.get_status(),
        "messages": [r["message"] for r in val_logger.results],
        "failures": {r["check_name"]: r["message"] 
                    for r in val_logger.results 
                    if r["status"] in ["HALT", "RED"]}
    }

def extract_adapter_content_from_fastqc(fastqc_zip_path):
    """
    Extract adapter content information from a FastQC zip file
    
    Args:
        fastqc_zip_path: Path to FastQC zip file
        
    Returns:
        Dictionary with adapter content information
    """
    adapter_stats = {
        "adapter_found": False,
        "adapter_percentage": 0.0,
        "adapter_details": {}
    }
    
    try:
        with zipfile.ZipFile(fastqc_zip_path, 'r') as zip_ref:
            # Find the adapter content data file
            data_file = None
            for file in zip_ref.namelist():
                if file.endswith('fastqc_data.txt'):
                    data_file = file
                    break
            
            if not data_file:
                logging.warning(f"Could not find fastqc_data.txt in {fastqc_zip_path}")
                return adapter_stats
            
            # Parse the file
            with zip_ref.open(data_file) as f:
                text_f = TextIOWrapper(f)
                in_adapter_section = False
                for line in text_f:
                    line = line.strip()
                    
                    # Find adapter content section
                    if line == ">>Adapter Content":
                        in_adapter_section = True
                        continue
                    
                    if in_adapter_section:
                        if line.startswith('>>END_MODULE'):
                            break
                        
                        if line.startswith('#'):
                            continue  # Skip header line
                        
                        # Parse adapter content data
                        parts = line.split('\t')
                        if len(parts) > 2:
                            position = parts[0]
                            # Get maximum adapter percentage from any adapter type
                            max_pct = max([float(x.strip()) for x in parts[1:]])
                            
                            adapter_stats["adapter_details"][position] = max_pct
                            
                            # If any position has >0.1% adapter, consider it found
                            if max_pct > 0.1:
                                adapter_stats["adapter_found"] = True
                                
                                # Use the maximum percentage found at any position
                                adapter_stats["adapter_percentage"] = max(adapter_stats["adapter_percentage"], max_pct)
    
    except Exception as e:
        logging.error(f"Error extracting adapter content from {fastqc_zip_path}: {e}")
    
    return adapter_stats

def validate_rseqc(outdir: Path,
                  samples_txt: Path,
                  paired_end: bool,
                  assay_suffix: str,
                  genebody_coverage_dir: Path = None,
                  infer_experiment_dir: Path = None,
                  inner_distance_dir: Path = None,
                  read_distribution_dir: Path = None) -> dict:
    """
    Validate RSeQC outputs and organize into the correct directory structure
    """
    # Create main RSeQC directory in the published directory structure
    published_rseqc_dir = outdir / "RSeQC_Analyses"
    if published_rseqc_dir.exists():
        logger.info(f"Removing existing RSeQC directory: {published_rseqc_dir}")
        shutil.rmtree(published_rseqc_dir)
    os.makedirs(published_rseqc_dir, exist_ok=True)
    
    # Also create working directory for Nextflow
    rseqc_dir = Path("RSeQC_Analyses")
    if rseqc_dir.exists():
        logger.info(f"Removing existing RSeQC directory in working dir: {rseqc_dir}")
        shutil.rmtree(rseqc_dir)
    os.makedirs(rseqc_dir, exist_ok=True)
    
    # Log validation parameters
    logger.info("Starting RSeQC validation:")
    logger.info(f"  Output directory: {outdir}")
    logger.info(f"  Samples file: {samples_txt}")
    logger.info(f"  Paired-end: {paired_end}")
    logger.info(f"  Assay suffix: {assay_suffix}")
    logger.info(f"  Genebody coverage dir: {genebody_coverage_dir}")
    logger.info(f"  Infer experiment dir: {infer_experiment_dir}")
    logger.info(f"  Inner distance dir: {inner_distance_dir if paired_end else 'N/A - single end data'}")
    logger.info(f"  Read distribution dir: {read_distribution_dir}")
    
    # Create section directories with correct names, but use actual RSeQC tool names for display
    # Define a mapping between directory names and the actual RSeQC script names for display in logs
    sections = {
        "02_geneBody_coverage": {
            "display_name": "geneBody_coverage.py",
            "dir": rseqc_dir / "02_geneBody_coverage",
            "pub_dir": published_rseqc_dir / "02_geneBody_coverage",
            "input": genebody_coverage_dir,
            "file_patterns": ["geneBodyCoverage"],
            "required": True
        },
        "03_infer_experiment": {
            "display_name": "infer_experiment.py",
            "dir": rseqc_dir / "03_infer_experiment",
            "pub_dir": published_rseqc_dir / "03_infer_experiment",
            "input": infer_experiment_dir,
            "file_patterns": ["infer_expt"],
            "required": True
        },
        "04_inner_distance": {
            "display_name": "inner_distance.py",
            "dir": rseqc_dir / "04_inner_distance" if paired_end else None,
            "pub_dir": published_rseqc_dir / "04_inner_distance" if paired_end else None,
            "input": inner_distance_dir if paired_end else None,
            "file_patterns": ["inner_distance"],
            "required": paired_end
        },
        "05_read_distribution": {
            "display_name": "read_distribution.py",
            "dir": rseqc_dir / "05_read_distribution",
            "pub_dir": published_rseqc_dir / "05_read_distribution",
            "input": read_distribution_dir,
            "file_patterns": ["read_dist"],
            "required": True
        }
    }
    
    # Create directories - both in work dir and published dir
    for name, section in sections.items():
        if section["dir"]:  # Skip inner_distance if single-end
            os.makedirs(section["dir"], exist_ok=True)
            os.makedirs(section["pub_dir"], exist_ok=True)
            logger.info(f"Created directory: {section['dir']}")
            logger.info(f"Created published directory: {section['pub_dir']}")
    
    # Read samples from CSV
    samples = []
    try:
        with open(samples_txt) as f:
            reader = csv.DictReader(f)
            for row in reader:
                if 'Sample Name' in row:
                    sample_name = row['Sample Name'].strip()
                    if sample_name:  # Only add non-empty sample names
                        samples.append(sample_name)
                        logger.debug(f"Added sample: {sample_name}")
                else:
                    logger.warning("Runsheet missing required 'Sample Name' column")
    except Exception as e:
        logger.error(f"Error reading samples from runsheet: {e}")
    
    logger.info(f"Found {len(samples)} samples to validate")
    logger.debug(f"Sample names: {samples}")
    
    # Create ValidationLogger
    val_logger = ValidationLogger()
    
    # Create sample directories in each section
    for name, section in sections.items():
        if section["dir"]:
            for sample in samples:
                # Create in working directory for Nextflow
                sample_dir = section["dir"] / sample
                os.makedirs(sample_dir, exist_ok=True)
                
                # Create in published directory structure
                pub_sample_dir = section["pub_dir"] / sample
                os.makedirs(pub_sample_dir, exist_ok=True)
                
                logger.debug(f"Created sample directories: {sample_dir} and {pub_sample_dir}")
    
    # Track which samples have files for each section
    # We'll use this to only report warnings if files are truly missing after all checks
    files_by_sample = {sample: set() for sample in samples}
    
    # Track which samples are in MultiQC reports
    multiqc_samples = {section_name: set() for section_name in sections}
    
    # Process each section
    for section_name, section in sections.items():
        if not section["dir"]:  # Skip inner_distance if single-end
            continue
            
        input_dir = section["input"]
        if not input_dir or not os.path.exists(input_dir):
            msg = f"Input directory not found: {input_dir}"
            logger.warning(msg)
            # Don't log warnings yet, wait until after rseqc-logs check
            continue
            
        # Check files exist for each sample
        for sample in samples:
            sample_dir = section["dir"] / sample
            pub_sample_dir = section["pub_dir"] / sample
            found_files = False
            
            # Process sample files
            for pattern in section["file_patterns"]:
                # Pattern match for this sample's files
                logger.info(f"Looking for files matching pattern '{pattern}' for sample '{sample}' in '{input_dir}'")
                
                for file in os.listdir(input_dir):
                    logger.debug(f"Checking file: {file}")
                    if sample in file and (pattern in file):
                        logger.info(f"MATCH! File {file} matched sample {sample} and pattern {pattern}")
                        src = os.path.realpath(os.path.join(input_dir, file))
                        
                        # Create symbolic links in both directories
                        work_dst = os.path.join(sample_dir, file)
                        pub_dst = os.path.join(pub_sample_dir, file)
                        
                        try:
                            # Link to work directory for Nextflow
                            if not os.path.exists(work_dst) and not os.path.islink(work_dst):
                                os.symlink(src, work_dst)
                                logger.info(f"Linked {src} to {work_dst}")
                            
                            # Link to published directory structure
                            if not os.path.exists(pub_dst) and not os.path.islink(pub_dst):
                                os.symlink(src, pub_dst)
                                logger.info(f"Linked {src} to {pub_dst}")
                                
                            found_files = True
                        except Exception as e:
                            msg = f"Failed to link {file}: {str(e)}"
                            logger.error(msg)
                            val_logger.log("rseqc", sample, f"{section['display_name']}_link", "RED", msg)
            
            # If files were found, update tracking
            if found_files:
                files_by_sample[sample].add(section_name)
                val_logger.log("rseqc", sample, f"{section['display_name']}_outputs_found", "GREEN", 
                              f"Found RSeQC outputs for {section['display_name']} for sample {sample}")
    
    # Check for MultiQC reports in each section and verify all samples are present
    for section_name, section in sections.items():
        if not section["dir"]:  # Skip inner_distance if single-end
            continue
            
        # Check for MultiQC report
        found_multiqc = False
        multiqc_path = None
        
        if section["input"] and os.path.exists(section["input"]):
            for file in os.listdir(section["input"]):
                if "multiqc" in file.lower() and (assay_suffix in file or file.endswith('.zip')):
                    src = os.path.realpath(os.path.join(section["input"], file))
                    multiqc_path = src
                    
                    # Create symbolic links in both directories
                    work_dst = os.path.join(section["dir"], file)
                    pub_dst = os.path.join(section["pub_dir"], file)
                    
                    try:
                        # Link to work directory for Nextflow
                        if not os.path.exists(work_dst) and not os.path.islink(work_dst):
                            os.symlink(src, work_dst)
                            logger.info(f"Linked MultiQC report: {src} to {work_dst}")
                        
                        # Link to published directory structure
                        if not os.path.exists(pub_dst) and not os.path.islink(pub_dst):
                            os.symlink(src, pub_dst)
                            logger.info(f"Linked MultiQC report: {src} to {pub_dst}")
                            
                        found_multiqc = True
                    except Exception as e:
                        msg = f"Failed to link MultiQC report {file}: {str(e)}"
                        logger.error(msg)
                        val_logger.log("rseqc", "ALL", f"{section['display_name']}_multiqc", "RED", msg)
        
        if found_multiqc:
            val_logger.log("rseqc", "ALL", f"{section['display_name']}_multiqc", "GREEN", 
                        f"MultiQC report found for {section['display_name']}")
            
            # Now verify all samples are present in the MultiQC report
            if multiqc_path and multiqc_path.endswith('.zip'):
                try:
                    # Create a temporary directory to extract the ZIP contents
                    with tempfile.TemporaryDirectory() as temp_dir:
                        logger.info(f"Extracting MultiQC zip to temporary directory: {temp_dir}")
                        
                        # Extract the ZIP file
                        with zipfile.ZipFile(multiqc_path, 'r') as zip_ref:
                            zip_ref.extractall(temp_dir)
                        
                        # Find the data JSON file
                        json_file = None
                        for root, dirs, files in os.walk(temp_dir):
                            for file in files:
                                if file == 'multiqc_data.json':
                                    json_file = os.path.join(root, file)
                                    break
                            if json_file:
                                break
                        
                        if json_file:
                            logger.info(f"Found MultiQC data JSON: {json_file}")
                            with open(json_file, 'r') as f:
                                multiqc_data = json.load(f)
                            
                            # Dump the top-level keys to understand the structure
                            logger.debug(f"MultiQC JSON top-level keys: {list(multiqc_data.keys())}")
                            
                            # Check if 'report_saved_raw_data' exists and what it contains
                            if 'report_saved_raw_data' in multiqc_data:
                                logger.debug(f"report_saved_raw_data keys: {list(multiqc_data['report_saved_raw_data'].keys())}")
                                
                                # Look for RSeQC modules in the raw data
                                rseqc_modules = [k for k in multiqc_data['report_saved_raw_data'].keys() 
                                                if 'rseqc' in k.lower()]
                                logger.debug(f"Found RSeQC modules in raw data: {rseqc_modules}")
                                
                                # Extract sample names from report_saved_raw_data if available
                                found_samples = set()
                                
                                for module in rseqc_modules:
                                    module_data = multiqc_data['report_saved_raw_data'][module]
                                    if isinstance(module_data, dict):
                                        logger.debug(f"Module {module} data keys: {list(module_data.keys())}")
                                        found_samples.update(module_data.keys())
                                    
                                    # If module data is a list, examine its structure
                                    elif isinstance(module_data, list) and len(module_data) > 0:
                                        logger.debug(f"Module {module} data is a list with {len(module_data)} items")
                                        for item in module_data:
                                            if isinstance(item, dict) and 's_name' in item:
                                                found_samples.add(item['s_name'])
                            
                            # Check the plot data if available
                            if 'report_plot_data' in multiqc_data:
                                logger.debug(f"report_plot_data keys: {list(multiqc_data['report_plot_data'].keys())}")
                                
                                # Look for RSeQC plots
                                rseqc_plots = [k for k in multiqc_data['report_plot_data'].keys() 
                                              if 'rseqc' in k.lower() or section_name.replace('_', '') in k.lower()]
                                logger.debug(f"Found RSeQC plots: {rseqc_plots}")
                                
                                # Extract sample names from plots
                                for plot in rseqc_plots:
                                    plot_data = multiqc_data['report_plot_data'][plot]
                                    if isinstance(plot_data, dict) and 'datasets' in plot_data:
                                        for dataset in plot_data['datasets']:
                                            if 'name' in dataset:
                                                found_samples.add(dataset['name'])
                            
                            # As a fallback, search for any key that might contain sample names
                            # by looking for dataset dictionaries throughout the structure
                            def search_datasets(data, path=""):
                                samples = set()
                                if isinstance(data, dict):
                                    # Look for datasets key with lists of samples
                                    if 'datasets' in data and isinstance(data['datasets'], list):
                                        for dataset in data['datasets']:
                                            if isinstance(dataset, dict) and 'name' in dataset:
                                                samples.add(dataset['name'])
                                                logger.debug(f"Found sample in {path}/datasets: {dataset['name']}")
                                    
                                    # Look for sample names as keys
                                    for k, v in data.items():
                                        # Check if key looks like a sample name (has some of our expected samples as substrings)
                                        if isinstance(k, str) and any(sample in k for sample in samples):
                                            samples.add(k)
                                            logger.debug(f"Found sample as key in {path}: {k}")
                                        
                                        # Recursively search nested structures
                                        found = search_datasets(v, f"{path}/{k}")
                                        samples.update(found)
                                
                                elif isinstance(data, list):
                                    for i, item in enumerate(data):
                                        found = search_datasets(item, f"{path}[{i}]")
                                        samples.update(found)
                                
                                return samples
                            
                            # Search the entire structure for datasets
                            extra_samples = search_datasets(multiqc_data)
                            if extra_samples:
                                logger.debug(f"Found additional samples through structure search: {extra_samples}")
                                found_samples.update(extra_samples)
                            
                            # Look for samples in file names within the ZIP
                            for item in zip_ref.namelist():
                                filename = os.path.basename(item)
                                for sample in samples:
                                    if sample in filename:
                                        found_samples.add(sample)
                                        logger.debug(f"Found sample in ZIP filename: {sample} in {filename}")
                            
                            # Log the raw sample names to help with debugging
                            logger.debug(f"Raw sample names from all sources in MultiQC: {sorted(list(found_samples))}")
                            
                            # Initialize the set to collect clean sample names
                            clean_found_samples = set()
                            
                            # For each sample we're expecting
                            for sample in samples:
                                # Direct match
                                if sample in found_samples:
                                    clean_found_samples.add(sample)
                                    logger.debug(f"Direct match: sample {sample}")
                                    continue
                                
                                # Check for exact matches with RSeQC-specific extensions
                                # These are the exact formats used in the MultiQC reports
                                exact_extensions = [
                                    f"{sample}.read_dist",
                                    f"{sample}.infer_expt",
                                    f"{sample}.geneBodyCoverage",
                                    f"{sample}.inner_distance",
                                    f"{sample}_sorted"
                                ]
                                
                                for ext_sample in exact_extensions:
                                    if ext_sample in found_samples:
                                        clean_found_samples.add(sample)
                                        logger.debug(f"Extension match: found {ext_sample} for sample {sample}")
                                        break
                                
                                # If still not found, try more flexible matching
                                if sample not in clean_found_samples:
                                    # Look for samples that start with the sample name followed by a delimiter
                                    for found_sample in found_samples:
                                        if found_sample.startswith(f"{sample}.") or found_sample.startswith(f"{sample}_"):
                                            clean_found_samples.add(sample)
                                            logger.debug(f"Prefix match: found {found_sample} for sample {sample}")
                                            break
                                    
                                    # Last resort - check if sample name is contained within the found sample
                                    if sample not in clean_found_samples:
                                        for found_sample in found_samples:
                                            if sample in found_sample:
                                                clean_found_samples.add(sample)
                                                logger.debug(f"Substring match: found {found_sample} for sample {sample}")
                                                break
                            
                            logger.info(f"Samples found in MultiQC after cleanup: {clean_found_samples}")
                            
                            # Record which samples were found
                            multiqc_samples[section_name] = clean_found_samples
                            
                            # Find missing samples
                            missing_samples = set(samples) - clean_found_samples
                            
                            # Log results
                            if missing_samples:
                                msg = f"Samples missing from MultiQC report for {section['display_name']}: {', '.join(missing_samples)}"
                                logger.warning(msg)
                                val_logger.log("rseqc", "ALL", f"{section['display_name']}_samples_in_multiqc", "YELLOW", msg)
                            else:
                                val_logger.log("rseqc", "ALL", f"{section['display_name']}_samples_in_multiqc", "GREEN",
                                            f"All samples found in MultiQC report for {section['display_name']}")
                        else:
                            logger.warning(f"Could not find multiqc_data.json in the extracted ZIP file")
                except Exception as e:
                    logger.error(f"Error checking samples in MultiQC report {multiqc_path}: {str(e)}")
                    logger.exception(e)  # Log the full exception traceback
        elif section["required"]:
            msg = f"No MultiQC report found for {section['display_name']}"
            logger.warning(msg)
            val_logger.log("rseqc", "ALL", f"{section['display_name']}_multiqc", "YELLOW", msg)
    
    # Look for the rseqc-logs directory in multiple locations
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
            "geneBodyCoverage": "02_geneBody_coverage",
            "infer_expt.out": "03_infer_experiment", 
            "inner_distance": "04_inner_distance",
            "read_dist.out": "05_read_distribution"
        }
        
        # Define mapping from directory name to display name for logs
        display_name_map = {
            "02_geneBody_coverage": "geneBody_coverage.py",
            "03_infer_experiment": "infer_experiment.py",
            "04_inner_distance": "inner_distance.py",
            "05_read_distribution": "read_distribution.py"
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
            for pattern, section in file_type_patterns.items():
                if pattern in file:
                    target_section = section
                    logger.info(f"File {file} matched to section {section} with pattern {pattern}")
                    break
                    
            if not target_section:
                logger.warning(f"Could not match file {file} to any RSeQC section - patterns: {file_type_patterns.keys()}")
                continue
                
            # Skip inner_distance section for single-end data
            if target_section == "04_inner_distance" and not paired_end:
                continue
                
            # Create symbolic links for the file
            src = os.path.realpath(os.path.join(rseqc_logs_dir, file))
            
            # Get the working and published directory paths
            work_dst_dir = os.path.join(rseqc_dir, target_section, sample_match)
            pub_dst_dir = os.path.join(published_rseqc_dir, target_section, sample_match)
            
            work_dst = os.path.join(work_dst_dir, file)
            pub_dst = os.path.join(pub_dst_dir, file)
            
            try:
                # Ensure destination directories exist
                os.makedirs(work_dst_dir, exist_ok=True)
                os.makedirs(pub_dst_dir, exist_ok=True)
                
                # Create symbolic links in both directories
                if not os.path.exists(work_dst) and not os.path.islink(work_dst):
                    os.symlink(src, work_dst)
                    logger.info(f"Linked {src} to {work_dst}")
                
                if not os.path.exists(pub_dst) and not os.path.islink(pub_dst):
                    os.symlink(src, pub_dst) 
                    logger.info(f"Linked {src} to {pub_dst}")
                
                # Update tracking
                files_by_sample[sample_match].add(target_section)
                
                # Log successful file linkage from rseqc-logs using the display name
                display_name = display_name_map.get(target_section, target_section)
                val_logger.log("rseqc", sample_match, f"{display_name}_outputs_found", "GREEN",
                             f"Found RSeQC outputs for {display_name}")
                
            except Exception as e:
                msg = f"Failed to link {file} from rseqc-logs: {str(e)}"
                logger.error(msg)
                display_name = display_name_map.get(target_section, target_section)
                val_logger.log("rseqc", sample_match, f"{display_name}_outputs_link", "RED", msg)
    else:
        val_logger.log("rseqc", "ALL", "rseqc_outputs_dir", "YELLOW", 
                     f"RSeQC outputs directory not found in any of the checked locations")
    
    # Now after checking both direct inputs and rseqc-logs, log warnings for any truly missing files
    for sample in samples:
        for section_name, section in sections.items():
            if section["dir"] and section["required"]:  # Skip non-required sections
                # Only log warning if no files were found for this sample and section
                if section_name not in files_by_sample[sample]:
                    msg = f"No {section['display_name']} outputs found for sample {sample} in any location"
                    logger.warning(msg)
                    val_logger.log("rseqc", sample, f"{section['display_name']}_outputs", "YELLOW", msg)
    
    # Return validation status
    return {
        "status": val_logger.get_status(),
        "messages": [r["message"] for r in val_logger.results],
        "failures": {r["check_name"]: r["message"] 
                    for r in val_logger.results 
                    if r["status"] in ["HALT", "RED"]}
    }

if __name__ == '__main__':
    vv()

# Component based testing examples:

# Raw reads example:
# stage_files('rnaseq', 'raw_reads', 
#            raw_fastq='path/to/fastq',
#            raw_fastqc='path/to/fastqc',
#            raw_multiqc='path/to/multiqc')

# Trimmed reads example:
# stage_files('rnaseq', 'trimmed_reads',
#            fastq='path/to/trimmed_fastq',
#            fastqc='path/to/trimmed_fastqc',
#            trimmed_multiqc='path/to/trimmed_multiqc',
#            trimming_reports='path/to/trimming_reports',
#            trimming_multiqc='path/to/trimming_multiqc')
