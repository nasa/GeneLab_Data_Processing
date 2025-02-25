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
def vv(assay_type, assay_suffix, runsheet_path, outdir, paired_end, mode, run_components, 
       raw_fastq, raw_fastqc, raw_multiqc,
       trimmed_fastq, trimmed_fastqc, trimmed_multiqc, trimming_reports, trimming_multiqc):
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
