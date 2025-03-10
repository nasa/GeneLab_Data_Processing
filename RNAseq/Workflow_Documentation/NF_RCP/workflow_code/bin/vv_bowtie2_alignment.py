#!/usr/bin/env python3
"""
Script to validate and verify raw reads based on runsheet information.
Expects to run from inside a output directory GLDS-##.

Parse the input runsheet to get:
Sample Name
Paired End
Has ERCC

Check that the expected output directories exist

Section-specific checks:

check_bowtie2_alignment_existence: Check if all expected bowtie2 alignment files exist for each sample

"""

import os
import sys
import argparse
import pandas as pd
import re
import subprocess
import gzip
import zipfile
import json
import tempfile
import numpy as np
from statistics import mean, median


def parse_runsheet(runsheet_path):
    """Parse the provided runsheet and extract relevant information."""
    if not os.path.exists(runsheet_path):
        print(f"Error: Runsheet not found at {runsheet_path}")
        sys.exit(1)
    
    try:
        # Try to read the runsheet using pandas
        df = pd.read_csv(runsheet_path)
        
        # Check for required columns
        required_columns = ['Sample Name', 'paired_end', 'has_ERCC', 'organism']
        missing_columns = [col for col in required_columns if col not in df.columns]
        
        if missing_columns:
            print(f"Error: Runsheet missing required columns: {', '.join(missing_columns)}")
            sys.exit(1)
            
        return df
    except Exception as e:
        print(f"Error parsing runsheet: {e}")
        sys.exit(1)

def check_directory_structure(outdir):
    """Check if the expected directory structure exists in the output directory."""
    raw_data_path = os.path.join(outdir, "00-RawData")
    fastq_path = os.path.join(raw_data_path, "Fastq")
    fastqc_path = os.path.join(raw_data_path, "FastQC_Reports")
    
    missing_dirs = []
    
    if not os.path.exists(raw_data_path):
        missing_dirs.append(str(raw_data_path))
    if not os.path.exists(fastq_path):
        missing_dirs.append(str(fastq_path))
    if not os.path.exists(fastqc_path):
        missing_dirs.append(str(fastqc_path))
    
    if missing_dirs:
        print(f"WARNING: The following directories are missing: {', '.join(missing_dirs)}")
        return False
    
    return True

def initialize_vv_log(outdir):
    """Initialize or append to the VV_log.csv file."""
    vv_log_path = os.path.join(outdir, "VV_log.csv")
    
    # Check if file exists
    if not os.path.exists(vv_log_path):
        # Create new file with header
        with open(vv_log_path, 'w') as f:
            f.write("component,sample_id,check_name,status,message,details\n")
    
    return vv_log_path

def log_check_result(log_path, component, sample_id, check_name, status, message="", details=""):
    """Log check result to the VV_log.csv file."""
    # Properly escape and format fields for CSV
    def escape_field(field, is_details=False):
        # Convert to string if not already
        field = str(field)
        
        # Replace newlines with semicolons to keep CSV valid
        field = field.replace('\n', '; ')
        
        # For details field, replace commas with semicolons to avoid CSV quoting
        if is_details and ',' in field:
            field = field.replace(', ', '; ')
        
        # Also replace commas in message fields that contain percentages
        if '%' in field and ',' in field:
            field = field.replace(', ', '; ')
        
        # If the field contains commas or quotes (but not just semicolons), wrap it in quotes
        if ',' in field or '"' in field:
            # Double any quotes within the field
            field = field.replace('"', '""')
            # Wrap in quotes
            field = f'"{field}"'
        return field
    
    # Format all fields
    component = escape_field(component)
    sample_id = escape_field(sample_id)
    check_name = escape_field(check_name)
    status = escape_field(status)
    message = escape_field(message)
    details = escape_field(details, True)  # Specify this is a details field
    
    # Write the formatted line
    with open(log_path, 'a') as f:
        f.write(f"{component},{sample_id},{check_name},{status},{message},{details}\n")

def check_gzip_integrity(outdir, samples, paired_end, log_path):
    """Check GZIP integrity for all FASTQ files using gzip -t."""
    fastq_dir = os.path.join(outdir, "00-RawData", "Fastq")
    failed_files = []
    all_files = []
    
    is_paired = paired_end[0]  # Assuming all samples have the same paired_end value
    
    for sample in samples:
        if is_paired:
            # Get R1 and R2 files for paired-end sequencing
            r1_file = os.path.join(fastq_dir, f"{sample}_R1_raw.fastq.gz")
            r2_file = os.path.join(fastq_dir, f"{sample}_R2_raw.fastq.gz")
            files_to_check = [r1_file, r2_file]
        else:
            # Get single file for single-end sequencing
            single_file = os.path.join(fastq_dir, f"{sample}_raw.fastq.gz")
            files_to_check = [single_file]
        
        for file_path in files_to_check:
            if os.path.exists(file_path):
                all_files.append(file_path)
                # Run gzip -t to test integrity
                result = subprocess.run(["gzip", "-t", file_path], capture_output=True)
                if result.returncode != 0:
                    failed_files.append(file_path)
                    print(f"GZIP integrity check failed for: {file_path}")
                    print(f"Error: {result.stderr.decode('utf-8')}")
    
    if failed_files:
        log_check_result(log_path, "raw_reads", "all", "check_gzip_integrity", "RED", 
                         f"Integrity check failed for {len(failed_files)} of {len(all_files)} files", 
                         ",".join(failed_files))
        return False
    
    if all_files:
        print(f"GZIP integrity check passed for all {len(all_files)} FASTQ files")
        log_check_result(log_path, "raw_reads", "all", "check_gzip_integrity", "GREEN", 
                         f"Integrity check passed for all {len(all_files)} files", "")
        return True
    else:
        print("No FASTQ files found to check for GZIP integrity")
        log_check_result(log_path, "raw_reads", "all", "check_gzip_integrity", "YELLOW", 
                         "No FASTQ files found to check", "")
        return False

def validate_fastq_format(outdir, samples, paired_end, log_path, max_lines=200000000):
    """Validate FASTQ format by checking that header lines start with @.
    Only checks the first 200 million lines for performance reasons."""
    fastq_dir = os.path.join(outdir, "00-RawData", "Fastq")
    invalid_files = []
    all_files = []
    
    is_paired = paired_end[0]  # Assuming all samples have the same paired_end value
    
    for sample in samples:
        if is_paired:
            # Get R1 and R2 files for paired-end sequencing
            r1_file = os.path.join(fastq_dir, f"{sample}_R1_raw.fastq.gz")
            r2_file = os.path.join(fastq_dir, f"{sample}_R2_raw.fastq.gz")
            files_to_check = [r1_file, r2_file]
        else:
            # Get single file for single-end sequencing
            single_file = os.path.join(fastq_dir, f"{sample}_raw.fastq.gz")
            files_to_check = [single_file]
        
        for file_path in files_to_check:
            if os.path.exists(file_path):
                all_files.append(file_path)
                print(f"Validating FASTQ format for: {file_path}")
                
                try:
                    with gzip.open(file_path, 'rt') as f:
                        line_count = 0
                        line_in_record = 0
                        
                        for line in f:
                            line_count += 1
                            line_in_record = (line_count - 1) % 4
                            
                            # Check if header line (every 4th line, starting with line 1)
                            if line_in_record == 0:
                                if not line.startswith('@'):
                                    invalid_files.append(file_path)
                                    print(f"Invalid FASTQ format in {file_path} at line {line_count}: Header line does not start with @")
                                    break
                            
                            # Stop after max_lines
                            if line_count >= max_lines:
                                print(f"Reached {max_lines} lines limit for {file_path}, stopping validation")
                                break
                                
                except Exception as e:
                    invalid_files.append(file_path)
                    print(f"Error validating {file_path}: {str(e)}")
    
    if invalid_files:
        log_check_result(log_path, "raw_reads", "all", "validate_fastq_format", "RED", 
                         f"Invalid FASTQ format in {len(invalid_files)} of {len(all_files)} files", 
                         ",".join(invalid_files))
        return False
    
    if all_files:
        print(f"FASTQ format validation passed for all {len(all_files)} files")
        log_check_result(log_path, "raw_reads", "all", "validate_fastq_format", "GREEN", 
                         f"Valid FASTQ format in all {len(all_files)} files", "")
        return True
    else:
        print("No FASTQ files found to validate format")
        log_check_result(log_path, "raw_reads", "all", "validate_fastq_format", "YELLOW", 
                         "No FASTQ files found to validate", "")
        return False

def check_raw_fastqc_existence(outdir, samples, paired_end, log_path):
    """Check if FastQC output files exist for each sample."""
    fastqc_dir = os.path.join(outdir, "00-RawData", "FastQC_Reports")
    missing_files = []
    
    is_paired = paired_end[0]  # Assuming all samples have the same paired_end value
    
    for sample in samples:
        if is_paired:
            # Check for R1 and R2 FastQC files for paired-end sequencing
            r1_html = os.path.join(fastqc_dir, f"{sample}_R1_raw_fastqc.html")
            r1_zip = os.path.join(fastqc_dir, f"{sample}_R1_raw_fastqc.zip")
            r2_html = os.path.join(fastqc_dir, f"{sample}_R2_raw_fastqc.html")
            r2_zip = os.path.join(fastqc_dir, f"{sample}_R2_raw_fastqc.zip")
            
            files_to_check = [r1_html, r1_zip, r2_html, r2_zip]
        else:
            # Check for single FastQC files for single-end sequencing
            html_file = os.path.join(fastqc_dir, f"{sample}_raw_fastqc.html")
            zip_file = os.path.join(fastqc_dir, f"{sample}_raw_fastqc.zip")
            
            files_to_check = [html_file, zip_file]
        
        for file_path in files_to_check:
            if not os.path.exists(file_path):
                missing_files.append(file_path)
    
    if missing_files:
        print(f"WARNING: The following FastQC output files are missing:")
        for file in missing_files:
            print(f"  - {file}")
        log_check_result(log_path, "raw_reads", "all", "check_raw_fastqc_existence", "RED", 
                         f"Missing {len(missing_files)} FastQC output files", ",".join(missing_files))
        return False
    
    print(f"All expected FastQC output files found for {len(samples)} samples")
    log_check_result(log_path, "raw_reads", "all", "check_raw_fastqc_existence", "GREEN", 
                     f"All FastQC output files found", "")
    return True

def check_samples_multiqc(outdir, samples, paired_end, log_path, assay_suffix="_GLbulkRNAseq"):
    """Check if all samples are present in the MultiQC report."""
    align_dir = os.path.join(outdir, "02-Bowtie2_Alignment")
    multiqc_zip = os.path.join(align_dir, f"align_multiqc{assay_suffix}_report.zip")
    
    if not os.path.exists(multiqc_zip):
        print(f"WARNING: MultiQC report not found at: {multiqc_zip}")
        log_check_result(log_path, "alignment", "all", "check_samples_multiqc", "RED", 
                        "MultiQC report not found", multiqc_zip)
        return False
    
    print(f"Checking samples in MultiQC report: {multiqc_zip}")
    
    # Create a temporary directory to extract files
    with tempfile.TemporaryDirectory() as temp_dir:
        try:
            # Extract the zip file
            with zipfile.ZipFile(multiqc_zip, 'r') as zip_ref:
                zip_ref.extractall(temp_dir)
            
            # Search for the MultiQC data JSON file
            json_files = []
            for root, dirs, files in os.walk(temp_dir):
                for file in files:
                    if file.endswith('.json') and 'multiqc_data' in file:
                        json_files.append(os.path.join(root, file))
            
            if not json_files:
                print(f"WARNING: No multiqc_data.json file found in the extracted zip")
                log_check_result(log_path, "alignment", "all", "check_samples_multiqc", "RED", 
                               "multiqc_data.json not found in zip", "")
                return False
            
            multiqc_data_file = json_files[0]
            
            # Parse the JSON file
            with open(multiqc_data_file) as f:
                multiqc_data = json.load(f)
            
            # Check for samples in report_data_sources
            if 'report_data_sources' not in multiqc_data:
                print(f"WARNING: No report_data_sources found in MultiQC data")
                log_check_result(log_path, "alignment", "all", "check_samples_multiqc", "RED", 
                               "No report_data_sources in MultiQC data", "")
                return False
            
            # Look for Bowtie2 section
            if 'Bowtie 2 / HiSAT2' not in multiqc_data['report_data_sources']:
                print(f"WARNING: No Bowtie2 data found in MultiQC report")
                log_check_result(log_path, "alignment", "all", "check_samples_multiqc", "RED", 
                               "No Bowtie2 data in MultiQC report", "")
                return False
            
            # Get samples from the Bowtie2 section
            bowtie_data = multiqc_data['report_data_sources']['Bowtie 2 / HiSAT2']['all_sections']
            multiqc_samples = set(bowtie_data.keys())
            
            # Check for missing samples
            missing_samples = []
            for sample in samples:
                if sample not in multiqc_samples:
                    missing_samples.append(sample)
            
            if missing_samples:
                print(f"WARNING: The following samples are missing from the MultiQC report:")
                for sample in missing_samples:
                    print(f"  - {sample}")
                log_check_result(log_path, "alignment", "all", "check_samples_multiqc", "RED", 
                               f"Missing {len(missing_samples)} samples in MultiQC report", 
                               "; ".join(missing_samples))
                return False
            
            print(f"All {len(samples)} samples found in MultiQC report")
            log_check_result(log_path, "alignment", "all", "check_samples_multiqc", "GREEN", 
                           f"All {len(samples)} samples found in report", 
                           f"Checked presence of {len(samples)} samples in Bowtie2 section of MultiQC report")
            return True
            
        except Exception as e:
            print(f"Error checking MultiQC report: {str(e)}")
            log_check_result(log_path, "alignment", "all", "check_samples_multiqc", "RED", 
                           f"Error checking MultiQC report: {str(e)}", "")
            return False

def get_bowtie2_multiqc_stats(outdir, samples, paired_end, log_path, assay_suffix="_GLbulkRNAseq"):
    """Get alignment statistics from MultiQC report."""
    align_dir = os.path.join(outdir, "02-Bowtie2_Alignment")
    multiqc_zip = os.path.join(align_dir, f"align_multiqc{assay_suffix}_report.zip")
    
    if not os.path.exists(multiqc_zip):
        print(f"WARNING: MultiQC report not found at: {multiqc_zip}")
        log_check_result(log_path, "alignment", "all", "get_bowtie2_multiqc_stats", "RED", 
                        "MultiQC report not found", multiqc_zip)
        return None
    
    # Create a temporary directory to extract files
    with tempfile.TemporaryDirectory() as temp_dir:
        try:
            # Extract the zip file
            with zipfile.ZipFile(multiqc_zip, 'r') as zip_ref:
                zip_ref.extractall(temp_dir)
            
            # Search for the MultiQC data JSON file
            json_files = []
            for root, dirs, files in os.walk(temp_dir):
                for file in files:
                    if file.endswith('.json') and 'multiqc_data' in file:
                        json_files.append(os.path.join(root, file))
            
            if not json_files:
                print(f"WARNING: No multiqc_data.json file found in the extracted zip")
                log_check_result(log_path, "alignment", "all", "get_bowtie2_multiqc_stats", "RED", 
                               "multiqc_data.json not found in zip", "")
                return None
            
            multiqc_data_file = json_files[0]
            
            # Parse the JSON file
            with open(multiqc_data_file) as f:
                data = json.load(f)
            
            samples_data = {}
            
            # Extract stats for each sample from general stats
            for section in data['report_general_stats_data']:
                for sample, stats in section.items():
                    # Clean sample name (remove read identifiers)
                    base_sample = re.sub(r'_R[12]$', '', sample)
                    
                    if base_sample not in samples_data:
                        samples_data[base_sample] = {}
                    
                    # Get the key metrics
                    if 'total_reads' in stats:
                        samples_data[base_sample]['total_reads'] = stats['total_reads']
                    
                    if 'overall_alignment_rate' in stats:
                        samples_data[base_sample]['overall_alignment_rate'] = stats['overall_alignment_rate']
                    
                    # Handle aligned reads stats for both paired and unpaired data
                    if 'paired_aligned_none' in stats:
                        samples_data[base_sample]['aligned_none'] = stats['paired_aligned_none']
                        samples_data[base_sample]['aligned_one'] = stats['paired_aligned_one']
                        samples_data[base_sample]['aligned_multi'] = stats['paired_aligned_multi']
                    elif 'unpaired_aligned_none' in stats:
                        samples_data[base_sample]['aligned_none'] = stats['unpaired_aligned_none']
                        samples_data[base_sample]['aligned_one'] = stats['unpaired_aligned_one']
                        samples_data[base_sample]['aligned_multi'] = stats['unpaired_aligned_multi']
            
            if not samples_data:
                print("No alignment statistics found in MultiQC report")
                log_check_result(log_path, "alignment", "all", "get_bowtie2_multiqc_stats", "RED", 
                               "No alignment statistics found", "")
                return None
                
            print(f"Retrieved alignment statistics for {len(samples_data)} samples")
            return samples_data
            
        except Exception as e:
            print(f"Error parsing MultiQC data: {str(e)}")
            log_check_result(log_path, "alignment", "all", "get_bowtie2_multiqc_stats", "RED", 
                           f"Error parsing MultiQC data: {str(e)}", "")
            return None

def report_multiqc_outliers(outdir, multiqc_data, log_path):
    """Report statistical outliers in MultiQC statistics using standard deviations."""
    metrics = ['total_reads', 'overall_alignment_rate', 'aligned_none', 'aligned_one', 'aligned_multi']
    outlier_samples = set()
    
    # First pass: collect values and convert to percentages
    metric_values = {metric: [] for metric in metrics}
    sample_values = {metric: {} for metric in metrics}
    
    for sample, stats in multiqc_data.items():
        if 'total_reads' not in stats:
            continue
            
        total = stats['total_reads']
        
        # Convert aligned counts to percentages
        if 'aligned_none' in stats:
            stats['aligned_none'] = (stats['aligned_none'] / total) * 100
        if 'aligned_one' in stats:
            stats['aligned_one'] = (stats['aligned_one'] / total) * 100
        if 'aligned_multi' in stats:
            stats['aligned_multi'] = (stats['aligned_multi'] / total) * 100
        
        # Collect values for each metric
        for metric in metrics:
            if metric in stats:
                metric_values[metric].append(stats[metric])
                sample_values[metric][sample] = stats[metric]
    
    # Calculate statistics and check for outliers
    for metric in metrics:
        if not metric_values[metric]:
            continue
            
        values = np.array(metric_values[metric])
        median = np.median(values)
        std = np.std(values)
        
        # Check each sample for this metric
        for sample, value in sample_values[metric].items():
            if std > 0:  # Avoid division by zero
                stdevs_from_median = abs(value - median) / std
                
                # Flag if more than 2 standard deviations from median
                if stdevs_from_median > 2:
                    outlier_samples.add(sample)
                    
                    # Format metric name for logging
                    metric_name = metric.replace('_', ' ')
                    
                    # Add detailed information for outliers
                    outlier_details = f"Median={median:.2f}, StdDev={std:.2f}, Threshold=2.0 standard deviations"
                    
                    # Log individual outlier
                    message = f"{sample}: {metric_name} ({value:.2f}) is {stdevs_from_median:.2f} standard deviations from the median ({median:.2f})"
                    log_check_result(log_path, "alignment", sample, f"outlier_{metric}", "YELLOW", message, outlier_details)
    
    # Log summary
    if outlier_samples:
        print(f"\nFound {len(outlier_samples)} samples with outlier metrics")
        log_check_result(log_path, "alignment", "all", "report_multiqc_outliers", "YELLOW",
                        f"Found {len(outlier_samples)} samples with outlier metrics",
                        "; ".join(sorted(outlier_samples)))
    else:
        print("\nNo statistical outliers found in alignment metrics")
        log_check_result(log_path, "alignment", "all", "report_multiqc_outliers", "GREEN",
                        "No outliers found", "All samples within 2 standard deviations of median values")
    
    # Print summary statistics
    print("\nAlignment Statistics Summary:")
    for metric in metrics:
        if metric_values[metric]:
            values = np.array(metric_values[metric])
            print(f"\n{metric}:")
            print(f"  Mean: {np.mean(values):.2f}")
            print(f"  Median: {np.median(values):.2f}")
            print(f"  Std Dev: {np.std(values):.2f}")
    
    return len(outlier_samples) == 0

def check_bowtie2_existence(outdir, samples, paired_end, log_path):
    """Check if all expected bowtie2 alignment files exist for each sample."""
    align_dir = os.path.join(outdir, "02-Bowtie2_Alignment")
    missing_files = []
    is_paired = paired_end[0]  # Assuming all samples have same paired_end value

    # First check for MultiQC report
    multiqc_file = os.path.join(align_dir, "align_multiqc_GLbulkRNAseq_report.zip")
    if not os.path.exists(multiqc_file):
        missing_files.append(multiqc_file)

    for sample in samples:
        # Check sample directory
        sample_dir = os.path.join(align_dir, sample)
        if not os.path.exists(sample_dir):
            missing_files.append(sample_dir)
            continue

        # Required files for both SE and PE
        required_files = [
            f"{sample}.bowtie2.log",
            f"{sample}_sorted.bam",
            f"{sample}_sorted.bam.bai"
        ]

        # Add unmapped files based on paired/single end
        if is_paired:
            required_files.extend([
                f"{sample}.unmapped.fastq.1.gz",
                f"{sample}.unmapped.fastq.2.gz"
            ])
        else:
            required_files.append(f"{sample}.unmapped.fastq.gz")

        # Check each required file
        for file_name in required_files:
            file_path = os.path.join(sample_dir, file_name)
            if not os.path.exists(file_path):
                missing_files.append(file_path)

    if missing_files:
        print(f"WARNING: Missing alignment files:")
        for file in missing_files:
            print(f"  - {file}")
        log_check_result(log_path, "alignment", "all", "check_bowtie2_existence", "RED", 
                        f"Missing {len(missing_files)} files", "; ".join(missing_files))
        return False

    print(f"All expected alignment files found")
    total_files = len(samples) * (4 if is_paired else 3) + 1  # +1 for multiqc report
    log_check_result(log_path, "alignment", "all", "check_bowtie2_existence", "GREEN", 
                    "All files found", f"Verified {total_files} files for {len(samples)} samples")
    return True

def validate_unmapped_fastq(outdir, samples, paired_end, log_path, max_lines=200000000):
    """Validate unmapped FASTQ files from bowtie2 alignment."""
    align_dir = os.path.join(outdir, "02-Bowtie2_Alignment")
    invalid_files = []
    all_files = []
    is_paired = paired_end[0]

    for sample in samples:
        sample_dir = os.path.join(align_dir, sample)
        
        if is_paired:
            files_to_check = [
                os.path.join(sample_dir, f"{sample}.unmapped.fastq.1.gz"),
                os.path.join(sample_dir, f"{sample}.unmapped.fastq.2.gz")
            ]
        else:
            files_to_check = [
                os.path.join(sample_dir, f"{sample}.unmapped.fastq.gz")
            ]
        
        for file_path in files_to_check:
            if os.path.exists(file_path):
                all_files.append(file_path)
                print(f"Validating unmapped FASTQ: {file_path}")
                
                try:
                    # Check GZIP integrity
                    result = subprocess.run(["gzip", "-t", file_path], capture_output=True)
                    if result.returncode != 0:
                        invalid_files.append(file_path)
                        print(f"GZIP integrity check failed for: {file_path}")
                        print(f"Error: {result.stderr.decode('utf-8')}")
                        continue

                    # Check FASTQ format and ensure file isn't empty
                    with gzip.open(file_path, 'rt') as f:
                        line_count = 0
                        line_in_record = 0
                        has_valid_record = False
                        
                        for line in f:
                            line_count += 1
                            line_in_record = (line_count - 1) % 4
                            
                            # Check if header line (every 4th line, starting with line 1)
                            if line_in_record == 0:
                                if not line.startswith('@'):
                                    invalid_files.append(file_path)
                                    print(f"Invalid FASTQ format in {file_path} at line {line_count}")
                                    break
                                has_valid_record = True
                            
                            # Stop after max_lines
                            if line_count >= max_lines:
                                print(f"Reached {max_lines} lines limit for {file_path}")
                                break
                        
                        # Check if file is empty or has incomplete records
                        if not has_valid_record:
                            invalid_files.append(file_path)
                            print(f"Empty or invalid FASTQ file: {file_path}")
                            continue
                        
                        if line_count % 4 != 0:
                            invalid_files.append(file_path)
                            print(f"Incomplete FASTQ record in {file_path}")
                                
                except Exception as e:
                    invalid_files.append(file_path)
                    print(f"Error validating {file_path}: {str(e)}")

    if invalid_files:
        print(f"WARNING: The following unmapped FASTQ files are invalid:")
        for file in invalid_files:
            print(f"  - {file}")
        log_check_result(log_path, "alignment", "all", "validate_unmapped_fastq", "RED", 
                        f"Invalid format in {len(invalid_files)} of {len(all_files)} files", 
                        "; ".join(invalid_files))
        return False

    print(f"All unmapped FASTQ files are valid")
    log_check_result(log_path, "alignment", "all", "validate_unmapped_fastq", "GREEN", 
                    "All files valid", f"Validated {len(all_files)} unmapped FASTQ files; Checked for valid FASTQ format and complete records")
    return True

def main():
    """Main function to process runsheet and validate raw reads."""
    parser = argparse.ArgumentParser(description='Validate raw reads based on runsheet information.')
    parser.add_argument('--runsheet', '-r', required=True, help='Path to the runsheet CSV file')
    parser.add_argument('--outdir', '-o', default=os.getcwd(), 
                        help='Output directory (GLDS-## folder), defaults to current directory')
    parser.add_argument('--assay-suffix', default="_GLbulkRNAseq", 
                        help='Assay suffix used in MultiQC report filenames (default: _GLbulkRNAseq)')
    args = parser.parse_args()

    # Initialize VV log
    vv_log_path = initialize_vv_log(args.outdir)
    
    # Check directory structure
    check_directory_structure(args.outdir)
    
    # Parse the runsheet
    runsheet_df = parse_runsheet(args.runsheet)
    
    # Extract sample names
    sample_names = runsheet_df['Sample Name'].tolist()
    
    # Check consistency of paired_end, has_ERCC, and organism values
    paired_end_values = runsheet_df['paired_end'].unique()
    has_ercc_values = runsheet_df['has_ERCC'].unique()
    organism_values = runsheet_df['organism'].unique()
    
    # Only print warnings for inconsistencies
    if len(paired_end_values) > 1:
        print(f"WARNING: Inconsistent paired_end values: {paired_end_values}")
    
    if len(has_ercc_values) > 1:
        print(f"WARNING: Inconsistent has_ERCC values: {has_ercc_values}")
        
    if len(organism_values) > 1:
        print(f"WARNING: Inconsistent organism values: {organism_values}")
    print(f'Got the following values:')
    print(f'paired_end: {paired_end_values}')
    print(f'has_ercc: {has_ercc_values}')
    print(f'organism: {organism_values}')

    # Initiate validation checks in logical order
    
    # 1. Check for bowtie2 alignment files existence
    check_bowtie2_existence(args.outdir, sample_names, paired_end_values, vv_log_path)

    # 2. Validate unmapped FASTQ files
    validate_unmapped_fastq(args.outdir, sample_names, paired_end_values, vv_log_path)

    # 3. Check that samples are present in MultiQC report
    check_samples_multiqc(args.outdir, sample_names, paired_end_values, vv_log_path, args.assay_suffix)

    # 4. Get MultiQC stats
    multiqc_data = get_bowtie2_multiqc_stats(args.outdir, sample_names, paired_end_values, vv_log_path, args.assay_suffix)

    # 5. Report MultiQC stats outliers
    if multiqc_data:
        report_multiqc_outliers(args.outdir, multiqc_data, vv_log_path)

if __name__ == "__main__":
    main()