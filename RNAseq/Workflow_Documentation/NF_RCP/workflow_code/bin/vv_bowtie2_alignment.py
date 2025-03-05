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
    """Check if all samples are included in the MultiQC report."""
    fastqc_dir = os.path.join(outdir, "00-RawData", "FastQC_Reports")
    multiqc_zip = os.path.join(fastqc_dir, f"raw_multiqc{assay_suffix}_report.zip")
    
    if not os.path.exists(multiqc_zip):
        print(f"WARNING: MultiQC report zip file not found: {multiqc_zip}")
        log_check_result(log_path, "raw_reads", "all", "check_samples_multiqc", "RED", 
                         "MultiQC report not found", multiqc_zip)
        return False
    
    print(f"Found MultiQC report: {multiqc_zip}")
    
    # Create a temporary directory to extract files
    with tempfile.TemporaryDirectory() as temp_dir:
        try:
            # Extract the zip file
            with zipfile.ZipFile(multiqc_zip, 'r') as zip_ref:
                zip_ref.extractall(temp_dir)
            
            # Check for the JSON file
            json_path = os.path.join(temp_dir, f"raw_multiqc{assay_suffix}_report", 
                                    f"raw_multiqc{assay_suffix}_data", "multiqc_data.json")
            
            if not os.path.exists(json_path):
                print(f"Could not find multiqc_data.json in the expected location")
                log_check_result(log_path, "raw_reads", "all", "check_samples_multiqc", "RED", 
                               "multiqc_data.json not found in zip", "")
                return False
                
            # Parse the MultiQC JSON
            with open(json_path, 'r') as f:
                multiqc_data = json.load(f)
            
            # Extract sample names from FastQC data
            mqc_samples = []
            is_paired = paired_end[0]  # Assuming all samples have the same paired_end value
            
            if ('report_data_sources' in multiqc_data and 
                'FastQC' in multiqc_data['report_data_sources'] and
                'all_sections' in multiqc_data['report_data_sources']['FastQC']):
                
                fastqc_sections = multiqc_data['report_data_sources']['FastQC']['all_sections']
                for mqc_sample in fastqc_sections.keys():
                    # For paired-end data, remove _R1 and _R2 suffixes
                    base_sample = mqc_sample.replace("_raw_fastqc", "").replace("_fastqc", "")
                    if is_paired:
                        base_sample = base_sample.replace("_R1", "").replace("_R2", "")
                    mqc_samples.append(base_sample)
            
            # Remove duplicates
            mqc_samples = list(set(mqc_samples))
            
            # Check if all runsheet samples are present in MultiQC report
            missing_samples = []
            for sample in samples:
                if sample not in mqc_samples:
                    missing_samples.append(sample)
            
            if missing_samples:
                print(f"WARNING: The following samples are missing from the MultiQC report:")
                for sample in missing_samples:
                    print(f"  - {sample}")
                log_check_result(log_path, "raw_reads", "all", "check_samples_multiqc", "RED", 
                                f"Missing {len(missing_samples)} samples in MultiQC report", 
                                ",".join(missing_samples))
                return False
            
            print(f"All {len(samples)} samples found in the MultiQC report")
            log_check_result(log_path, "raw_reads", "all", "check_samples_multiqc", "GREEN", 
                            "All samples found in MultiQC report", "")
            return True
            
        except Exception as e:
            print(f"Error processing MultiQC report: {str(e)}")
            log_check_result(log_path, "raw_reads", "all", "check_samples_multiqc", "RED", 
                           f"Error processing MultiQC report: {str(e)}", "")
            return False

def report_multiqc_outliers(outdir, multiqc_data, log_path):
    """Identify and report outliers in MultiQC statistics."""
    if not multiqc_data:
        print("No MultiQC data to analyze for outliers")
        return False

    # Keys to check for outliers as specified in protocol.py
    metrics_to_check = {
        "raw_percent_gc_f": "percent GC content (forward)",
        "raw_avg_sequence_length_f": "average sequence length (forward)",
        "raw_total_sequences_f": "total sequences (forward)",
        "raw_percent_duplicates_f": "percent duplicates (forward)",
    }
    
    # For paired-end data, also check the reverse reads
    if any("_r" in key for sample in multiqc_data.values() for key in sample):
        metrics_to_check.update({
            "raw_percent_gc_r": "percent GC content (reverse)",
            "raw_avg_sequence_length_r": "average sequence length (reverse)",
            "raw_total_sequences_r": "total sequences (reverse)",
            "raw_percent_duplicates_r": "percent duplicates (reverse)",
        })
    
    # Thresholds for outlier detection
    thresholds = [
        {"code": "YELLOW", "stdev_threshold": 2, "middle_fcn": "median"},
        {"code": "RED", "stdev_threshold": 4, "middle_fcn": "median"},
    ]
    
    # Set to keep track of outlier samples
    outlier_samples = set()
    outlier_details = []
    
    # Check each metric for outliers
    for metric_key, metric_description in metrics_to_check.items():
        # Skip if metric doesn't exist in any sample
        if not any(metric_key in sample_data for sample_data in multiqc_data.values()):
            continue
            
        # Collect values for this metric across all samples
        values = []
        for sample, sample_data in multiqc_data.items():
            if metric_key in sample_data:
                try:
                    values.append((sample, float(sample_data[metric_key])))
                except (ValueError, TypeError):
                    print(f"Warning: Non-numeric value for {metric_key} in sample {sample}")
        
        if not values:
            continue
            
        # Calculate median and standard deviation
        all_values = [v[1] for v in values]
        median_value = np.median(all_values)
        stdev_value = np.std(all_values)
        
        if stdev_value == 0:
            # Skip if there's no variation
            continue
            
        # Check each sample against thresholds
        for sample, value in values:
            # Calculate deviation in terms of standard deviations
            deviation = abs(value - median_value) / stdev_value
            
            # Check against thresholds (largest threshold first)
            for threshold in sorted(thresholds, key=lambda x: x["stdev_threshold"], reverse=True):
                if deviation >= threshold["stdev_threshold"]:
                    outlier_samples.add(sample)
                    outlier_message = (
                        f"{sample}: {metric_description} ({value}) is {deviation:.2f} "
                        f"standard deviations from the median ({median_value:.2f})"
                    )
                    outlier_details.append((sample, metric_key, threshold["code"], outlier_message))
                    break  # Stop after finding the first matching threshold
    
    # Report results
    if outlier_details:
        print(f"Found {len(outlier_samples)} samples with outlier metrics:")
        for sample, metric, code, message in outlier_details:
            print(f"  - {message} [{code}]")
            log_check_result(log_path, "raw_reads", sample, f"outlier_{metric}", code, message, "")
        
        # Keep the list of outlier samples in the summary log
        log_check_result(log_path, "raw_reads", "all", "check_for_outliers", "YELLOW", 
                        f"Found {len(outlier_samples)} samples with outlier metrics", 
                        ",".join(outlier_samples))  # Include outlier samples list
    else:
        print("No outliers found in MultiQC metrics")
        log_check_result(log_path, "raw_reads", "all", "check_for_outliers", "GREEN", 
                        "No outliers found", "")
    
    return len(outlier_details) > 0

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
                        f"Missing {len(missing_files)} files", ",".join(missing_files))
        return False

    print(f"All expected alignment files found")
    log_check_result(log_path, "alignment", "all", "check_bowtie2_existence", "GREEN", 
                    "All files found", "")
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
                        ",".join(invalid_files))
        return False

    print(f"All unmapped FASTQ files are valid")
    log_check_result(log_path, "alignment", "all", "validate_unmapped_fastq", "GREEN", 
                    "All files valid", "")
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

if __name__ == "__main__":
    main()