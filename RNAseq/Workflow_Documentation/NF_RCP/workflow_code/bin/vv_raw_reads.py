#!/usr/bin/env python
"""
Script to validate and verify raw reads based on runsheet information.
Expects to run from inside a output directory GLDS-##.

Parse the input runsheet to get:
Sample Name
Paired End
Has ERCC

Check that the expected output directories exist

Section-specific checks:

check_raw_fastq_existence: Check if all expected raw FASTQ files exist for each sample
check_gzip_integrity: Verify GZIP integrity of FASTQ files using gzip -t
validate_fastq_format: Ensure FASTQ files have proper formatting (headers start with @)
check_raw_fastqc_existence: Verify FastQC output files exist (HTML and ZIP files)
check_samples_multiqc: Confirm all samples are included in the MultiQC report
get_raw_multiqc_stats: Extract FastQC metrics from MultiQC report
report_multiqc_outliers: Identify and report outliers in FastQC metrics
check_paired_read_counts: Ensure paired-end reads have matching read counts
report_read_depth_stats: Log statistics about sequencing depth across all samples
report_duplication_rate_stats: Log statistics about duplication rates across all samples

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
    multiqc_path = os.path.join(raw_data_path, "MultiQC_Reports")
    
    missing_dirs = []
    
    if not os.path.exists(raw_data_path):
        missing_dirs.append(str(raw_data_path))
    if not os.path.exists(fastq_path):
        missing_dirs.append(str(fastq_path))
    if not os.path.exists(fastqc_path):
        missing_dirs.append(str(fastqc_path))
    if not os.path.exists(multiqc_path):
        missing_dirs.append(str(multiqc_path))
    
    if missing_dirs:
        print(f"WARNING: The following directories are missing: {', '.join(missing_dirs)}")
        return False
    
    return True

def initialize_vv_log(outdir):
    """Initialize or append to the VV_log.csv file."""
    vv_log_path = os.path.join(outdir, "VV_log.csv")
    
    # Check if file exists
    if not os.path.exists(vv_log_path):
        # Create new file with header - both status (color) and flag_code (number)
        with open(vv_log_path, 'w') as f:
            f.write("component,sample_id,check_name,status,flag_code,message,details\n")
    
    return vv_log_path

def log_check_result(log_path, component, sample_id, check_name, status, message="", details=""):
    """Log check result to the VV_log.csv file."""
    def escape_field(field, is_details=False):
        # Convert to string if not already
        if not isinstance(field, str):
            field = str(field)
        
        # Replace newlines with semicolons to keep CSV valid
        field = field.replace('\n', '; ')
        
        # For details field, replace commas with semicolons to avoid CSV quoting
        if is_details and ',' in field:
            field = field.replace(', ', '; ')
        
        # If the field contains commas or quotes (but not just semicolons), wrap it in quotes
        if ',' in field or '"' in field:
            # Double any quotes within the field
            field = field.replace('"', '""')
            # Wrap in quotes
            field = f'"{field}"'
        return field

    # Map status (color) to flag_code (number)
    flag_codes = {
        "GREEN": "20",   # Using strings for consistency in CSV
        "YELLOW": "30",
        "RED": "50",
        "HALT": "80"
    }

    # Get numeric flag code based on status color
    flag_code = flag_codes.get(status, "80")  # Default to HALT if unknown status
    
    # Format all fields
    component = escape_field(component)
    sample_id = escape_field(sample_id)
    check_name = escape_field(check_name)
    status = escape_field(status)  # The color (GREEN/YELLOW/RED/HALT)
    message = escape_field(message)
    details = escape_field(details, True)
    
    # Write both status (color) and flag_code (number)
    with open(log_path, 'a') as f:
        f.write(f"{component},{sample_id},{check_name},{status},{flag_code},{message},{details}\n")

def check_raw_fastq_existence(outdir, samples, paired_end, log_path):
    """Check if the expected FASTQ files exist for each sample."""
    fastq_dir = os.path.join(outdir, "00-RawData", "Fastq")
    missing_files = []
    
    is_paired = paired_end[0]  # Assuming all samples have the same paired_end value
    
    for sample in samples:
        if is_paired:
            # Check for R1 and R2 files for paired-end sequencing
            r1_file = os.path.join(fastq_dir, f"{sample}_R1_raw.fastq.gz")
            r2_file = os.path.join(fastq_dir, f"{sample}_R2_raw.fastq.gz")
            
            if not os.path.exists(r1_file):
                missing_files.append(r1_file)
            if not os.path.exists(r2_file):
                missing_files.append(r2_file)
        else:
            # Check for single file for single-end sequencing
            single_file = os.path.join(fastq_dir, f"{sample}_raw.fastq.gz")
            
            if not os.path.exists(single_file):
                missing_files.append(single_file)
    
    if missing_files:
        print(f"WARNING: The following FASTQ files are missing:")
        for file in missing_files:
            print(f"  - {file}")
        log_check_result(log_path, "raw_reads", "all", "check_raw_fastq_existence", "HALT", 
                         f"Missing {len(missing_files)} FASTQ files", ",".join(missing_files))
        return False
    
    print(f"All expected FASTQ files found for {len(samples)} samples")
    log_check_result(log_path, "raw_reads", "all", "check_raw_fastq_existence", "GREEN", 
                     f"All FASTQ files found", "")
    return True

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
        log_check_result(log_path, "raw_reads", "all", "check_gzip_integrity", "HALT", 
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
        log_check_result(log_path, "raw_reads", "all", "check_gzip_integrity", "HALT", 
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
        print(f"WARNING: The following FASTQ files are invalid:")
        for file in invalid_files:
            print(f"  - {file}")
        log_check_result(log_path, "raw_reads", "all", "validate_fastq_format", "HALT", 
                        f"Invalid format in {len(invalid_files)} of {len(all_files)} files", 
                        ",".join(invalid_files))
        return False
    
    if all_files:
        print(f"All {len(all_files)} FASTQ files are valid")
        log_check_result(log_path, "raw_reads", "all", "validate_fastq_format", "GREEN", 
                        "All files valid", "")
        return True
    else:
        print("No FASTQ files found to validate")
        log_check_result(log_path, "raw_reads", "all", "validate_fastq_format", "HALT", 
                        "No files found to validate", "")
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
        log_check_result(log_path, "raw_reads", "all", "check_raw_fastqc_existence", "HALT", 
                         f"Missing {len(missing_files)} FastQC output files", ",".join(missing_files))
        return False
    
    print(f"All expected FastQC output files found for {len(samples)} samples")
    log_check_result(log_path, "raw_reads", "all", "check_raw_fastqc_existence", "GREEN", 
                     f"All FastQC output files found", "")
    return True

def check_samples_multiqc(outdir, samples, paired_end, log_path, assay_suffix="_GLbulkRNAseq"):
    """Check if all samples are included in the MultiQC report."""
    multiqc_dir = os.path.join(outdir, "00-RawData", "MultiQC_Reports")
    multiqc_data_zip = os.path.join(multiqc_dir, f"raw_multiqc{assay_suffix}_data.zip")
    
    if not os.path.exists(multiqc_data_zip):
        print(f"WARNING: MultiQC data zip file not found: {multiqc_data_zip}")
        log_check_result(log_path, "raw_reads", "all", "check_samples_multiqc", "RED", 
                         "MultiQC data zip not found", multiqc_data_zip)
        return False
    
    print(f"Found MultiQC data zip: {multiqc_data_zip}")
    
    # Create a temporary directory to extract files
    with tempfile.TemporaryDirectory() as temp_dir:
        try:
            # Extract the zip file
            with zipfile.ZipFile(multiqc_data_zip, 'r') as zip_ref:
                zip_ref.extractall(temp_dir)
            
            # Check for the JSON file (new path structure)
            json_path = os.path.join(temp_dir, f"raw_multiqc{assay_suffix}_data", "multiqc_data.json")
            
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

def get_raw_multiqc_stats(outdir, samples, paired_end, log_path, assay_suffix="_GLbulkRNAseq"):
    """Extract raw MultiQC stats for all samples."""
    multiqc_dir = os.path.join(outdir, "00-RawData", "MultiQC_Reports")
    multiqc_data_zip = os.path.join(multiqc_dir, f"raw_multiqc{assay_suffix}_data.zip")
    multiqc_html = os.path.join(multiqc_dir, f"raw_multiqc{assay_suffix}.html")
    
    if not os.path.exists(multiqc_data_zip):
        print(f"WARNING: MultiQC data zip file not found: {multiqc_data_zip}")
        log_check_result(log_path, "raw_reads", "all", "get_raw_multiqc_stats", "RED", 
                         "MultiQC data zip not found", "")
        return False
    
    print(f"Extracting stats from MultiQC data: {multiqc_data_zip}")
    
    # Create a temporary directory to extract files
    with tempfile.TemporaryDirectory() as temp_dir:
        try:
            # Extract the zip file
            with zipfile.ZipFile(multiqc_data_zip, 'r') as zip_ref:
                zip_ref.extractall(temp_dir)
            
            # Path directly to the extracted data directory
            multiqc_data_dir = os.path.join(temp_dir, f"raw_multiqc{assay_suffix}_data")
            
            # Parse FastQC data using the top-level temp directory as the base path
            fastqc_data = parse_fastqc(os.path.join(temp_dir, "raw"), assay_suffix)
            
            # Collect all possible column names across all samples
            all_columns = set()
            for sample_data in fastqc_data.values():
                all_columns.update(sample_data.keys())
            
            # Sort column names for consistent output and clean up paths
            sorted_columns = sorted(list(all_columns))
            
            # Clean the column names - strip the temporary directory path
            clean_columns = []
            for col in sorted_columns:
                # If the column name contains a path, extract just the metric part
                if '/' in col:
                    # Get the last part of the path (the metric name)
                    metric_part = col.split('/')[-1]
                    clean_columns.append(metric_part)
                else:
                    clean_columns.append(col)
            
            # Create a mapping from original columns to clean columns
            column_mapping = dict(zip(sorted_columns, clean_columns))
            
            log_check_result(log_path, "raw_reads", "all", "get_raw_multiqc_stats", "GREEN", 
                            f"Extracted MultiQC stats for {len(fastqc_data)} samples", "")
            
            # Create a new version of the data with clean column names
            clean_data = {}
            for sample, sample_data in fastqc_data.items():
                clean_data[sample] = {}
                for col, value in sample_data.items():
                    clean_data[sample][column_mapping[col]] = value
            
            # Check for missing samples in MultiQC stats
            missing_samples = []
            for sample in samples:
                if sample not in fastqc_data:
                    missing_samples.append(sample)
            
            if missing_samples:
                print(f"WARNING: {len(missing_samples)} samples missing from MultiQC stats:")
                for sample in missing_samples[:10]:
                    print(f"  - {sample}")
                log_check_result(log_path, "raw_reads", "all", "get_raw_multiqc_stats", "RED", 
                               f"Missing {len(missing_samples)} samples in MultiQC stats", 
                               ",".join(missing_samples[:20]))
                return clean_data  # Still return data for partial analysis
            
            return clean_data
            
        except Exception as e:
            print(f"Error extracting MultiQC stats: {str(e)}")
            log_check_result(log_path, "raw_reads", "all", "get_raw_multiqc_stats", "RED", 
                           f"Error extracting MultiQC stats: {str(e)}", "")
            return False

def parse_fastqc(prefix, assay_suffix):
    """Parse MultiQC JSON data to extract FastQC metrics."""
    # Updated path to directly access the multiqc_data.json in the data directory
    with open(f'{prefix}_multiqc{assay_suffix}_data/multiqc_data.json') as f:
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
                    code = threshold["code"]
                    details = (f"Median={median_value:.2f}, StDev={stdev_value:.2f}, " 
                            f"Threshold={threshold['stdev_threshold']} standard deviations")
                    outlier_details.append((sample, metric_key, code, outlier_message, details))
                    break
    
    # Report results
    if outlier_details:
        print(f"Found {len(outlier_samples)} samples with outlier metrics:")
        for sample, metric, code, message, details in outlier_details:
            print(f"  - {message} [{code}]")
            log_check_result(log_path, "raw_reads", sample, f"outlier_{metric}", code, message, details)
        
        # Keep the list of outlier samples in the summary log
        log_check_result(log_path, "raw_reads", "all", "check_for_outliers", "YELLOW", 
                        f"Found {len(outlier_samples)} samples with outlier metrics", 
                        "; ".join(outlier_samples))  # Include outlier samples list
    else:
        print("No outliers found in MultiQC metrics")
        log_check_result(log_path, "raw_reads", "all", "check_for_outliers", "GREEN", 
                        "No outliers found", "")
    
    return len(outlier_details) > 0

def check_paired_read_counts(multiqc_data, log_path, paired_end):
    """Check if R1 and R2 read counts match for paired-end samples."""
    if not multiqc_data:
        print("No MultiQC data to analyze read count comparison")
        log_check_result(log_path, "raw_reads", "all", "check_paired_read_counts", "RED", 
                        "No MultiQC data available", "Cannot perform read count comparison")
        return False
    
    # Check if we have paired-end data by looking for _r keys
    has_paired_data = any("raw_total_sequences_r" in sample_data for sample_data in multiqc_data.values())
    
    # Check if there's a mismatch between runsheet and actual data
    is_paired_in_runsheet = paired_end[0]  # From runsheet
    if is_paired_in_runsheet and not has_paired_data:
        print("WARNING: Runsheet specifies paired-end but data appears to be single-end")
        log_check_result(log_path, "raw_reads", "all", "check_paired_read_counts", "RED", 
                        "Paired-end/single-end mismatch", 
                        "Runsheet specifies paired-end but data appears to be single-end")
        return False
    elif not is_paired_in_runsheet and has_paired_data:
        print("WARNING: Runsheet specifies single-end but data appears to be paired-end")
        log_check_result(log_path, "raw_reads", "all", "check_paired_read_counts", "RED", 
                        "Paired-end/single-end mismatch", 
                        "Runsheet specifies single-end but data appears to be paired-end")
        return False
    
    # Only report as GREEN for single-end data
    if not is_paired_in_runsheet and not has_paired_data:
        print("Single-end data detected (as specified in runsheet)")
        log_check_result(log_path, "raw_reads", "all", "check_paired_read_counts", "GREEN", 
                        "Single-end data confirmed", "Data matches runsheet specification")
        return True

    # Rest of the function remains the same for paired-end checks
    mismatched_samples = []
    
    # Compare forward and reverse read counts for each sample
    for sample, sample_data in multiqc_data.items():
        if "raw_total_sequences_f" in sample_data and "raw_total_sequences_r" in sample_data:
            # Get the read counts
            forward_reads = float(sample_data["raw_total_sequences_f"])
            reverse_reads = float(sample_data["raw_total_sequences_r"])
            
            # Calculate the difference percentage
            if forward_reads == 0 and reverse_reads == 0:
                continue  # Skip if both are zero
                
            max_reads = max(forward_reads, reverse_reads)
            difference_pct = abs(forward_reads - reverse_reads) / max_reads * 100
            
            # Flag if difference is more than 0.1%
            if difference_pct > 0.1:
                message = (f"Read count mismatch: forward={int(forward_reads)}, "
                          f"reverse={int(reverse_reads)}, difference={difference_pct:.2f}%")
                mismatched_samples.append((sample, message))
                log_check_result(log_path, "raw_reads", sample, "check_paired_read_counts", "RED", 
                                message, "")
    
    # Report results
    if mismatched_samples:
        print(f"Found {len(mismatched_samples)} samples with mismatched read counts:")
        for sample, message in mismatched_samples:
            print(f"  - {sample}: {message}")
        log_check_result(log_path, "raw_reads", "all", "check_paired_read_counts", "RED", 
                        f"Found {len(mismatched_samples)} samples with mismatched read counts", 
                        ",".join([s[0] for s in mismatched_samples]))
        return False
    else:
        print("All paired-end samples have matching read counts")
        log_check_result(log_path, "raw_reads", "all", "check_paired_read_counts", "GREEN", 
                        "All paired-end read counts match", "")
        return True

def report_read_depth_stats(multiqc_data, log_path):
    """Report min, max, and median read depths from MultiQC data."""
    if not multiqc_data:
        print("No MultiQC data to analyze read depths")
        log_check_result(log_path, "raw_reads", "all", "read_depth_stats", "YELLOW", 
                        "No data available", "Cannot report read depth statistics")
        return

    # Extract read counts for all samples
    read_counts = []
    for sample, data in multiqc_data.items():
        # The field name in the processed data has a prefix and suffix
        if 'raw_total_sequences_f' in data:
            read_counts.append((sample, data['raw_total_sequences_f']))
    
    if not read_counts:
        print("No read count data found in MultiQC stats")
        # Print all available keys for debugging
        if multiqc_data:
            first_sample = next(iter(multiqc_data.keys()))
            print(f"Available fields for sample '{first_sample}': {sorted(multiqc_data[first_sample].keys())}")
        log_check_result(log_path, "raw_reads", "all", "read_depth_stats", "YELLOW", 
                        "No read count data found", "")
        return
    
    # Calculate statistics
    counts = [count for _, count in read_counts]
    min_count = min(counts)
    max_count = max(counts)
    median_count = median(counts)
    
    # Just convert everything to strings without any special formatting
    min_str = str(int(min_count))
    max_str = str(int(max_count))
    median_str = str(int(median_count))
    unit = "reads"
    
    # Report results
    print(f"Read depth statistics:")
    print(f"  Minimum: {min_str}")
    print(f"  Maximum: {max_str}")
    print(f"  Median: {median_str}")
    
    # Log the result
    message = f"Read depth range: {min_str} - {max_str} {unit}"
    details = f"Min: {min_str}; Max: {max_str}; Median: {median_str}; From {len(counts)} samples"
    
    log_check_result(log_path, "raw_reads", "all", "read_depth_stats", "GREEN", message, details)

def report_duplication_rate_stats(multiqc_data, log_path):
    """Report min, max, and median duplication rates from MultiQC data."""
    if not multiqc_data:
        print("No MultiQC data to analyze duplication rates")
        log_check_result(log_path, "raw_reads", "all", "duplication_rate_stats", "YELLOW", 
                        "No data available", "Cannot report duplication rate statistics")
        return

    # Extract duplication rates for all samples
    duplication_rates = []
    for sample, data in multiqc_data.items():
        if 'raw_percent_duplicates_f' in data:
            duplication_rates.append((sample, data['raw_percent_duplicates_f']))
    
    if not duplication_rates:
        print("No duplication rate data found in MultiQC stats")
        # Print all available keys for debugging
        if multiqc_data:
            first_sample = next(iter(multiqc_data.keys()))
            print(f"Available fields for sample '{first_sample}': {sorted(multiqc_data[first_sample].keys())}")
        log_check_result(log_path, "raw_reads", "all", "duplication_rate_stats", "YELLOW", 
                        "No duplication rate data found", "")
        return
    
    # Calculate statistics
    rates = [rate for _, rate in duplication_rates]
    min_rate = min(rates)
    max_rate = max(rates)
    median_rate = median(rates)
    
    # Report results
    print(f"Duplication rate statistics:")
    print(f"  Minimum: {min_rate:.2f}%")
    print(f"  Maximum: {max_rate:.2f}%")
    print(f"  Median: {median_rate:.2f}%")
    
    # Log the result
    message = f"Duplication rate range: {min_rate:.2f}% - {max_rate:.2f}%"
    details = f"Min: {min_rate:.2f}%; Max: {max_rate:.2f}%; Median: {median_rate:.2f}%; From {len(rates)} samples"
    
    log_check_result(log_path, "raw_reads", "all", "duplication_rate_stats", "GREEN", message, details)

def report_gc_content_stats(multiqc_data, log_path):
    """Report GC content statistics from MultiQC data, focusing on peak values and outliers."""
    if not multiqc_data:
        print("No MultiQC data to analyze GC content")
        log_check_result(log_path, "raw_reads", "all", "gc_content_stats", "YELLOW", 
                        "No data available", "Cannot report GC content statistics")
        return

    # Extract GC content percentages for all samples
    gc_contents = []
    search_fields = ['raw_percent_gc_f', 'percent_gc_f', 'percent_gc']
    field_used = None
    
    for sample, data in multiqc_data.items():
        # Try different possible field names
        for field in search_fields:
            if field in data:
                gc_contents.append((sample, data[field]))
                if not field_used:
                    field_used = field
                break
    
    if not gc_contents:
        print(f"No GC content data found. Looked for fields: {search_fields}")
        log_check_result(log_path, "raw_reads", "all", "gc_content_stats", "YELLOW", 
                        f"No GC content data found. Tried fields: {', '.join(search_fields)}", "")
        return
    
    # Extract just the GC percentages
    gc_values = [gc for _, gc in gc_contents]
    
    # Find the most common GC percentage (peak)
    # Round to nearest integer to identify the peak
    gc_counts = {}
    for gc in gc_values:
        gc_int = int(round(gc))
        gc_counts[gc_int] = gc_counts.get(gc_int, 0) + 1
    
    # Find the peak GC percentage (mode)
    peak_gc = max(gc_counts.items(), key=lambda x: x[1])[0]
    
    # Calculate statistics
    min_gc = min(gc_values)
    max_gc = max(gc_values)
    median_gc = median(gc_values)
    
    # Identify outliers (samples more than 2 standard deviations from median)
    stddev_gc = np.std(gc_values) if len(gc_values) > 1 else 0
    outliers = []
    for sample, gc in gc_contents:
        if abs(gc - median_gc) > 2 * stddev_gc:
            deviation = abs(gc - median_gc) / stddev_gc
            outliers.append(f"{sample} ({gc:.1f}%; {deviation:.2f} SD)")
    
    # Report results
    print(f"GC content statistics (using field '{field_used}'):")
    print(f"  Mode: {peak_gc}%")
    print(f"  Range: {min_gc:.1f}% - {max_gc:.1f}%")
    print(f"  Median: {median_gc:.1f}%")
    
    if outliers:
        print(f"  Outliers ({len(outliers)}):")
        for outlier in outliers:
            print(f"    - {outlier}")
    
    # Log the result
    message = f"GC content mode: {peak_gc}%"
    details = f"Range: {min_gc:.1f}% - {max_gc:.1f}%; Median: {median_gc:.1f}%; From {len(gc_values)} samples"
    if outliers:
        details += f"; Outliers: {'; '.join(outliers)}"
    
    log_check_result(log_path, "raw_reads", "all", "gc_content_stats", "GREEN", message, details)

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
    
    # 1. Check for raw FASTQ files existence
    check_raw_fastq_existence(args.outdir, sample_names, paired_end_values, vv_log_path)
    
    # 2. Check GZIP integrity of FASTQ files
    check_gzip_integrity(args.outdir, sample_names, paired_end_values, vv_log_path)
    
    # 3. Validate FASTQ format
    validate_fastq_format(args.outdir, sample_names, paired_end_values, vv_log_path)
    
    # 4. Check FastQC outputs existence
    check_raw_fastqc_existence(args.outdir, sample_names, paired_end_values, vv_log_path)
    
    # 5. Check all samples are in MultiQC report
    check_samples_multiqc(args.outdir, sample_names, paired_end_values, vv_log_path, args.assay_suffix)

    # 6. Get MultiQC stats
    multiqc_data = get_raw_multiqc_stats(args.outdir, sample_names, paired_end_values, vv_log_path, args.assay_suffix)
    
    # 7. Report MultiQC stats outliers
    if multiqc_data:
        report_multiqc_outliers(args.outdir, multiqc_data, vv_log_path)
        
        # 8. Check paired read counts match (for paired-end data)
        check_paired_read_counts(multiqc_data, vv_log_path, paired_end_values)

    # 9. Report read depth stats
    if multiqc_data:
        report_read_depth_stats(multiqc_data, vv_log_path)
        
    # 10. Report duplication rate stats
    if multiqc_data:
        report_duplication_rate_stats(multiqc_data, vv_log_path)
        
    # 11. Report GC content stats (new)
    if multiqc_data:
        report_gc_content_stats(multiqc_data, vv_log_path)

if __name__ == "__main__":
    main()