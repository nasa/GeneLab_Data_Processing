#!/usr/bin/env python
"""
Script to validate and verify trimmed reads based on runsheet information.
Expects to run from inside a output directory GLDS-##.

Parse the input runsheet to get:
Sample Name
Paired End
Has ERCC

Check that the expected output directories exist

Section-specific checks:

check_trimmed_fastq_existence: Check if all expected trimmed FASTQ files exist for each sample
check_gzip_integrity: Verify GZIP integrity of FASTQ files using gzip -t
validate_fastq_format: Ensure FASTQ files have proper formatting (headers start with @)
check_trimmed_fastqc_existence: Verify FastQC output files exist (HTML and ZIP files)
check_samples_multiqc: Confirm all samples are included in the MultiQC report
get_trimmed_multiqc_stats: Extract FastQC metrics from MultiQC report
report_multiqc_outliers: Identify and report outliers in FastQC metrics
check_paired_read_counts: Ensure paired-end reads have matching read counts
check_trimming_report_existence: Verify trimming report files exist
check_trimming_multiqc_samples: Confirm all samples are included in the trimming MultiQC report
check_adapters_presence: Check for adapter contamination in trimmed reads
report_read_depth_stats: Log statistics about sequencing depth across all samples
report_duplication_rate_stats: Log statistics about duplication rates across all samples
report_gc_content_stats: Log statistics about GC content across all samples

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
    trimmed_data_path = os.path.join(outdir, "01-TG_Preproc")
    fastq_path = os.path.join(trimmed_data_path, "Fastq")
    fastqc_path = os.path.join(trimmed_data_path, "FastQC_Reports")
    multiqc_path = os.path.join(trimmed_data_path, "MultiQC_Reports")
    
    missing_dirs = []
    
    if not os.path.exists(trimmed_data_path):
        missing_dirs.append(str(trimmed_data_path))
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
    # Properly escape and format fields for CSV
    def escape_field(field, is_details=False):
        # Convert to string if not already
        if not isinstance(field, str):
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
    status = escape_field(status)
    message = escape_field(message)
    details = escape_field(details, True)  # Specify this is a details field
    
    # Write the formatted line
    with open(log_path, 'a') as f:
        f.write(f"{component},{sample_id},{check_name},{status},{flag_code},{message},{details}\n")

def check_trimmed_fastq_existence(outdir, samples, paired_end, log_path):
    """Check if the expected trimmed FASTQ files exist for each sample."""
    fastq_dir = os.path.join(outdir, "01-TG_Preproc", "Fastq")
    missing_files = []
    
    is_paired = paired_end[0]  # Assuming all samples have the same paired_end value
    
    for sample in samples:
        if is_paired:
            # Check for R1 and R2 files for paired-end sequencing
            r1_file = os.path.join(fastq_dir, f"{sample}_R1_trimmed.fastq.gz")
            r2_file = os.path.join(fastq_dir, f"{sample}_R2_trimmed.fastq.gz")
            
            if not os.path.exists(r1_file):
                missing_files.append(r1_file)
            if not os.path.exists(r2_file):
                missing_files.append(r2_file)
        else:
            # Check for single file for single-end sequencing
            single_file = os.path.join(fastq_dir, f"{sample}_trimmed.fastq.gz")
            
            if not os.path.exists(single_file):
                missing_files.append(single_file)
    
    if missing_files:
        print(f"WARNING: Missing {len(missing_files)} trimmed FASTQ files")
        log_check_result(log_path, "trimmed_reads", "all", "check_trimmed_fastq_existence", "HALT", 
                         f"Missing {len(missing_files)} trimmed FASTQ files", 
                         ",".join(missing_files))
        return False
    
    print(f"Found all expected trimmed FASTQ files ({len(samples)} samples)")
    log_check_result(log_path, "trimmed_reads", "all", "check_trimmed_fastq_existence", "GREEN", 
                     "All trimmed FASTQ files found", "")
    return True

def check_gzip_integrity(outdir, samples, paired_end, log_path):
    """Check GZIP integrity for trimmed FASTQ files."""
    fastq_dir = os.path.join(outdir, "01-TG_Preproc", "Fastq")
    failed_files = []
    all_files = []
    
    is_paired = paired_end[0]  # Assuming all samples have the same paired_end value
    
    for sample in samples:
        if is_paired:
            # Get R1 and R2 files for paired-end sequencing
            r1_file = os.path.join(fastq_dir, f"{sample}_R1_trimmed.fastq.gz")
            r2_file = os.path.join(fastq_dir, f"{sample}_R2_trimmed.fastq.gz")
            files_to_check = [r1_file, r2_file]
        else:
            # Get single file for single-end sequencing
            single_file = os.path.join(fastq_dir, f"{sample}_trimmed.fastq.gz")
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
        log_check_result(log_path, "trimmed_reads", "all", "check_gzip_integrity", "HALT", 
                         f"Integrity check failed for {len(failed_files)} of {len(all_files)} files", 
                         ",".join(failed_files))
        return False
    
    if all_files:
        log_check_result(log_path, "trimmed_reads", "all", "check_gzip_integrity", "GREEN", 
                         f"Integrity check passed for all {len(all_files)} files", "")
        return True
    else:
        log_check_result(log_path, "trimmed_reads", "all", "check_gzip_integrity", "HALT", 
                         "No trimmed FASTQ files found to check", "")
        return False

def validate_fastq_format(outdir, samples, paired_end, log_path, max_lines=200000000):
    """Validate FASTQ format by checking that header lines start with @.
    Only checks the first 200 million lines for performance reasons."""
    trimmed_dir = os.path.join(outdir, "01-TG_Preproc", "Fastq")
    invalid_files = []
    all_files = []
    
    is_paired = paired_end[0]  # Assuming all samples have the same paired_end value
    
    for sample in samples:
        if is_paired:
            # Get R1 and R2 files for paired-end sequencing
            r1_file = os.path.join(trimmed_dir, f"{sample}_R1_trimmed.fastq.gz")
            r2_file = os.path.join(trimmed_dir, f"{sample}_R2_trimmed.fastq.gz")
            files_to_check = [r1_file, r2_file]
        else:
            # Get single file for single-end sequencing
            single_file = os.path.join(trimmed_dir, f"{sample}_trimmed.fastq.gz")
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
        log_check_result(log_path, "trimmed_reads", "all", "validate_fastq_format", "HALT", 
                        f"Invalid format in {len(invalid_files)} of {len(all_files)} files", 
                        ",".join(invalid_files))
        return False
    
    if not all_files:
        print("No FASTQ files found to validate")
        log_check_result(log_path, "trimmed_reads", "all", "validate_fastq_format", "HALT", 
                        "No trimmed FASTQ files found to validate", "")
        return False

    print(f"All {len(all_files)} FASTQ files are valid")
    log_check_result(log_path, "trimmed_reads", "all", "validate_fastq_format", "GREEN", 
                    "All files valid", "")
    return True

def check_trimmed_fastqc_existence(outdir, samples, paired_end, log_path):
    """Check if FastQC output files exist for each sample."""
    fastqc_dir = os.path.join(outdir, "01-TG_Preproc", "FastQC_Reports")
    missing_files = []
    
    is_paired = paired_end[0]  # Assuming all samples have the same paired_end value
    
    for sample in samples:
        if is_paired:
            # Check for R1 and R2 FastQC files for paired-end sequencing
            r1_html = os.path.join(fastqc_dir, f"{sample}_R1_trimmed_fastqc.html")
            r1_zip = os.path.join(fastqc_dir, f"{sample}_R1_trimmed_fastqc.zip")
            r2_html = os.path.join(fastqc_dir, f"{sample}_R2_trimmed_fastqc.html")
            r2_zip = os.path.join(fastqc_dir, f"{sample}_R2_trimmed_fastqc.zip")
            
            files_to_check = [r1_html, r1_zip, r2_html, r2_zip]
        else:
            # Check for single FastQC files for single-end sequencing
            html_file = os.path.join(fastqc_dir, f"{sample}_trimmed_fastqc.html")
            zip_file = os.path.join(fastqc_dir, f"{sample}_trimmed_fastqc.zip")
            
            files_to_check = [html_file, zip_file]
        
        for file_path in files_to_check:
            if not os.path.exists(file_path):
                missing_files.append(file_path)
    
    if missing_files:
        print(f"WARNING: The following FastQC output files are missing:")
        for file in missing_files:
            print(f"  - {file}")
        log_check_result(log_path, "trimmed_reads", "all", "check_trimmed_fastqc_existence", "HALT", 
                         f"Missing {len(missing_files)} FastQC output files", ",".join(missing_files))
        return False
    
    print(f"All expected FastQC output files found for {len(samples)} samples")
    log_check_result(log_path, "trimmed_reads", "all", "check_trimmed_fastqc_existence", "GREEN", 
                     f"All FastQC output files found", "")
    return True

def check_samples_multiqc(outdir, samples, paired_end, log_path, assay_suffix="_GLbulkRNAseq"):
    """Check if all samples are included in the MultiQC report."""
    multiqc_dir = os.path.join(outdir, "01-TG_Preproc", "MultiQC_Reports")
    multiqc_data_zip = os.path.join(multiqc_dir, f"trimmed_multiqc{assay_suffix}_data.zip")
    
    if not os.path.exists(multiqc_data_zip):
        print(f"WARNING: MultiQC data zip file not found: {multiqc_data_zip}")
        log_check_result(log_path, "trimmed_reads", "all", "check_samples_multiqc", "RED", 
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
            json_path = os.path.join(temp_dir, f"trimmed_multiqc{assay_suffix}_data", "multiqc_data.json")
            
            if not os.path.exists(json_path):
                print(f"Could not find multiqc_data.json in the expected location")
                log_check_result(log_path, "trimmed_reads", "all", "check_samples_multiqc", "RED", 
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
                    base_sample = mqc_sample.replace("_trimmed_fastqc", "").replace("_fastqc", "")
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
                log_check_result(log_path, "trimmed_reads", "all", "check_samples_multiqc", "RED", 
                                f"Missing {len(missing_samples)} samples in MultiQC report", 
                                ",".join(missing_samples))
                return False
            
            print(f"All {len(samples)} samples found in the MultiQC report")
            log_check_result(log_path, "trimmed_reads", "all", "check_samples_multiqc", "GREEN", 
                            "All samples found in MultiQC report", "")
            return True
            
        except Exception as e:
            print(f"Error processing MultiQC report: {str(e)}")
            log_check_result(log_path, "trimmed_reads", "all", "check_samples_multiqc", "RED", 
                           f"Error processing MultiQC report: {str(e)}", "")
            return False

def get_trimmed_multiqc_stats(outdir, samples, paired_end, log_path, assay_suffix="_GLbulkRNAseq"):
    """Extract trimmed MultiQC stats for all samples and write to a stats file for analysis."""
    multiqc_dir = os.path.join(outdir, "01-TG_Preproc", "MultiQC_Reports")
    multiqc_data_zip = os.path.join(multiqc_dir, f"trimmed_multiqc{assay_suffix}_data.zip")
    
    if not os.path.exists(multiqc_data_zip):
        print(f"WARNING: MultiQC data zip file not found: {multiqc_data_zip}")
        log_check_result(log_path, "trimmed_reads", "all", "get_trimmed_multiqc_stats", "RED", 
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
            multiqc_data_dir = os.path.join(temp_dir, f"trimmed_multiqc{assay_suffix}_data")
            
            # Parse FastQC data using the top-level temp directory as the base path
            fastqc_data = parse_fastqc(os.path.join(temp_dir, "trimmed"), assay_suffix)
            
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
            
            log_check_result(log_path, "trimmed_reads", "all", "get_trimmed_multiqc_stats", "GREEN", 
                            f"Extracted MultiQC stats for {len(fastqc_data)} samples", "")
            
            # Create a new version of the data with clean column names
            clean_data = {}
            for sample, sample_data in fastqc_data.items():
                clean_data[sample] = {}
                for col, value in sample_data.items():
                    clean_data[sample][column_mapping[col]] = value
            
            return clean_data
            
        except Exception as e:
            print(f"Error extracting MultiQC stats: {str(e)}")
            log_check_result(log_path, "trimmed_reads", "all", "get_trimmed_multiqc_stats", "RED", 
                           f"Error extracting MultiQC stats: {str(e)}", "")
            return False

def parse_fastqc(prefix, assay_suffix):
    """Parse MultiQC JSON data to extract FastQC metrics with more flexible field name handling."""
    with open(f'{prefix}_multiqc{assay_suffix}_data/multiqc_data.json') as f:
        j = json.loads(f.read())

    # Print debug info about the data structure
    print(f"\nDEBUG: Parsing FastQC data from {prefix}_multiqc{assay_suffix}_data/multiqc_data.json")
    
    # Check if we have what we expect in the data structure
    if 'report_general_stats_data' not in j or not j['report_general_stats_data']:
        print("WARNING: No report_general_stats_data found in MultiQC JSON")
        return {}
    
    # Debug info about stats modules
    stats_modules = j['report_general_stats_data']
    print(f"Found {len(stats_modules)} stats modules in general_stats_data")
    
    # Find which module contains the FastQC data
    fastqc_module_index = -1
    for i, module in enumerate(stats_modules):
        if module:  # Check if module has any data
            sample = next(iter(module.keys()))
            sample_data = module[sample]
            # Look for typical FastQC fields
            if any(field in sample_data for field in ['total_sequences', 'percent_gc']):
                fastqc_module_index = i
                print(f"Found FastQC data in stats module {i}")
                break
    
    if fastqc_module_index == -1:
        print("WARNING: Could not find FastQC data in any stats module")
        return {}
    
    # Group the samples by base name for paired end data
    sample_groups = {}
    for sample in j['report_general_stats_data'][fastqc_module_index].keys():
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
        
        # Map between expected field names and actual field names in the MultiQC
        field_mapping = {
            'total_sequences': ['total_sequences'],
            'percent_gc': ['percent_gc'],
            'avg_sequence_length': ['avg_sequence_length'],
            'median_sequence_length': ['median_sequence_length'],
            'percent_duplicates': ['percent_duplicates']
        }
        
        # Process forward read
        if reads['f']:
            sample_fields = j['report_general_stats_data'][fastqc_module_index][reads['f']]
            
            # Copy data with appropriate field name translation
            for expected_field, possible_actual_fields in field_mapping.items():
                for actual_field in possible_actual_fields:
                    if actual_field in sample_fields:
                        # Add field with prefix and _f suffix for consistency
                        data[base_name][prefix + '_' + expected_field + '_f'] = sample_fields[actual_field]
                        break
        
        # Process reverse read
        if reads['r']:
            sample_fields = j['report_general_stats_data'][fastqc_module_index][reads['r']]
            
            # Copy data with appropriate field name translation
            for expected_field, possible_actual_fields in field_mapping.items():
                for actual_field in possible_actual_fields:
                    if actual_field in sample_fields:
                        # Add field with prefix and _r suffix for consistency
                        data[base_name][prefix + '_' + expected_field + '_r'] = sample_fields[actual_field]
                        break

    # Process other stats sections (quality, GC, etc) - keep this logic mostly the same
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

    # Print summary of what was found
    print(f"Extracted data for {len(data)} samples")
    if data:
        first_sample = next(iter(data.keys()))
        print(f"Sample '{first_sample}' has fields: {sorted(data[first_sample].keys())}")
    
    return data

def report_multiqc_outliers(outdir, multiqc_data, log_path):
    """Identify and report outliers in MultiQC statistics."""
    if not multiqc_data:
        print("No MultiQC data to analyze for outliers")
        return False

    # Keys to check for outliers as specified in protocol.py
    metrics_to_check = {
        "trimmed_percent_gc_f": "percent GC content (forward)",
        "trimmed_avg_sequence_length_f": "average sequence length (forward)",
        "trimmed_total_sequences_f": "total sequences (forward)",
        "trimmed_percent_duplicates_f": "percent duplicates (forward)",
    }
    
    # For paired-end data, also check the reverse reads
    if any("_r" in key for sample in multiqc_data.values() for key in sample):
        metrics_to_check.update({
            "trimmed_percent_gc_r": "percent GC content (reverse)",
            "trimmed_avg_sequence_length_r": "average sequence length (reverse)",
            "trimmed_total_sequences_r": "total sequences (reverse)",
            "trimmed_percent_duplicates_r": "percent duplicates (reverse)",
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
                    break  # Stop after finding the first matching threshold
    
    # Report results
    if outlier_details:
        print(f"Found {len(outlier_samples)} samples with outlier metrics:")
        for sample, metric, code, message, details in outlier_details:
            print(f"  - {message} [{code}]")
            log_check_result(log_path, "trimmed_reads", sample, f"outlier_{metric}", code, message, details)
        
        # Keep the list of outlier samples in the summary log
        log_check_result(log_path, "trimmed_reads", "all", "check_for_outliers", "YELLOW", 
                        f"Found {len(outlier_samples)} samples with outlier metrics", 
                        "; ".join(outlier_samples))  # Include outlier samples list
    else:
        print("No outliers found in MultiQC metrics")
        log_check_result(log_path, "trimmed_reads", "all", "check_for_outliers", "GREEN", 
                        "No outliers found", "")
    
    return len(outlier_details) > 0

def check_paired_read_counts(multiqc_data, log_path, paired_end):
    """Check if R1 and R2 read counts match for paired-end samples."""
    if not multiqc_data:
        print("No MultiQC data to analyze read count comparison")
        log_check_result(log_path, "trimmed_reads", "all", "check_paired_read_counts", "RED", 
                        "No MultiQC data available", "Cannot perform read count comparison")
        return False
    
    # Print debug info about the data structure
    print("\nDEBUG: Examining MultiQC data structure for paired-end detection")
    first_sample = next(iter(multiqc_data.keys()))
    print(f"Sample name example: {first_sample}")
    print(f"Available fields: {sorted(multiqc_data[first_sample].keys())}")
    
    # Look for paired-end indicators with multiple possible field name patterns
    paired_indicators = ["_r", "_R2", "Read 2", "read2"]
    paired_data_detected = False
    
    # Check all samples and fields for any indication of paired data
    for sample, data in multiqc_data.items():
        # Check if sample name itself indicates paired-end (contains _R1, etc.)
        for indicator in ["_R1", "_r1", " Read 1"]:
            if indicator in sample and sample.replace(indicator, indicator.replace("1", "2")) in multiqc_data:
                print(f"DEBUG: Paired-end data detected from sample names: {sample}")
                paired_data_detected = True
                break
        
        # Check field names for paired indicators
        for field in data.keys():
            for indicator in paired_indicators:
                if indicator in field:
                    print(f"DEBUG: Paired-end data detected from field name: {field}")
                    paired_data_detected = True
                    break
            if paired_data_detected:
                break
        
        if paired_data_detected:
            break
    
    # The special "find _f/_r pairs" logic to detect paired data
    sequence_fields = []
    for sample, data in multiqc_data.items():
        for field in data.keys():
            if "sequence" in field.lower() and field.endswith("_f"):
                base_field = field[:-2]  # Remove the _f
                if base_field + "_r" in data:
                    print(f"DEBUG: Paired-end data detected from field pair: {field} and {base_field}_r")
                    paired_data_detected = True
                    sequence_fields = [field, base_field + "_r"]
                    break
        if paired_data_detected:
            break
    
    # Check if there's a mismatch between runsheet and actual data
    is_paired_in_runsheet = paired_end[0]  # From runsheet
    print(f"DEBUG: Runsheet indicates paired-end: {is_paired_in_runsheet}")
    print(f"DEBUG: Data indicates paired-end: {paired_data_detected}")
    
    if is_paired_in_runsheet and not paired_data_detected:
        print("WARNING: Runsheet specifies paired-end but data appears to be single-end")
        log_check_result(log_path, "trimmed_reads", "all", "check_paired_read_counts", "RED", 
                        "Paired-end/single-end mismatch", 
                        "Runsheet specifies paired-end but data appears to be single-end")
        return False
    elif not is_paired_in_runsheet and paired_data_detected:
        print("WARNING: Runsheet specifies single-end but data appears to be paired-end")
        log_check_result(log_path, "trimmed_reads", "all", "check_paired_read_counts", "RED", 
                        "Paired-end/single-end mismatch", 
                        "Runsheet specifies single-end but data appears to be paired-end")
        return False
    
    # Handle single-end case properly like in the raw reads script
    if not is_paired_in_runsheet and not paired_data_detected:
        print("Single-end data detected (as specified in runsheet)")
        log_check_result(log_path, "trimmed_reads", "all", "check_paired_read_counts", "GREEN", 
                        "Single-end data confirmed", "Data matches runsheet specification")
        return True

    # If we get here, we're dealing with paired-end data as expected
    # We need to check if the read counts match for all samples
    # But first we need to find the actual field names used for read counts
    if sequence_fields and len(sequence_fields) == 2:
        # We already found the fields above
        forward_field, reverse_field = sequence_fields
    else:
        # Try to find the fields by common naming patterns
        possible_pairs = [
            ("trimmed_total_sequences_f", "trimmed_total_sequences_r"),
            ("total_sequences_f", "total_sequences_r"),
            ("total_sequences", "total_sequences")  # Same field for single-end
        ]
        
        forward_field = None
        reverse_field = None
        
        for f_field, r_field in possible_pairs:
            if f_field in multiqc_data[first_sample]:
                forward_field = f_field
                reverse_field = r_field
                break
        
        if not forward_field:
            print("ERROR: Could not determine field names for read counts")
            log_check_result(log_path, "trimmed_reads", "all", "check_paired_read_counts", "RED", 
                            "Could not determine field names for read counts", "")
            return False
    
    print(f"DEBUG: Using fields '{forward_field}' and '{reverse_field}' for read count comparison")
    
    # Now check the read counts
    mismatched_samples = []
    
    # Compare forward and reverse read counts for each sample
    for sample, sample_data in multiqc_data.items():
        if forward_field in sample_data and reverse_field in sample_data:
            # Get the read counts
            forward_reads = float(sample_data[forward_field])
            reverse_reads = float(sample_data[reverse_field])
            
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
                log_check_result(log_path, "trimmed_reads", sample, "check_paired_read_counts", "RED", 
                                message, "")
                print(f"DEBUG: Mismatch found for {sample}: {message}")
    
    # Report results
    if mismatched_samples:
        print(f"Found {len(mismatched_samples)} samples with mismatched read counts:")
        for sample, message in mismatched_samples:
            print(f"  - {sample}: {message}")
        log_check_result(log_path, "trimmed_reads", "all", "check_paired_read_counts", "RED", 
                        f"Found {len(mismatched_samples)} samples with mismatched read counts", 
                        ",".join([s[0] for s in mismatched_samples]))
        return False
    else:
        print("All paired-end samples have matching read counts")
        log_check_result(log_path, "trimmed_reads", "all", "check_paired_read_counts", "GREEN", 
                        "All paired-end read counts match", "")
        return True

def check_trimming_report_existence(outdir, samples, paired_end, log_path):
    """Check if the expected trimming report files exist for each sample."""
    trimming_dir = os.path.join(outdir, "01-TG_Preproc", "Trimming_Reports")
    missing_files = []
    
    is_paired = paired_end[0]  # Assuming all samples have the same paired_end value
    
    for sample in samples:
        if is_paired:
            # For paired-end data, check for both R1 and R2 reports
            r1_report = os.path.join(trimming_dir, f"{sample}_R1_raw.fastq.gz_trimming_report.txt")
            r2_report = os.path.join(trimming_dir, f"{sample}_R2_raw.fastq.gz_trimming_report.txt")
            
            if not os.path.exists(r1_report):
                missing_files.append(r1_report)
            if not os.path.exists(r2_report):
                missing_files.append(r2_report)
        else:
            # For single-end data, check for just one report per sample
            report_file = os.path.join(trimming_dir, f"{sample}_raw.fastq.gz_trimming_report.txt")
            
            if not os.path.exists(report_file):
                missing_files.append(report_file)
    
    if missing_files:
        print(f"WARNING: The following trimming report files are missing:")
        for file in missing_files:
            print(f"  - {file}")
        log_check_result(log_path, "trimmed_reads", "all", "check_trimming_report_existence", "HALT", 
                         f"Missing {len(missing_files)} trimming report files", ",".join(missing_files))
        return False
    
    print(f"All expected trimming report files found for {len(samples)} samples")
    log_check_result(log_path, "trimmed_reads", "all", "check_trimming_report_existence", "GREEN", 
                     f"All trimming report files found", "")
    return True

def check_trimming_multiqc_samples(outdir, samples, log_path, assay_suffix="_GLbulkRNAseq"):
    """Check if all samples are included in the trimming MultiQC report (now combined in the trimmed FastQC MultiQC)."""
    # Use MultiQC_Reports directory instead of FastQC_Reports
    multiqc_dir = os.path.join(outdir, "01-TG_Preproc", "MultiQC_Reports")
    multiqc_data_zip = os.path.join(multiqc_dir, f"trimmed_multiqc{assay_suffix}_data.zip")
    
    if not os.path.exists(multiqc_data_zip):
        print(f"WARNING: Trimmed MultiQC data zip file not found: {multiqc_data_zip}")
        log_check_result(log_path, "trimmed_reads", "all", "check_trimming_multiqc_samples", "RED", 
                         "Trimmed MultiQC data zip not found", "")
        return False
    
    print(f"Found Trimmed MultiQC data zip: {multiqc_data_zip}")
    
    # Create a temporary directory to extract files
    with tempfile.TemporaryDirectory() as temp_dir:
        try:
            # Extract the zip file
            with zipfile.ZipFile(multiqc_data_zip, 'r') as zip_ref:
                zip_ref.extractall(temp_dir)
            
            # Check for the JSON file (new path structure)
            json_path = os.path.join(temp_dir, f"trimmed_multiqc{assay_suffix}_data", "multiqc_data.json")
            
            if not os.path.exists(json_path):
                print(f"WARNING: No multiqc_data.json file found in the expected location")
                log_check_result(log_path, "trimmed_reads", "all", "check_trimming_multiqc_samples", "RED", 
                               "multiqc_data.json not found in zip", "")
                return False
            
            print(f"Using MultiQC data file: {json_path}")
            
            # Parse the JSON file
            with open(json_path) as f:
                multiqc_data = json.load(f)
            
            # Extract sample names from the Cutadapt section
            mqc_samples = set()
            
            # Find Cutadapt data in the combined MultiQC report
            if ('report_data_sources' in multiqc_data and 
                'Cutadapt' in multiqc_data['report_data_sources'] and
                'all_sections' in multiqc_data['report_data_sources']['Cutadapt']):
                
                cutadapt_sections = multiqc_data['report_data_sources']['Cutadapt']['all_sections']
                for mqc_sample in cutadapt_sections.keys():
                    # Extract sample name from the path or filename
                    # First get filename from path
                    raw_name = os.path.basename(cutadapt_sections[mqc_sample])
                    
                    # Handle both paired-end and single-end naming conventions
                    # For paired-end, filenames have _R1 or _R2 before _raw
                    # For single-end, no _R1 or _R2 is present
                    if "_R1_raw" in raw_name:
                        sample_name = raw_name.split("_R1_raw")[0]
                        mqc_samples.add(sample_name)
                    elif "_R2_raw" in raw_name:
                        sample_name = raw_name.split("_R2_raw")[0]
                        mqc_samples.add(sample_name)
                    else:
                        # Handle single-end case where there's just _raw
                        sample_name = raw_name.split("_raw")[0]
                        mqc_samples.add(sample_name)
            
            # Check if all runsheet samples are present in MultiQC report
            missing_samples = []
            for sample in samples:
                if sample not in mqc_samples:
                    missing_samples.append(sample)
            
            if missing_samples:
                print(f"WARNING: The following samples are missing from the Trimming section in the combined MultiQC report:")
                for sample in missing_samples:
                    print(f"  - {sample}")
                log_check_result(log_path, "trimmed_reads", "all", "check_trimming_multiqc_samples", "RED", 
                                f"Missing {len(missing_samples)} samples in Trimming section of MultiQC report", 
                                ",".join(missing_samples))
                return False
            
            print(f"All {len(samples)} samples found in the Trimming section of the MultiQC report")
            log_check_result(log_path, "trimmed_reads", "all", "check_trimming_multiqc_samples", "GREEN", 
                            "All samples found in Trimming section of MultiQC report", "")
            return True
            
        except Exception as e:
            print(f"Error processing Trimming section in MultiQC report: {str(e)}")
            log_check_result(log_path, "trimmed_reads", "all", "check_trimming_multiqc_samples", "RED", 
                           f"Error processing Trimming section in MultiQC report: {str(e)}", "")
            return False

def check_adapters_presence(outdir, samples, paired_end, log_path, threshold=0.001, assay_suffix="_GLbulkRNAseq"):
    """Check if adapter sequences were present in raw reads and properly removed during trimming."""
    # Use the combined trimmed FastQC MultiQC report
    multiqc_dir = os.path.join(outdir, "01-TG_Preproc", "MultiQC_Reports")
    multiqc_data_zip = os.path.join(multiqc_dir, f"trimmed_multiqc{assay_suffix}_data.zip")
    
    if not os.path.exists(multiqc_data_zip):
        print(f"WARNING: Trimmed MultiQC data zip file not found: {multiqc_data_zip}")
        log_check_result(log_path, "trimmed_reads", "all", "check_adapters_presence", "RED", 
                         "Trimmed MultiQC data zip not found", "")
        return False
    
    print(f"Checking adapter presence in trimming report data within MultiQC: {multiqc_data_zip}")

    # Create a temporary directory to extract files
    with tempfile.TemporaryDirectory() as temp_dir:
        try:
            # Extract the zip file
            with zipfile.ZipFile(multiqc_data_zip, 'r') as zip_ref:
                zip_ref.extractall(temp_dir)
            
            # Check for the JSON file (new path structure)
            json_path = os.path.join(temp_dir, f"trimmed_multiqc{assay_suffix}_data", "multiqc_data.json")
            
            if not os.path.exists(json_path):
                print(f"WARNING: No multiqc_data.json file found in the expected location")
                log_check_result(log_path, "trimmed_reads", "all", "check_adapters_presence", "RED", 
                               "multiqc_data.json not found in zip", "")
                return False
            
            # Parse the JSON file
            with open(json_path) as f:
                multiqc_data = json.load(f)
            
            # Get the general stats data from the Cutadapt section in the combined MultiQC
            if 'report_general_stats_data' not in multiqc_data or not multiqc_data['report_general_stats_data']:
                print(f"WARNING: No report_general_stats_data found in MultiQC data")
                log_check_result(log_path, "trimmed_reads", "all", "check_adapters_presence", "RED", 
                               "No stats data found in MultiQC report", "")
                return False
            
            # Find the Cutadapt data in the general stats
            stats_data = None
            for module_data in multiqc_data['report_general_stats_data']:
                # Check if this module contains Cutadapt data
                for sample_key in module_data:
                    if 'r_processed' in module_data[sample_key] and 'r_with_adapters' in module_data[sample_key]:
                        stats_data = module_data
                        break
                if stats_data:
                    break
            
            if not stats_data:
                print(f"WARNING: No Cutadapt data found in MultiQC report")
                log_check_result(log_path, "trimmed_reads", "all", "check_adapters_presence", "RED", 
                               "No Cutadapt data found in combined MultiQC report", "")
                return False
            
            # Check if this is paired-end data
            is_paired = paired_end[0]
            
            # Keep track of samples with low adapter presence
            low_adapter_samples = []
            sample_adapter_ratios = {}  # Store adapter ratios by sample
            
            for sample in samples:
                sample_ratios = []  # Collect ratios for both R1 and R2 if paired
                
                if is_paired:
                    # For paired end data, look for "Read 1" and "Read 2" entries
                    r1_key = f"{sample} Read 1"
                    r2_key = f"{sample} Read 2"
                    
                    # Skip if sample not found in stats
                    if r1_key not in stats_data or r2_key not in stats_data:
                        # Try with _R1 and _R2 format instead
                        r1_key = f"{sample}_R1"
                        r2_key = f"{sample}_R2"
                        
                        if r1_key not in stats_data or r2_key not in stats_data:
                            print(f"WARNING: Sample {sample} not found in paired-end trimming stats")
                            continue
                    
                    # Check R1
                    if 'r_processed' in stats_data[r1_key] and 'r_with_adapters' in stats_data[r1_key]:
                        r_processed = stats_data[r1_key]['r_processed']
                        r_with_adapters = stats_data[r1_key]['r_with_adapters']
                        
                        if r_processed > 0:
                            adapter_ratio_r1 = r_with_adapters / r_processed
                            sample_ratios.append(adapter_ratio_r1)
                    
                    # Check R2
                    if 'r_processed' in stats_data[r2_key] and 'r_with_adapters' in stats_data[r2_key]:
                        r_processed = stats_data[r2_key]['r_processed']
                        r_with_adapters = stats_data[r2_key]['r_with_adapters']
                        
                        if r_processed > 0:
                            adapter_ratio_r2 = r_with_adapters / r_processed
                            sample_ratios.append(adapter_ratio_r2)
                else:
                    # For single end data
                    if sample not in stats_data:
                        print(f"WARNING: Sample {sample} not found in single-end trimming stats")
                        continue
                    
                    if 'r_processed' in stats_data[sample] and 'r_with_adapters' in stats_data[sample]:
                        r_processed = stats_data[sample]['r_processed']
                        r_with_adapters = stats_data[sample]['r_with_adapters']
                        
                        if r_processed > 0:
                            adapter_ratio = r_with_adapters / r_processed
                            sample_ratios.append(adapter_ratio)
                
                # Calculate the average ratio for this sample (combining R1 and R2 for paired-end)
                if sample_ratios:
                    avg_ratio = sum(sample_ratios) / len(sample_ratios)
                    sample_adapter_ratios[sample] = avg_ratio
                    
                    # Check if below threshold
                    if avg_ratio < threshold:
                        low_adapter_samples.append(f"{sample} ({avg_ratio:.2%})")
            
            # Get all adapter ratios as a flat list for statistics
            all_adapter_ratios = list(sample_adapter_ratios.values())
            
            # Calculate statistics if we have data
            if all_adapter_ratios:
                mean_adapter_ratio = sum(all_adapter_ratios) / len(all_adapter_ratios)
                median_adapter_ratio = np.median(all_adapter_ratios)
                min_adapter_ratio = min(all_adapter_ratios)
                max_adapter_ratio = max(all_adapter_ratios)
                stdev_adapter_ratio = np.std(all_adapter_ratios) if len(all_adapter_ratios) > 1 else 0
                stats_message = f"Mean: {mean_adapter_ratio:.2%}; Median: {median_adapter_ratio:.2%}; Range: {min_adapter_ratio:.2%} - {max_adapter_ratio:.2%}; StdDev: {stdev_adapter_ratio:.2%}"
            else:
                stats_message = "No adapter ratio data available"
            
            # Report results
            if low_adapter_samples:
                print(f"WARNING: Found {len(low_adapter_samples)} samples with low adapter presence (<{threshold:.2%}):")
                for sample in low_adapter_samples:
                    print(f"  - {sample}")
                print(f"Adapter presence stats: {stats_message}")
                if all_adapter_ratios:
                    log_check_result(
                        log_path, "trimmed_reads", "all", "check_adapters_presence", "RED",
                        f"Adapter content: {min_adapter_ratio:.2%} - {max_adapter_ratio:.2%} (below threshold)",
                        f"{stats_message}; From {len(sample_adapter_ratios)} samples"
                    )
                else:
                    log_check_result(
                        log_path, "trimmed_reads", "all", "check_adapters_presence", "RED",
                        "No adapter ratio data available",
                        f"{stats_message}; From {len(sample_adapter_ratios)} samples"
                    )
                return False
            
            print(f"All {len(sample_adapter_ratios)} samples had adapter content above threshold (≥{threshold:.2%}) in raw reads")
            print(f"Adapter presence stats: {stats_message}")
            if all_adapter_ratios:
                log_check_result(
                    log_path, "trimmed_reads", "all", "check_adapters_presence", "GREEN",
                    f"Adapter content: {min_adapter_ratio:.2%} - {max_adapter_ratio:.2%}",
                    f"{stats_message}; From {len(sample_adapter_ratios)} samples"
                )
            else:
                log_check_result(
                    log_path, "trimmed_reads", "all", "check_adapters_presence", "GREEN",
                    "No adapter ratio data available",
                    f"{stats_message}; From {len(sample_adapter_ratios)} samples"
                )
            return True
            
        except Exception as e:
            print(f"ERROR checking adapter presence: {str(e)}")
            log_check_result(log_path, "trimmed_reads", "all", "check_adapters_presence", "RED", 
                           f"Failed to analyze adapter presence: {str(e)}", 
                           "Check trimming reports and MultiQC data for completeness")
            return False

def report_read_depth_stats(multiqc_data, log_path):
    """Report min, max, and median read depths from MultiQC data."""
    if not multiqc_data:
        print("No MultiQC data to analyze read depths")
        log_check_result(log_path, "trimmed_reads", "all", "read_depth_stats", "YELLOW", 
                        "No data available", "Cannot report read depth statistics")
        return

    # Extract read counts for all samples
    read_counts = []
    for sample, data in multiqc_data.items():
        # The field name in the processed data has a prefix and suffix
        if 'trimmed_total_sequences_f' in data:
            read_counts.append((sample, data['trimmed_total_sequences_f']))
    
    if not read_counts:
        print("No read count data found in MultiQC stats")
        # Print all available keys for debugging
        if multiqc_data:
            first_sample = next(iter(multiqc_data.keys()))
            print(f"Available fields for sample '{first_sample}': {sorted(multiqc_data[first_sample].keys())}")
        log_check_result(log_path, "trimmed_reads", "all", "read_depth_stats", "YELLOW", 
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
    
    log_check_result(log_path, "trimmed_reads", "all", "read_depth_stats", "GREEN", message, details)

def report_duplication_rate_stats(multiqc_data, log_path):
    """Report min, max, and median duplication rates from MultiQC data."""
    if not multiqc_data:
        print("No MultiQC data to analyze duplication rates")
        log_check_result(log_path, "trimmed_reads", "all", "duplication_rate_stats", "YELLOW", 
                        "No data available", "Cannot report duplication rate statistics")
        return

    # Extract duplication rates for all samples
    duplication_rates = []
    for sample, data in multiqc_data.items():
        if 'trimmed_percent_duplicates_f' in data:
            duplication_rates.append((sample, data['trimmed_percent_duplicates_f']))
    
    if not duplication_rates:
        print("No duplication rate data found in MultiQC stats")
        # Print all available keys for debugging
        if multiqc_data:
            first_sample = next(iter(multiqc_data.keys()))
            print(f"Available fields for sample '{first_sample}': {sorted(multiqc_data[first_sample].keys())}")
        log_check_result(log_path, "trimmed_reads", "all", "duplication_rate_stats", "YELLOW", 
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
    
    log_check_result(log_path, "trimmed_reads", "all", "duplication_rate_stats", "GREEN", message, details)

def report_gc_content_stats(multiqc_data, log_path):
    """Report GC content statistics from MultiQC data, focusing on peak values and outliers."""
    if not multiqc_data:
        print("No MultiQC data to analyze GC content")
        log_check_result(log_path, "trimmed_reads", "all", "gc_content_stats", "YELLOW", 
                        "No data available", "Cannot report GC content statistics")
        return

    # Extract GC content percentages for all samples
    gc_contents = []
    for sample, data in multiqc_data.items():
        if 'trimmed_percent_gc_f' in data:
            gc_contents.append((sample, data['trimmed_percent_gc_f']))
    
    if not gc_contents:
        print("No GC content data found in MultiQC stats")
        # Print all available keys for debugging
        if multiqc_data:
            first_sample = next(iter(multiqc_data.keys()))
            print(f"Available fields for sample '{first_sample}': {sorted(multiqc_data[first_sample].keys())}")
        log_check_result(log_path, "trimmed_reads", "all", "gc_content_stats", "YELLOW", 
                        "No GC content data found", "")
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
    print(f"GC content statistics:")
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
    
    log_check_result(log_path, "trimmed_reads", "all", "gc_content_stats", "GREEN", message, details)

def main():
    """Main function to process runsheet and validate trimmed reads."""
    parser = argparse.ArgumentParser(description='Validate trimmed reads based on runsheet information.')
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
    
    # 1. Check for trimmed FASTQ files existence
    check_trimmed_fastq_existence(args.outdir, sample_names, paired_end_values, vv_log_path)
    
    # 2. Check GZIP integrity of FASTQ files
    check_gzip_integrity(args.outdir, sample_names, paired_end_values, vv_log_path)
    
    # 3. Validate FASTQ format
    validate_fastq_format(args.outdir, sample_names, paired_end_values, vv_log_path)
    
    # 4. Check FastQC outputs existence
    check_trimmed_fastqc_existence(args.outdir, sample_names, paired_end_values, vv_log_path)
    
    # 5. Check all samples are in MultiQC report
    check_samples_multiqc(args.outdir, sample_names, paired_end_values, vv_log_path, args.assay_suffix)

    # 6. Get MultiQC stats
    multiqc_data = get_trimmed_multiqc_stats(args.outdir, sample_names, paired_end_values, vv_log_path, args.assay_suffix)
    
    # 7. Report MultiQC stats outliers
    if multiqc_data:
        report_multiqc_outliers(args.outdir, multiqc_data, vv_log_path)
        
        # 8. Check paired read counts match (for paired-end data)
        check_paired_read_counts(multiqc_data, vv_log_path, paired_end_values)
    
    # 9. Check for trimming report files existence
    # Trimming reports are now incorporated into the combined MultiQC
    # Rather than calling check_trimming_report_existence, add a dummy passing check
    # print(f"Skipping separate trimming report checks - now incorporated into trimmed MultiQC")
    # log_check_result(vv_log_path, "trimmed_reads", "all", "check_trimming_report_existence", "GREEN", 
    #                 "Trimming reports verified", "")

    # 10. Skip check_trimming_multiqc_samples - also combined into main MultiQC
    # print(f"Skipping trimming MultiQC samples check - now part of main trimmed MultiQC")
    # log_check_result(vv_log_path, "trimmed_reads", "all", "check_trimming_multiqc_samples", "GREEN", 
    #                 "All samples found in MultiQC", "")

    # 11. Check >0.1% of reads had adapters present
    check_adapters_presence(args.outdir, sample_names, paired_end_values, vv_log_path, assay_suffix=args.assay_suffix)

    # 12. Report read depth stats
    if multiqc_data:
        report_read_depth_stats(multiqc_data, vv_log_path)
        
    # 13. Report duplication rate stats
    if multiqc_data:
        report_duplication_rate_stats(multiqc_data, vv_log_path)
        
    # 14. Report GC content stats
    if multiqc_data:
        report_gc_content_stats(multiqc_data, vv_log_path)

if __name__ == "__main__":
    main()