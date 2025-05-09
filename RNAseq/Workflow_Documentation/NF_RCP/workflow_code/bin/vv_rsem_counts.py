#!/usr/bin/env python

"""
Script to validate and verify RSEM count outputs based on runsheet information.
Expects to run from inside a output directory GLDS-##.

Parse the input runsheet to get:
Sample Name
Paired End
Has ERCC

Check that the expected output directories exist

Section-specific checks:
- check_rsem_output_existence: Check if all expected RSEM output files exist for each sample
- check_rsem_multiqc_stats: Extract and validate RSEM metrics from MultiQC
- report_rsem_outliers: Detect and report outliers in RSEM metrics
- check_all_samples_in_multiqc: Verify all samples are included in the MultiQC report
- add_rsem_group_stats: Log groupwise statistics for RSEM metrics
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
        df = pd.read_csv(runsheet_path)
        required_columns = ['Sample Name', 'paired_end', 'has_ERCC', 'organism']
        for col in required_columns:
            if col not in df.columns:
                print(f"Error: Required column '{col}' not found in runsheet")
                sys.exit(1)
        
        # Convert boolean-like columns to actual booleans
        if 'paired_end' in df.columns:
            df['paired_end'] = df['paired_end'].astype(bool)
        if 'has_ERCC' in df.columns:
            df['has_ERCC'] = df['has_ERCC'].astype(bool)
            
        return df
    except Exception as e:
        print(f"Error parsing runsheet: {str(e)}")
        sys.exit(1)


def check_directory_structure(outdir):
    """Check if the expected directory structure exists."""
    expected_dirs = [
        '03-RSEM_Counts'
    ]
    
    missing_dirs = []
    for dir_path in expected_dirs:
        full_path = os.path.join(outdir, dir_path)
        if not os.path.exists(full_path):
            missing_dirs.append(dir_path)
    
    if missing_dirs:
        print(f"Warning: The following expected directories are missing:")
        for dir_path in missing_dirs:
            print(f"  - {dir_path}")
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


def check_rsem_output_existence(outdir, samples, log_path, assay_suffix="_GLbulkRNAseq"):
    """Check if all expected RSEM output files exist for each sample."""
    rsem_dir = os.path.join(outdir, '03-RSEM_Counts')
    
    # Expected file patterns for each sample in sample-specific subdirectories
    expected_patterns = [
        "{sample}/{sample}.genes.results",
        "{sample}/{sample}.isoforms.results",
        "{sample}/{sample}.stat/{sample}.cnt",
        "{sample}/{sample}.stat/{sample}.model",
        "{sample}/{sample}.stat/{sample}.theta"
    ]
    
    # Dataset-level files (directly in the RSEM directory)
    dataset_files = [
        "RSEM_NumNonZeroGenes{assay_suffix}.csv",
        "RSEM_Unnormalized_Counts{assay_suffix}.csv"
    ]
    
    # MultiQC files in the MultiQC_Reports subdirectory
    multiqc_files = [
        "RSEM_count_multiqc{assay_suffix}_data.zip",
        "RSEM_count_multiqc{assay_suffix}.html"
    ]
    
    missing_files_by_sample = {}
    
    # Check sample-specific files
    for sample in samples:
        missing_files = []
        for pattern in expected_patterns:
            file_path = os.path.join(rsem_dir, pattern.format(sample=sample))
            if not os.path.exists(file_path):
                missing_files.append(os.path.basename(file_path))
        
        if missing_files:
            missing_files_by_sample[sample] = missing_files
    
    # Check dataset-level files
    missing_dataset_files = []
    for file_name in dataset_files:
        # Format with assay_suffix if present in the pattern
        if "{assay_suffix}" in file_name:
            file_path = os.path.join(rsem_dir, file_name.format(assay_suffix=assay_suffix))
        else:
            file_path = os.path.join(rsem_dir, file_name)
        
        if not os.path.exists(file_path):
            missing_dataset_files.append(os.path.basename(file_path))

    # Check MultiQC files
    multiqc_dir = os.path.join(rsem_dir, "MultiQC_Reports")
    missing_multiqc_files = []
    for file_name in multiqc_files:
        file_path = os.path.join(multiqc_dir, file_name.format(assay_suffix=assay_suffix))
        if not os.path.exists(file_path):
            missing_multiqc_files.append(os.path.basename(file_path))

    # Log results
    if missing_files_by_sample or missing_dataset_files or missing_multiqc_files:
        # Log sample-specific missing files
        for sample, missing_files in missing_files_by_sample.items():
            log_check_result(
                log_path, 
                "RSEM_counts", 
                sample, 
                "check_rsem_output_existence", 
                "HALT", 
                f"Missing {len(missing_files)} expected RSEM output files", 
                ",".join(missing_files)
            )
        
        # Log dataset-level missing files
        if missing_dataset_files:
            log_check_result(
                log_path, 
                "RSEM_counts", 
                "all", 
                "check_rsem_output_existence", 
                "HALT", 
                f"Missing {len(missing_dataset_files)} expected dataset-level RSEM output files", 
                ",".join(missing_dataset_files)
            )
        
        # Log MultiQC files
        if missing_multiqc_files:
            log_check_result(
                log_path,
                "RSEM_counts",
                "all",
                "check_rsem_output_existence",
                "HALT",
                f"Missing {len(missing_multiqc_files)} expected MultiQC files",
                ",".join(missing_multiqc_files)
            )
        
        print(f"WARNING: Some expected RSEM output files are missing")
        return False
    else:
        log_check_result(
            log_path, 
            "RSEM_counts", 
            "all", 
            "check_rsem_output_existence", 
            "GREEN", 
            "All expected RSEM output files exist", 
            ""
        )
        print("All expected RSEM output files exist")
        return True


def parse_rsem(multiqc_data_dir, assay_suffix="_GLbulkRNAseq"):
    """Parse RSEM data from MultiQC data directory."""
    multiqc_data_json = os.path.join(multiqc_data_dir, "multiqc_data.json")
    
    if not os.path.exists(multiqc_data_json):
        print(f"WARNING: MultiQC data file not found: {multiqc_data_json}")
        return {}
    
    try:
        with open(multiqc_data_json) as f:
            j = json.loads(f.read())
        
        data = {}
        
        if 'report_saved_raw_data' not in j or 'multiqc_rsem' not in j['report_saved_raw_data']:
            print("WARNING: No RSEM data found in MultiQC report")
            return {}
        
        for sample, count_data in j['report_saved_raw_data']['multiqc_rsem'].items():
            sample_name = sample
            # Clean up sample name if needed
            if "/" in sample:
                sample_name = sample.split("/")[-1]
            
            total_reads = count_data['Unique'] + count_data['Multi'] + count_data['Filtered'] + count_data['Unalignable']
            
            data[sample_name] = {
                'num_uniquely_aligned': count_data['Unique'],
                'pct_uniquely_aligned': count_data['Unique'] / total_reads * 100 if total_reads > 0 else 0,
                'pct_multi_aligned': count_data['Multi'] / total_reads * 100 if total_reads > 0 else 0,
                'pct_filtered': count_data['Filtered'] / total_reads * 100 if total_reads > 0 else 0,
                'pct_unalignable': count_data['Unalignable'] / total_reads * 100 if total_reads > 0 else 0,
                'total_reads': total_reads
            }
        
        return data
    
    except Exception as e:
        print(f"ERROR parsing RSEM data: {str(e)}")
        return {}


def get_rsem_multiqc_stats(outdir, samples, log_path, assay_suffix="_GLbulkRNAseq"):
    """Extract RSEM metrics from MultiQC data."""
    rsem_dir = os.path.join(outdir, '03-RSEM_Counts')
    multiqc_dir = os.path.join(rsem_dir, "MultiQC_Reports")
    multiqc_zip = os.path.join(multiqc_dir, f"RSEM_count_multiqc{assay_suffix}_data.zip")
    
    if not os.path.exists(multiqc_zip):
        print(f"WARNING: MultiQC data zip file not found: {multiqc_zip}")
        log_check_result(
            log_path, 
            "RSEM_counts", 
            "all", 
            "get_rsem_multiqc_stats", 
            "HALT", 
            "RSEM MultiQC data not found", 
            f"Expected at {multiqc_zip}"
        )
        return False
    
    print(f"Extracting RSEM stats from MultiQC data: {multiqc_zip}")
    
    # Create a temporary directory to extract files
    with tempfile.TemporaryDirectory() as temp_dir:
        try:
            # Extract the zip file
            with zipfile.ZipFile(multiqc_zip, 'r') as zip_ref:
                zip_ref.extractall(temp_dir)
            
            # Path to the MultiQC data JSON file (new structure)
            json_path = os.path.join(temp_dir, f"RSEM_count_multiqc{assay_suffix}_data", "multiqc_data.json")
            
            if not os.path.exists(json_path):
                print("WARNING: multiqc_data.json not found in the expected location")
                log_check_result(
                    log_path, 
                    "RSEM_counts", 
                    "all", 
                    "get_rsem_multiqc_stats", 
                    "HALT", 
                    "multiqc_data.json not found in MultiQC zip", 
                    f"Expected at {json_path}"
                )
                return False
            
            # Use the directory containing the JSON file for parse_rsem
            data_dir = os.path.dirname(json_path)
            
            # Parse RSEM data
            rsem_data = parse_rsem(data_dir, assay_suffix)
            
            if not rsem_data:
                print("WARNING: No RSEM data found in MultiQC data")
                log_check_result(
                    log_path, 
                    "RSEM_counts", 
                    "all", 
                    "get_rsem_multiqc_stats", 
                    "RED", 
                    "No RSEM data found in MultiQC data", 
                    ""
                )
                return False
            
            # Check if all samples are in the stats
            missing_samples = []
            for sample in samples:
                if sample not in rsem_data:
                    missing_samples.append(sample)
            
            if missing_samples:
                print(f"WARNING: {len(missing_samples)} samples missing from RSEM stats:")
                for sample in missing_samples[:10]:  # Show first 10 only
                    print(f"  - {sample}")
                if len(missing_samples) > 10:
                    print(f"  ... and {len(missing_samples) - 10} more")
                    
                log_check_result(
                    log_path, 
                    "RSEM_counts", 
                    "all", 
                    "get_rsem_multiqc_stats", 
                    "YELLOW", 
                    f"Missing {len(missing_samples)} samples in RSEM stats", 
                    ",".join(missing_samples[:20])  # Limit to 20 sample names
                )
            else:
                print("All samples found in RSEM stats")
                log_check_result(
                    log_path, 
                    "RSEM_counts", 
                    "all", 
                    "get_rsem_multiqc_stats", 
                    "GREEN", 
                    "All samples found in RSEM stats", 
                    ""
                )
            
            return rsem_data
            
        except Exception as e:
            print(f"ERROR extracting RSEM stats: {str(e)}")
            log_check_result(
                log_path, 
                "RSEM_counts", 
                "all", 
                "get_rsem_multiqc_stats", 
                "RED", 
                f"Error extracting RSEM stats: {str(e)}", 
                ""
            )
            return False


def report_rsem_outliers(outdir, rsem_data, log_path):
    """Identify and report outliers in RSEM metrics."""
    if not rsem_data:
        print("No RSEM data to analyze for outliers")
        log_check_result(
            log_path, 
            "RSEM_counts", 
            "all", 
            "report_rsem_outliers", 
            "RED", 
            "No RSEM data to analyze for outliers", 
            ""
        )
        return False
    
    # Metrics to check for outliers
    metrics_to_check = {
        "pct_uniquely_aligned": "uniquely aligned reads percentage",
        "pct_multi_aligned": "multi-mapped reads percentage",
        "pct_filtered": "filtered reads percentage",
        "pct_unalignable": "unalignable reads percentage"
    }
    
    # Thresholds for outlier detection, based on dp_tools config
    thresholds = [
        {"code": "YELLOW", "stdev_threshold": 2, "middle_fcn": "median"},
        {"code": "RED", "stdev_threshold": 4, "middle_fcn": "median"}
    ]
    
    # Set to keep track of outlier samples
    outlier_samples = set()
    outlier_details = []
    
    # Check each metric for outliers
    for metric_key, metric_description in metrics_to_check.items():
        # Collect values for this metric across all samples
        values = []
        for sample, sample_data in rsem_data.items():
            if metric_key in sample_data:
                values.append((sample, sample_data[metric_key]))
        
        if not values:
            continue
        
        # Sort by sample name for consistent reporting
        values.sort(key=lambda x: x[0])
        
        # Calculate median and standard deviation
        values_only = [v[1] for v in values]
        median_value = np.median(values_only)
        stdev = np.std(values_only)
        
        # Print metric summary
        print(f"\nAnalyzing {metric_description} ({metric_key}):")
        print(f"  Median: {median_value:.2f}%")
        print(f"  Standard Deviation: {stdev:.2f}")
        
        # Check for outliers based on thresholds
        for sample, value in values:
            for threshold in thresholds:
                stdev_threshold = threshold["stdev_threshold"]
                code = threshold["code"]
                
                # Skip if standard deviation is zero (all values are the same)
                if stdev == 0:
                    continue
                
                # Calculate deviation from the median in terms of standard deviations
                deviation = abs(value - median_value)
                stdev_multiples = deviation / stdev
                
                # Check if sample is an outlier based on threshold
                if stdev_multiples >= stdev_threshold:
                    outlier_samples.add(sample)
                    # Format detail with semicolons instead of commas to avoid quotes in CSV
                    detail = f"{sample}: {metric_description}={value:.2f}%; {stdev_multiples:.2f} stdevs from median"
                    outlier_details.append(detail)
                    
                    print(f"  {code} outlier: {detail}")
                    log_check_result(
                        log_path, 
                        "RSEM_counts", 
                        sample, 
                        f"outlier_{metric_key}", 
                        code, 
                        f"{stdev_multiples:.2f} stdevs from median", 
                        f"value={value:.2f}%; median={median_value:.2f}%; stdev={stdev:.2f}"
                    )
                    break  # Once we've identified the highest threshold, we can stop
    
    # Summarize results
    if outlier_samples:
        print(f"\nFound {len(outlier_samples)} samples with outlier RSEM metrics:")
        for sample in sorted(outlier_samples)[:10]:  # Show first 10 only
            print(f"  - {sample}")
        if len(outlier_samples) > 10:
            print(f"  ... and {len(outlier_samples) - 10} more")
            
        log_check_result(
            log_path, 
            "RSEM_counts", 
            "all", 
            "report_rsem_outliers", 
            "YELLOW", 
            f"Found {len(outlier_samples)} samples with outlier metrics", 
            ";".join(outlier_details[:20])  # Limit to 20 outlier details
        )
        return True
    else:
        print("\nNo outliers detected in RSEM metrics")
        log_check_result(
            log_path, 
            "RSEM_counts", 
            "all", 
            "report_rsem_outliers", 
            "GREEN", 
            "No outliers detected in RSEM metrics", 
            ""
        )
        return True


def check_all_samples_in_multiqc(outdir, samples, log_path, assay_suffix="_GLbulkRNAseq"):
    """Check if all samples are present in the RSEM MultiQC data."""
    rsem_dir = os.path.join(outdir, '03-RSEM_Counts')
    multiqc_dir = os.path.join(rsem_dir, "MultiQC_Reports")
    multiqc_zip = os.path.join(multiqc_dir, f"RSEM_count_multiqc{assay_suffix}_data.zip")
    
    if not os.path.exists(multiqc_zip):
        print(f"WARNING: MultiQC data zip file not found: {multiqc_zip}")
        log_check_result(
            log_path, 
            "RSEM_counts", 
            "all", 
            "check_all_samples_in_multiqc", 
            "HALT", 
            "RSEM MultiQC data not found", 
            f"Expected at {multiqc_zip}"
        )
        return False
    
    # Create a temporary directory to extract files
    with tempfile.TemporaryDirectory() as temp_dir:
        try:
            # Extract the zip file
            with zipfile.ZipFile(multiqc_zip, 'r') as zip_ref:
                zip_ref.extractall(temp_dir)
            
            # Path to the expected data directory
            data_dir = os.path.join(temp_dir, f"RSEM_count_multiqc{assay_suffix}_data")
            
            # Check for the multiqc_sources.txt file in the expected location
            sources_file = os.path.join(data_dir, "multiqc_sources.txt")
            
            if not os.path.exists(sources_file):
                print("WARNING: multiqc_sources.txt not found in the expected location")
                log_check_result(
                    log_path, 
                    "RSEM_counts", 
                    "all", 
                    "check_all_samples_in_multiqc", 
                    "HALT", 
                    "multiqc_sources.txt not found in MultiQC zip", 
                    f"Expected at {sources_file}"
                )
                return False
            
            # Path to the MultiQC data JSON file in the expected location
            data_json_file = os.path.join(data_dir, "multiqc_data.json")
            
            # Get the sample names used in the RSEM module if possible
            rsem_module_samples = []
            if os.path.exists(data_json_file):
                try:
                    with open(data_json_file, 'r') as f:
                        mqc_data = json.load(f)
                        if 'report_saved_raw_data' in mqc_data and 'multiqc_rsem' in mqc_data['report_saved_raw_data']:
                            rsem_module_samples = list(mqc_data['report_saved_raw_data']['multiqc_rsem'].keys())
                            print(f"Found {len(rsem_module_samples)} samples in RSEM module data")
                except Exception as e:
                    print(f"Warning: Could not parse MultiQC data.json: {e}")
            
            # Read the sources file to get the samples included in the MultiQC report
            multiqc_samples = []
            with open(sources_file, 'r') as f:
                # Skip the header line
                next(f)
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        sample_name = parts[0]
                        multiqc_samples.append(sample_name)
            
            # Combine both sources of sample names
            all_multiqc_samples = list(set(multiqc_samples + rsem_module_samples))
            
            # Debug output to show what's in the MultiQC report            
            print(f"Found {len(all_multiqc_samples)} total sample entries in MultiQC:")
            for sample in all_multiqc_samples[:5]:  # Show just a few examples
                print(f"  - {sample}")
            if len(all_multiqc_samples) > 5:
                print(f"  ... and {len(all_multiqc_samples) - 5} more")
            
            # Helper function to extract core sample name from various formats
            def extract_core_sample_name(name):
                # Strategy 1: Extract just the sample ID by removing extensions and paths
                base_name = os.path.basename(name)
                # Remove common extensions
                for ext in ['.genes.results', '.isoforms.results', '.stat', '.bam', '.fastq.gz', '.fastq']:
                    base_name = base_name.replace(ext, '')
                
                # Strategy 2: Try to match sample IDs with common prefixes/patterns
                # Look for patterns like "RR23_LVR_FLT_F1" in the string
                sample_pattern = re.search(r'(RR\d+_[A-Za-z]+_[A-Za-z]+_[A-Za-z]\d+(?:_techrep\d+)?)', name)
                if sample_pattern:
                    return sample_pattern.group(1)
                
                return base_name
            
            # More flexible matching using multiple strategies
            missing_samples = []
            matched_samples = {}  # sample -> matching multiqc name
            
            for sample in samples:
                found = False
                # Try to find the sample name as a substring in any of the multiqc sample names
                for multiqc_sample in all_multiqc_samples:
                    # Try different ways to match:
                    # 1. Direct match after path stripping
                    base_multiqc_sample = os.path.basename(multiqc_sample)
                    
                    # 2. Sample name is a substring of multiqc sample name
                    if (sample == base_multiqc_sample or 
                        sample in multiqc_sample or
                        extract_core_sample_name(multiqc_sample) == sample):
                        found = True
                        matched_samples[sample] = multiqc_sample
                        break
                
                if not found:
                    missing_samples.append(sample)
            
            if missing_samples:
                print(f"WARNING: {len(missing_samples)} samples missing from RSEM MultiQC report:")
                for sample in missing_samples[:10]:  # Show first 10 only
                    print(f"  - {sample}")
                if len(missing_samples) > 10:
                    print(f"  ... and {len(missing_samples) - 10} more")
                
                # Show examples of what's in the MultiQC report for debugging
                print(f"Examples of sample names in MultiQC report:")
                for sample in list(all_multiqc_samples)[:5]:
                    print(f"  - {sample}")
                    if sample in multiqc_samples and sample in rsem_module_samples:
                        print(f"    (found in both sources.txt and RSEM module)")
                    elif sample in multiqc_samples:
                        print(f"    (found in sources.txt only)")
                    else:
                        print(f"    (found in RSEM module only)")
                
                log_check_result(
                    log_path, 
                    "RSEM_counts", 
                    "all", 
                    "check_all_samples_in_multiqc", 
                    "YELLOW",  # Changed to YELLOW since we can proceed even with this warning
                    f"Missing {len(missing_samples)} samples in RSEM MultiQC report", 
                    f"Missing: {';'.join(missing_samples[:20])}"  # Limit to 20 sample names
                )
                # Return True so we can continue with validation
                return True
            else:
                print(f"All samples found in RSEM MultiQC report")
                print(f"Successfully matched {len(matched_samples)} runsheet samples to MultiQC entries")
                
                log_check_result(
                    log_path, 
                    "RSEM_counts", 
                    "all", 
                    "check_all_samples_in_multiqc", 
                    "GREEN", 
                    "All samples found in RSEM MultiQC report", 
                    ""
                )
                return True
                
        except Exception as e:
            print(f"ERROR checking samples in RSEM MultiQC report: {str(e)}")
            log_check_result(
                log_path, 
                "RSEM_counts", 
                "all", 
                "check_all_samples_in_multiqc", 
                "RED", 
                f"Error checking samples in RSEM MultiQC report: {str(e)}", 
                ""
            )
            return False


def add_rsem_group_stats(outdir, rsem_data, log_path):
    """Calculate and log summary statistics for RSEM metrics across all samples."""
    if not rsem_data or len(rsem_data) == 0:
        print("No RSEM data to analyze for group statistics")
        log_check_result(
            log_path, 
            "RSEM_counts", 
            "all", 
            "add_rsem_group_stats", 
            "RED", 
            "No RSEM data to analyze for group statistics", 
            ""
        )
        return False
    
    # Metrics to analyze
    metrics = {
        "num_uniquely_aligned": "uniquely aligned reads",
        "pct_uniquely_aligned": "uniquely aligned reads percentage",
        "pct_multi_aligned": "multi-mapped reads percentage",
        "pct_filtered": "filtered reads percentage",
        "pct_unalignable": "unalignable reads percentage",
        "total_reads": "total reads processed"
    }
    
    # Collect values for each metric
    metric_values = {metric: [] for metric in metrics.keys()}
    
    for sample, data in rsem_data.items():
        for metric in metrics.keys():
            if metric in data:
                metric_values[metric].append(data[metric])
    
    # Calculate statistics for each metric
    stats_summary = {}
    
    for metric, values in metric_values.items():
        if not values:
            continue
            
        values = np.array(values)
        stats_summary[metric] = {
            "min": np.min(values),
            "max": np.max(values),
            "mean": np.mean(values),
            "median": np.median(values),
            "stddev": np.std(values),
            "count": len(values)
        }
    
    # Print and log the statistics
    print("\nRSEM Metrics Summary Statistics:")
    
    # No summary header line - each metric will be logged individually
    
    for metric, description in metrics.items():
        if metric not in stats_summary:
            continue
            
        stats = stats_summary[metric]
        
        # Format appropriately based on metric type
        if metric.startswith("pct_") or metric.endswith("_percent"):
            # Percentage metrics - use semicolons instead of commas to avoid CSV quoting
            detail_str = (f"Range: {stats['min']:.2f}% - {stats['max']:.2f}%; "
                  f"Median: {stats['median']:.2f}%; "
                  f"Mean: {stats['mean']:.2f}%; "
                  f"StdDev: {stats['stddev']:.2f}")
            
            print(f"  {description}: {detail_str}")
            
            # Log each metric as a separate entry
            log_check_result(
                log_path, 
                "RSEM_counts", 
                "all", 
                f"rsem_stats_{metric}", 
                "GREEN", 
                f"{description}", 
                detail_str
            )
        else:
            # Count metrics - use semicolons instead of commas to avoid CSV quoting
            detail_str = (f"Range: {stats['min']:.0f} - {stats['max']:.0f}; "
                  f"Median: {stats['median']:.0f}; "
                  f"Mean: {stats['mean']:.1f}; "
                  f"StdDev: {stats['stddev']:.1f}")
            
            print(f"  {description}: {detail_str}")
            
            # Log each metric as a separate entry
            log_check_result(
                log_path, 
                "RSEM_counts", 
                "all", 
                f"rsem_stats_{metric}", 
                "GREEN", 
                f"{description}", 
                detail_str
            )
    
    return True


def main():
    """Main function to process runsheet and validate RSEM count outputs."""
    parser = argparse.ArgumentParser(description='Validate RSEM count outputs based on runsheet information.')
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
    
    # Extract paired_end status
    paired_end = runsheet_df['paired_end'].iloc[0] if len(runsheet_df) > 0 else False
    
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
    
    # Run validation checks
    
    # Check if RSEM output files exist
    check_rsem_output_existence(args.outdir, sample_names, vv_log_path, args.assay_suffix)
    
    # Check if all samples are in the MultiQC report
    check_all_samples_in_multiqc(args.outdir, sample_names, vv_log_path, args.assay_suffix)
    
    # Get RSEM MultiQC stats and look for outliers
    rsem_data = get_rsem_multiqc_stats(args.outdir, sample_names, vv_log_path, args.assay_suffix)
    
    # Add group statistics for RSEM metrics
    if rsem_data:
        add_rsem_group_stats(args.outdir, rsem_data, vv_log_path)
        
        # Report outliers in RSEM metrics
        report_rsem_outliers(args.outdir, rsem_data, vv_log_path)
    
    print("RSEM counts validation complete")


if __name__ == "__main__":
    main() 