#!/usr/bin/env python

import os
import sys
import argparse
import pandas as pd
import re
import gzip
import json
import numpy as np
from collections import defaultdict
import tempfile
import zipfile
import statistics
import glob
import csv
import datetime

def parse_runsheet(runsheet_path):
    """
    Parse the runsheet to extract sample information and dataset metadata.
    
    Args:
        runsheet_path: Path to the runsheet CSV file
        
    Returns:
        tuple: (sample_names, metadata_dict)
            - sample_names: List of sample names
            - metadata_dict: Dictionary of metadata values
    """
    if not os.path.exists(runsheet_path):
        print(f"Error: Runsheet file not found: {runsheet_path}")
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
        
        # Extract sample names
        sample_names = df['Sample Name'].tolist()
        
        # Extract metadata - store raw values, handle type conversion when using
        metadata = {
            'paired_end': df['paired_end'].iloc[0] if len(df['paired_end'].unique()) == 1 else "Mixed",
            'has_ercc': df['has_ERCC'].iloc[0] if len(df['has_ERCC'].unique()) == 1 else "Mixed",
            'organism': df['organism'].iloc[0] if len(df['organism'].unique()) == 1 else "Mixed"
        }
        
        # Check for inconsistencies
        if len(df['paired_end'].unique()) > 1:
            print(f"WARNING: Inconsistent paired_end values: {df['paired_end'].unique()}")
        
        if len(df['has_ERCC'].unique()) > 1:
            print(f"WARNING: Inconsistent has_ERCC values: {df['has_ERCC'].unique()}")
            
        if len(df['organism'].unique()) > 1:
            print(f"WARNING: Inconsistent organism values: {df['organism'].unique()}")
            
        return sample_names, metadata
    except Exception as e:
        print(f"Error parsing runsheet: {e}")
        sys.exit(1)

def check_directory_structure(outdir):
    """Check if the expected directory structure exists in the output directory."""
    rseqc_path = os.path.join(outdir, "RSeQC_Analyses")
    
    # Check for expected subdirectories
    expected_subdirs = [
        os.path.join(rseqc_path, "02_geneBody_coverage"),
        os.path.join(rseqc_path, "03_infer_experiment"),
        os.path.join(rseqc_path, "05_read_distribution")
    ]
    
    missing_dirs = []
    
    if not os.path.exists(rseqc_path):
        missing_dirs.append(str(rseqc_path))
    
    for subdir in expected_subdirs:
        if not os.path.exists(subdir):
            missing_dirs.append(str(subdir))
    
    if missing_dirs:
        print(f"WARNING: The following directories are missing: {', '.join(missing_dirs)}")
        return False
    
    return True

def initialize_vv_log(outdir):
    """Initialize the validation log file."""
    # Create the validation log file directly in the output directory
    log_path = os.path.join(outdir, "VV_log.csv")
    
    # Check if the log file already exists
    if not os.path.exists(log_path):
        # Initialize with header
        with open(log_path, 'w') as f:
            f.write("component,sample_id,check_name,status,message,details\n")
    
    print(f"Validation log initialized at: {log_path}")
    return log_path

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

def check_gene_body_coverage_existence(outdir, samples, log_path):
    """Check if gene body coverage files exist for all samples."""
    component = "rseqc"
    check_name = "check_geneBody_coverage_files"
    print(f"Checking for gene body coverage files...")
    
    # Get the expected directory path
    rseqc_dir = os.path.join(outdir, "RSeQC_Analyses", "02_geneBody_coverage")
    
    # Check if the directory exists
    if not os.path.exists(rseqc_dir):
        print(f"ERROR: Gene body coverage directory not found: {rseqc_dir}")
        log_check_result(log_path, component, "all", check_name, "RED", 
                         "Gene body coverage directory not found", f"Directory: {rseqc_dir}")
        return "RED"
    
    # Check MultiQC report existence
    multiqc_glob = os.path.join(rseqc_dir, "geneBody_cov_multiqc*_report.zip")
    multiqc_files = glob.glob(multiqc_glob)
    
    if not multiqc_files:
        print(f"WARNING: No gene body coverage MultiQC report found in {rseqc_dir}")
        log_check_result(log_path, component, "all", "check_genebody_coverage_multiqc", "RED", 
                         "geneBody coverage MultiQC report not found", f"Searched: {multiqc_glob}")
        status = "RED"
    else:
        print(f"Found gene body coverage MultiQC report: {os.path.basename(multiqc_files[0])}")
        log_check_result(log_path, component, "all", "check_genebody_coverage_multiqc", "GREEN", 
                         "geneBody coverage MultiQC report exists")
    
    # Track missing files
    missing_files = []
    status = "GREEN"
    
    # Check for individual files for each sample
    for sample in samples:
        # Check for .geneBodyCoverage.* files in sample subdirectory
        sample_dir = os.path.join(rseqc_dir, sample)
        coverage_txt = os.path.join(sample_dir, f"{sample}.geneBodyCoverage.txt")
        coverage_r = os.path.join(sample_dir, f"{sample}.geneBodyCoverage.r")
        coverage_curves = os.path.join(sample_dir, f"{sample}.geneBodyCoverage.curves.pdf")
        
        if not os.path.exists(coverage_txt):
            missing_files.append(f"{sample}.geneBodyCoverage.txt")
            
        if not os.path.exists(coverage_r):
            missing_files.append(f"{sample}.geneBodyCoverage.r")
            
        if not os.path.exists(coverage_curves):
            missing_files.append(f"{sample}.geneBodyCoverage.curves.pdf")
    
    # Log result based on missing files
    if missing_files:
        print(f"WARNING: Missing {len(missing_files)} gene body coverage file(s)")
        print("\n".join(missing_files[:5]) + ("..." if len(missing_files) > 5 else ""))
        
        log_check_result(log_path, component, "all", check_name, "RED", 
                         f"Missing {len(missing_files)} gene body coverage file(s)", 
                         "\n".join(missing_files[:20]) + ("..." if len(missing_files) > 20 else ""))
        status = "RED"
    else:
        print(f"All gene body coverage files exist for {len(samples)} samples")
        log_check_result(log_path, component, "all", check_name, "GREEN", 
                         "All geneBody coverage files exist")
        status = "GREEN"
    
    return status

def check_infer_experiment_existence(outdir, samples, log_path):
    """Check if infer experiment files exist for all samples."""
    component = "rseqc"
    check_name = "check_infer_experiment_files"
    print(f"Checking for infer experiment files...")
    
    # Get the expected directory path
    rseqc_dir = os.path.join(outdir, "RSeQC_Analyses", "03_infer_experiment")
    
    # Check if the directory exists
    if not os.path.exists(rseqc_dir):
        print(f"ERROR: Infer experiment directory not found: {rseqc_dir}")
        log_check_result(log_path, component, "all", check_name, "RED", 
                         "Infer experiment directory not found", f"Directory: {rseqc_dir}")
        return "RED"
    
    # Check MultiQC report existence
    multiqc_glob = os.path.join(rseqc_dir, "infer_exp_multiqc*_report.zip")
    multiqc_files = glob.glob(multiqc_glob)
    
    if not multiqc_files:
        print(f"WARNING: No infer experiment MultiQC report found in {rseqc_dir}")
        log_check_result(log_path, component, "all", "check_infer_experiment_multiqc", "RED", 
                         "Infer experiment MultiQC report not found", f"Searched: {multiqc_glob}")
        status = "RED"
    else:
        print(f"Found infer experiment MultiQC report: {os.path.basename(multiqc_files[0])}")
        log_check_result(log_path, component, "all", "check_infer_experiment_multiqc", "GREEN", 
                         "Infer experiment MultiQC report exists")
    
    # Track missing files
    missing_files = []
    status = "GREEN"
    
    # Check for individual files for each sample
    for sample in samples:
        # Check for infer experiment output file with correct pattern
        infer_out = os.path.join(rseqc_dir, f"{sample}.infer_expt.out")
        
        if not os.path.exists(infer_out):
            missing_files.append(f"{sample}.infer_experiment.txt")
    
    # Log result based on missing files
    if missing_files:
        print(f"WARNING: Missing {len(missing_files)} infer experiment file(s)")
        print("\n".join(missing_files[:5]) + ("..." if len(missing_files) > 5 else ""))
        
        log_check_result(log_path, component, "all", check_name, "RED", 
                         f"Missing {len(missing_files)} infer experiment file(s)", 
                         "\n".join(missing_files[:20]) + ("..." if len(missing_files) > 20 else ""))
        status = "RED"
    else:
        print(f"All infer experiment files exist for {len(samples)} samples")
        log_check_result(log_path, component, "all", check_name, "GREEN", 
                         "All infer experiment files exist")
        status = "GREEN"
    
    return status

def check_read_distribution_existence(outdir, samples, log_path):
    """Check if read distribution files exist for all samples."""
    component = "rseqc"
    check_name = "check_read_distribution_files"
    print(f"Checking for read distribution files...")
    
    # Get the expected directory path
    rseqc_dir = os.path.join(outdir, "RSeQC_Analyses", "05_read_distribution")
    
    # Check if the directory exists
    if not os.path.exists(rseqc_dir):
        print(f"ERROR: Read distribution directory not found: {rseqc_dir}")
        log_check_result(log_path, component, "all", check_name, "RED", 
                         "Read distribution directory not found", f"Directory: {rseqc_dir}")
        return "RED"
    
    # Check MultiQC report existence
    multiqc_glob = os.path.join(rseqc_dir, "read_dist_multiqc*_report.zip")
    multiqc_files = glob.glob(multiqc_glob)
    
    if not multiqc_files:
        print(f"WARNING: No read distribution MultiQC report found in {rseqc_dir}")
        log_check_result(log_path, component, "all", "check_read_distribution_multiqc", "RED", 
                         "Read distribution MultiQC report not found", f"Searched: {multiqc_glob}")
        status = "RED"
    else:
        print(f"Found read distribution MultiQC report: {os.path.basename(multiqc_files[0])}")
        log_check_result(log_path, component, "all", "check_read_distribution_multiqc", "GREEN", 
                         "Read distribution MultiQC report exists")
    
    # Track missing files
    missing_files = []
    status = "GREEN"
    
    # Check for individual files for each sample
    for sample in samples:
        # Check for read distribution output file with correct pattern
        read_dist = os.path.join(rseqc_dir, f"{sample}.read_dist.out")
        
        if not os.path.exists(read_dist):
            missing_files.append(f"{sample}.read_distribution.txt")
    
    # Log result based on missing files
    if missing_files:
        print(f"WARNING: Missing {len(missing_files)} read distribution file(s)")
        print("\n".join(missing_files[:5]) + ("..." if len(missing_files) > 5 else ""))
        
        log_check_result(log_path, component, "all", check_name, "RED", 
                         f"Missing {len(missing_files)} read distribution file(s)", 
                         "\n".join(missing_files[:20]) + ("..." if len(missing_files) > 20 else ""))
        status = "RED"
    else:
        print(f"All read distribution files exist for {len(samples)} samples")
        log_check_result(log_path, component, "all", check_name, "GREEN", 
                         "All read distribution files exist")
        status = "GREEN"
    
    return status

def check_inner_distance_existence(outdir, samples, log_path):
    """Check if inner distance files exist for all samples (paired-end only)."""
    component = "rseqc"
    check_name = "check_inner_distance_files"
    print(f"Checking for inner distance files...")
    
    # Get the expected directory path
    rseqc_dir = os.path.join(outdir, "RSeQC_Analyses", "04_inner_distance")
    
    # Check if the directory exists
    if not os.path.exists(rseqc_dir):
        print(f"ERROR: Inner distance directory not found: {rseqc_dir}")
        log_check_result(log_path, component, "all", check_name, "RED", 
                         "Inner distance directory not found", f"Directory: {rseqc_dir}")
        return "RED"
    
    # Check MultiQC report existence
    multiqc_glob = os.path.join(rseqc_dir, "inner_dist_multiqc*_report.zip")
    multiqc_files = glob.glob(multiqc_glob)
    
    if not multiqc_files:
        print(f"WARNING: No inner distance MultiQC report found in {rseqc_dir}")
        log_check_result(log_path, component, "all", "check_inner_distance_multiqc", "RED", 
                         "Inner distance MultiQC report not found", f"Searched: {multiqc_glob}")
        status = "RED"
    else:
        print(f"Found inner distance MultiQC report: {os.path.basename(multiqc_files[0])}")
        log_check_result(log_path, component, "all", "check_inner_distance_multiqc", "GREEN", 
                         "Inner distance MultiQC report exists")
    
    # Track missing files
    missing_files = []
    status = "GREEN"
    
    # Check for individual files for each sample
    for sample in samples:
        # Check for inner distance output files in sample subdirectory with correct patterns
        sample_dir = os.path.join(rseqc_dir, sample)
        inner_dist_txt = os.path.join(sample_dir, f"{sample}.inner_distance.txt")
        inner_dist_freq = os.path.join(sample_dir, f"{sample}.inner_distance_freq.txt")
        inner_dist_plot = os.path.join(sample_dir, f"{sample}.inner_distance_plot.pdf")
        
        if not os.path.exists(inner_dist_txt):
            missing_files.append(f"{sample}.inner_distance.txt")
            
        if not os.path.exists(inner_dist_freq):
            missing_files.append(f"{sample}.inner_distance_freq.txt")
            
        if not os.path.exists(inner_dist_plot):
            missing_files.append(f"{sample}.inner_distance.pdf")
    
    # Log result based on missing files
    if missing_files:
        print(f"WARNING: Missing {len(missing_files)} inner distance file(s)")
        print("\n".join(missing_files[:5]) + ("..." if len(missing_files) > 5 else ""))
        
        log_check_result(log_path, component, "all", check_name, "RED", 
                         f"Missing {len(missing_files)} inner distance file(s)", 
                         "\n".join(missing_files[:20]) + ("..." if len(missing_files) > 20 else ""))
        status = "RED"
    else:
        print(f"All inner distance files exist for {len(samples)} samples")
        log_check_result(log_path, component, "all", check_name, "GREEN", 
                         "All inner distance files exist")
        status = "GREEN"
    
    return status

def get_genebody_coverage_multiqc_stats(outdir, samples, log_path, assay_suffix="_GLbulkRNAseq"):
    """Extract gene body coverage MultiQC stats for all samples and write to a stats file for analysis."""
    rseqc_dir = os.path.join(outdir, "RSeQC_Analyses", "02_geneBody_coverage")
    multiqc_zip = os.path.join(rseqc_dir, f"geneBody_cov_multiqc{assay_suffix}_report.zip")
    
    if not os.path.exists(multiqc_zip):
        print(f"WARNING: Gene body coverage MultiQC report zip file not found: {multiqc_zip}")
        log_check_result(log_path, "rseqc", "all", "get_genebody_coverage_multiqc_stats", "RED", 
                         "Gene body coverage MultiQC report not found", "")
        return None
    
    print(f"Extracting stats from Gene body coverage MultiQC report: {multiqc_zip}")
    
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
                log_check_result(log_path, "rseqc", "all", "get_genebody_coverage_multiqc_stats", "RED", 
                               "multiqc_data.json not found in zip", "")
                return None
            
            multiqc_data_file = json_files[0]
            
            # Parse the MultiQC data
            with open(multiqc_data_file) as f:
                multiqc_data = json.load(f)
            
            genebody_data = {}
            
            # Check if genebody coverage data exists in the MultiQC report
            if ('report_saved_raw_data' not in multiqc_data or 
                'rseqc_gene_body_cov' not in multiqc_data['report_saved_raw_data']):
                print("WARNING: No gene body coverage data found in MultiQC report")
                log_check_result(log_path, "rseqc", "all", "get_genebody_coverage_multiqc_stats", "RED", 
                               "No gene body coverage data in MultiQC report", "")
                return None
            
            # Extract the gene body coverage data from the raw data section
            genebody_raw_data = multiqc_data['report_saved_raw_data']['rseqc_gene_body_cov']
            
            # Process each sample's data
            for sample_name, sample_data in genebody_raw_data.items():
                if isinstance(sample_data, dict):
                    # Convert the data to percentages from 0-100
                    genebody_data[sample_name] = {}
                    total_positions = len(sample_data)
                    
                    # Check if positions are already 1-100 or need conversion
                    positions = [int(pos) for pos in sample_data.keys() if pos.isdigit()]
                    max_position = max(positions) if positions else 0
                    
                    for pos, value in sample_data.items():
                        try:
                            # Handle positions that are already integers
                            pos_int = int(pos)
                            
                            # Calculate position as percentage (0-100)
                            if max_position <= 100:  # Already a percentage format (1-100)
                                pos_pct = pos_int
                            else:  # Convert to percentage based on position in total range
                                pos_pct = int((pos_int / max_position) * 100)
                                
                            # Ensure we have values at each percentage point for consistency
                            genebody_data[sample_name][str(pos_pct)] = value
                        except (ValueError, TypeError):
                            # Skip any non-numeric positions
                            continue
            
            # Check if we found data for all samples
            missing_samples = [s for s in samples if s not in genebody_data]
            
            if missing_samples:
                print(f"WARNING: Missing gene body coverage data for {len(missing_samples)} samples")
                log_check_result(log_path, "rseqc", "all", "get_genebody_coverage_multiqc_stats", "YELLOW", 
                               f"Missing gene body coverage data for {len(missing_samples)} samples", 
                               "\n".join(missing_samples[:20]) + ("..." if len(missing_samples) > 20 else ""))
            else:
                print(f"Found gene body coverage data for all {len(samples)} samples")
                log_check_result(log_path, "rseqc", "all", "get_genebody_coverage_multiqc_stats", "GREEN", 
                               "All samples found in gene body coverage stats")
            
            return genebody_data
            
        except Exception as e:
            print(f"ERROR: Failed to extract gene body coverage stats: {str(e)}")
            log_check_result(log_path, "rseqc", "all", "get_genebody_coverage_multiqc_stats", "RED", 
                           f"Failed to extract gene body coverage stats: {str(e)}", "")
            return None

def report_genebody_coverage_issues(outdir, multiqc_data, log_path):
    """
    Analyze gene body coverage data for issues and report them.
    
    This function looks for:
    1. Unusual 5' to 3' bias
    2. Severe coverage drops
    3. Outliers in overall coverage 
    """
    if multiqc_data is None:
        print("No gene body coverage data to analyze")
        return "RED"
    
    component = "rseqc"
    check_name = "coverage_metrics"
    
    # Print header
    print("\n=== Gene Body Coverage Analysis ===")
    print(f"Analyzing {len(multiqc_data)} samples for coverage issues")
    
    # Calculate metrics for each sample
    sample_metrics = {}
    
    for sample, data in multiqc_data.items():
        # Get coverage at specific gene body positions
        positions = {}
        for pos_str in data.keys():
            if pos_str.isdigit():
                pos = int(pos_str)
                positions[pos] = data[pos_str]
        
        # Skip samples with insufficient data
        if not positions:
            print(f"WARNING: Sample {sample} has no position data, skipping")
            continue
        
        # Calculate metrics
        sample_metrics[sample] = {}
        
        # Get 5' coverage (average of 0-33% positions)
        five_prime_pos = [p for p in positions.keys() if 0 <= p <= 33]
        if five_prime_pos:
            five_prime_coverage = sum(positions[p] for p in five_prime_pos) / len(five_prime_pos)
            sample_metrics[sample]['5_prime'] = five_prime_coverage
        
        # Get middle coverage (average of 34-66% positions)
        middle_pos = [p for p in positions.keys() if 34 <= p <= 66]
        if middle_pos:
            middle_coverage = sum(positions[p] for p in middle_pos) / len(middle_pos)
            sample_metrics[sample]['middle'] = middle_coverage
        
        # Get 3' coverage (average of 67-100% positions)
        three_prime_pos = [p for p in positions.keys() if 67 <= p <= 100]
        if three_prime_pos:
            three_prime_coverage = sum(positions[p] for p in three_prime_pos) / len(three_prime_pos)
            sample_metrics[sample]['3_prime'] = three_prime_coverage
        
        # Calculate 3'/5' ratio - best practice is to use median values at each end
        if five_prime_pos and three_prime_pos:
            # Get values at extreme ends (first 10% and last 10%)
            five_prime_extreme = [positions[p] for p in five_prime_pos if 0 <= p <= 10]
            three_prime_extreme = [positions[p] for p in three_prime_pos if 90 <= p <= 100]
            
            if five_prime_extreme and three_prime_extreme:
                # Use median values at each end for more robust ratio calculation
                five_prime_median = statistics.median(five_prime_extreme)
                three_prime_median = statistics.median(three_prime_extreme)
                
                if five_prime_median > 0:
                    ratio = three_prime_median / five_prime_median
                    sample_metrics[sample]['3_to_5_ratio'] = ratio
                else:
                    sample_metrics[sample]['3_to_5_ratio'] = 0
    
    # Calculate overall statistics for metrics
    all_five_prime = [m['5_prime'] for m in sample_metrics.values() if '5_prime' in m]
    all_middle = [m['middle'] for m in sample_metrics.values() if 'middle' in m]
    all_three_prime = [m['3_prime'] for m in sample_metrics.values() if '3_prime' in m]
    all_ratios = [m['3_to_5_ratio'] for m in sample_metrics.values() if '3_to_5_ratio' in m]
    
    # Calculate statistics if there are enough samples
    if len(all_five_prime) >= 3:
        # Calculate mean and standard deviation for each metric
        mean_five_prime = sum(all_five_prime) / len(all_five_prime)
        stdev_five_prime = statistics.stdev(all_five_prime) if len(all_five_prime) > 1 else 0
        
        mean_middle = sum(all_middle) / len(all_middle) if all_middle else 0
        stdev_middle = statistics.stdev(all_middle) if len(all_middle) > 1 else 0
        
        mean_three_prime = sum(all_three_prime) / len(all_three_prime)
        stdev_three_prime = statistics.stdev(all_three_prime) if len(all_three_prime) > 1 else 0
        
        mean_ratio = sum(all_ratios) / len(all_ratios)
        stdev_ratio = statistics.stdev(all_ratios) if len(all_ratios) > 1 else 0
        
        # Log the overall statistics
        stats_message = (f"5' region (0-33%): {mean_five_prime:.2f} ±{stdev_five_prime:.2f}; "
                        f"Middle (33-67%): {mean_middle:.2f} ±{stdev_middle:.2f}; "
                        f"3' region (67-100%): {mean_three_prime:.2f} ±{stdev_three_prime:.2f}; "
                        f"3'/5' ratio: {mean_ratio:.2f} ±{stdev_ratio:.2f}")
        
        log_check_result(log_path, component, "all", "genebody_coverage_stats", "GREEN", 
                       "Gene body coverage statistics summary (mean ± stdev)", stats_message)
        
        # Print summary statistics
        print("\nGene Body Coverage Summary (mean ± stdev):")
        print(f"5' region (0-33%): {mean_five_prime:.2f} ±{stdev_five_prime:.2f}")
        print(f"Middle (33-67%): {mean_middle:.2f} ±{stdev_middle:.2f}")
        print(f"3' region (67-100%): {mean_three_prime:.2f} ±{stdev_three_prime:.2f}")
        print(f"3'/5' ratio: {mean_ratio:.2f} ±{stdev_ratio:.2f}")
        
        # Identify outliers
        outliers = []
        severe_outliers = []
        
        for sample, metrics in sample_metrics.items():
            issues = []
            
            # Check 5' coverage
            if '5_prime' in metrics and stdev_five_prime > 0:
                zscore = abs(metrics['5_prime'] - mean_five_prime) / stdev_five_prime
                if zscore > 2:
                    if zscore > 3:
                        issues.append(f"Severe 5' coverage outlier: {metrics['5_prime']:.2f} ({zscore:.2f} stdev from median)")
                        severe_outliers.append(sample)
                    else:
                        issues.append(f"5' coverage outlier: {metrics['5_prime']:.2f} ({zscore:.2f} stdev from median)")
                        outliers.append(sample)
            
            # Check 3' coverage
            if '3_prime' in metrics and stdev_three_prime > 0:
                zscore = abs(metrics['3_prime'] - mean_three_prime) / stdev_three_prime
                if zscore > 2:
                    if zscore > 3:
                        issues.append(f"Severe 3' coverage outlier: {metrics['3_prime']:.2f} ({zscore:.2f} stdev from median)")
                        severe_outliers.append(sample)
                    else:
                        issues.append(f"3' coverage outlier: {metrics['3_prime']:.2f} ({zscore:.2f} stdev from median)")
                        outliers.append(sample)
            
            # Check 3'/5' ratio
            if '3_to_5_ratio' in metrics and stdev_ratio > 0:
                zscore = abs(metrics['3_to_5_ratio'] - mean_ratio) / stdev_ratio
                if zscore > 2:
                    if zscore > 3:
                        issues.append(f"Severe 3'/5' ratio outlier: {metrics['3_to_5_ratio']:.2f} ({zscore:.2f} stdev from median)")
                        severe_outliers.append(sample)
                    else:
                        issues.append(f"3'/5' ratio outlier: {metrics['3_to_5_ratio']:.2f} ({zscore:.2f} stdev from median)")
                        outliers.append(sample)
            
            # Log issues for this sample
            if issues:
                status = "RED" if sample in severe_outliers else "YELLOW"
                log_check_result(log_path, component, sample, check_name, status, 
                               f"Unusual gene body coverage profile", "; ".join(issues))
        
        # Log overall result
        if severe_outliers:
            log_check_result(log_path, component, "all", check_name, "RED", 
                           f"Major outliers found in {len(severe_outliers)} samples", 
                           "; ".join(severe_outliers))
            status = "RED"
        elif outliers:
            log_check_result(log_path, component, "all", check_name, "YELLOW", 
                           f"Minor outliers found in {len(outliers)} samples", 
                           "; ".join(outliers))
            status = "YELLOW"
        else:
            log_check_result(log_path, component, "all", check_name, "GREEN", 
                           "No outliers found in gene body coverage metrics")
            status = "GREEN"
    else:
        print("WARNING: Not enough samples to calculate reliable statistics")
        log_check_result(log_path, component, "all", check_name, "YELLOW", 
                       "Not enough samples for reliable statistics", 
                       f"Only {len(all_five_prime)} samples with data")
        status = "YELLOW"
    
    return status

def get_infer_experiment_multiqc_stats(outdir, samples, log_path, assay_suffix="_GLbulkRNAseq"):
    """
    Extract infer experiment stats from MultiQC report.
    Retrieves: pct_sense, pct_antisense, pct_undetermined
    """
    component = "rseqc"
    check_name = "get_infer_experiment_multiqc_stats"
    
    # Path to MultiQC report
    multiqc_zip = os.path.join(outdir, "RSeQC_Analyses", "03_infer_experiment", f"infer_exp_multiqc{assay_suffix}_report.zip")
    
    if not os.path.exists(multiqc_zip):
        print(f"WARNING: Infer experiment MultiQC report not found: {multiqc_zip}")
        log_check_result(log_path, component, "all", check_name, "RED", 
                         "Infer experiment MultiQC report not found", "")
        return None
    
    # Create a temporary directory to extract the zip file
    with tempfile.TemporaryDirectory() as tmpdirname:
        # Extract the zip file
        with zipfile.ZipFile(multiqc_zip, 'r') as zip_ref:
            zip_ref.extractall(tmpdirname)
        
        # Read the JSON data
        json_path = os.path.join(tmpdirname, f"infer_exp_multiqc{assay_suffix}_report", 
                                 f"infer_exp_multiqc{assay_suffix}_data", "multiqc_data.json")
        
        if not os.path.exists(json_path):
            print(f"WARNING: MultiQC data file not found in ZIP: {json_path}")
            log_check_result(log_path, component, "all", check_name, "RED", 
                             "MultiQC data file not found in ZIP", "")
            return None
        
        try:
            with open(json_path, 'r') as f:
                multiqc_data = json.load(f)
                
            # Extract infer experiment data
            if 'report_saved_raw_data' not in multiqc_data or 'multiqc_rseqc_infer_experiment' not in multiqc_data['report_saved_raw_data']:
                print("WARNING: Infer experiment data not found in MultiQC report")
                log_check_result(log_path, component, "all", check_name, "RED", 
                                "Infer experiment data not found in MultiQC report", "")
                return None
            
            # Use the same approach as parse_multiqc.py
            key_dict = {
                'se_sense': 'pct_sense',
                'se_antisense': 'pct_antisense',
                'pe_sense': 'pct_sense',
                'pe_antisense': 'pct_antisense',
                'failed': 'pct_undetermined'
            }
            
            infer_exp_data = multiqc_data['report_saved_raw_data']['multiqc_rseqc_infer_experiment']
            processed_data = {s.replace('_infer_expt', ''): {key_dict[k]: v * 100 for k, v in d.items()} 
                             for s, d in infer_exp_data.items()}
            
            # Check if all samples are present
            missing_samples = [sample for sample in samples if sample not in processed_data]
            
            if missing_samples:
                print(f"WARNING: {len(missing_samples)} samples not found in infer experiment data")
                log_check_result(log_path, component, "all", check_name, "YELLOW", 
                                f"{len(missing_samples)} samples not found in infer experiment data", 
                                ', '.join(missing_samples))
            else:
                log_check_result(log_path, component, "all", check_name, "GREEN", 
                                "All samples found in infer experiment stats", "")
            
            return processed_data
            
        except Exception as e:
            print(f"ERROR: Failed to parse infer experiment MultiQC data: {str(e)}")
            log_check_result(log_path, component, "all", check_name, "RED", 
                             f"Failed to parse infer experiment MultiQC data: {str(e)}", "")
            return None

def report_infer_experiment_issues(outdir, infer_exp_data, log_path):
    """
    Report issues with infer experiment metrics.
    Checks for outliers in:
    - pct_sense (percentage of reads mapped to sense strand)
    - pct_antisense (percentage of reads mapped to antisense strand)
    - pct_undetermined (percentage of reads with undetermined strandedness)
    """
    if infer_exp_data is None or len(infer_exp_data) == 0:
        print("No infer experiment data available, skipping outlier detection")
        return "GREEN"
    
    component = "rseqc"
    check_name = "infer_experiment_metrics"
    
    # Collect values for each metric
    metrics = ['pct_sense', 'pct_antisense', 'pct_undetermined']
    metric_values = {metric: [] for metric in metrics}
    
    for sample, data in infer_exp_data.items():
        for metric in metrics:
            if metric in data:
                metric_values[metric].append((sample, data[metric]))
    
    # Calculate summary statistics for each metric
    summary_stats = {}
    summary_texts = []
    
    for metric in metrics:
        if not metric_values[metric]:
            print(f"No data for {metric}, skipping")
            continue
        
        # Extract values for all samples
        values = [value for _, value in metric_values[metric]]
        
        # Calculate mean and standard deviation
        mean = statistics.mean(values)
        stdev = statistics.stdev(values) if len(values) > 1 else 0
        
        # Store stats
        summary_stats[metric] = {
            'mean': mean,
            'stdev': stdev,
            'median': statistics.median(values)
        }
        
        # Create readable metric name
        readable_metric = metric.replace('pct_', '').capitalize()
        summary_texts.append(f"{readable_metric}: {mean:.2f}% ±{stdev:.2f}")
    
    # Log summary statistics
    if summary_texts:
        log_check_result(log_path, component, "all", "infer_experiment_summary", "GREEN",
                        "Infer experiment metrics summary (mean ± stdev)",
                        "; ".join(summary_texts))
        
        print("\nInfer Experiment Summary (mean ± stdev):")
        for text in summary_texts:
            print(f"{text}")
        print()
    
    # Check each metric for outliers
    any_outliers = False
    
    for metric in metrics:
        if not metric_values[metric] or metric not in summary_stats:
            continue
        
        median = summary_stats[metric]['median']
        stdev = summary_stats[metric]['stdev']
        
        # Don't attempt to calculate outliers if stdev is zero (all values are the same)
        if stdev == 0:
            continue
        
        # Check each sample for outliers
        for sample, value in metric_values[metric]:
            # Calculate how many standard deviations from the median
            deviation = abs(value - median) / stdev
            
            # Flag as outlier if more than 2 or 4 standard deviations from median
            if deviation > 4.0:
                any_outliers = True
                log_check_result(log_path, component, sample, metric, "RED", 
                                f"Outlier {metric}: {value:.2f} ({deviation:.2f} stdev from median {median:.2f})", "")
            elif deviation > 2.0:
                any_outliers = True
                log_check_result(log_path, component, sample, metric, "YELLOW", 
                                f"Possible outlier {metric}: {value:.2f} ({deviation:.2f} stdev from median {median:.2f})", "")
    
    # Log overall result
    if any_outliers:
        log_check_result(log_path, component, "all", check_name, "YELLOW", 
                        "Outliers found in infer experiment metrics", "")
        print("Outliers found in infer experiment metrics, see validation log for details")
        return "YELLOW"
    else:
        log_check_result(log_path, component, "all", check_name, "GREEN", 
                        "No outliers found in infer experiment metrics", "")
        print("No outliers found in infer experiment metrics")
        return "GREEN"

def detect_coverage_bin_outliers(outdir, genebody_data, log_path):
    """
    Detect outliers in the coverage bin data for each sample.
    
    Args:
        outdir: Output directory
        genebody_data: Dictionary containing gene body coverage data
        log_path: Path to the validation log
    
    Returns:
        Status: "RED", "YELLOW", "GREEN", or "SKIPPED" based on outlier detection
    """
    if not genebody_data:
        log_check_result(log_path, "rseqc", "all", "genebody_coverage_bins", "SKIPPED", "No gene body coverage data available")
        return "SKIPPED"
    
    print("\n=== Gene Body Coverage Bin Analysis (5% bins) ===")
    print(f"Analyzing {len(genebody_data)} samples for outliers in 5% coverage bins")
    
    bin_positions = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100]
    
    
    # Collect bin values for all samples
    sample_bin_values = {}
    for sample, data in genebody_data.items():
        sample_bin_values[sample] = {}
        bin_values = []
        for key in data:
            if key.isdigit() or (isinstance(key, int)):
                bin_values.append(data[key])
                
        # Map bin values to positions
        for i, pos in enumerate(bin_positions):
            if i < len(bin_values):
                sample_bin_values[sample][pos] = bin_values[i]
    
    # Calculate statistics for each bin position
    bin_stats = {}
    for pos in bin_positions:
        values = [data[pos] for sample, data in sample_bin_values.items() if pos in data]
        if values:
            bin_stats[pos] = {
                'median': statistics.median(values),
                'mean': statistics.mean(values),
                'stdev': statistics.stdev(values) if len(values) > 1 else 0
            }
    
    # Calculate overall coverage statistics
    all_mean_coverages = []
    for sample, bin_data in sample_bin_values.items():
        if bin_data:
            mean_coverage = statistics.mean(list(bin_data.values()))
            all_mean_coverages.append(mean_coverage)
    
    overall_median = statistics.median(all_mean_coverages) if all_mean_coverages else 0
    overall_stdev = statistics.stdev(all_mean_coverages) if len(all_mean_coverages) > 1 else 0
    
    # Check each sample for outliers
    outlier_samples = set()
    outlier_count = defaultdict(int)
    sample_anomalies = defaultdict(list)
    
    for sample, bin_data in sample_bin_values.items():
        if not bin_data:
            continue
            
        # Check overall coverage level
        mean_coverage = statistics.mean(list(bin_data.values()))
        if overall_stdev > 0:
            # Calculate zscore with sign to indicate direction
            mean_zscore = (mean_coverage - overall_median) / overall_stdev
            
            if abs(mean_zscore) > 2:
                outlier_samples.add(sample)
                outlier_count[sample] += 1
                direction = "higher than" if mean_zscore > 0 else "lower than"
                anomaly_msg = f"Unusual overall coverage: {mean_coverage:.2f} ({mean_zscore:.2f} stdev {direction} median)"
                sample_anomalies[sample].append(("genebody_coverage_bins_mean_coverage", "YELLOW", anomaly_msg))
        
        # Check individual bins
        for i, pos in enumerate(bin_positions):
            if pos in bin_stats and bin_stats[pos]['stdev'] > 0:
                value = bin_data[pos]
                median = bin_stats[pos]['median']
                stdev = bin_stats[pos]['stdev']
                
                # Calculate zscore with sign to indicate direction
                zscore = (value - median) / stdev
                
                if abs(zscore) > 2:
                    outlier_samples.add(sample)
                    outlier_count[sample] += 1
                    direction = "higher than" if zscore > 0 else "lower than"
                    anomaly_msg = f"Possible outlier at {pos}% position: {value:.2f} coverage ({zscore:.2f} stdev {direction} median)"
                    # Include all outliers in sample_anomalies for detailed reporting
                    sample_anomalies[sample].append((f"genebody_coverage_bins_{pos}pct", "YELLOW", anomaly_msg))
            
            # Check for unusual coverage changes between adjacent bins
            if i < len(bin_positions) - 1:
                next_pos = bin_positions[i + 1]
                value_diff = bin_data[next_pos] - bin_data[pos]
                
                # Collect all differences between these positions
                all_diffs = []
                for s, s_data in sample_bin_values.items():
                    if pos in s_data and next_pos in s_data:
                        all_diffs.append(s_data[next_pos] - s_data[pos])
                
                if all_diffs:
                    median_diff = statistics.median(all_diffs)
                    if len(all_diffs) > 1:
                        stdev_diff = statistics.stdev(all_diffs)
                        if stdev_diff > 0:
                            # Calculate zscore with sign to indicate direction
                            diff_zscore = (value_diff - median_diff) / stdev_diff
                            
                            if abs(diff_zscore) > 2:
                                outlier_samples.add(sample)
                                outlier_count[sample] += 1
                                direction = "higher than" if diff_zscore > 0 else "lower than"
                                anomaly_msg = f"Unusual coverage change from {pos}% to {next_pos}%: {value_diff:.2f} ({diff_zscore:.2f} stdev {direction} median)"
                                sample_anomalies[sample].append((f"genebody_coverage_bins_{pos}to{next_pos}pct", "YELLOW", anomaly_msg))
    
    # Log individual anomalies
    for sample, anomalies in sample_anomalies.items():
        # Skip logging individual anomalies since we'll include them in the summary
        # Keep the loop to collect anomalies but don't log each one
        pass
        # for check_name, status, message in anomalies:
        #    log_check_result(log_path, "rseqc", sample, check_name, status, message)
    
    # Determine overall status and log summary
    if outlier_samples:
        # Create CSV-compatible details string with pipe separator instead of newlines
        sample_severity = {}
        for sample in outlier_samples:
            sample_detail_issues = []
            sample_severity[sample] = "YELLOW"  # Default to YELLOW for all outlier samples
            
            for _, severity, message in sample_anomalies[sample]:
                sample_detail_issues.append(message)
                if severity == "RED":
                    sample_severity[sample] = "RED"
            
            # Limit to first 5 issues for brevity if there are many
            if len(sample_detail_issues) > 5:
                detail_message = f"{len(sample_detail_issues)} issues - First 5: {' | '.join(sample_detail_issues[:5])}"
            else:
                detail_message = f"{len(sample_detail_issues)} issues: {' | '.join(sample_detail_issues)}"
                
            log_check_result(log_path, "rseqc", sample, "genebody_coverage_bins_summary", sample_severity[sample], 
                           f"Coverage anomalies detected", detail_message)
        
        # Create a CSV-safe summary for all samples
        summary_details = []
        for sample in sorted(outlier_samples):
            summary_details.append(f"{sample}: {outlier_count[sample]} issues")
        
        # Set overall status to the maximum severity level of any individual sample
        status = "YELLOW"
        if any(severity == "RED" for severity in sample_severity.values()):
            status = "RED"
            
        message = f"{'Major' if status == 'RED' else 'Minor'} coverage anomalies detected in {len(outlier_samples)} samples"
        log_check_result(log_path, "rseqc", "all", "genebody_coverage_bins", status, message, "; ".join(summary_details))
        return status
    else:
        log_check_result(log_path, "rseqc", "all", "genebody_coverage_bins", "GREEN", "No coverage anomalies detected")
        return "GREEN"

def get_inner_distance_multiqc_stats(outdir, samples, log_path, assay_suffix="_GLbulkRNAseq"):
    """Extract inner distance metrics from MultiQC report for each sample.
    
    Args:
        outdir: The output directory path
        samples: List of sample names
        log_path: Path to the validation log file
        assay_suffix: Suffix for the assay in MultiQC report names
        
    Returns:
        Dictionary of inner distance data by sample or None if data not found
    """
    component = "rseqc"
    check_name = "inner_distance_multiqc_stats"
    print(f"Extracting inner distance stats from MultiQC...")
    
    # Get the RSeQC inner distance directory
    inner_dist_dir = os.path.join(outdir, "RSeQC_Analyses", "04_inner_distance")
    
    # Find the MultiQC report zip file
    multiqc_glob = os.path.join(inner_dist_dir, f"inner_dist_multiqc{assay_suffix}_report.zip")
    multiqc_files = glob.glob(multiqc_glob)
    
    if not multiqc_files:
        print(f"ERROR: Inner distance MultiQC report not found: {multiqc_glob}")
        log_check_result(log_path, component, "all", check_name, "RED", 
                         "Inner distance MultiQC report not found", 
                         f"Expected: {multiqc_glob}")
        return None
    
    # Get the most recent MultiQC report (should only be one, but just in case)
    multiqc_report = sorted(multiqc_files)[-1]
    print(f"Using MultiQC report: {multiqc_report}")
    
    # Extract the zip to a temporary directory
    try:
        with tempfile.TemporaryDirectory() as tmpdirname:
            with zipfile.ZipFile(multiqc_report, 'r') as zip_ref:
                zip_ref.extractall(tmpdirname)
            
            # Search for the MultiQC data JSON file recursively
            json_files = []
            for root, dirs, files in os.walk(tmpdirname):
                for file in files:
                    if file.endswith('.json') and 'multiqc_data' in file:
                        json_files.append(os.path.join(root, file))
            
            if not json_files:
                print(f"ERROR: No multiqc_data.json file found in the extracted zip")
                log_check_result(log_path, component, "all", check_name, "RED", 
                                "No multiqc_data.json found in extracted zip", "")
                return None
                
            multiqc_data_file = json_files[0]
            
            # Parse the MultiQC data file
            try:
                with open(multiqc_data_file, 'r') as f:
                    multiqc_data = json.load(f)
                
                # Extract inner distance data
                inner_distance_data = {}
                
                # Check if the expected data structure exists
                if ('report_plot_data' in multiqc_data and 
                    'rseqc_inner_distance_plot' in multiqc_data['report_plot_data'] and
                    'datasets' in multiqc_data['report_plot_data']['rseqc_inner_distance_plot'] and
                    len(multiqc_data['report_plot_data']['rseqc_inner_distance_plot']['datasets']) > 1 and
                    'lines' in multiqc_data['report_plot_data']['rseqc_inner_distance_plot']['datasets'][1]):
                    
                    # Get the plot data for each sample
                    plot_data = multiqc_data['report_plot_data']['rseqc_inner_distance_plot']['datasets'][1]['lines']
                    
                    # Process each sample's data
                    for dist_data in plot_data:
                        sample_name = dist_data['name']
                        
                        # Remove any file extension or path from the sample name
                        sample_name = os.path.basename(sample_name)
                        if '.' in sample_name:
                            sample_name = sample_name.split('.')[0]
                        
                        # Find the peak inner distance (the one with the highest percentage of reads)
                        if len(dist_data.get('pairs', [])) > 0:
                            # Sort the pairs by read percentage (highest first)
                            sorted_pairs = sorted(dist_data['pairs'], key=lambda i: i[1], reverse=True)
                            max_dist = sorted_pairs[0]
                            
                            inner_distance_data[sample_name] = {
                                'peak_inner_dist': max_dist[0],
                                'peak_inner_dist_pct_reads': max_dist[1]
                            }
                else:
                    print(f"ERROR: MultiQC data does not contain expected inner distance plot data")
                    log_check_result(log_path, component, "all", check_name, "RED", 
                                "MultiQC data doesn't contain inner distance plot data", 
                                "Missing rseqc_inner_distance_plot section in MultiQC data")
                    return None
                
                # Check if we got data for all expected samples
                missing_samples = [s for s in samples if s not in inner_distance_data]
                if missing_samples:
                    print(f"WARNING: Missing inner distance data for {len(missing_samples)} sample(s)")
                    log_check_result(log_path, component, "all", check_name, "YELLOW", 
                                    f"Missing inner distance data for {len(missing_samples)} sample(s)", 
                                    f"Samples: {', '.join(missing_samples[:20])}" + 
                                    ("..." if len(missing_samples) > 20 else ""))
                else:
                    print(f"Found inner distance data for all {len(samples)} samples")
                    log_check_result(log_path, component, "all", check_name, "GREEN", 
                                    "Found inner distance data for all samples")
                
                return inner_distance_data
                
            except Exception as e:
                print(f"ERROR: Failed to parse inner distance MultiQC data: {str(e)}")
                log_check_result(log_path, component, "all", check_name, "RED", 
                            "Failed to parse inner distance MultiQC data", 
                            f"Error: {str(e)}")
                return None
    except Exception as e:
        print(f"ERROR: Failed to extract MultiQC report: {str(e)}")
        log_check_result(log_path, component, "all", check_name, "RED", 
                      "Failed to extract MultiQC report", 
                      f"Error: {str(e)}")
        return None

def report_inner_distance_issues(outdir, inner_dist_data, log_path):
    """Report potential issues with inner distance metrics.
    
    Args:
        outdir: The output directory path
        inner_dist_data: Dictionary of inner distance data by sample
        log_path: Path to the validation log file
        
    Returns:
        Status string ("GREEN", "YELLOW", or "RED")
    """
    component = "rseqc"
    check_name = "inner_distance_metrics"
    print(f"Analyzing inner distance metrics for {len(inner_dist_data)} samples...")
    
    # Extract metrics for summary statistics
    peak_distances = [metrics['peak_inner_dist'] for metrics in inner_dist_data.values()]
    peak_percentages = [metrics['peak_inner_dist_pct_reads'] for metrics in inner_dist_data.values()]
    
    # Calculate summary statistics
    mean_distance = statistics.mean(peak_distances)
    stdev_distance = statistics.stdev(peak_distances) if len(peak_distances) > 1 else 0
    mean_percentage = statistics.mean(peak_percentages)
    stdev_percentage = statistics.stdev(peak_percentages) if len(peak_percentages) > 1 else 0
    
    # Calculate min and max
    min_distance = min(peak_distances)
    max_distance = max(peak_distances)
    min_percentage = min(peak_percentages)
    max_percentage = max(peak_percentages)
    
    # Initialize density threshold crossings
    min_density_crossing = float('inf')
    max_density_crossing = float('-inf')
    density_threshold = 1.0  # 1% density threshold
    
    # Get the inner distance MultiQC report zip file
    inner_dist_dir = os.path.join(outdir, "RSeQC_Analyses", "04_inner_distance")
    multiqc_report = os.path.join(inner_dist_dir, f"inner_dist_multiqc_GLbulkRNAseq_report.zip")
    
    try:
        if os.path.exists(multiqc_report):
            # Extract the zip to a temporary directory
            with tempfile.TemporaryDirectory() as tmpdirname:
                with zipfile.ZipFile(multiqc_report, 'r') as zip_ref:
                    zip_ref.extractall(tmpdirname)
                
                # Look for the percentages file in the extracted directory
                percentages_file = None
                for root, dirs, files in os.walk(tmpdirname):
                    for file in files:
                        if file == "rseqc_inner_distance_plot_Percentages.txt":
                            percentages_file = os.path.join(root, file)
                            break
                    if percentages_file:
                        break
                
                if percentages_file and os.path.exists(percentages_file):
                    # Parse the distribution data
                    with open(percentages_file, 'r') as f:
                        lines = f.readlines()
                        
                    if len(lines) > 1:  # Ensure there's data after the header
                        # First line is the header with bin positions
                        header = lines[0].strip().split('\t')
                        
                        # Process each sample's distribution
                        for i in range(1, len(lines)):
                            sample_data = lines[i].strip().split('\t')
                            sample_name = sample_data[0]
                            
                            # Skip if we don't have enough data points
                            if len(sample_data) < 2:
                                continue
                            
                            # Parse each bin's distance and percentage
                            bins_data = []
                            for j in range(1, len(sample_data)):
                                try:
                                    # Data is stored as (distance, percentage)
                                    bin_data = sample_data[j].strip('()').split(',')
                                    distance = float(bin_data[0])
                                    percentage = float(bin_data[1])
                                    bins_data.append((distance, percentage))
                                except (IndexError, ValueError):
                                    continue
                            
                            # Find first crossing from left (ascending)
                            for j in range(len(bins_data) - 1):
                                if bins_data[j][1] < density_threshold and bins_data[j+1][1] >= density_threshold:
                                    # Linear interpolation to find exact crossing
                                    d1, p1 = bins_data[j]
                                    d2, p2 = bins_data[j+1]
                                    if p2 - p1 > 0:  # Avoid division by zero
                                        crossing = d1 + (density_threshold - p1) * (d2 - d1) / (p2 - p1)
                                        min_density_crossing = min(min_density_crossing, crossing)
                                    break
                            
                            # Find last crossing from right (descending)
                            for j in range(len(bins_data) - 1, 0, -1):
                                if bins_data[j][1] < density_threshold and bins_data[j-1][1] >= density_threshold:
                                    # Linear interpolation to find exact crossing
                                    d1, p1 = bins_data[j]
                                    d2, p2 = bins_data[j-1]
                                    if p2 - p1 > 0:  # Avoid division by zero
                                        crossing = d1 + (density_threshold - p1) * (d2 - d1) / (p2 - p1)
                                        max_density_crossing = max(max_density_crossing, crossing)
                                    break
    except Exception as e:
        print(f"WARNING: Could not calculate density crossings: {str(e)}")
        # Fall back to using simple min/max if we can't calculate density crossings
        min_density_crossing = min_distance - stdev_distance
        max_density_crossing = max_distance + stdev_distance
    
    # Use the density crossings if we found them, otherwise fall back to min/max
    if min_density_crossing == float('inf') or max_density_crossing == float('-inf'):
        min_density_crossing = min_distance - stdev_distance
        max_density_crossing = max_distance + stdev_distance
    
    # Format the summary message
    summary_message = "Inner distance metrics summary (mean ± stdev)"
    summary_message += f", Range (1% PDF thresholds): {min_density_crossing:.2f}-{max_density_crossing:.2f} bp"
    summary_message += f", Peak distance: {mean_distance:.2f} ±{stdev_distance:.2f} bp"
    summary_message += f", Peak read %: {mean_percentage:.4f} ±{stdev_percentage:.4f}"
    summary_message += f", Min/Max peak %: {min_percentage:.4f}-{max_percentage:.4f}"
    
    # Log the summary statistics
    log_check_result(log_path, component, "all", "inner_distance_summary", "GREEN", 
                     summary_message, "")
    print(f"SUMMARY: {summary_message}")
    
    # Prepare to track issues
    distance_outlier_samples = []
    percentage_outlier_samples = []
    
    # Create data-driven thresholds (using standard deviation outlier detection)
    # Flag samples that are more than 3 standard deviations from the mean
    # or that fall outside the density crossing range
    STDEV_THRESHOLD = 3  # More than 3 stdev from mean is an outlier
    
    # Analyze each sample
    for sample, metrics in inner_dist_data.items():
        peak_dist = metrics['peak_inner_dist']
        peak_pct = metrics['peak_inner_dist_pct_reads']
        
        # Check for distance outliers using Z-score (standard deviations from mean)
        z_score_dist = abs(peak_dist - mean_distance) / stdev_distance if stdev_distance > 0 else 0
        is_distance_outlier = (z_score_dist > STDEV_THRESHOLD) or (peak_dist < min_density_crossing) or (peak_dist > max_density_crossing)
            
        # Check for read percentage outliers
        z_score_pct = abs(peak_pct - mean_percentage) / stdev_percentage if stdev_percentage > 0 else 0
        is_percentage_outlier = z_score_pct > STDEV_THRESHOLD
        
        # Track samples with issues
        if is_distance_outlier:
            distance_outlier_samples.append((sample, peak_dist, z_score_dist))
        
        if is_percentage_outlier:
            percentage_outlier_samples.append((sample, peak_pct, z_score_pct))
    
    # Determine overall status
    status = "GREEN"
    
    if distance_outlier_samples or percentage_outlier_samples:
        status = "YELLOW"
    
    # Create a summary entry only if issues exist
    if status != "GREEN":
        # Create a summarized message without the verbose description
        message = ""
        if distance_outlier_samples:
            message += f"{len(distance_outlier_samples)} samples with outlier peak inner distances"
        if percentage_outlier_samples:
            if message:
                message += "; "
            message += f"{len(percentage_outlier_samples)} samples with outlier peak read percentages"
            
        # Log the summary with sample ids in the details
        details = ""
        if distance_outlier_samples:
            distance_sample_ids = [s[0] for s in distance_outlier_samples]
            if details:
                details += "; "
            details += f"Peak distance outliers: {', '.join(distance_sample_ids)}"
        if percentage_outlier_samples:
            percentage_sample_ids = [s[0] for s in percentage_outlier_samples]
            if details:
                details += "; "
            details += f"Peak percentage outliers: {', '.join(percentage_sample_ids)}"
            
        log_check_result(log_path, component, "all", check_name, status, message, details)
        print(f"WARNING: Inner distance issues found: {message}")
        
        # Log individual sample entries for each outlier
        for sample, peak_dist, z_score in distance_outlier_samples:
            sample_message = f"{z_score:.2f} stdevs from mean"
            sample_details = f"value={peak_dist:.2f} bp, mean={mean_distance:.2f} bp, stdev={stdev_distance:.2f} bp"
            log_check_result(log_path, component, sample, "outlier_peak_inner_distance", status, sample_message, sample_details)
            
        for sample, peak_pct, z_score in percentage_outlier_samples:
            sample_message = f"{z_score:.2f} stdevs from mean"
            sample_details = f"value={peak_pct:.4f}%, mean={mean_percentage:.4f}%, stdev={stdev_percentage:.4f}%"
            log_check_result(log_path, component, sample, "outlier_peak_read_percentage", status, sample_message, sample_details)
    else:
        message = "All samples have consistent inner distance metrics"
        print(f"PASSED: {message}")
        log_check_result(log_path, component, "all", check_name, status, message, "")
    
    return status

def get_read_distribution_multiqc_stats(outdir, samples, log_path, assay_suffix="_GLbulkRNAseq"):
    """Extract read distribution metrics from MultiQC report for each sample.
    
    Args:
        outdir: The output directory path
        samples: List of sample names
        log_path: Path to the validation log file
        assay_suffix: Suffix for the assay in MultiQC report names
        
    Returns:
        Dictionary of read distribution data by sample or None if data not found
    """
    component = "rseqc"
    check_name = "read_distribution_multiqc_stats"
    print(f"Extracting read distribution stats from MultiQC...")
    
    # Get the RSeQC read distribution directory
    read_dist_dir = os.path.join(outdir, "RSeQC_Analyses", "05_read_distribution")
    
    # Find the MultiQC report zip file
    multiqc_glob = os.path.join(read_dist_dir, f"read_dist_multiqc{assay_suffix}_report.zip")
    multiqc_files = glob.glob(multiqc_glob)
    
    if not multiqc_files:
        print(f"ERROR: Read distribution MultiQC report not found: {multiqc_glob}")
        log_check_result(log_path, component, "all", check_name, "RED", 
                         "Read distribution MultiQC report not found", 
                         f"Expected: {multiqc_glob}")
        return None
    
    # Get the most recent MultiQC report (should only be one, but just in case)
    multiqc_report = sorted(multiqc_files)[-1]
    print(f"Using MultiQC report: {multiqc_report}")
    
    # Extract the zip to a temporary directory
    try:
        with tempfile.TemporaryDirectory() as tmpdirname:
            with zipfile.ZipFile(multiqc_report, 'r') as zip_ref:
                zip_ref.extractall(tmpdirname)
            
            # Search for the MultiQC data JSON file recursively
            json_files = []
            for root, dirs, files in os.walk(tmpdirname):
                for file in files:
                    if file.endswith('.json') and 'multiqc_data' in file:
                        json_files.append(os.path.join(root, file))
            
            if not json_files:
                print(f"ERROR: No multiqc_data.json file found in the extracted zip")
                log_check_result(log_path, component, "all", check_name, "RED", 
                                "No multiqc_data.json found in extracted zip", "")
                return None
                
            multiqc_data_file = json_files[0]
            
            # Parse the MultiQC data
            try:
                # Read and parse the JSON data
                with open(multiqc_data_file, 'r') as f:
                    multiqc_data = json.load(f)
                
                # Extract read distribution data for each sample
                read_dist_data = {}
                
                if 'report_saved_raw_data' not in multiqc_data or 'multiqc_rseqc_read_distribution' not in multiqc_data['report_saved_raw_data']:
                    print(f"ERROR: Read distribution data not found in MultiQC report")
                    log_check_result(log_path, component, "all", check_name, "RED", 
                                    "Read distribution data not found in MultiQC report")
                    return None
                
                # Extract the read distribution data for each sample
                raw_data = multiqc_data['report_saved_raw_data']['multiqc_rseqc_read_distribution']
                
                for sample_name, read_data in raw_data.items():
                    # Clean up sample name by removing '_read_dist' suffix
                    clean_sample = sample_name.replace('_read_dist', '')
                    
                    # Extract the core metrics
                    read_dist_data[clean_sample] = {}
                    
                    # Extract basic fields that match directly
                    for field in ['cds_exons_tag_pct', '5_utr_exons_tag_pct', '3_utr_exons_tag_pct', 
                                'introns_tag_pct', 'tss_up_1kb_tag_pct', 'tes_down_1kb_tag_pct',
                                'other_intergenic_tag_pct']:
                        if field in read_data:
                            # Remove '_tag' from the field name
                            clean_field = field.replace('_tag', '')
                            read_dist_data[clean_sample][clean_field] = read_data[field]
                    
                    # Calculate derived fields as in parse_multiqc.py
                    if 'tss_up_5kb_tag_pct' in read_data and 'tss_up_1kb_tag_pct' in read_data:
                        read_dist_data[clean_sample]['tss_up_1kb_5kb_pct'] = read_data['tss_up_5kb_tag_pct'] - read_data['tss_up_1kb_tag_pct']
                    
                    if 'tss_up_10kb_tag_pct' in read_data and 'tss_up_5kb_tag_pct' in read_data:
                        read_dist_data[clean_sample]['tss_up_5kb_10kb_pct'] = read_data['tss_up_10kb_tag_pct'] - read_data['tss_up_5kb_tag_pct']
                    
                    if 'tes_down_5kb_tag_pct' in read_data and 'tes_down_1kb_tag_pct' in read_data:
                        read_dist_data[clean_sample]['tss_down_1kb_5kb_pct'] = read_data['tes_down_5kb_tag_pct'] - read_data['tes_down_1kb_tag_pct']
                    
                    if 'tes_down_10kb_tag_pct' in read_data and 'tes_down_5kb_tag_pct' in read_data:
                        read_dist_data[clean_sample]['tss_down_5kb_10kb_pct'] = read_data['tes_down_10kb_tag_pct'] - read_data['tes_down_5kb_tag_pct']
                
                # Check if we got data for all expected samples
                missing_samples = [s for s in samples if s not in read_dist_data]
                if missing_samples:
                    print(f"WARNING: Missing read distribution data for {len(missing_samples)} sample(s)")
                    log_check_result(log_path, component, "all", check_name, "YELLOW", 
                                    f"Missing read distribution data for {len(missing_samples)} sample(s)", 
                                    f"Samples: {', '.join(missing_samples[:20])}" + 
                                    ("..." if len(missing_samples) > 20 else ""))
                else:
                    print(f"Found read distribution data for all {len(samples)} samples")
                    log_check_result(log_path, component, "all", check_name, "GREEN", 
                                    "Found read distribution data for all samples")
                
                return read_dist_data
                
            except Exception as e:
                print(f"ERROR: Failed to parse read distribution MultiQC data: {str(e)}")
                log_check_result(log_path, component, "all", check_name, "RED", 
                            "Failed to parse read distribution MultiQC data", 
                            f"Error: {str(e)}")
                return None
    except Exception as e:
        print(f"ERROR: Failed to extract MultiQC report: {str(e)}")
        log_check_result(log_path, component, "all", check_name, "RED", 
                      "Failed to extract MultiQC report", 
                      f"Error: {str(e)}")
        return None

def report_read_distribution_issues(outdir, read_dist_data, log_path):
    """Report potential issues with read distribution metrics.
    
    Args:
        outdir: The output directory path
        read_dist_data: Dictionary of read distribution data by sample
        log_path: Path to the validation log file
        
    Returns:
        Status string ("GREEN", "YELLOW", or "RED")
    """
    component = "rseqc"
    check_name = "read_distribution_metrics"
    print(f"Analyzing read distribution metrics for {len(read_dist_data)} samples...")
    
    # Extract metrics for summary statistics
    cds_exons_pct = [data.get('cds_exons_pct', 0) for data in read_dist_data.values()]
    five_utr_pct = [data.get('5_utr_exons_pct', 0) for data in read_dist_data.values()]
    three_utr_pct = [data.get('3_utr_exons_pct', 0) for data in read_dist_data.values()]
    introns_pct = [data.get('introns_pct', 0) for data in read_dist_data.values()]
    intergenic_pct = [data.get('other_intergenic_pct', 0) for data in read_dist_data.values()]
    
    # Calculate summary statistics
    mean_cds = statistics.mean(cds_exons_pct) if cds_exons_pct else 0
    stdev_cds = statistics.stdev(cds_exons_pct) if len(cds_exons_pct) > 1 else 0
    
    mean_five_utr = statistics.mean(five_utr_pct) if five_utr_pct else 0
    stdev_five_utr = statistics.stdev(five_utr_pct) if len(five_utr_pct) > 1 else 0
    
    mean_three_utr = statistics.mean(three_utr_pct) if three_utr_pct else 0
    stdev_three_utr = statistics.stdev(three_utr_pct) if len(three_utr_pct) > 1 else 0
    
    mean_introns = statistics.mean(introns_pct) if introns_pct else 0
    stdev_introns = statistics.stdev(introns_pct) if len(introns_pct) > 1 else 0
    
    mean_intergenic = statistics.mean(intergenic_pct) if intergenic_pct else 0
    stdev_intergenic = statistics.stdev(intergenic_pct) if len(intergenic_pct) > 1 else 0
    
    # Create a list of regions and their mean percentages
    regions = [
        ("CDS exons", mean_cds, stdev_cds),
        ("5' UTR", mean_five_utr, stdev_five_utr),
        ("3' UTR", mean_three_utr, stdev_three_utr),
        ("Introns", mean_introns, stdev_introns),
        ("Intergenic", mean_intergenic, stdev_intergenic)
    ]
    
    # Sort by percentage (highest first)
    regions.sort(key=lambda x: x[1], reverse=True)
    
    # Format the summary message
    summary_message = "Read distribution summary (mean ± stdev)"
    
    # Add significant regions (>0.1%) sorted by abundance
    significant_regions = []
    for name, mean_value, stdev_value in regions:
        if mean_value > 0.1:  # Only include regions with >0.1% reads
            significant_regions.append(f"{name}: {mean_value:.2f}% ±{stdev_value:.2f}")
    
    # Create details string with regions
    details = "; ".join(significant_regions)
    
    # Log the summary statistics
    log_check_result(log_path, component, "all", "read_distribution_summary", "GREEN", 
                     summary_message, details)
    print(f"SUMMARY: {summary_message} - {details}")
    
    # Prepare to track issues
    cds_outlier_samples = []
    introns_outlier_samples = []
    intergenic_outlier_samples = []
    
    # Create data-driven thresholds using standard deviation outlier detection
    YELLOW_THRESHOLD = 2  # More than 2 stdev from mean is a minor outlier (YELLOW)
    RED_THRESHOLD = 4     # More than 4 stdev from mean is a major outlier (RED)
    
    # Track sample severities
    sample_severities = {}
    
    # Analyze each sample
    for sample, metrics in read_dist_data.items():
        cds_pct = metrics.get('cds_exons_pct', 0)
        introns_pct = metrics.get('introns_pct', 0)
        intergenic_pct = metrics.get('other_intergenic_pct', 0)
        
        # Check for CDS outliers - only if stdev is not too small
        if stdev_cds > 0.01:  # Avoid division by very small numbers
            z_score_cds = abs(cds_pct - mean_cds) / stdev_cds
            
            if z_score_cds > RED_THRESHOLD:
                severity = "RED"
                cds_outlier_samples.append((sample, cds_pct, z_score_cds, severity))
                sample_severities[sample] = max(sample_severities.get(sample, "GREEN"), severity)
            elif z_score_cds > YELLOW_THRESHOLD:
                severity = "YELLOW"
                cds_outlier_samples.append((sample, cds_pct, z_score_cds, severity))
                sample_severities[sample] = max(sample_severities.get(sample, "GREEN"), severity)
        
        # Check for introns outliers - only if stdev is not too small
        if stdev_introns > 0.01:  # Avoid division by very small numbers
            z_score_introns = abs(introns_pct - mean_introns) / stdev_introns
            
            if z_score_introns > RED_THRESHOLD:
                severity = "RED"
                introns_outlier_samples.append((sample, introns_pct, z_score_introns, severity))
                sample_severities[sample] = max(sample_severities.get(sample, "GREEN"), severity)
            elif z_score_introns > YELLOW_THRESHOLD:
                severity = "YELLOW"
                introns_outlier_samples.append((sample, introns_pct, z_score_introns, severity))
                sample_severities[sample] = max(sample_severities.get(sample, "GREEN"), severity)
        
        # Check for intergenic outliers - only if stdev is not too small
        if stdev_intergenic > 0.01:  # Avoid division by very small numbers
            z_score_intergenic = abs(intergenic_pct - mean_intergenic) / stdev_intergenic
            
            if z_score_intergenic > RED_THRESHOLD:
                severity = "RED"
                intergenic_outlier_samples.append((sample, intergenic_pct, z_score_intergenic, severity))
                sample_severities[sample] = max(sample_severities.get(sample, "GREEN"), severity)
            elif z_score_intergenic > YELLOW_THRESHOLD:
                severity = "YELLOW"
                intergenic_outlier_samples.append((sample, intergenic_pct, z_score_intergenic, severity))
                sample_severities[sample] = max(sample_severities.get(sample, "GREEN"), severity)
    
    # Determine overall status
    status = "GREEN"
    if any(severity == "RED" for severity in sample_severities.values()):
        status = "RED"
    elif any(severity == "YELLOW" for severity in sample_severities.values()):
        status = "YELLOW"
    
    # Create a summary entry only if issues exist
    if status != "GREEN":
        # Create a summarized message 
        red_count = sum(1 for severity in sample_severities.values() if severity == "RED")
        yellow_count = sum(1 for severity in sample_severities.values() if severity == "YELLOW")
        
        message = ""
        if red_count > 0:
            message += f"{red_count} samples with major read distribution anomalies"
        if yellow_count > 0:
            if message:
                message += "; "
            message += f"{yellow_count} samples with minor read distribution anomalies"
            
        # Log the summary with sample ids in the details
        details = ""
        if cds_outlier_samples:
            red_samples = [s[0] for s in cds_outlier_samples if s[3] == "RED"]
            yellow_samples = [s[0] for s in cds_outlier_samples if s[3] == "YELLOW"]
            
            if red_samples:
                details += f"CDS major outliers: {', '.join(red_samples)}"
            if yellow_samples:
                if details:
                    details += "; "
                details += f"CDS minor outliers: {', '.join(yellow_samples)}"
        
        if introns_outlier_samples:
            red_samples = [s[0] for s in introns_outlier_samples if s[3] == "RED"]
            yellow_samples = [s[0] for s in introns_outlier_samples if s[3] == "YELLOW"]
            
            if red_samples:
                if details:
                    details += "; "
                details += f"Intron major outliers: {', '.join(red_samples)}"
            if yellow_samples:
                if details:
                    details += "; "
                details += f"Intron minor outliers: {', '.join(yellow_samples)}"
        
        if intergenic_outlier_samples:
            red_samples = [s[0] for s in intergenic_outlier_samples if s[3] == "RED"]
            yellow_samples = [s[0] for s in intergenic_outlier_samples if s[3] == "YELLOW"]
            
            if red_samples:
                if details:
                    details += "; "
                details += f"Intergenic major outliers: {', '.join(red_samples)}"
            if yellow_samples:
                if details:
                    details += "; "
                details += f"Intergenic minor outliers: {', '.join(yellow_samples)}"
            
        log_check_result(log_path, component, "all", check_name, status, message, details)
        print(f"WARNING: Read distribution issues found: {message}")
        
        # Log individual sample entries for each outlier
        for sample, cds_pct, z_score, severity in cds_outlier_samples:
            sample_message = f"{z_score:.2f} stdevs from mean"
            sample_details = f"value={cds_pct:.2f}%, mean={mean_cds:.2f}%, stdev={stdev_cds:.2f}%"
            log_check_result(log_path, component, sample, "outlier_cds_exons_pct", severity, sample_message, sample_details)
            
        for sample, introns_pct, z_score, severity in introns_outlier_samples:
            sample_message = f"{z_score:.2f} stdevs from mean"
            sample_details = f"value={introns_pct:.2f}%, mean={mean_introns:.2f}%, stdev={stdev_introns:.2f}%"
            log_check_result(log_path, component, sample, "outlier_introns_pct", severity, sample_message, sample_details)
            
        for sample, intergenic_pct, z_score, severity in intergenic_outlier_samples:
            sample_message = f"{z_score:.2f} stdevs from mean"
            sample_details = f"value={intergenic_pct:.2f}%, mean={mean_intergenic:.2f}%, stdev={stdev_intergenic:.2f}%"
            log_check_result(log_path, component, sample, "outlier_intergenic_pct", severity, sample_message, sample_details)
    else:
        message = "All samples have consistent read distribution metrics"
        details = f"Using thresholds: {YELLOW_THRESHOLD} stdev (warning), {RED_THRESHOLD} stdev (error); Regions analyzed: CDS exons, introns, intergenic regions"
        print(f"PASSED: {message}")
        log_check_result(log_path, component, "all", check_name, status, message, details)
    
    return status

def main():
    """Main function to execute the validation checks."""
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Validate RSeQC output files")
    parser.add_argument("--runsheet", required=True, help="Path to the runsheet CSV file")
    parser.add_argument("--outdir", required=True, help="Path to the output directory")
    parser.add_argument("--assay_suffix", default="_GLbulkRNAseq", help="Assay suffix for MultiQC reports")
    
    args = parser.parse_args()
    
    # Validate runsheet existence
    if not os.path.exists(args.runsheet):
        print(f"ERROR: Runsheet not found: {args.runsheet}")
        sys.exit(1)
    
    # Validate output directory existence
    if not os.path.isdir(args.outdir):
        print(f"ERROR: Output directory not found: {args.outdir}")
        sys.exit(1)
    
    # Initialize validation log
    vv_log_path = initialize_vv_log(args.outdir)
    
    # Parse the runsheet to get samples
    try:
        samples, metadata = parse_runsheet(args.runsheet)
    except Exception as e:
        print(f"ERROR: Failed to parse runsheet: {str(e)}")
        sys.exit(1)
    
    # Determine if data is paired-end - handle cases when it's a numpy.bool
    paired_end_value = metadata.get('paired_end', False)
    # For string values, check for truthy strings
    if isinstance(paired_end_value, str):
        is_paired_end = paired_end_value.lower() in ['true', 't', 'yes', 'y', '1']
    else:
        # For boolean or numpy.bool values, convert directly to standard Python bool
        is_paired_end = bool(paired_end_value)
    
    print("\n=== RSeQC Validation ===")
    print(f"Analyzing {len(samples)} samples")
    print(f"Data type: {'Paired-end' if is_paired_end else 'Single-end'}")
    print("Dataset properties:")
    for key, value in metadata.items():
        print(f"  - {key}: {value}")
    
    # Check result status tracking
    check_results = {}
    multiqc_cache = {}  # Cache for MultiQC data to avoid reparsing
    overall_status = "GREEN"  # Track overall validation status
    
    # Check for existence of gene body coverage files
    print("\n--- Checking Gene Body Coverage Files ---")
    result = check_gene_body_coverage_existence(args.outdir, samples, vv_log_path)
    check_results["Gene Body Coverage Files"] = result
    overall_status = max(overall_status, result)
    
    # Check for existence of infer experiment files
    print("\n--- Checking Infer Experiment Files ---")
    result = check_infer_experiment_existence(args.outdir, samples, vv_log_path)
    check_results["Infer Experiment Files"] = result
    overall_status = max(overall_status, result)
    
    # Check for existence of read distribution files
    print("\n--- Checking Read Distribution Files ---")
    result = check_read_distribution_existence(args.outdir, samples, vv_log_path)
    check_results["Read Distribution Files"] = result
    overall_status = max(overall_status, result)
    
    # Check for existence of inner distance files (paired-end only)
    if is_paired_end:
        print("\n--- Checking Inner Distance Files ---")
        result = check_inner_distance_existence(args.outdir, samples, vv_log_path)
        check_results["Inner Distance Files"] = result
        overall_status = max(overall_status, result)
    else:
        print("\n--- Skipping Inner Distance Files (single-end data) ---")
        check_results["Inner Distance Files"] = "SKIPPED"
    
    # Analyze gene body coverage if files exist
    print("\n--- Analyzing Gene Body Coverage Stats ---")
    genebody_data = get_genebody_coverage_multiqc_stats(args.outdir, samples, vv_log_path, args.assay_suffix)
    multiqc_cache['genebody'] = genebody_data
    
    # Report genebody coverage issues
    if genebody_data:
        result = report_genebody_coverage_issues(args.outdir, genebody_data, vv_log_path)
        check_results["Gene Body Coverage Metrics"] = result
        overall_status = max(overall_status, result)
    
    # Get infer experiment stats from its MultiQC report
    print("\n--- Analyzing Infer Experiment Stats ---")
    infer_exp_data = get_infer_experiment_multiqc_stats(args.outdir, samples, vv_log_path, args.assay_suffix)
    multiqc_cache['infer_exp'] = infer_exp_data
    
    # Report infer experiment issues
    if infer_exp_data:
        result = report_infer_experiment_issues(args.outdir, infer_exp_data, vv_log_path)
        check_results["Infer Experiment Metrics"] = result
        overall_status = max(overall_status, result)
    else:
        print("No infer experiment data found, skipping analysis")
        check_results["Infer Experiment Metrics"] = "SKIPPED"
    
    # Analyze coverage bin outliers if data is available
    print("\n--- Analyzing Coverage Bins for Outliers ---")
    if genebody_data:
        result = detect_coverage_bin_outliers(args.outdir, genebody_data, vv_log_path)
        check_results["Coverage Bin Analysis"] = result
        overall_status = max(overall_status, result)
    else:
        check_results["Coverage Bin Analysis"] = "SKIPPED"
    
    # Analyze inner distance metrics if data is available (paired-end only)
    if is_paired_end:
        print("\n--- Analyzing Inner Distance Metrics ---")
        inner_dist_data = get_inner_distance_multiqc_stats(args.outdir, samples, vv_log_path, args.assay_suffix)
        multiqc_cache['inner_dist'] = inner_dist_data
        
        # Report inner distance issues
        if inner_dist_data:
            result = report_inner_distance_issues(args.outdir, inner_dist_data, vv_log_path)
            check_results["Inner Distance Metrics"] = result
            overall_status = max(overall_status, result)
        else:
            print("No inner distance data found, skipping analysis")
            check_results["Inner Distance Metrics"] = "SKIPPED"
    else:
        print("\n--- Skipping Inner Distance Metrics (single-end data) ---")
        check_results["Inner Distance Metrics"] = "SKIPPED"
    
    # Get read distribution stats from MultiQC
    print("\n--- Analyzing Read Distribution Metrics ---")
    read_dist_data = get_read_distribution_multiqc_stats(args.outdir, samples, vv_log_path, args.assay_suffix)
    multiqc_cache['read_dist'] = read_dist_data
    
    # Analyze read distribution metrics if data is available
    if read_dist_data:
        result = report_read_distribution_issues(args.outdir, read_dist_data, vv_log_path)
        check_results["Read Distribution Metrics"] = result
        overall_status = max(overall_status, result)
    else:
        print("No read distribution data found, skipping analysis")
        check_results["Read Distribution Metrics"] = "SKIPPED"
    
    print_summary(check_results, vv_log_path, overall_status)

def print_summary(check_results, vv_log_path, overall_status="GREEN"):
    """Print a summary table of check results."""
    # Print summary table
    print("\n=== RSeQC Validation Summary ===")
    print(f"{'Check':<40}{'Status'}")
    print("-------------------------------------------------------")
    for check, status in check_results.items():
        print(f"{check:<40}{status}")
    
    # Overall status
    print("\n=== Overall Status ===")
    if overall_status == "GREEN":
        print("All checks passed!")
    elif overall_status == "YELLOW":
        print("Some checks have warnings")
    else:
        print("Some checks failed")
    
    # Final message
    print(f"\nRSeQC validation complete. Full details available in the validation log:")
    print(f"  {vv_log_path}")

if __name__ == "__main__":
    main()
