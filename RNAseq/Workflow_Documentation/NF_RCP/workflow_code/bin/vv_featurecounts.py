#!/usr/bin/env python

"""
Script to validate featureCounts results based on runsheet information.
Expects to run from inside an output directory GLDS-##.

Parse the input runsheet to get:
Sample Name
Paired End
Has ERCC

Check that the expected output directories exist

Section-specific checks:
- check_featurecounts_files_existence: Check if all expected featureCounts output files exist
- get_featurecounts_multiqc_stats: Extract featureCounts metrics from MultiQC report
- report_multiqc_outliers: Identify and report outliers in featureCounts metrics
- check_rrna_removal: Verify that rRNA genes were properly removed in the rRNArm file

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
        
        # Convert Sample Name to string to ensure compatibility with numeric sample names
        df['Sample Name'] = df['Sample Name'].astype(str)
        
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
    # Define expected directories for featureCounts validation
    expected_dirs = [
        '03-FeatureCounts'
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
    
    if not os.path.exists(vv_log_path):
        with open(vv_log_path, 'w') as f:
            f.write("component,sample_id,check_name,status,flag_code,message,details\n")
    
    return vv_log_path


def log_check_result(log_path, component, sample_id, check_name, status, message="", details=""):
    """Log check result to the VV_log.csv file."""
    def escape_field(field, is_details=False):
        if not isinstance(field, str):
            field = str(field)
        field = field.replace('\n', '; ')
        if is_details and ',' in field:
            field = field.replace(', ', '; ')
        if ',' in field or '"' in field:
            field = field.replace('"', '""')
            field = f'"{field}"'
        return field

    # Map status strings to flag codes
    flag_codes = {
        "GREEN": "20",
        "YELLOW": "30",
        "RED": "50",
        "HALT": "80"
    }

    # Get flag code based on status
    flag_code = flag_codes.get(status, "80")
    
    component = escape_field(component)
    sample_id = escape_field(sample_id)
    check_name = escape_field(check_name)
    status = escape_field(status)
    message = escape_field(message)
    details = escape_field(details, True)
    
    with open(log_path, 'a') as f:
        f.write(f"{component},{sample_id},{check_name},{status},{flag_code},{message},{details}\n")


def check_featurecounts_files_existence(outdir, log_path, assay_suffix="_GLbulkRNAseq"):
    """Check if the expected featureCounts output files exist."""
    featurecounts_dir = os.path.join(outdir, '03-FeatureCounts')
    
    # Define expected files
    expected_files = [
        f"FeatureCounts{assay_suffix}.tsv",
        f"FeatureCounts{assay_suffix}.tsv.summary"
        # f"FeatureCounts_rRNArm{assay_suffix}.tsv"
    ]
    
    # Define MultiQC files in the MultiQC_Reports subdirectory
    multiqc_files = [
        f"MultiQC_Reports/FeatureCounts_multiqc{assay_suffix}_data.zip",
        f"MultiQC_Reports/FeatureCounts_multiqc{assay_suffix}.html"
    ]
    
    missing_files = []
    for file_name in expected_files:
        file_path = os.path.join(featurecounts_dir, file_name)
        if not os.path.exists(file_path):
            missing_files.append(file_name)
    
    if missing_files:
        print(f"WARNING: The following expected featureCounts files are missing:")
        for file_name in missing_files:
            print(f"  - {file_name}")
        log_check_result(log_path, "featurecounts", "all", "check_featurecounts_files_existence", "HALT", 
                       f"Missing {len(missing_files)} expected featureCounts files", 
                       ",".join(missing_files))
        return False
    else:
        print("All expected featureCounts output files exist")
        log_check_result(log_path, "featurecounts", "all", "check_featurecounts_files_existence", "GREEN", 
                       "All expected featureCounts files exist", "")
        return True


def get_featurecounts_multiqc_stats(outdir, samples, log_path, assay_suffix="_GLbulkRNAseq"):
    """Extract featureCounts MultiQC stats for all samples and write to a stats file for analysis."""
    featurecounts_dir = os.path.join(outdir, "03-FeatureCounts")
    multiqc_dir = os.path.join(featurecounts_dir, "MultiQC_Reports")
    multiqc_zip = os.path.join(multiqc_dir, f"FeatureCounts_multiqc{assay_suffix}_data.zip")
    
    if not os.path.exists(multiqc_zip):
        print(f"WARNING: MultiQC data zip file not found: {multiqc_zip}")
        log_check_result(log_path, "featurecounts", "all", "get_featurecounts_multiqc_stats", "RED", 
                         "MultiQC data not found", "")
        return False
    
    print(f"Extracting stats from MultiQC data: {multiqc_zip}")
    
    # Convert all samples to strings for consistent comparison
    samples = [str(sample) for sample in samples]
    
    # Create a temporary directory to extract files
    with tempfile.TemporaryDirectory() as temp_dir:
        try:
            # Extract the zip file
            with zipfile.ZipFile(multiqc_zip, 'r') as zip_ref:
                zip_ref.extractall(temp_dir)
            
            # Path to the MultiQC data JSON file (new structure)
            json_path = os.path.join(temp_dir, f"FeatureCounts_multiqc{assay_suffix}_data", "multiqc_data.json")
            
            if not os.path.exists(json_path):
                print(f"WARNING: No multiqc_data.json file found in the expected location")
                log_check_result(log_path, "featurecounts", "all", "get_featurecounts_multiqc_stats", "RED", 
                               "multiqc_data.json not found in zip", "")
                return False
            
            # Parse the MultiQC data
            with open(json_path) as f:
                multiqc_data = json.load(f)
            
            fc_data = {}
            
            # Check if featureCounts data exists in the MultiQC report
            if ('report_saved_raw_data' not in multiqc_data or 
                'multiqc_featurecounts' not in multiqc_data['report_saved_raw_data']):
                print("WARNING: No featureCounts data found in MultiQC data")
                log_check_result(log_path, "featurecounts", "all", "get_featurecounts_multiqc_stats", "RED", 
                               "No featureCounts data in MultiQC data", "")
                return False
            
            # Extract featureCounts data
            for sample, count_data in multiqc_data['report_saved_raw_data']['multiqc_featurecounts'].items():
                # Clean sample name (remove path if present)
                sample_name = os.path.basename(sample) if '/' in sample else sample
                
                # Convert to string to ensure compatibility with numeric sample names
                sample_name = str(sample_name)
                
                # Extract key metrics
                fc_data[sample_name] = {
                    'total_count': count_data['Total'],
                    'num_assigned': count_data['Assigned'],
                    'pct_assigned': count_data['percent_assigned'],
                    'num_unassigned_nofeatures': count_data['Unassigned_NoFeatures'],
                    'num_unassigned_ambiguity': count_data['Unassigned_Ambiguity'],
                    'pct_unassigned_nofeatures': count_data['Unassigned_NoFeatures'] / count_data['Total'] * 100 if count_data['Total'] > 0 else 0,
                    'pct_unassigned_ambiguity': count_data['Unassigned_Ambiguity'] / count_data['Total'] * 100 if count_data['Total'] > 0 else 0
                }
            
            if not fc_data:
                print("WARNING: No valid featureCounts data found in MultiQC data")
                log_check_result(log_path, "featurecounts", "all", "get_featurecounts_multiqc_stats", "RED", 
                               "No valid featureCounts data found", "")
                return False
            
            # Check if all samples are in the stats
            missing_samples = []
            for sample in samples:
                if sample not in fc_data:
                    missing_samples.append(sample)
            
            if missing_samples:
                print(f"WARNING: {len(missing_samples)} samples missing from featureCounts stats:")
                for sample in missing_samples:
                    print(f"  - {sample}")
                log_check_result(log_path, "featurecounts", "all", "get_featurecounts_multiqc_stats", "RED", 
                               f"Missing {len(missing_samples)} samples in featureCounts stats", 
                               ",".join(missing_samples[:20]))
                return fc_data
            
            print("All samples found in featureCounts stats")
            log_check_result(log_path, "featurecounts", "all", "get_featurecounts_multiqc_stats", "GREEN", 
                           "All samples found in featureCounts stats", "")
            return fc_data
            
        except Exception as e:
            print(f"ERROR extracting featureCounts stats: {str(e)}")
            log_check_result(log_path, "featurecounts", "all", "get_featurecounts_multiqc_stats", "RED", 
                           f"Error extracting featureCounts stats: {str(e)}", "")
            return False


def parse_featurecounts(multiqc_data_dir, assay_suffix="_GLbulkRNAseq"):
    """Parse featureCounts data from MultiQC data directory."""
    multiqc_data_json = os.path.join(multiqc_data_dir, "multiqc_data.json")
    
    if not os.path.exists(multiqc_data_json):
        print(f"WARNING: MultiQC data file not found: {multiqc_data_json}")
        return None
    
    try:
        with open(multiqc_data_json) as f:
            j = json.load(f)
        
        fc_data = {}
        
        if 'report_saved_raw_data' not in j or 'multiqc_featurecounts' not in j['report_saved_raw_data']:
            print("WARNING: No featureCounts data found in MultiQC data")
            return None
        
        for sample, count_data in j['report_saved_raw_data']['multiqc_featurecounts'].items():
            sample_name = sample
            # Clean up sample name if needed
            if "/" in sample:
                sample_name = sample.split("/")[-1]
            
            # Convert to string to ensure compatibility with numeric sample names
            sample_name = str(sample_name)
            
            fc_data[sample_name] = {
                'total_count': count_data['Total'],
                'num_assigned': count_data['Assigned'],
                'pct_assigned': count_data['percent_assigned'],
                'num_unassigned_nofeatures': count_data['Unassigned_NoFeatures'],
                'num_unassigned_ambiguity': count_data['Unassigned_Ambiguity'],
                'pct_unassigned_nofeatures': count_data['Unassigned_NoFeatures'] / count_data['Total'] * 100 if count_data['Total'] > 0 else 0,
                'pct_unassigned_ambiguity': count_data['Unassigned_Ambiguity'] / count_data['Total'] * 100 if count_data['Total'] > 0 else 0
            }
        
        return fc_data
    
    except Exception as e:
        print(f"ERROR parsing featureCounts data: {str(e)}")
        return None


def report_multiqc_outliers(outdir, multiqc_data, log_path):
    """Identify and report outliers in featureCounts MultiQC statistics."""
    if not multiqc_data:
        print("No MultiQC data to analyze for outliers")
        log_check_result(log_path, "featurecounts", "all", "report_multiqc_outliers", "RED", 
                        "No MultiQC data to analyze", "")
        return False

    # Keys to check for outliers in featureCounts data
    metrics_to_check = {
        "pct_assigned": "percentage of reads assigned to features",
        "pct_unassigned_ambiguity": "percentage of reads with ambiguous features",
        "pct_unassigned_nofeatures": "percentage of reads with no features",
        "total_count": "total number of reads counted"
    }
    
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
            print(f"Metric {metric_key} not found in any sample, skipping")
            continue
            
        # Collect values for this metric across all samples
        values = []
        for sample, sample_data in multiqc_data.items():
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
        print(f"  Median: {median_value:.2f}")
        print(f"  Standard Deviation: {stdev:.2f}")
        
        # Check for outliers based on thresholds
        for sample, value in values:
            for threshold in thresholds:
                stdev_threshold = threshold["stdev_threshold"]
                code = threshold["code"]
                
                # Calculate deviation from the median in terms of standard deviations
                deviation = abs(value - median_value)
                stdev_multiples = deviation / stdev if stdev > 0 else 0
                
                # Check if sample is an outlier based on threshold
                if stdev_multiples >= stdev_threshold:
                    outlier_samples.add(sample)
                    detail = f"{sample}: {metric_description}={value:.2f}, {stdev_multiples:.2f} stdevs from median"
                    outlier_details.append(detail)
                    
                    print(f"  {code} outlier: {detail}")
                    log_check_result(log_path, "featurecounts", sample, f"outlier_{metric_key}", code, 
                                   f"{stdev_multiples:.2f} stdevs from median", 
                                   f"value={value:.2f}; median={median_value:.2f}; stdev={stdev:.2f}")
                    break  # Once we've identified the highest threshold, we can stop
    
    # Summarize results
    if outlier_samples:
        print(f"\nFound {len(outlier_samples)} samples with outlier metrics:")
        for sample in sorted(outlier_samples):
            print(f"  - {sample}")
        
        log_check_result(log_path, "featurecounts", "all", "report_multiqc_outliers", "YELLOW", 
                       f"Found {len(outlier_samples)} samples with outlier metrics", 
                       ";".join(outlier_details))
        return True
    else:
        print("\nNo outliers detected in featureCounts metrics")
        log_check_result(log_path, "featurecounts", "all", "report_multiqc_outliers", "GREEN", 
                       "No outliers detected in featureCounts metrics", "")
        return True


def check_rrna_removal(outdir, log_path, assay_suffix="_GLbulkRNAseq"):
    """Check if rRNA genes were properly removed in the rRNArm file."""
    featurecounts_dir = os.path.join(outdir, "03-FeatureCounts")
    
    # Define paths to both files
    regular_counts_file = os.path.join(featurecounts_dir, f"FeatureCounts{assay_suffix}.tsv")
    rrna_removed_file = os.path.join(featurecounts_dir, f"FeatureCounts_rRNArm{assay_suffix}.tsv")
    
    # Check if both files exist
    if not os.path.exists(regular_counts_file):
        print(f"WARNING: Regular FeatureCounts file not found: {regular_counts_file}")
        log_check_result(log_path, "featurecounts", "all", "check_rrna_removal", "RED", 
                        "Regular FeatureCounts file not found", "")
        return False
        
    if not os.path.exists(rrna_removed_file):
        print(f"WARNING: rRNA-removed FeatureCounts file not found: {rrna_removed_file}")
        log_check_result(log_path, "featurecounts", "all", "check_rrna_removal", "RED", 
                        "rRNA-removed FeatureCounts file not found", "")
        return False
    
    try:
        # Read both files
        regular_df = pd.read_csv(regular_counts_file, sep='\t', comment='#')
        rrna_removed_df = pd.read_csv(rrna_removed_file, sep='\t', comment='#')
        
        # Check the gene IDs column name (usually first column)
        id_column = regular_df.columns[0]
        
        # Get the sets of gene IDs from both files
        regular_genes = set(regular_df[id_column])
        rrna_removed_genes = set(rrna_removed_df[id_column])
        
        # Find the removed genes
        removed_genes = regular_genes - rrna_removed_genes
        
        # Check if any genes were removed
        if len(removed_genes) > 0:
            # Sort the removed genes for consistent output
            removed_list = sorted(list(removed_genes))
            
            # Calculate total counts in both files (sum across all sample columns)
            # Identify sample count columns (usually all columns except gene identifiers and metadata)
            sample_columns = [col for col in regular_df.columns if col != id_column and not col.startswith('Chr') and not col.startswith('Start') and not col.startswith('End') and not col.startswith('Strand') and not col.startswith('Length')]
            
            total_original_counts = regular_df[sample_columns].sum().sum()
            total_removed_counts = rrna_removed_df[sample_columns].sum().sum()
            counts_difference = total_original_counts - total_removed_counts
            percent_removed = (counts_difference / total_original_counts) * 100 if total_original_counts > 0 else 0
            
            # Get counts for each removed gene
            removed_genes_counts = {}
            for gene in removed_list:
                gene_row = regular_df[regular_df[id_column] == gene]
                if not gene_row.empty:
                    gene_counts = gene_row[sample_columns].sum().sum()
                    removed_genes_counts[gene] = gene_counts
            
            # Format removed genes with counts
            removed_genes_with_counts = [f"{gene}:{removed_genes_counts.get(gene, 0)}" for gene in removed_list]
            
            print(f"Found {len(removed_genes)} rRNA genes removed:")
            print(f"Total removed counts: {counts_difference:,} ({percent_removed:.2f}% of original)")
            for i, gene in enumerate(removed_list[:10]):  # Show first 10 only
                count = removed_genes_counts.get(gene, 0)
                pct = (count / total_original_counts) * 100 if total_original_counts > 0 else 0
                print(f"  - {gene}: {count:,} counts ({pct:.2f}% of total)")
            
            if len(removed_list) > 10:
                print(f"  ... and {len(removed_list) - 10} more")
                
            # Prepare details for logging
            removed_str = ",".join(removed_genes_with_counts)
            log_message = f"Removed {len(removed_genes)} rRNA genes ({counts_difference:,} counts, {percent_removed:.2f}% of total)"
            
            log_check_result(log_path, "featurecounts", "all", "check_rrna_removal", "GREEN", 
                           log_message, 
                           removed_str[:1000] if len(removed_str) > 1000 else removed_str)  # Limit details length
            return True
        else:
            print("WARNING: No rRNA genes were removed between the two files")
            log_check_result(log_path, "featurecounts", "all", "check_rrna_removal", "RED", 
                           "No rRNA genes were removed", "")
            return False
            
    except Exception as e:
        print(f"ERROR checking rRNA removal: {str(e)}")
        log_check_result(log_path, "featurecounts", "all", "check_rrna_removal", "RED", 
                       f"Error checking rRNA removal: {str(e)}", "")
        return False


def add_featurecounts_group_stats(outdir, multiqc_data, log_path):
    """Calculate and log summary statistics for featureCounts metrics across all samples."""
    if not multiqc_data or len(multiqc_data) == 0:
        print("No featureCounts data to analyze for group statistics")
        log_check_result(
            log_path, 
            "featurecounts", 
            "all", 
            "add_featurecounts_group_stats", 
            "RED", 
            "No featureCounts data to analyze for group statistics", 
            ""
        )
        return False
    
    # Metrics to analyze
    metrics = {
        "total_count": "total number of reads counted",
        "num_assigned": "number of reads assigned to features",
        "pct_assigned": "percentage of reads assigned to features",
        "num_unassigned_nofeatures": "number of reads with no features",
        "num_unassigned_ambiguity": "number of reads with ambiguous features",
        "pct_unassigned_nofeatures": "percentage of reads with no features",
        "pct_unassigned_ambiguity": "percentage of reads with ambiguous features"
    }
    
    # Collect values for each metric
    metric_values = {metric: [] for metric in metrics.keys()}
    
    for sample, data in multiqc_data.items():
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
    print("\nfeatureCounts Metrics Summary Statistics:")
    
    for metric, description in metrics.items():
        if metric not in stats_summary:
            continue
            
        stats = stats_summary[metric]
        
        # Format appropriately based on metric type
        if metric.startswith("pct_") or metric.endswith("_percent"):
            # Percentage metrics
            detail_str = (f"Range: {stats['min']:.2f}% - {stats['max']:.2f}%; "
                         f"Median: {stats['median']:.2f}%; "
                         f"Mean: {stats['mean']:.2f}%; "
                         f"StdDev: {stats['stddev']:.2f}")
            
            print(f"  {description}: {detail_str}")
            
            # Log each metric as a separate entry
            log_check_result(
                log_path, 
                "featurecounts", 
                "all", 
                f"featurecounts_stats_{metric}", 
                "GREEN", 
                f"{description}", 
                detail_str
            )
        else:
            # Count metrics
            detail_str = (f"Range: {stats['min']:.0f} - {stats['max']:.0f}; "
                         f"Median: {stats['median']:.0f}; "
                         f"Mean: {stats['mean']:.1f}; "
                         f"StdDev: {stats['stddev']:.1f}")
            
            print(f"  {description}: {detail_str}")
            
            # Log each metric as a separate entry
            log_check_result(
                log_path, 
                "featurecounts", 
                "all", 
                f"featurecounts_stats_{metric}", 
                "GREEN", 
                f"{description}", 
                detail_str
            )
    
    return True


def main():
    """Main function to process runsheet and validate featureCounts output."""
    parser = argparse.ArgumentParser(description='Validate featureCounts output based on runsheet information.')
    parser.add_argument('--runsheet', '-r', required=True, help='Path to the runsheet CSV file')
    parser.add_argument('--outdir', '-o', default=os.getcwd(), 
                        help='Output directory (GLDS-## folder), defaults to current directory')
    parser.add_argument('--assay-suffix', default="_GLbulkRNAseq", 
                        help='Assay suffix used in MultiQC data filenames (default: _GLbulkRNAseq)')
    args = parser.parse_args()

    # Initialize VV log
    vv_log_path = initialize_vv_log(args.outdir)
    
    # Check directory structure
    check_directory_structure(args.outdir)
    
    # Parse the runsheet
    runsheet_df = parse_runsheet(args.runsheet)
    
    # Extract sample names and convert to strings to handle numeric sample names
    sample_names = [str(sample) for sample in runsheet_df['Sample Name'].tolist()]
    
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

    # Validation checks will be implemented here
    
    # Check if featureCounts output files exist
    check_featurecounts_files_existence(args.outdir, vv_log_path, args.assay_suffix)
    
    # Get featureCounts MultiQC stats
    multiqc_data = get_featurecounts_multiqc_stats(args.outdir, sample_names, vv_log_path, args.assay_suffix)
    
    # Add group statistics for featureCounts metrics
    if multiqc_data:
        add_featurecounts_group_stats(args.outdir, multiqc_data, vv_log_path)
        
        # Report outliers in featureCounts metrics
        report_multiqc_outliers(args.outdir, multiqc_data, vv_log_path)
    
    # Remove rRNA check since the file is no longer provided in the pipeline
    # check_rrna_removal(args.outdir, vv_log_path, args.assay_suffix)
    
    print("FeatureCounts validation complete")


if __name__ == "__main__":
    main()
