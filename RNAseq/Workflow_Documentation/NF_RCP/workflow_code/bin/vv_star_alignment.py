#!/usr/bin/env python

"""
Script to validate and verify STAR alignment output based on runsheet information.
Expects to run from inside a output directory GLDS-##.

Parse the input runsheet to get:
Sample Name
Paired End
Has ERCC

Check that the expected output directories exist

Section-specific checks:
- check_star_output_existence: Check if all expected STAR output files exist for each sample
- check_bam_file_integrity: Verify BAM file integrity using samtools quickcheck
- check_star_alignment_multiqc: Extract and validate STAR alignment metrics from MultiQC
- report_star_alignment_outliers: Detect and report outliers in alignment metrics
- check_mapping_rates: Validate mapping rates (uniquely mapped, multimapped) meet thresholds
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
        '02-STAR_Alignment'
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
        # Create new file with header
        with open(vv_log_path, 'w') as f:
            f.write("component,sample_id,check_name,status,flag_code,message,details\n")
    
    return vv_log_path


def log_check_result(log_path, component, sample_id, check_name, status, message="", details=""):
    """Log check result to the VV_log.csv file."""
    
    def escape_field(field, is_details=False):
        # Convert to string if not already
        if not isinstance(field, str):
            field = str(field)
        
        # Escape commas, quotes and newlines
        field = field.replace('"', '""')  # Double quotes to escape them
        
        # If it's the details field, truncate to 1000 chars if too long
        if is_details and len(field) > 1000:
            field = field[:997] + "..."
        
        # If field contains commas, quotes or newlines, enclose in quotes
        if ',' in field or '"' in field or '\n' in field:
            field = f'"{field}"'
            
        return field
    
    # Map status strings to flag codes
    flag_codes = {
        "GREEN": "20",   # Using strings for consistency in CSV
        "YELLOW": "30",
        "RED": "50",
        "HALT": "80"
    }

    # Get numeric flag code based on status color
    flag_code = flag_codes.get(status, "80")  # Default to HALT if unknown status
    
    with open(log_path, 'a') as f:
        component = escape_field(component)
        sample_id = escape_field(sample_id)
        check_name = escape_field(check_name)
        status = escape_field(status)
        message = escape_field(message)
        details = escape_field(details, is_details=True)
        
        f.write(f"{component},{sample_id},{check_name},{status},{flag_code},{message},{details}\n")


def check_star_output_existence(outdir, samples, paired_end, log_path, assay_suffix="_GLbulkRNAseq"):
    """Check if all expected STAR alignment output files exist for each sample."""
    alignment_dir = os.path.join(outdir, '02-STAR_Alignment')
    
    # Expected file patterns for each sample in sample-specific subdirectories
    expected_patterns = [
        "{sample}/{sample}_Aligned.sortedByCoord.out.bam",
        "{sample}/{sample}_Aligned.toTranscriptome.out.bam",
        "{sample}/{sample}_Log.final.out",
        "{sample}/{sample}_Log.progress.out",
        "{sample}/{sample}_Log.out",
        "{sample}/{sample}_ReadsPerGene.out.tab",
        "{sample}/{sample}_SJ.out.tab"
    ]
    
    # Add mate-specific files if paired-end
    if paired_end:
        expected_patterns.extend([
            "{sample}/{sample}_R1_unmapped.fastq.gz", 
            "{sample}/{sample}_R2_unmapped.fastq.gz"
        ])
    else:
        expected_patterns.append("{sample}/{sample}_R1_unmapped.fastq.gz")
    
    # Dataset-level files (directly in the alignment directory)
    dataset_files = [
        "STAR_NumNonZeroGenes{assay_suffix}.csv",
        "STAR_Unnormalized_Counts{assay_suffix}.csv"
    ]
    
    # MultiQC files in the MultiQC_Reports subdirectory
    multiqc_files = [
        "MultiQC_Reports/align_multiqc{assay_suffix}_data.zip"
    ]
    
    missing_files_by_sample = {}
    
    # Check sample-specific files
    for sample in samples:
        missing_files = []
        for pattern in expected_patterns:
            file_path = os.path.join(alignment_dir, pattern.format(sample=sample))
            if not os.path.exists(file_path):
                missing_files.append(os.path.basename(file_path))
        
        if missing_files:
            missing_files_by_sample[sample] = missing_files
    
    # Check dataset-level files
    missing_dataset_files = []
    for file_name in dataset_files:
        # Format with assay_suffix if present in the pattern
        if "{assay_suffix}" in file_name:
            file_path = os.path.join(alignment_dir, file_name.format(assay_suffix=assay_suffix))
        else:
            file_path = os.path.join(alignment_dir, file_name)
            
        if not os.path.exists(file_path):
            missing_dataset_files.append(os.path.basename(file_path))
    
    # Log results
    if missing_files_by_sample or missing_dataset_files:
        # Log sample-specific missing files
        for sample, missing_files in missing_files_by_sample.items():
            log_check_result(
                log_path, 
                "STAR_alignment", 
                sample, 
                "check_star_output_existence", 
                "HALT", 
                f"Missing {len(missing_files)} expected STAR output files", 
                ",".join(missing_files)
            )
        
        # Log dataset-level missing files
        if missing_dataset_files:
            log_check_result(
                log_path, 
                "STAR_alignment", 
                "all", 
                "check_star_output_existence", 
                "HALT", 
                f"Missing {len(missing_dataset_files)} expected dataset-level STAR output files", 
                ",".join(missing_dataset_files)
            )
            
        print(f"WARNING: Some expected STAR output files are missing")
        return False
    else:
        log_check_result(
            log_path, 
            "STAR_alignment", 
            "all", 
            "check_star_output_existence", 
            "GREEN", 
            "All expected STAR output files exist", 
            ""
        )
        print("All expected STAR output files exist")
        return True


def check_bam_file_integrity(outdir, samples, log_path):
    """Verify BAM file integrity using samtools quickcheck."""
    alignment_dir = os.path.join(outdir, '02-STAR_Alignment')
    
    # BAM file patterns to check for each sample
    bam_patterns = [
        "{sample}/{sample}_Aligned.toTranscriptome.out.bam",
        "{sample}/{sample}_Aligned.sortedByCoord.out.bam"
    ]
    
    failed_samples = {}
    
    for sample in samples:
        failed_bams = []
        for pattern in bam_patterns:
            bam_path = os.path.join(alignment_dir, pattern.format(sample=sample))
            
            if not os.path.exists(bam_path):
                # Skip if file doesn't exist (already logged in check_star_output_existence)
                continue
                
            try:
                # Run samtools quickcheck
                result = subprocess.run(
                    ["samtools", "quickcheck", bam_path],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    check=False
                )
                
                if result.returncode != 0:
                    failed_bams.append(os.path.basename(bam_path))
                    stderr = result.stderr.decode('utf-8', errors='replace').strip()
                    print(f"BAM integrity check failed for {bam_path}: {stderr}")
            except Exception as e:
                failed_bams.append(os.path.basename(bam_path))
                print(f"Error checking BAM integrity for {bam_path}: {str(e)}")
        
        if failed_bams:
            failed_samples[sample] = failed_bams
    
    # Log results
    if failed_samples:
        for sample, failed_bams in failed_samples.items():
            log_check_result(
                log_path, 
                "STAR_alignment", 
                sample, 
                "check_bam_file_integrity", 
                "HALT", 
                f"{len(failed_bams)} BAM files failed integrity check", 
                ",".join(failed_bams)
            )
        print(f"WARNING: {len(failed_samples)} samples have BAM files that failed integrity check")
        return False
    else:
        log_check_result(
            log_path, 
            "STAR_alignment", 
            "all", 
            "check_bam_file_integrity", 
            "GREEN", 
            "All BAM files passed integrity check", 
            ""
        )
        print("All BAM files passed integrity check")
        return True


def get_star_multiqc_stats(outdir, samples, log_path, assay_suffix="_GLbulkRNAseq"):
    """Extract STAR alignment metrics from MultiQC data."""
    alignment_dir = os.path.join(outdir, '02-STAR_Alignment')
    multiqc_dir = os.path.join(alignment_dir, "MultiQC_Reports")
    multiqc_zip = os.path.join(multiqc_dir, f"align_multiqc{assay_suffix}_data.zip")
    
    if not os.path.exists(multiqc_zip):
        print(f"WARNING: MultiQC data zip file not found: {multiqc_zip}")
        log_check_result(
            log_path, 
            "STAR_alignment", 
            "all", 
            "get_star_multiqc_stats", 
            "HALT", 
            "STAR MultiQC data not found", 
            f"Expected at {multiqc_zip}"
        )
        return False
    
    print(f"Extracting STAR stats from MultiQC data: {multiqc_zip}")
    
    # Create a temporary directory to extract files
    with tempfile.TemporaryDirectory() as temp_dir:
        try:
            # Extract the zip file
            with zipfile.ZipFile(multiqc_zip, 'r') as zip_ref:
                zip_ref.extractall(temp_dir)
            
            # Path to the MultiQC data JSON file (new structure)
            json_path = os.path.join(temp_dir, f"align_multiqc{assay_suffix}_data", "multiqc_data.json")
            
            if not os.path.exists(json_path):
                print("WARNING: multiqc_data.json not found in the expected location")
                log_check_result(
                    log_path, 
                    "STAR_alignment", 
                    "all", 
                    "get_star_multiqc_stats", 
                    "HALT", 
                    "multiqc_data.json not found in MultiQC zip", 
                    f"Expected at {json_path}"
                )
                return False
            
            # Parse the MultiQC data
            with open(json_path) as f:
                multiqc_data = json.load(f)
            
            # Extract STAR metrics from MultiQC data
            star_data = {}
            
            # Check if STAR data exists in the MultiQC report
            if 'report_general_stats_data' not in multiqc_data:
                print("WARNING: No general stats data found in MultiQC data")
                log_check_result(
                    log_path, 
                    "STAR_alignment", 
                    "all", 
                    "get_star_multiqc_stats", 
                    "HALT", 
                    "No general stats data in MultiQC data", 
                    ""
                )
                return False
            
            # Look for STAR data in general_stats_data
            star_metrics_found = False
            for stats_section in multiqc_data['report_general_stats_data']:
                for sample, stats in stats_section.items():
                    if 'uniquely_mapped_percent' in stats and sample not in star_data:
                        # This is STAR data
                        star_metrics_found = True
                        
                        # Clean sample name (remove path prefix if present)
                        sample_name = os.path.basename(sample)
                        
                        # Extract key STAR metrics
                        star_data[sample_name] = {
                            'uniquely_mapped_percent': stats.get('uniquely_mapped_percent', 0),
                            'multimapped_percent': stats.get('multimapped_percent', 0),
                            'multimapped_toomany_percent': stats.get('multimapped_toomany_percent', 0),
                            'unmapped_tooshort_percent': stats.get('unmapped_tooshort_percent', 0),
                            'unmapped_other_percent': stats.get('unmapped_other_percent', 0)
                        }
            
            if not star_metrics_found:
                print("WARNING: No STAR metrics found in MultiQC data")
                log_check_result(
                    log_path, 
                    "STAR_alignment", 
                    "all", 
                    "get_star_multiqc_stats", 
                    "HALT", 
                    "No STAR metrics found in MultiQC data", 
                    ""
                )
                return False
            
            # Check if all samples are in the stats
            missing_samples = []
            for sample in samples:
                if sample not in star_data:
                    missing_samples.append(sample)
            
            if missing_samples:
                print(f"WARNING: {len(missing_samples)} samples missing from STAR stats:")
                for sample in missing_samples[:10]:  # Show first 10 only
                    print(f"  - {sample}")
                if len(missing_samples) > 10:
                    print(f"  ... and {len(missing_samples) - 10} more")
                    
                log_check_result(
                    log_path, 
                    "STAR_alignment", 
                    "all", 
                    "get_star_multiqc_stats", 
                    "YELLOW", 
                    f"Missing {len(missing_samples)} samples in STAR stats", 
                    ",".join(missing_samples[:20])  # Limit to 20 sample names
                )
            else:
                print("All samples found in STAR alignment stats")
                log_check_result(
                    log_path, 
                    "STAR_alignment", 
                    "all", 
                    "get_star_multiqc_stats", 
                    "GREEN", 
                    "All samples found in STAR alignment stats", 
                    ""
                )
            
            return star_data
            
        except Exception as e:
            print(f"ERROR extracting STAR stats: {str(e)}")
            log_check_result(
                log_path, 
                "STAR_alignment", 
                "all", 
                "get_star_multiqc_stats", 
                "HALT", 
                f"Error extracting STAR stats: {str(e)}", 
                ""
            )
            return False


def report_star_alignment_outliers(outdir, star_data, log_path):
    """Identify and report outliers in STAR alignment metrics."""
    if not star_data:
        print("No STAR data to analyze for outliers")
        log_check_result(
            log_path, 
            "STAR_alignment", 
            "all", 
            "report_star_alignment_outliers", 
            "HALT", 
            "No STAR data to analyze for outliers", 
            ""
        )
        return False
    
    # Metrics to check for outliers
    metrics_to_check = {
        "uniquely_mapped_percent": "uniquely mapped reads percentage",
        "multimapped_percent": "multi-mapped reads percentage",
        "multimapped_toomany_percent": "too many multi-mapped reads percentage",
        "unmapped_tooshort_percent": "unmapped reads (too short) percentage",
        "unmapped_other_percent": "unmapped reads (other) percentage"
    }
    
    # Thresholds for outlier detection
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
        for sample, sample_data in star_data.items():
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
                
                # Calculate deviation from the median in terms of standard deviations
                deviation = abs(value - median_value)
                stdev_multiples = deviation / stdev if stdev > 0 else 0
                
                # Check if sample is an outlier based on threshold
                if stdev_multiples >= stdev_threshold:
                    outlier_samples.add(sample)
                    detail = f"{sample}: {metric_description}={value:.2f}%, {stdev_multiples:.2f} stdevs from median"
                    outlier_details.append(detail)
                    
                    print(f"  {code} outlier: {detail}")
                    log_check_result(
                        log_path, 
                        "STAR_alignment", 
                        sample, 
                        f"outlier_{metric_key}", 
                        code, 
                        f"{stdev_multiples:.2f} stdevs from median", 
                        f"value={value:.2f}%, median={median_value:.2f}%, stdev={stdev:.2f}"
                    )
                    break  # Once we've identified the highest threshold, we can stop
    
    # Summarize results
    if outlier_samples:
        print(f"\nFound {len(outlier_samples)} samples with outlier metrics:")
        for sample in sorted(outlier_samples)[:10]:  # Show first 10 only
            print(f"  - {sample}")
        if len(outlier_samples) > 10:
            print(f"  ... and {len(outlier_samples) - 10} more")
            
        log_check_result(
            log_path, 
            "STAR_alignment", 
            "all", 
            "report_star_alignment_outliers", 
            "YELLOW", 
            f"Found {len(outlier_samples)} samples with outlier metrics", 
            ";".join(outlier_details[:20])  # Limit to 20 outlier details
        )
        return True
    else:
        print("\nNo outliers detected in STAR alignment metrics")
        log_check_result(
            log_path, 
            "STAR_alignment", 
            "all", 
            "report_star_alignment_outliers", 
            "GREEN", 
            "No outliers detected in STAR alignment metrics", 
            ""
        )
        return True


def check_mapping_rates(outdir, star_data, log_path):
    """Check if mapping rates meet expected thresholds."""
    if not star_data:
        print("No STAR data to analyze mapping rates")
        log_check_result(
            log_path, 
            "STAR_alignment", 
            "all", 
            "check_mapping_rates", 
            "HALT", 
            "No STAR data to analyze mapping rates", 
            ""
        )
        return False
    
    # Define thresholds based on protocol.py configuration
    thresholds = {
        "total_mapped": [
            {"code": "YELLOW", "type": "lower", "value": 70},
            {"code": "RED", "type": "lower", "value": 50}
        ],
        "multi_mapped": [
            {"code": "YELLOW", "type": "lower", "value": 30},
            {"code": "RED", "type": "lower", "value": 15}
        ]
    }
    
    # Track samples that fail thresholds
    failed_total_mapped = []
    failed_multi_mapped = []
    
    # Check each sample against thresholds
    for sample, metrics in star_data.items():
        # Calculate total mapped percentage
        total_mapped_pct = metrics.get('uniquely_mapped_percent', 0) + metrics.get('multimapped_percent', 0)
        
        # Calculate multi-mapped percentage
        multi_mapped_pct = metrics.get('multimapped_percent', 0) + metrics.get('multimapped_toomany_percent', 0)
        
        # Check total mapped percentage against thresholds
        for threshold in thresholds["total_mapped"]:
            if threshold["type"] == "lower" and total_mapped_pct < threshold["value"]:
                failed_total_mapped.append({
                    "sample": sample,
                    "value": total_mapped_pct,
                    "threshold": threshold["value"],
                    "code": threshold["code"]
                })
                print(f"{threshold['code']} flag: {sample} has low mapping rate ({total_mapped_pct:.2f}% < {threshold['value']}%)")
                log_check_result(
                    log_path, 
                    "STAR_alignment", 
                    sample, 
                    "check_total_mapping_rate", 
                    threshold["code"], 
                    f"Low mapping rate: {total_mapped_pct:.2f}% < {threshold['value']}%", 
                    f"uniquely_mapped={metrics.get('uniquely_mapped_percent', 0):.2f}%, multimapped={metrics.get('multimapped_percent', 0):.2f}%"
                )
                break  # Stop after finding the first threshold that's violated
        
        # Check multi-mapped percentage against thresholds
        for threshold in thresholds["multi_mapped"]:
            if threshold["type"] == "lower" and multi_mapped_pct < threshold["value"]:
                failed_multi_mapped.append({
                    "sample": sample,
                    "value": multi_mapped_pct,
                    "threshold": threshold["value"],
                    "code": threshold["code"]
                })
                print(f"{threshold['code']} flag: {sample} has low multi-mapping rate ({multi_mapped_pct:.2f}% < {threshold['value']}%)")
                log_check_result(
                    log_path, 
                    "STAR_alignment", 
                    sample, 
                    "check_multi_mapping_rate", 
                    threshold["code"], 
                    f"Low multi-mapping rate: {multi_mapped_pct:.2f}% < {threshold['value']}%", 
                    f"multimapped={metrics.get('multimapped_percent', 0):.2f}%, multimapped_toomany={metrics.get('multimapped_toomany_percent', 0):.2f}%"
                )
                break  # Stop after finding the first threshold that's violated
    
    # Summarize results
    if failed_total_mapped or failed_multi_mapped:
        print(f"\nFound mapping rate issues:")
        if failed_total_mapped:
            print(f"  - {len(failed_total_mapped)} samples with low total mapping rate")
        if failed_multi_mapped:
            print(f"  - {len(failed_multi_mapped)} samples with low multi-mapping rate")
            
        log_check_result(
            log_path, 
            "STAR_alignment", 
            "all", 
            "check_mapping_rates_summary", 
            "YELLOW", 
            f"Found mapping rate issues in {len(set([f['sample'] for f in failed_total_mapped + failed_multi_mapped]))} samples", 
            f"total_mapped_failures={len(failed_total_mapped)}, multi_mapped_failures={len(failed_multi_mapped)}"
        )
        return False
    else:
        print("\nAll samples meet mapping rate thresholds")
        log_check_result(
            log_path, 
            "STAR_alignment", 
            "all", 
            "check_mapping_rates_summary", 
            "GREEN", 
            "All samples meet mapping rate thresholds", 
            ""
        )
        return True


def main():
    """Main function to process runsheet and validate STAR alignment output."""
    parser = argparse.ArgumentParser(description='Validate STAR alignment output based on runsheet information.')
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
    
    # Check if STAR output files exist
    check_star_output_existence(args.outdir, sample_names, paired_end, vv_log_path, args.assay_suffix)
    
    # Check BAM file integrity
    check_bam_file_integrity(args.outdir, sample_names, vv_log_path)
    
    # Get STAR MultiQC stats
    star_data = get_star_multiqc_stats(args.outdir, sample_names, vv_log_path, args.assay_suffix)
    
    # Report outliers in STAR alignment metrics
    if star_data:
        report_star_alignment_outliers(args.outdir, star_data, vv_log_path)
        
        # Check mapping rates against thresholds
        check_mapping_rates(args.outdir, star_data, vv_log_path)
    
    print("STAR alignment validation complete")


if __name__ == "__main__":
    main() 