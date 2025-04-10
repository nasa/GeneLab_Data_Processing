#!/usr/bin/env python3

import os
import sys
import hashlib
import argparse
import re

def calculate_md5(filepath):
    """Calculate MD5 hash for a file."""
    md5_hash = hashlib.md5()
    
    # Follow symlinks to get the actual file
    actual_path = os.path.realpath(filepath) if os.path.islink(filepath) else filepath
    
    try:
        with open(actual_path, "rb") as f:
            # Read in chunks in case of large files
            for chunk in iter(lambda: f.read(4096), b""):
                md5_hash.update(chunk)
        return md5_hash.hexdigest()
    except Exception as e:
        sys.stderr.write(f"Error calculating MD5 for {filepath}: {str(e)}\n")
        return "ERROR"

def is_raw_file(filepath):
    """Check if file is a raw FASTQ or raw multiqc file."""
    # Match raw fastq files but not trimming reports
    if "/Fastq/" in filepath and filepath.endswith("raw.fastq.gz") and "_raw.fastq.gz_" not in filepath:
        return True
    # Match raw multiqc reports
    if "raw_multiqc" in filepath and filepath.endswith(".zip"):
        return True
    return False

def should_include(filepath, outdir):
    """Check if file should be included in MD5 calculation."""
    # Skip files in VV_Logs
    if "/VV_Logs/" in filepath:
        return False
    
    # Skip files in GeneLab except for qc_metrics
    if "/GeneLab/" in filepath and not filepath.endswith("qc_metrics" + args.assay_suffix + ".csv"):
        return False
        
    return True

def main():
    parser = argparse.ArgumentParser(description='Generate MD5 sum files for GeneLab data.')
    parser.add_argument('--outdir', required=True, help='Output directory containing files to process')
    parser.add_argument('--assay_suffix', required=True, help='Suffix for assay type (e.g., _GLbulkRNAseq)')
    
    global args
    args = parser.parse_args()
    
    # Create output files
    raw_md5_file = f"raw_md5sum{args.assay_suffix}.tsv"
    processed_md5_file = f"processed_md5sum{args.assay_suffix}.tsv"
    
    # Initialize files without headers
    with open(raw_md5_file, 'w') as f:
        pass  # Create empty file
    
    with open(processed_md5_file, 'w') as f:
        pass  # Create empty file
    
    # Make sure outdir is absolute path
    outdir = os.path.abspath(args.outdir)
    
    # Track processed files for reporting
    raw_count = 0
    processed_count = 0
    
    # Walk through all files recursively
    print(f"Scanning directory: {outdir}")
    for root, _, files in os.walk(outdir):
        for filename in files:
            filepath = os.path.join(root, filename)
            
            # Skip files that shouldn't be included
            if not should_include(filepath, outdir):
                continue
            
            # Get just the filename (basename)
            basename = os.path.basename(filepath)
            
            # Calculate MD5 sum
            md5sum = calculate_md5(filepath)
            
            # Add to appropriate output file
            if is_raw_file(filepath):
                with open(raw_md5_file, 'a') as f:
                    f.write(f"{basename}\t{md5sum}\n")
                raw_count += 1
            else:
                with open(processed_md5_file, 'a') as f:
                    f.write(f"{basename}\t{md5sum}\n")
                processed_count += 1
    
    print(f"Added {raw_count} files to {raw_md5_file}")
    print(f"Added {processed_count} files to {processed_md5_file}")

if __name__ == "__main__":
    main()
