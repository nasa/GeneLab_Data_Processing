#!/usr/bin/env python3

import os
import sys
import hashlib
import argparse

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

def should_include(filepath):
    """Check if file should be included in MD5 calculation."""
    # Skip files in GeneLab except for HTML report, software_versions, and purged processing_info.txt
    allowed_files = [
        "NF_MAAffymetrix_v" + args.workflow_version + "_GLmicroarray.html",
        "software_versions_GLmicroarray.md",
        "nextflow_processing_info_GLmicroarray.txt"
    ]
    if "/GeneLab/" in filepath and not any(filepath.endswith(f) for f in allowed_files):
        return False
    
    # Skip ISA.zip
    if filepath.endswith("ISA.zip"):
        return False

    # Skip VV logs
    if "/VV_Logs/" in filepath:
        return False
    
    return True

def main():
    parser = argparse.ArgumentParser(description='Generate MD5 sum files for GeneLab data.')
    parser.add_argument('--outdir', required=True, help='Output directory containing files to process')
    parser.add_argument('--workflow_version', default='', help='Version of the NF_MAAffymetrix workflow manifest')
    
    global args
    args = parser.parse_args()
    
    # Make sure outdir is absolute path
    outdir = os.path.abspath(args.outdir)
    
    # Create output files and initialize them without headers
    processed_md5_file = f"processed_md5sum_GLmicroarray.tsv"
    with open(processed_md5_file, 'w') as f:
        f.write("File Name\tmd5sum\n")
    
    
    # Track processed files for reporting
    processed_count = 0
    
    # Walk through all files recursively
    print(f"Scanning directory: {outdir}")
    for root, _, files in os.walk(outdir):
        for filename in files:
            filepath = os.path.join(root, filename)
            
            # Skip files that shouldn't be included
            if not should_include(filepath):
                continue
            
            # Get just the filename (basename)
            basename = os.path.basename(filepath)
                      
            md5sum = calculate_md5(filepath)
            with open(processed_md5_file, 'a') as f:
                f.write(f"{basename}\t{md5sum}\n")
            processed_count += 1
    
    print(f"Added {processed_count} files to {processed_md5_file}")

    def dedup_file(filename):
        seen = set()
        lines = []
        with open(filename, 'r') as f:
            for line in f:
                key = line.split('\t', 1)[0]  # dedup by basename
                if key not in seen:
                    seen.add(key)
                    lines.append(line)
        with open(filename, 'w') as f:
            f.writelines(lines)

    dedup_file(processed_md5_file)

if __name__ == "__main__":
    main()