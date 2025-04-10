#! /usr/bin/env python

from pathlib import Path
import shutil
import sys
import argparse

import pandas as pd

# Parse command line arguments
parser = argparse.ArgumentParser(description='Sort files into subdirectories based on sample names from a runsheet.')
parser.add_argument('-f', '--from', dest='from_dir', required=True, help='Source directory containing files')
parser.add_argument('-t', '--to', dest='to_dir', required=True, help='Destination directory for sorted files')
parser.add_argument('-r', '--runsheet', dest='runsheet_path', required=True, help='Path to the CSV runsheet file')
parser.add_argument('-g', '--glob', dest='glob_suffix', default='*', help='File pattern suffix to match (default: matches any file starting with sample name)')

args = parser.parse_args()

# Determine samples list
df = pd.read_csv(args.runsheet_path)
samples = list(df["Sample Name"])
print(samples)

# For a given directory, sort all files into {sample: str, [files: str]}
files_by_sample = dict()
for sample in samples:
    pattern = f"{sample}{args.glob_suffix}"
    print(f"Looking for files matching: {pattern}")
    files_for_this_sample = list(Path(args.from_dir).glob(pattern))
    
    # Move files
    for file in files_for_this_sample:
        dest = Path(args.to_dir) / sample / file.name
        print(f"Moving {file} to {dest}")
        dest.parent.mkdir( parents=True, exist_ok=True )
        shutil.move(file, dest)