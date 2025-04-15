#! /usr/bin/env python

from pathlib import Path
import shutil
import sys

import pandas as pd

FROM = sys.argv[1]
TO = sys.argv[2]
# Removed GLOB_SUFFIX - we'll always use the "*" pattern to match any files starting with the sample name

# Determine samples list
[runsheetPath] = (f for f in Path.cwd().glob("Metadata/*runsheet*.csv"))
df = pd.read_csv(runsheetPath)
samples = list(df["Sample Name"])
print(f"Found {len(samples)} samples: {samples}")

# For a given directory, sort all files into sample subdirectories
for sample in samples:
    # Use simple pattern to match any file starting with the sample name
    files_for_this_sample = list(Path(FROM).glob(f"{sample}*"))
    
    print(f"Found {len(files_for_this_sample)} files for sample {sample}")

    # Move files
    for file in files_for_this_sample:
        dest = Path(TO) / sample / file.name
        print(f"Moving {file} to {dest}")
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(file, dest)