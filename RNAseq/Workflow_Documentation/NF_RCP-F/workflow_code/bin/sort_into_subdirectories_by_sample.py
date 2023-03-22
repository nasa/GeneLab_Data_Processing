#! /usr/bin/env python

from pathlib import Path
import shutil
import sys

import pandas as pd

FROM = sys.argv[1]
TO = sys.argv[2]
GLOB_SUFFIX = sys.argv[3]

# Determine samples list
[runsheetPath] = (f for f in Path.cwd().glob("Metadata/*runsheet*.csv"))
df = pd.read_csv(runsheetPath)
samples = list(df["Sample Name"])
print(samples)

# For a given directory, sort all files into {sample: str, [files: str]}
files_by_sample = dict()
for sample in samples:
    files_for_this_sample = list(Path(FROM).glob(f"{sample}{GLOB_SUFFIX}"))

    # Move files
    for file in files_for_this_sample:
        dest = Path(TO) / sample / file.name
        print(f"Moving {file} to {dest}")
        dest.parent.mkdir( parents=True, exist_ok=True )
        shutil.move(file, dest)

