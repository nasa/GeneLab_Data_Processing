#! /usr/bin/env python
""" Script that extracts max readlength from fastQC report
"""
import sys
from pathlib import Path

def main(fastqc_path: Path) -> int:
    with fastqc_path.open() as f:
        for line in f.readlines():
            # find line like below
            # Sequence length	151
            # or like
            # Sequence length	146-151
            if line.startswith("Sequence length"):
                read_length_token = line.split()[-1]
                if "-" in read_length_token:
                    max_read_length = read_length_token.split("-")[-1]
                else:
                    max_read_length = read_length_token
                return max_read_length

if __name__ == "__main__":
    max_read_length = main(Path(sys.argv[1]))
    print(max_read_length)
