#!/usr/bin/env python

import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='This script swaps the MAG IDs back to what they were prior to running KEGGDecoder.')

required = parser.add_argument_group('required arguments')

required.add_argument("-i", "--input-tsv", help='Output table from KEGGDecoder', action="store", required=True)
required.add_argument("-m", "--map-tsv", help='Tab-delimited map with 1st column holding original name, and 2nd column holding modified name', action="store", required=True)

parser.add_argument("-o", "--output-tsv", help='Output table with adjusted MAG IDs (Default: "out.tsv")', action="store", default="out.tsv")

args = parser.parse_args()

# reading in mapping file into dictionary
map_dict = {}
with open(args.map_tsv) as mapping:
    for line in mapping:
        line = line.strip().split("\t")
        map_dict[line[1]] = line[0]

# reading in output table from KEGGDecoder
in_tab = pd.read_csv(args.input_tsv, sep = "\t", index_col = 0)

# renaming back to what they were before modifying to be compliant with KEGGDecoder
mod_tab = in_tab.rename(index = map_dict)

# writing out modified file
mod_tab.to_csv(args.output_tsv, sep = "\t")
