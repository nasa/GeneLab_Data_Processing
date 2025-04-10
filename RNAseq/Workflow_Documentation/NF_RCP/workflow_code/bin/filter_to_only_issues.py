#! /usr/bin/env python
import pandas as pd
import argparse

# Set up command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('--assay_suffix', type=str, required=True, 
                    help='Suffix for input/output files (e.g. "_GLbulkRNAseq")')
args = parser.parse_args()

INPUT_FN = f"VV_log_final{args.assay_suffix}.csv"
OUTPUT_FN = f"VV_log_final_only_issues{args.assay_suffix}.csv"

df = pd.read_csv(INPUT_FN, sep=",")

# Filter out GREEN status (code 20) to keep only issues
# Handle both string and integer flag_codes
df_filtered = df.loc[~((df["flag_code"] == "20") | (df["flag_code"] == 20))]
df_filtered.to_csv(OUTPUT_FN, sep=",", index=False)