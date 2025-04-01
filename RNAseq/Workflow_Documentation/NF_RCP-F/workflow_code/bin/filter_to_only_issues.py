#! /usr/bin/env python
import pandas as pd

INPUT_FN="VV_log_final_GLbulkRNAseq.tsv"
OUTPUT_FN="VV_log_final_only_issues_GLbulkRNAseq.tsv"

df = pd.read_csv(INPUT_FN, sep="\t")

df_filtered = df.loc[~df["code"].isin(["FlagCode.GREEN"])]
df_filtered.to_csv(OUTPUT_FN, sep="\t", index=False)