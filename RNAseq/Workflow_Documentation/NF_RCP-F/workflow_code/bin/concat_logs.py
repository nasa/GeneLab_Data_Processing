#! /usr/bin/env python
# Concatenate V&V logs

from pathlib import Path

logs = sorted(list(Path.cwd().glob("VV_in.tsv*")))

OUTPUT_FN = Path("VV_log_final_GLbulkRNAseq.tsv")

for i, log in enumerate(logs):
    with open(log, "r") as in_f:
        if i == 0:  # first file
            contents = in_f.read()
        else:
            contents += "".join(in_f.readlines()[1:])  # skip header

with open(OUTPUT_FN, "w") as out_f:
    out_f.write(contents)
