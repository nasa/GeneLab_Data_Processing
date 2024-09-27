#!/usr/bin/env python
import subprocess
import sys

jobid = sys.argv[1]

# if wanting to use, this should be added to the snakemake call from the root workflow dir: `--cluster-status scripts/slurm-status.py`

output = str(subprocess.check_output("sacct -j %s --format State --noheader | head -1 | awk '{print $1}'" % jobid, shell=True).strip())

running_status=["PENDING", "CONFIGURING", "COMPLETING", "RUNNING", "SUSPENDED"]
if "COMPLETED" in output:
    print("success")
elif any(r in output for r in running_status):
    print("running")
else:
    print("failed")
