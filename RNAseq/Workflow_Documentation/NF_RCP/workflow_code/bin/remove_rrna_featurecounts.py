#!/usr/bin/env python
import sys
import pandas as pd

def main(counts_file, rrna_ids_file, filtered_output, summary_output="rRNA_counts.txt"):
    # Load rRNA IDs
    with open(rrna_ids_file) as f:
        rrna_ids = set(line.strip() for line in f if line.strip())
    
    # Extract header lines
    header_lines = []
    with open(counts_file) as f:
        for line in f:
            if line.startswith('#'):
                header_lines.append(line)
            else:
                break
    
    # Load data (skip header lines)
    df = pd.read_csv(counts_file, sep='\t', comment='#')
    
    # Filter and summarize
    gene_id_col = df.columns[0]
    rna_mask = df[gene_id_col].isin(rrna_ids)
    sample_stats = {col.replace('.bam', ''): 
                    (df[rna_mask][col] >= 1).sum() 
                    for col in df.columns[6:]}
    
    # Write outputs
    with open(summary_output, 'w') as f:
        f.write(f"Total entries listed: {len(df)}\n")
        f.write(f"Total rRNA entries listed: {rna_mask.sum()}\n\n")
        for sample, count in sample_stats.items():
            f.write(f"{sample}: {count} rRNA entries removed\n")
    
    # Write filtered output with original header
    with open(filtered_output, 'w') as f:
        # Write header lines first
        for line in header_lines:
            f.write(line)
        # Then write the filtered data
        df[~rna_mask].to_csv(f, sep='\t', index=False)

if __name__ == '__main__':
    if len(sys.argv) not in [4, 5]:
        sys.exit("Usage: python remove_rrna_featurecounts.py <counts_file> <rrna_ids_file> <filtered_output> [summary_output]")
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4] if len(sys.argv) == 5 else "rRNA_counts.txt")
