#!/usr/bin/env python3

import sys
import re

def gtf_to_bed(input_file, output_file):
    """
    Convert a prokaryotic GTF file to BED format.
    
    Parameters:
        input_file (str): Path to input GTF file.
        output_file (str): Path to output BED file.
    """
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            # Skip comment lines
            if line.startswith('#'):
                continue
            
            # Parse GTF line
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            chrom = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attributes = fields[8]

            # Extract gene_id
            gene_id_match = re.search(r'gene_id "([^"]+)"', attributes)
            gene_id = gene_id_match.group(1) if gene_id_match else "unknown"

            # Compute block size
            block_size = end - start
            block_start = 0  # Single exon, so start at 0

            # Write to BED format
            bed_line = f"{chrom}\t{start}\t{end}\t{gene_id}\t0\t{strand}\t{start}\t{end}\t0\t1\t{block_size},\t{block_start},\n"
            f_out.write(bed_line)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python gtf_to_bed.py input.gtf output.bed")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    gtf_to_bed(input_file, output_file)