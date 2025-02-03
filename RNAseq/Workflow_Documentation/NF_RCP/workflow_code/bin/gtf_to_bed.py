#!/usr/bin/env python3

import sys
import re

def gtf_to_bed(input_file, output_file):
    """
    Convert GTF file to BED format
    Parameters:
        input_file (str): Path to input GTF file
        output_file (str): Path to output BED file
    """
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            # Skip comment lines
            if line.startswith('#'):
                continue
            
            try:
                # Parse GTF line
                fields = line.strip().split('\t')
                if len(fields) < 9:  # GTF must have 9 fields
                    continue
                
                chrom = fields[0]
                start = int(fields[3]) - 1  # Convert to 0-based
                end = int(fields[4])
                strand = fields[6]
                
                # Extract gene_id from attributes
                attributes = fields[8]
                gene_id_match = re.search(r'gene_id "([^"]+)"', attributes)
                gene_id = gene_id_match.group(1) if gene_id_match else "unknown"
                
                # Write BED line
                bed_line = f"{chrom}\t{start}\t{end}\t{gene_id}\t0\t{strand}\n"
                f_out.write(bed_line)
                
            except (IndexError, ValueError) as e:
                # Skip malformed lines
                continue

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python gtf_to_bed.py input.gtf output.bed")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    gtf_to_bed(input_file, output_file) 