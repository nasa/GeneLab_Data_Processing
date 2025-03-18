#!/usr/bin/env python

"""
gtf_to_bed.py - Direct conversion of GTF to BED format

This script directly converts a GTF file to BED12 format without
requiring intermediate formats. It properly handles transcript structures
and extracts identifiers in a general way that works with both NCBI 
and Ensembl-style GTF files.

Usage:
    gtf_to_bed.py ${genome_gtf} ${genome_gtf.baseName}.bed
    gtf_to_bed.py input.gtf output.bed --chromosome-prefix "Chromosome"

Conversion Logic:
----------------
1. Parse the GTF file and group features by transcript
   - Identify exons, CDS, start/stop codons for each transcript
   - Convert 1-based GTF coordinates to 0-based BED coordinates
   - Skip gene-level features (only process transcript-level)

2. Handle transcript structures
   - Create exons from explicit exon features
   - If no exons defined, create a single exon spanning the transcript
   - Sort exons by position to maintain correct structure
   - Calculate block sizes and starts for BED12 blocks

3. Determine CDS boundaries
   - For protein-coding genes: use actual CDS boundaries
   - For pseudogenes/non-coding: set both cds_start and cds_end to transcript end
   - This follows Ensembl's convention, even when CDS features exist in pseudogenes

4. Identifier selection (prioritized)
   - Use gene_id as the primary identifier (e.g., PA14_RS00005)
   - If gene_id not available, fall back to locus_tag
   - Last resort, use transcript_id
   - For transcripts without explicit transcript_id, derive from locus_tag or gene_id

5. Gene Type Detection
   - Identify pseudogenes via:
     - gene_biotype/gene_type = 'pseudogene'
     - pseudo="true" attribute
   - Identify non-coding RNAs via non-protein_coding biotypes
   - These identifications determine CDS boundary handling

Features:
--------
- Maintains proper exon structure in BED12 format
- Uses consistent chromosome naming (configurable)
- Handles both coding and non-coding transcripts
- Correctly processes pseudogenes with CDS features that are present in NCBI GTF files
- Works with both NCBI and Ensembl GTF formats
"""

import argparse
import re
import sys
from collections import defaultdict

def parse_attributes(attr_string):
    """Parse the GTF attribute string into a dictionary.
    
    Extracts key-value pairs from GTF attribute field (9th column).
    GTF attributes have the format: key "value"; key "value";
    
    Args:
        attr_string (str): GTF attribute string from 9th column
        
    Returns:
        dict: Dictionary mapping attribute names to values
    """
    attributes = {}
    # GTF format has attributes like: gene_id "ENSG00000139618"; gene_name "BRCA2";
    pattern = re.compile(r'(\w+)\s+"([^"]+)"')
    for match in pattern.finditer(attr_string):
        key, value = match.groups()
        attributes[key] = value
    return attributes

class Transcript:
    """Represents a transcript with features and coordinates.
    
    This class stores all information about a transcript including its identifiers,
    genomic coordinates, exons, CDS regions, and feature types. It processes
    GTF features (exons, CDS, etc.) and converts them to the BED12 format.
    
    Attributes:
        gene_id (str): Gene identifier from GTF
        transcript_id (str): Transcript identifier
        protein_id (str): Protein identifier (if available)
        locus_tag (str): Locus tag (if available)
        name (str): Name used in BED output (based on prioritization)
        chrom (str): Chromosome/sequence name
        start (int): 0-based transcript start position (for BED)
        end (int): Transcript end position
        strand (str): Strand ("+" or "-")
        cds_start (int): 0-based CDS start position
        cds_end (int): CDS end position
        features (list): All GTF features associated with this transcript
        exons (list): List of exon coordinates (start, end)
        gene_type (str): Gene biotype/type
        transcript_type (str): Transcript biotype/type
        has_cds (bool): Whether the transcript has CDS features
        is_noncoding (bool): Whether the transcript is non-coding
        is_pseudogene (bool): Whether the transcript is a pseudogene
    """
    def __init__(self, gene_id, transcript_id):
        self.gene_id = gene_id
        self.transcript_id = transcript_id
        self.protein_id = None
        self.locus_tag = None
        self.name = None  # Will be set based on available IDs
        self.chrom = None
        self.start = None  # 0-based for BED
        self.end = None
        self.strand = None
        self.cds_start = None
        self.cds_end = None
        self.features = []  # List of features (exon, CDS, etc.)
        self.exons = []  # List of exon coordinates (start, end)
        self.gene_type = None
        self.transcript_type = None
        self.has_cds = False
        self.is_noncoding = False
        self.is_pseudogene = False
        
    def add_feature(self, feature_type, start, end, strand, attributes):
        """Add a GTF feature to the transcript and update transcript properties.
        
        This method processes a GTF feature line and:
        1. Tracks transcript boundaries
        2. Records exon coordinates
        3. Tracks CDS boundaries
        4. Extracts identifiers and type information
        5. Detects pseudogenes and non-coding transcripts
        
        Args:
            feature_type (str): GTF feature type (exon, CDS, etc.)
            start (int): 1-based feature start position from GTF
            end (int): Feature end position
            strand (str): Strand ("+" or "-")
            attributes (dict): Parsed attributes from GTF line
        """
        self.features.append((feature_type, start, end, strand, attributes))
        
        # Track transcript boundaries (converting 1-based GTF to 0-based BED)
        start_0based = start - 1
        if self.start is None or start_0based < self.start:
            self.start = start_0based
        if self.end is None or end > self.end:
            self.end = end
            
        # Set chromosome and strand if not set
        if self.chrom is None and 'seq' in attributes:
            self.chrom = attributes['seq']
        if self.strand is None:
            self.strand = strand
            
        # Store exons separately
        if feature_type == 'exon':
            self.exons.append((start_0based, end))
            
        # Track CDS boundaries for thickStart/thickEnd
        if feature_type in ('CDS', 'start_codon', 'stop_codon'):
            self.has_cds = True
            if self.cds_start is None or start_0based < self.cds_start:
                self.cds_start = start_0based
            if self.cds_end is None or end > self.cds_end:
                self.cds_end = end
                
        # Store identifiers
        if 'protein_id' in attributes and not self.protein_id:
            self.protein_id = attributes['protein_id'].split('.')[0]  # Remove version
        if 'locus_tag' in attributes and not self.locus_tag:
            self.locus_tag = attributes['locus_tag']
            
        # Check for pseudogene or non-coding status from gene_biotype or gene_type
        if 'gene_biotype' in attributes:
            self.gene_type = attributes['gene_biotype']
            if attributes['gene_biotype'] == 'pseudogene':
                self.is_pseudogene = True
                self.is_noncoding = True
            elif attributes['gene_biotype'] != 'protein_coding':
                self.is_noncoding = True
        if 'gene_type' in attributes and not self.gene_type:
            self.gene_type = attributes['gene_type']
            if attributes['gene_type'] == 'pseudogene':
                self.is_pseudogene = True
                self.is_noncoding = True
            elif attributes['gene_type'] != 'protein_coding':
                self.is_noncoding = True
                
        # Also check for pseudo attribute
        if 'pseudo' in attributes and attributes['pseudo'] == 'true':
            self.is_pseudogene = True
            self.is_noncoding = True
                
        if 'transcript_type' in attributes and not self.transcript_type:
            self.transcript_type = attributes['transcript_type']
        if 'transcript_biotype' in attributes and not self.transcript_type:
            self.transcript_type = attributes['transcript_biotype']
            
    def finalize(self):
        """Finalize transcript for BED output.
        
        This method:
        1. Sets the final name based on ID priority
        2. Ensures exons are defined (creates if missing)
        3. Sorts exons by position
        4. Handles CDS boundaries based on transcript type
           - For pseudogenes/non-coding: sets to transcript end
           - For coding genes: uses actual CDS boundaries
        """
        # Set name based on available identifiers with prioritization
        if self.gene_id:
            self.name = self.gene_id
        elif self.locus_tag:
            self.name = self.locus_tag
        else:
            self.name = self.transcript_id
        
        # If no exons defined, create one using the transcript boundaries
        if not self.exons:
            self.exons = [(self.start, self.end)]
            
        # Sort exons by position
        self.exons.sort()
        
        # Handle CDS boundaries
        # For pseudogenes and non-coding RNAs, always set both cds_start and cds_end to transcript end
        # This matches Ensembl's convention, even if CDS features exist
        if self.is_noncoding or self.is_pseudogene:
            self.cds_start = self.end
            self.cds_end = self.end
        # For transcripts without CDS features, also set to transcript end
        elif not self.has_cds:
            self.cds_start = self.end
            self.cds_end = self.end
    
    def to_bed12(self, chromosome_prefix=""):
        """Convert transcript to BED12 format.
        
        Generates a BED12 line for the transcript with:
        - Proper exon blocks
        - CDS boundaries (thickStart/thickEnd)
        - Consistent chromosome naming
        
        Args:
            chromosome_prefix (str): Optional prefix to override chromosome name
            
        Returns:
            str: Tab-delimited BED12 line
        """
        self.finalize()
        
        # Calculate block sizes and starts
        block_count = len(self.exons)
        block_sizes = []
        block_starts = []
        
        for exon_start, exon_end in self.exons:
            block_sizes.append(exon_end - exon_start)
            block_starts.append(exon_start - self.start)  # Relative to start
            
        # Format BED12 line
        chrom = self.chrom
        if chromosome_prefix:
            # Use specified chromosome prefix if provided
            chrom = chromosome_prefix
            
        bed_line = [
            chrom,                                  # chrom
            str(self.start),                        # chromStart (0-based)
            str(self.end),                          # chromEnd
            self.name,                              # name
            "0",                                    # score (always 0)
            self.strand if self.strand else ".",    # strand
            str(self.cds_start),                    # thickStart (CDS start)
            str(self.cds_end),                      # thickEnd (CDS end)
            "0",                                    # itemRgb (always 0)
            str(block_count),                       # blockCount
            ",".join(map(str, block_sizes)) + ",",  # blockSizes (comma-terminated)
            ",".join(map(str, block_starts)) + ","  # blockStarts (comma-terminated)
        ]
        
        return "\t".join(bed_line)

def convert_gtf_to_bed(in_file, out_file, chromosome_prefix=""):
    """Convert GTF to BED format in a general way that works with both NCBI and Ensembl.
    
    This function:
    1. Reads the GTF file line by line
    2. Parses features and groups them by transcript
    3. Creates Transcript objects to store feature information
    4. Processes each transcript into BED12 format
    5. Writes the output BED file
    
    Args:
        in_file (str): Path to input GTF file
        out_file (str): Path to output BED file
        chromosome_prefix (str): Optional prefix to override chromosome names
    """
    # Dictionary to store transcripts by ID
    transcripts = {}
    
    # Process GTF file
    print(f"Reading GTF file: {in_file}")
    with open(in_file, 'r') as f:
        for line_num, line in enumerate(f):
            if line_num % 10000 == 0 and line_num > 0:
                print(f"Processed {line_num} lines...")
                
            line = line.strip()
            if not line or line.startswith('#'):
                continue
                
            parts = line.split('\t')
            if len(parts) < 9:
                continue
                
            # Extract fields
            chrom = parts[0]  # Use original chromosome name
            source = parts[1]
            feature_type = parts[2]
            start = int(parts[3])  # 1-based in GTF
            end = int(parts[4])
            score = parts[5]
            strand = parts[6]
            frame = parts[7]
            attr_string = parts[8]
            
            # Parse attributes
            attrs = parse_attributes(attr_string)
            attrs['seq'] = chrom  # Store the chromosome
            
            # Skip if no gene_id
            if 'gene_id' not in attrs:
                continue
                
            gene_id = attrs['gene_id']
            
            # Skip gene features for now, focus on transcript-level features
            if feature_type == 'gene':
                continue
                
            # Get/create transcript_id
            transcript_id = attrs.get('transcript_id', '')
            if transcript_id == '' or transcript_id.startswith('unassigned_transcript'):
                # For genes without proper transcript_id, use protein_id or locus_tag
                if 'protein_id' in attrs:
                    transcript_id = attrs['protein_id'].split('.')[0]
                elif 'locus_tag' in attrs:
                    transcript_id = attrs['locus_tag']
                else:
                    # Last resort, use gene_id
                    transcript_id = gene_id
            
            # Get or create transcript
            key = (gene_id, transcript_id)
            if key not in transcripts:
                transcripts[key] = Transcript(gene_id, transcript_id)
                
            # Add feature to transcript
            transcript = transcripts[key]
            transcript.add_feature(feature_type, start, end, strand, attrs)
    
    # Finalize and write all transcripts
    print(f"Writing BED file: {out_file}")
    with open(out_file, 'w') as f:
        for key, transcript in transcripts.items():
            bed_line = transcript.to_bed12(chromosome_prefix)
            f.write(bed_line + "\n")
            
    print(f"Converted {len(transcripts)} transcripts to BED format")

def main():
    """Parse command line arguments and run the GTF to BED conversion.
    
    Command-line interface for the script, handling argument parsing
    and executing the conversion process.
    """
    parser = argparse.ArgumentParser(description='Convert GTF to BED format')
    parser.add_argument('input_gtf', help='GTF input file')
    parser.add_argument('output_bed', help='BED output file')
    parser.add_argument('--chromosome-prefix', default="", 
                      help='Override chromosome name with this prefix (e.g., "Chromosome")')
    
    args = parser.parse_args()
    convert_gtf_to_bed(args.input_gtf, args.output_bed, args.chromosome_prefix)

if __name__ == "__main__":
    main() 