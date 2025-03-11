#!/usr/bin/env python

import os
import sys
import argparse
import zipfile
import glob
import pandas as pd
import shutil
import tempfile

def parse_args():
    parser = argparse.ArgumentParser(description='Extract RNA-Seq assay table from ISA.zip')
    parser.add_argument('--outdir', required=True, help='Directory containing Metadata folder with ISA.zip')
    parser.add_argument('--mode', required=True, help='Processing mode')
    parser.add_argument('--assay_suffix', required=True, help='Suffix for output file')
    parser.add_argument('--glds_accession', required=True, help='GLDS accession number (e.g. GLDS-123)')
    return parser.parse_args()

def find_file(directory, pattern, error_msg=None):
    """Generic function to find files matching a pattern in a directory."""
    matches = glob.glob(os.path.join(directory, pattern))
    if matches:
        return matches[0]
    if error_msg:
        print(error_msg)
        sys.exit(1)
    return None

def find_runsheet(outdir, glds_accession):
    """Find the runsheet file and load it as a dataframe."""
    metadata_dir = os.path.join(outdir, 'Metadata')
    if not os.path.exists(metadata_dir):
        print(f"Warning: Metadata directory {metadata_dir} not found")
        return None
    
    # Try different patterns to find the runsheet
    patterns = [
        f"*{glds_accession}*runsheet*.csv",
        f"*{glds_accession.replace('GLDS-', '')}*runsheet*.csv",
        "*runsheet*.csv"
    ]
    
    runsheet_file = None
    for pattern in patterns:
        runsheet_file = find_file(metadata_dir, pattern)
        if runsheet_file:
            break
    
    if not runsheet_file:
        print(f"Warning: Could not find runsheet in {metadata_dir}")
        return None
    
    try:
        print(f"Found runsheet: {runsheet_file}")
        runsheet_df = pd.read_csv(runsheet_file)
        print(f"Runsheet has {len(runsheet_df)} rows and {len(runsheet_df.columns)} columns")
        return runsheet_df
    except Exception as e:
        print(f"Error reading runsheet: {e}")
        return None

def is_paired_end_data(runsheet_df):
    """Determine if the data is paired-end based on the runsheet."""
    if runsheet_df is None:
        return False
    
    # Check if the paired_end column exists
    if 'paired_end' in runsheet_df.columns:
        # Convert to string to handle different formats, then check for truthy values
        paired_end_values = runsheet_df['paired_end'].astype(str).str.lower()
        return any(val in ['true', 'yes', '1', 't', 'y'] for val in paired_end_values)
    
    # If no explicit paired_end column, check for file path patterns
    if 'read1_path' in runsheet_df.columns and 'read2_path' in runsheet_df.columns:
        # If read2_path has values, it's paired-end
        return not runsheet_df['read2_path'].isna().all()
    
    return False

def extract_and_find_assay(outdir, glds_accession):
    """Extract ISA.zip and find the RNA-Seq assay file."""
    # Setup paths
    metadata_dir = os.path.join(outdir, 'Metadata')
    if not os.path.exists(metadata_dir):
        print(f"Error: Metadata directory {metadata_dir} not found")
        sys.exit(1)

    # Find and extract ISA zip
    isa_zip = find_file(
        metadata_dir, 
        '*ISA*.zip', 
        error_msg=f"Error: No ISA zip found in {metadata_dir}"
    ) or find_file(metadata_dir, '*.zip')
    
    print(f"Found ISA zip file: {isa_zip}")
    
    # Create and use a temporary directory for extraction
    with tempfile.TemporaryDirectory() as extraction_dir:
        print(f"Extracting to temporary directory: {extraction_dir}")
        
        with zipfile.ZipFile(isa_zip, 'r') as zip_ref:
            zip_ref.extractall(extraction_dir)
        
        # Find investigation file
        i_file = find_file(
            extraction_dir, 
            'i_*.txt', 
            error_msg=f"Error: No investigation file found in {extraction_dir}"
        )
        
        # Try to find RNA-Seq assay file directly by filename pattern
        assay_file = (
            find_file(extraction_dir, 'a_*rna-seq*.txt') or 
            find_file(extraction_dir, 'a_*transcription*.txt')
        )
        
        if assay_file:
            print(f"Found RNA-Seq assay file by filename: {assay_file}")
            # Need to read the file here since the temp dir will be gone
            return pd.read_csv(assay_file, sep='\t')
        
        # Parse investigation file to find assay file
        with open(i_file, 'r') as f:
            for line in f:
                if line.startswith('Study Assay File Name'):
                    # Extract assay filenames
                    parts = line.strip().split('\t')
                    if len(parts) > 1:
                        assay_files = [p.strip() for p in parts[1:] if p.strip()]
                        
                        # First check for RNA-Seq in filename
                        for file in assay_files:
                            if 'rna-seq' in file.lower() or 'transcription' in file.lower():
                                full_path = os.path.join(extraction_dir, file)
                                if os.path.exists(full_path):
                                    print(f"Found RNA-Seq assay file: {full_path}")
                                    return pd.read_csv(full_path, sep='\t')
                        
                        # If not found, use first assay file
                        if assay_files:
                            first_assay = os.path.join(extraction_dir, assay_files[0])
                            if os.path.exists(first_assay):
                                print(f"Using first assay file: {first_assay}")
                                return pd.read_csv(first_assay, sep='\t')
                        
                        # Try pattern matching
                        assay_file = find_file(extraction_dir, 'a_*.txt')
                        if assay_file:
                            print(f"Using first available assay file: {assay_file}")
                            return pd.read_csv(assay_file, sep='\t')
        
        # Last resort - find any assay file
        assay_file = find_file(
            extraction_dir, 
            'a_*.txt',
            error_msg=f"Error: No assay file found in {extraction_dir}"
        )
        
        if assay_file:
            return pd.read_csv(assay_file, sep='\t')
        else:
            print(f"Error: No valid assay file found in {extraction_dir}")
            print(f"Files in directory: {os.listdir(extraction_dir)}")
            sys.exit(1)

def add_parameter_column(df, column_name, value, prefix=None):
    """Add a parameter column to the dataframe if it doesn't exist already.
    
    Args:
        df: The assay table dataframe
        column_name: The parameter column name to add (e.g., "Parameter Value[Entry]")
        value: The value to set for all rows
        prefix: Optional prefix to add to the value (e.g., "GLDS-123_rna_seq_")
        
    Returns:
        The modified dataframe
    """
    # Apply prefix to value if provided
    if prefix and isinstance(value, str):
        # Check if value already has the prefix
        if not value.startswith(prefix):
            prefixed_value = f"{prefix}{value}"
        else:
            prefixed_value = value
    else:
        prefixed_value = value
    
    if column_name not in df.columns:
        print(f"Adding new column: {column_name}")
        df[column_name] = prefixed_value
    else:
        print(f"Column {column_name} already exists, updating values")
        df[column_name] = prefixed_value
    
    return df

def add_unmapped_reads_column(df, glds_prefix, runsheet_df=None):
    """Add the Unmapped Reads column to the dataframe.
    
    Args:
        df: The assay table dataframe
        glds_prefix: The GLDS prefix to add to filenames
        runsheet_df: Optional runsheet dataframe with sample information
        
    Returns:
        The modified dataframe
    """
    column_name = "Parameter Value[Unmapped Reads]"
    
    # Determine if paired-end from runsheet
    is_paired_end = is_paired_end_data(runsheet_df)
    print(f"Data is {'paired-end' if is_paired_end else 'single-end'} based on runsheet")
    
    # Find the sample name column in the assay table
    sample_col = next((col for col in df.columns if 'Sample Name' in col), None)
    
    if not sample_col:
        print("Warning: Could not find Sample Name column in assay table")
        # If no sample column, just use placeholder values
        if is_paired_end:
            values = [f"{glds_prefix}sample{i+1}.unmapped.fastq.1.gz,{glds_prefix}sample{i+1}.unmapped.fastq.2.gz" for i in range(len(df))]
        else:
            values = [f"{glds_prefix}sample{i+1}.unmapped.fastq.gz" for i in range(len(df))]
    else:
        # Get sample names from assay table
        assay_sample_names = df[sample_col].tolist()
        
        # Create a mapping from assay table sample names to runsheet sample names if available
        sample_name_map = {}
        if runsheet_df is not None and 'Sample Name' in runsheet_df.columns:
            # Check for 'Original Sample Name' column to map between assay table and runsheet
            if 'Original Sample Name' in runsheet_df.columns:
                for _, row in runsheet_df.iterrows():
                    orig_name = row['Original Sample Name']
                    rs_name = row['Sample Name']
                    if orig_name in assay_sample_names:
                        sample_name_map[orig_name] = rs_name
        
        # Generate file paths using the appropriate sample names
        values = []
        for assay_sample in assay_sample_names:
            # Use mapped name if available, otherwise use the assay table name
            sample = sample_name_map.get(assay_sample, assay_sample)
            
            if is_paired_end:
                # For paired-end data, create entries with both .1 and .2 files
                values.append(f"{glds_prefix}{sample}.unmapped.fastq.1.gz,{glds_prefix}{sample}.unmapped.fastq.2.gz")
            else:
                # For single-end data
                values.append(f"{glds_prefix}{sample}.unmapped.fastq.gz")
    
    # Add the column to the dataframe
    if column_name not in df.columns:
        print(f"Adding column: {column_name}")
        df[column_name] = values
    else:
        print(f"Updating column: {column_name}")
        df[column_name] = values
    
    return df

def add_trimmed_data_column(df, glds_prefix, runsheet_df=None):
    """Add the Trimmed Sequence Data column to the dataframe.
    
    Args:
        df: The assay table dataframe
        glds_prefix: The GLDS prefix to add to filenames
        runsheet_df: Optional runsheet dataframe with sample information
        
    Returns:
        The modified dataframe
    """
    column_name = "Parameter Value[Trimmed Sequence Data]"
    
    # Determine if paired-end from runsheet
    is_paired_end = is_paired_end_data(runsheet_df)
    print(f"Data is {'paired-end' if is_paired_end else 'single-end'} based on runsheet")
    
    # Find the sample name column in the assay table
    sample_col = next((col for col in df.columns if 'Sample Name' in col), None)
    
    if not sample_col:
        print("Warning: Could not find Sample Name column in assay table")
        # If no sample column, just use placeholder values
        if is_paired_end:
            values = [f"{glds_prefix}sample{i+1}_R1_trimmed.fastq.gz,{glds_prefix}sample{i+1}_R2_trimmed.fastq.gz" for i in range(len(df))]
        else:
            values = [f"{glds_prefix}sample{i+1}_trimmed.fastq.gz" for i in range(len(df))]
    else:
        # Get sample names from assay table
        assay_sample_names = df[sample_col].tolist()
        
        # Create a mapping from assay table sample names to runsheet sample names if available
        sample_name_map = {}
        if runsheet_df is not None and 'Sample Name' in runsheet_df.columns:
            # Check for 'Original Sample Name' column to map between assay table and runsheet
            if 'Original Sample Name' in runsheet_df.columns:
                for _, row in runsheet_df.iterrows():
                    orig_name = row['Original Sample Name']
                    rs_name = row['Sample Name']
                    if orig_name in assay_sample_names:
                        sample_name_map[orig_name] = rs_name
        
        # Generate file paths using the appropriate sample names
        values = []
        for assay_sample in assay_sample_names:
            # Use mapped name if available, otherwise use the assay table name
            sample = sample_name_map.get(assay_sample, assay_sample)
            
            if is_paired_end:
                # For paired-end data, create entries with both R1 and R2 files, comma-separated without spaces
                values.append(f"{glds_prefix}{sample}_R1_trimmed.fastq.gz,{glds_prefix}{sample}_R2_trimmed.fastq.gz")
            else:
                # For single-end data
                values.append(f"{glds_prefix}{sample}_trimmed.fastq.gz")
    
    # Add the column to the dataframe
    if column_name not in df.columns:
        print(f"Adding column: {column_name}")
        df[column_name] = values
    else:
        print(f"Updating column: {column_name}")
        df[column_name] = values
    
    return df

def add_trimming_reports_column(df, glds_prefix, runsheet_df=None):
    """Add the Trimming Reports column to the dataframe.
    
    Args:
        df: The assay table dataframe
        glds_prefix: The GLDS prefix to add to filenames
        runsheet_df: Optional runsheet dataframe with sample information
        
    Returns:
        The modified dataframe
    """
    column_name = "Parameter Value[Trimmed Sequence Data/Trimming Reports]"
    
    # Determine if paired-end from runsheet
    is_paired_end = is_paired_end_data(runsheet_df)
    print(f"Data is {'paired-end' if is_paired_end else 'single-end'} based on runsheet")
    
    # Find the sample name column in the assay table
    sample_col = next((col for col in df.columns if 'Sample Name' in col), None)
    
    if not sample_col:
        print("Warning: Could not find Sample Name column in assay table")
        # If no sample column, just use placeholder values
        if is_paired_end:
            values = [f"{glds_prefix}sample{i+1}_R1_raw.fastq.gz_trimming_report.txt,{glds_prefix}sample{i+1}_R2_raw.fastq.gz_trimming_report.txt" for i in range(len(df))]
        else:
            values = [f"{glds_prefix}sample{i+1}_raw.fastq.gz_trimming_report.txt" for i in range(len(df))]
    else:
        # Get sample names from assay table
        assay_sample_names = df[sample_col].tolist()
        
        # Create a mapping from assay table sample names to runsheet sample names if available
        sample_name_map = {}
        if runsheet_df is not None and 'Sample Name' in runsheet_df.columns:
            # Check for 'Original Sample Name' column to map between assay table and runsheet
            if 'Original Sample Name' in runsheet_df.columns:
                for _, row in runsheet_df.iterrows():
                    orig_name = row['Original Sample Name']
                    rs_name = row['Sample Name']
                    if orig_name in assay_sample_names:
                        sample_name_map[orig_name] = rs_name
        
        # Generate file paths using the appropriate sample names
        values = []
        for assay_sample in assay_sample_names:
            # Use mapped name if available, otherwise use the assay table name
            sample = sample_name_map.get(assay_sample, assay_sample)
            
            if is_paired_end:
                # For paired-end data, create report entries for both R1 and R2 files
                values.append(f"{glds_prefix}{sample}_R1_raw.fastq.gz_trimming_report.txt,{glds_prefix}{sample}_R2_raw.fastq.gz_trimming_report.txt")
            else:
                # For single-end data
                values.append(f"{glds_prefix}{sample}_raw.fastq.gz_trimming_report.txt")
    
    # Add the column to the dataframe
    if column_name not in df.columns:
        print(f"Adding column: {column_name}")
        df[column_name] = values
    else:
        print(f"Updating column: {column_name}")
        df[column_name] = values
    
    return df

def add_multiqc_reports_column(df, glds_prefix, assay_suffix):
    """Add the Multiqc Reports column to the dataframe.
    
    Args:
        df: The assay table dataframe
        glds_prefix: The GLDS prefix to add to filenames
        assay_suffix: The suffix for the multiqc report
        
    Returns:
        The modified dataframe
    """
    column_name = "Parameter Value[Trimmed Sequence Data/Multiqc Reports]"
    
    # The multiqc report is the same for all samples
    report_value = f"{glds_prefix}trimmed_multiqc{assay_suffix}_report.zip"
    
    # Add the column to the dataframe with the same value for all rows
    if column_name not in df.columns:
        print(f"Adding column: {column_name}")
        df[column_name] = report_value
    else:
        print(f"Updating column: {column_name}")
        df[column_name] = report_value
    
    return df

def add_aligned_sequence_data_column(df, glds_prefix, runsheet_df=None):
    """Add the aligned sequence data column to the dataframe."""
    # Title Case column name
    column_name = "Parameter Value[Aligned Sequence Data]"
    
    # Find the sample name column in the assay table
    sample_col = next((col for col in df.columns if 'Sample Name' in col), None)
    
    if not sample_col:
        print("Warning: Could not find Sample Name column in assay table")
        # If no sample column, just use placeholder values
        values = [f"{glds_prefix}sample{i+1}_sorted.bam,{glds_prefix}sample{i+1}_sorted.bam.bai" for i in range(len(df))]
    else:
        # Get sample names from assay table
        assay_sample_names = df[sample_col].tolist()
        
        # Create a mapping from assay table sample names to runsheet sample names if available
        sample_name_map = {}
        if runsheet_df is not None and 'Sample Name' in runsheet_df.columns:
            # Check for 'Original Sample Name' column to map between assay table and runsheet
            if 'Original Sample Name' in runsheet_df.columns:
                for _, row in runsheet_df.iterrows():
                    orig_name = row['Original Sample Name']
                    rs_name = row['Sample Name']
                    if orig_name in assay_sample_names:
                        sample_name_map[orig_name] = rs_name
        
        # Generate file paths using the appropriate sample names
        values = []
        for assay_sample in assay_sample_names:
            # Use mapped name if available, otherwise use the assay table name
            sample = sample_name_map.get(assay_sample, assay_sample)
            
            # Create only the BAM and BAI files, not the log
            files = [
                f"{glds_prefix}{sample}_sorted.bam",
                f"{glds_prefix}{sample}_sorted.bam.bai"
            ]
            values.append(",".join(files))
    
    # Add the column to the dataframe
    if column_name not in df.columns:
        print(f"Adding column: {column_name}")
        df[column_name] = values
    else:
        print(f"Updating column: {column_name}")
        df[column_name] = values
    
    return df

def add_alignment_logs_column(df, glds_prefix, runsheet_df=None):
    """Add the Alignment Logs column to the dataframe."""
    column_name = "Parameter Value[Aligned Sequence Data/Alignment Logs]"
    
    # Find the sample name column in the assay table
    sample_col = next((col for col in df.columns if 'Sample Name' in col), None)
    
    if not sample_col:
        print("Warning: Could not find Sample Name column in assay table")
        # If no sample column, just use placeholder values
        values = [f"{glds_prefix}sample{i+1}.bowtie2.log" for i in range(len(df))]
    else:
        # Get sample names from assay table
        assay_sample_names = df[sample_col].tolist()
        
        # Create a mapping from assay table sample names to runsheet sample names if available
        sample_name_map = {}
        if runsheet_df is not None and 'Sample Name' in runsheet_df.columns:
            # Check for 'Original Sample Name' column to map between assay table and runsheet
            if 'Original Sample Name' in runsheet_df.columns:
                for _, row in runsheet_df.iterrows():
                    orig_name = row['Original Sample Name']
                    rs_name = row['Sample Name']
                    if orig_name in assay_sample_names:
                        sample_name_map[orig_name] = rs_name
        
        # Generate file paths using the appropriate sample names
        values = []
        for assay_sample in assay_sample_names:
            # Use mapped name if available, otherwise use the assay table name
            sample = sample_name_map.get(assay_sample, assay_sample)
            
            # Create just the log file
            values.append(f"{glds_prefix}{sample}.bowtie2.log")
    
    # Add the column to the dataframe
    if column_name not in df.columns:
        print(f"Adding column: {column_name}")
        df[column_name] = values
    else:
        print(f"Updating column: {column_name}")
        df[column_name] = values
    
    return df

def main():
    args = parse_args()
    
    # Find and parse the runsheet
    runsheet_df = find_runsheet(args.outdir, args.glds_accession)
    
    # Find RNA-Seq assay file and get its contents
    assay_df = extract_and_find_assay(args.outdir, args.glds_accession)
    print(f"Original assay table has {len(assay_df)} rows and {len(assay_df.columns)} columns")
    
    # Process and save assay file
    try:
        # Create GLDS prefix for filenames
        glds_prefix = f"{args.glds_accession.upper()}_rna_seq_"
        
        # Add Trimmed Sequence Data column using information from the runsheet
        assay_df = add_trimmed_data_column(assay_df, glds_prefix, runsheet_df=runsheet_df)
        
        # Add Multiqc Reports column first
        assay_df = add_multiqc_reports_column(assay_df, glds_prefix, args.assay_suffix)
        
        # Add Trimming Reports column second
        assay_df = add_trimming_reports_column(assay_df, glds_prefix, runsheet_df=runsheet_df)
        
        # Add aligned sequence data column
        assay_df = add_aligned_sequence_data_column(assay_df, glds_prefix, runsheet_df=runsheet_df)
        
        # Add alignment logs column
        assay_df = add_alignment_logs_column(assay_df, glds_prefix, runsheet_df=runsheet_df)
        
        # Add unmapped reads column
        assay_df = add_unmapped_reads_column(assay_df, glds_prefix, runsheet_df=runsheet_df)
        
        # Get the original assay name from the ISA zip for output filename
        metadata_dir = os.path.join(args.outdir, 'Metadata')
        isa_zip = find_file(metadata_dir, '*ISA*.zip') or find_file(metadata_dir, '*.zip')
        
        # Default filename if we can't find anything else
        orig_filename = f"a_{args.glds_accession}_assay.txt"
        
        if isa_zip:
            with zipfile.ZipFile(isa_zip, 'r') as zip_ref:
                file_list = zip_ref.namelist()
                assay_files = [f for f in file_list if os.path.basename(f).startswith('a_')]
                if assay_files:
                    orig_filename = os.path.basename(assay_files[0])
        
        # Create both original and modified output files
        # Original file (preserving the original name)
        assay_df.to_csv(orig_filename, sep='\t', index=False)
        print(f"Original assay table saved as: {orig_filename}")
        
        # Modified file with GLDS prefix
        if not orig_filename.startswith(glds_prefix):
            mod_filename = f"{glds_prefix}{orig_filename}"
        else:
            mod_filename = orig_filename
            
        assay_df.to_csv(mod_filename, sep='\t', index=False)
        print(f"Modified assay table saved as: {mod_filename}")
        
    except Exception as e:
        print(f"Error processing assay file: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
