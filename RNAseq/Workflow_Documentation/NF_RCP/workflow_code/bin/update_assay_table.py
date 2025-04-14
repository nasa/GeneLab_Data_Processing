#!/usr/bin/env python
"""
RNA-Seq Assay Table Updater for the NASA GeneLab Data Processing Pipeline

This script processes and updates RNA-Seq assay tables from ISA.zip files within a GLDS dataset.
It adds or updates various parameter columns required for the GeneLab RNA-Seq processing pipeline,
including paths to processed files like trimmed sequence data, aligned sequence data, raw counts,
normalized counts, and differential expression analysis results.

The script supports both default mode (STAR/RSEM) and microbes mode (Bowtie2/FeatureCounts).
It also handles both paired-end and single-end sequencing data and can detect the presence
of ERCC spike-ins to conditionally add related columns.

Usage:
    python update_assay_table.py --outdir <directory> --assay_suffix <suffix> --glds_accession <GLDS-XXX> [--mode <mode>]

Parameters:
    --outdir: Directory containing the Metadata folder with ISA.zip
    --assay_suffix: Suffix to append to output filenames (e.g., "_GLbulkRNAseq")
    --glds_accession: GLDS accession number (e.g., "GLDS-123")
    --mode: Processing mode - "microbes" for microbial datasets or empty for default

The script automatically:
1. Extracts and locates the RNA-Seq assay table from the ISA.zip file
2. Loads the runsheet if available for additional metadata
3. Detects if data is paired-end or single-end
4. Detects if ERCC spike-ins are used
5. Updates the assay table with appropriate file paths for all processing outputs
6. Saves the updated assay table using the original filename
"""

import os
import sys
import argparse
import zipfile
import glob
import pandas as pd
import shutil
import tempfile
import json
import re

def parse_args():
    parser = argparse.ArgumentParser(description='Extract RNA-Seq assay table from ISA.zip')
    parser.add_argument('--outdir', required=True, help='Directory containing Metadata folder with ISA.zip')
    parser.add_argument('--mode', required=False, default="", help='Processing mode (microbes or empty for default)')
    parser.add_argument('--assay_suffix', required=True, help='Suffix for output file')
    parser.add_argument('--glds_accession', required=True, help='GLDS accession number (e.g. GLDS-123)')
    return parser.parse_args()

# Create a global variable to track changes
column_changes = []

def find_column_case_insensitive(df, column_name):
    """Find a column name in the dataframe case-insensitively.
    
    Args:
        df: The DataFrame to search in
        column_name: The column name to find
        
    Returns:
        The actual column name if found, None otherwise
    """
    # Convert column name to lowercase for comparison
    column_lower = column_name.lower()
    
    # Check if any existing column matches (case-insensitive)
    for col in df.columns:
        if col.lower() == column_lower:
            return col
    
    # If no match found
    return None

def update_column(df, column_name, values):
    """Add or update a column in the dataframe, with case-insensitive matching.
    
    Args:
        df: The DataFrame to update
        column_name: The column name to add or update
        values: The values to set
        
    Returns:
        The updated DataFrame
    """
    # Find if column exists (case-insensitive)
    existing_col = find_column_case_insensitive(df, column_name)
    
    if existing_col:
        # Column exists (might be different case)
        print(f"Updating column: {existing_col}")
        df[existing_col] = values
        column_changes.append(f"Updated: {existing_col}")
    else:
        # Column doesn't exist
        print(f"Adding new column: {column_name}")
        df[column_name] = values
        column_changes.append(f"Added: {column_name}")
    
    return df

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
    """Determine if this is paired-end data from the runsheet.
    
    Args:
        runsheet_df: DataFrame containing the runsheet
        
    Returns:
        bool: True if paired-end, False if single-end
    """
    if runsheet_df is None:
        print("Warning: No runsheet provided, assuming single-end data")
        return False
    
    # Check if there's a paired_end column
    if 'paired_end' in runsheet_df.columns:
        paired_end_values = runsheet_df['paired_end'].unique()
        if len(paired_end_values) == 1:
            # If all values are the same, use that
            value = paired_end_values[0]
            # Handle different types of values (string or boolean)
            if isinstance(value, bool):
                return value
            elif isinstance(value, str):
                return value.lower() == 'true'
            else:
                # Try to convert to string if it's not a boolean or string
                return str(value).lower() == 'true'
    
    # Check based on R1/R2 file presence
    if any(col for col in runsheet_df.columns if col.endswith('_R2.fastq.gz')):
        print("Detected paired-end data based on _R2 files in runsheet")
        return True
    else:
        print("Assuming single-end data (no _R2 files in runsheet)")
        return False

def has_ercc_spikes(runsheet_df):
    """Determine if ERCC spike-ins are used based on the runsheet."""
    if runsheet_df is None:
        print("Warning: No runsheet provided, assuming no ERCC spike-ins")
        return False
    
    # Check for either has_ercc or has_ERCC column (case-insensitive)
    ercc_column = None
    for col in runsheet_df.columns:
        if col.lower() == 'has_ercc':
            ercc_column = col
            break
            
    if ercc_column:
        has_ercc_values = runsheet_df[ercc_column].unique()
        if len(has_ercc_values) == 1:
            # If all values are the same, use that
            value = has_ercc_values[0]
            # Handle different types of values (string or boolean)
            if isinstance(value, bool):
                ercc_used = value
            elif isinstance(value, str):
                ercc_used = value.lower() == 'true'
            else:
                # Try to convert to string if it's not a boolean or string
                ercc_used = str(value).lower() == 'true'
            
            print(f"ERCC spike-ins {'are' if ercc_used else 'are not'} used based on runsheet")
            return ercc_used
        
    # Check for other indicators of ERCC usage
    # Look for ERCC in sample names as a fallback
    if 'Sample Name' in runsheet_df.columns and any('ercc' in str(name).lower() for name in runsheet_df['Sample Name']):
        print("Detected ERCC spike-ins based on sample names")
        return True
    
    print("Assuming no ERCC spike-ins (not indicated in runsheet)")
    return False

def extract_and_find_assay(outdir, glds_accession):
    """Extract ISA.zip and find the RNA-Seq assay file using ONLY investigation file logic."""
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
        
        # Define valid combination for RNA-seq
        valid_measurement = "transcription profiling"
        valid_technology = "RNA Sequencing (RNA-Seq)"
        
        # Parse investigation file to build STUDY ASSAYS table
        study_assays_table = {}
        study_assays_section = False
        
        with open(i_file, 'r') as f:
            for line in f:
                line = line.strip()
                
                # Track STUDY ASSAYS section
                if line == "STUDY ASSAYS":
                    study_assays_section = True
                    continue
                elif study_assays_section and not line:
                    study_assays_section = False
                    continue
                
                # Extract data from section
                if study_assays_section and line:
                    parts = line.split('\t')
                    if parts and parts[0]:
                        key = parts[0]
                        values = [v.strip() for v in parts[1:] if v.strip()]
                        study_assays_table[key] = values
        
        # Check if we have all required keys
        required_keys = ["Study Assay Measurement Type", "Study Assay Technology Type", "Study Assay File Name"]
        if not all(key in study_assays_table for key in required_keys):
            print(f"Error: Missing required keys in STUDY ASSAYS section")
            sys.exit(1)
        
        # Get the values from the table
        measurement_types = study_assays_table["Study Assay Measurement Type"]
        technology_types = study_assays_table["Study Assay Technology Type"]
        file_names = study_assays_table["Study Assay File Name"]
        
        # Ensure all lists have equal length
        if not (len(measurement_types) == len(technology_types) == len(file_names)):
            print(f"Error: Measurement types, technology types, and file names have different lengths")
            sys.exit(1)
        
        # Find matching assay file
        matched_file = None
        for i in range(len(measurement_types)):
            if (measurement_types[i].lower() == valid_measurement.lower() and 
                technology_types[i].lower() == valid_technology.lower()):
                matched_file = file_names[i]
                break
        
        if not matched_file:
            print(f"Error: No assay file matched for RNA-Seq. Measurement types: {measurement_types}, Technology types: {technology_types}")
            sys.exit(1)
        
        # Load the matched assay file
        assay_path = os.path.join(extraction_dir, matched_file)
        if not os.path.exists(assay_path):
            print(f"Error: Matched assay file doesn't exist: {assay_path}")
            sys.exit(1)
        
        print(f"Using RNA-Seq assay file: {assay_path}")
        # Return both the dataframe and the filename
        return pd.read_csv(assay_path, sep='\t'), matched_file

def add_read_counts(df, outdir, glds_accession, assay_suffix, runsheet_df=None):
    """Add the read counts column to the dataframe.
    
    Args:
        df: The assay table dataframe
        outdir: The output directory
        glds_accession: The GLDS accession number
        assay_suffix: The assay suffix for MultiQC report files
        runsheet_df: Optional runsheet dataframe with sample information
        
    Returns:
        The modified dataframe
    """
    column_name = "Parameter Value[Read Count]"
    
    # Determine if paired-end from runsheet
    is_paired_end = is_paired_end_data(runsheet_df)
    print(f"Data is {'paired-end' if is_paired_end else 'single-end'} based on runsheet")
    
    # Path to raw MultiQC data zip
    fastqc_dir = os.path.join(outdir, "00-RawData", "FastQC_Reports")
    multiqc_data_zip = os.path.join(fastqc_dir, f"raw_multiqc{assay_suffix}_data.zip")
    
    if not os.path.exists(multiqc_data_zip):
        print(f"WARNING: MultiQC data zip file not found at {multiqc_data_zip}")
        # If zip not found, just add placeholder values
        df[column_name] = "N/A"
        column_changes.append(f"Added: {column_name} (with placeholder values)")
        return df
    
    print(f"Found MultiQC data zip: {multiqc_data_zip}")
    
    # Find the sample name column in the assay table
    sample_col = next((col for col in df.columns if 'Sample Name' in col), None)
    if not sample_col:
        print("Warning: Could not find Sample Name column in assay table")
        # If no sample column, just add placeholder values
        df[column_name] = "N/A"
        column_changes.append(f"Added: {column_name} (with placeholder values)")
        return df
    
    # Get sample names from assay table
    assay_sample_names = df[sample_col].tolist()
    
    # Get read counts from MultiQC data
    read_counts = {}
    
    # Create a temporary directory to extract files
    with tempfile.TemporaryDirectory() as temp_dir:
        try:
            # Extract the zip file
            with zipfile.ZipFile(multiqc_data_zip, 'r') as zip_ref:
                zip_ref.extractall(temp_dir)
            
            # Expected path to the JSON file in the extracted data directory
            expected_dir = f"raw_multiqc{assay_suffix}_data"
            json_path = os.path.join(temp_dir, expected_dir, "multiqc_data.json")
            
            # Check if JSON file exists in the expected location
            if not os.path.exists(json_path):
                # Fallback: try to find it elsewhere
                print(f"JSON not found at expected path: {json_path}, searching elsewhere...")
                
                # Try directly in the temp directory
                direct_path = os.path.join(temp_dir, "multiqc_data.json")
                if os.path.exists(direct_path):
                    json_path = direct_path
                else:
                    # Search all subdirectories
                    for root, dirs, files in os.walk(temp_dir):
                        if "multiqc_data.json" in files:
                            json_path = os.path.join(root, "multiqc_data.json")
                            print(f"Found JSON at: {json_path}")
                            break
            
            if not os.path.exists(json_path):
                print(f"ERROR: Could not find multiqc_data.json in the extracted zip")
                # Add placeholder values
                df[column_name] = "N/A"
                column_changes.append(f"Added: {column_name} (with placeholder values)")
                return df
            
            print(f"Using MultiQC data from: {json_path}")
            
            # Parse the MultiQC JSON
            with open(json_path, 'r') as f:
                multiqc_data = json.load(f)
            
            # First, try to extract directly from report_general_stats_data
            if 'report_general_stats_data' in multiqc_data and multiqc_data['report_general_stats_data']:
                stats_module = multiqc_data['report_general_stats_data'][0]  # Use first module
                print("Extracting read counts directly from report_general_stats_data")
                
                for sample_name, sample_data in stats_module.items():
                    if 'total_sequences' in sample_data:
                        read_count = int(sample_data['total_sequences'])
                        print(f"Found count for {sample_name}: {read_count}")
                        read_counts[sample_name] = read_count
            
            # Fallback to FastQC module specific extraction if needed
            elif ('report_data_sources' in multiqc_data and 
                'FastQC' in multiqc_data['report_data_sources']):
                
                # Find the index for FastQC in the general stats data
                fastqc_index = None
                for i, module_data in enumerate(multiqc_data.get('report_general_stats_data', [])):
                    if module_data and any('total_sequences' in sample_data for sample_data in module_data.values()):
                        fastqc_index = i
                        break
                
                if fastqc_index is not None:
                    fastqc_stats = multiqc_data['report_general_stats_data'][fastqc_index]
                    
                    # Process each sample to extract read counts
                    for sample_name, sample_data in fastqc_stats.items():
                        if 'total_sequences' in sample_data:
                            read_count = int(sample_data['total_sequences'])
                            print(f"Found count for {sample_name}: {read_count}")
                            read_counts[sample_name] = read_count
            
            if not read_counts:
                print("WARNING: Could not extract any read counts from MultiQC data")
                
            # Debug output to show read counts found
            print(f"Successfully extracted {len(read_counts)} read counts:")
            for sample_name, count in read_counts.items():
                print(f"  - {sample_name}: {count}")
                
        except Exception as e:
            print(f"Error extracting read counts from MultiQC data: {str(e)}")
            # If error occurs, add placeholder values
            df[column_name] = "N/A"
            column_changes.append(f"Added: {column_name} (with placeholder values)")
            return df
    
    # Debug output to show sample names in assay table
    print("Sample names in assay table:")
    for sample_name in assay_sample_names:
        print(f"  - {sample_name}")
    
    # Generate read count values for each sample in the assay table
    values = []
    for assay_sample in assay_sample_names:
        print(f"Looking for read count for sample: {assay_sample}")
        
        # Try direct match first
        if assay_sample in read_counts:
            print(f"Direct match found for {assay_sample}")
            values.append(str(read_counts[assay_sample]))
        else:
            # Try a more flexible match if direct match fails
            found_match = False
            for mqc_sample, count in read_counts.items():
                # Check if assay sample name is contained in MultiQC sample name or vice versa
                if assay_sample in mqc_sample or mqc_sample in assay_sample:
                    values.append(str(count))
                    found_match = True
                    print(f"Flexible match found: {assay_sample} -> {mqc_sample} = {count}")
                    break
            
            if not found_match:
                print(f"WARNING: No read count found for sample {assay_sample}")
                values.append("N/A")
    
    # Add the column to the dataframe
    df = update_column(df, column_name, values)
    
    return df

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

def add_unmapped_reads_column(df, glds_prefix, runsheet_df=None, mode=""):
    """Add the Unmapped Reads column to the dataframe."""
    column_name = "Parameter Value[Unmapped Reads]"
    
    # Determine if paired-end from runsheet
    is_paired_end = is_paired_end_data(runsheet_df)
    print(f"Data is {'paired-end' if is_paired_end else 'single-end'} based on runsheet")
    
    # Find the sample name column in the assay table
    sample_col = next((col for col in df.columns if 'Sample Name' in col), None)
    
    if not sample_col:
        print("Warning: Could not find Sample Name column in assay table")
        # If no sample column, just use placeholder values
        if mode == "microbes":
            # Microbes mode (Bowtie2)
            if is_paired_end:
                values = [f"{glds_prefix}sample{i+1}.unmapped.fastq.1.gz,{glds_prefix}sample{i+1}.unmapped.fastq.2.gz" for i in range(len(df))]
            else:
                values = [f"{glds_prefix}sample{i+1}.unmapped.fastq.gz" for i in range(len(df))]
        else:
            # Default mode (STAR)
            if is_paired_end:
                values = [f"{glds_prefix}sample{i+1}_Unmapped.out.mate1,{glds_prefix}sample{i+1}_Unmapped.out.mate2" for i in range(len(df))]
            else:
                values = [f"{glds_prefix}sample{i+1}_Unmapped.out.mate1" for i in range(len(df))]
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
            
            if mode == "microbes":
                # Microbes mode (Bowtie2)
                if is_paired_end:
                    # For paired-end data, create entries with both .1 and .2 files
                    values.append(f"{glds_prefix}{sample}.unmapped.fastq.1.gz,{glds_prefix}{sample}.unmapped.fastq.2.gz")
                else:
                    # For single-end data
                    values.append(f"{glds_prefix}{sample}.unmapped.fastq.gz")
            else:
                # Default mode (STAR)
                if is_paired_end:
                    # For paired-end data with STAR
                    values.append(f"{glds_prefix}{sample}_Unmapped.out.mate1,{glds_prefix}{sample}_Unmapped.out.mate2")
                else:
                    # For single-end data with STAR
                    values.append(f"{glds_prefix}{sample}_Unmapped.out.mate1")
    
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

def add_trimmed_multiqc_reports_column(df, glds_prefix, assay_suffix):
    """Add the trimmed sequence data MultiQC reports column to the dataframe."""
    column_name = "Parameter Value[Trimmed Sequence Data/MultiQC Reports]"
    
    # Create the multiqc report filenames - both data zip and html
    multiqc_html = f"trimmed_multiqc{assay_suffix}.html"
    multiqc_data = f"trimmed_multiqc{assay_suffix}_data.zip"
    
    # Join the files with commas
    combined_files = f"{multiqc_data},{multiqc_html}"
    
    # Add the column to the dataframe with the same value for all rows
    df = update_column(df, column_name, combined_files)
    
    return df

def add_align_multiqc_reports_column(df, glds_prefix, assay_suffix):
    """Add the aligned sequence data MultiQC reports column to the dataframe."""
    column_name = "Parameter Value[Aligned Sequence Data/MultiQC Reports]"
    
    # Create the alignment multiqc report filenames - both data zip and html
    multiqc_html = f"align_multiqc{assay_suffix}.html"
    multiqc_data = f"align_multiqc{assay_suffix}_data.zip"
    
    # Join the files with commas
    combined_files = f"{multiqc_data},{multiqc_html}"
    
    # Add the column to the dataframe with the same value for all rows
    df = update_column(df, column_name, combined_files)
    
    return df

def add_rseqc_multiqc_reports_column(df, glds_prefix, assay_suffix):
    """Add the RSeQC MultiQC reports column to the dataframe."""
    column_name = "Parameter Value[RSeQC/MultiQC Reports]"
    
    # Create the four RSeQC multiqc report filenames - data zip and html for each
    multiqc_reports = []
    for report_type in ["geneBody_cov", "infer_exp", "inner_dist", "read_dist"]:
        data_zip = f"{report_type}_multiqc{assay_suffix}_data.zip"
        html = f"{report_type}_multiqc{assay_suffix}.html"
        multiqc_reports.extend([data_zip, html])
    
    # Join the reports with commas
    combined_reports = ",".join(multiqc_reports)
    
    # Add the column to the dataframe with the same value for all rows
    df = update_column(df, column_name, combined_reports)
    
    return df

def add_raw_counts_data_column(df, glds_prefix, assay_suffix, mode=""):
    """Add the Raw Counts Data column to the dataframe."""
    column_name = "Parameter Value[Raw Counts Data]"
    
    # Find the sample name column in the assay table
    sample_col = next((col for col in df.columns if 'Sample Name' in col), None)
    
    if not sample_col:
        print("Warning: Could not find Sample Name column in assay table")
        
        if mode == "microbes":
            # Microbes mode (FeatureCounts)
            # Create the feature counts file names - same for all samples
            feature_counts_files = [
                f"{glds_prefix}FeatureCounts{assay_suffix}.tsv",
                f"{glds_prefix}FeatureCounts{assay_suffix}.tsv.summary"
            ]
            # Join the files with commas
            combined_files = ",".join(feature_counts_files)
            values = [combined_files] * len(df)
        else:
            # Default mode (RSEM) - placeholder values
            values = [f"{glds_prefix}sample{i+1}.genes.results,{glds_prefix}sample{i+1}.isoforms.results" for i in range(len(df))]
    else:
        # Get sample names from assay table
        assay_sample_names = df[sample_col].tolist()
        
        # Create a mapping from assay table sample names to runsheet sample names if available
        sample_name_map = {}
        runsheet_df = None  # Define outside the if to avoid UnboundLocalError
        if 'runsheet_df' in locals() and runsheet_df is not None and 'Sample Name' in runsheet_df.columns:
            # Check for 'Original Sample Name' column to map between assay table and runsheet
            if 'Original Sample Name' in runsheet_df.columns:
                for _, row in runsheet_df.iterrows():
                    orig_name = row['Original Sample Name']
                    rs_name = row['Sample Name']
                    if orig_name in assay_sample_names:
                        sample_name_map[orig_name] = rs_name
        
        values = []
        if mode == "microbes":
            # Microbes mode (FeatureCounts) - same for all samples
            feature_counts_files = [
                f"{glds_prefix}FeatureCounts{assay_suffix}.tsv",
                f"{glds_prefix}FeatureCounts{assay_suffix}.tsv.summary"
            ]
            # Join the files with commas
            combined_files = ",".join(feature_counts_files)
            values = [combined_files] * len(df)
        else:
            # Default mode (RSEM) - sample-specific files
            for assay_sample in assay_sample_names:
                # Use mapped name if available, otherwise use the assay table name
                sample = sample_name_map.get(assay_sample, assay_sample)
                
                # For RSEM, each sample has genes and isoforms result files
                rsem_files = [
                    f"{glds_prefix}{sample}.genes.results",
                    f"{glds_prefix}{sample}.isoforms.results"
                ]
                values.append(",".join(rsem_files))
    
    # Add the column to the dataframe
    if column_name not in df.columns:
        print(f"Adding column: {column_name}")
        df[column_name] = values
    else:
        print(f"Updating column: {column_name}")
        df[column_name] = values
    
    return df

def add_raw_counts_multiqc_column(df, glds_prefix, assay_suffix, mode=""):
    """Add the Raw Counts Data MultiQC reports column to the dataframe."""
    column_name = "Parameter Value[Raw Counts Data/MultiQC Reports]"
    
    # Create the appropriate multiqc report filenames based on mode - both data zip and html
    if mode == "microbes":
        # Microbes mode (FeatureCounts)
        multiqc_html = f"featureCounts_multiqc{assay_suffix}.html"
        multiqc_data = f"featureCounts_multiqc{assay_suffix}_data.zip"
    else:
        # Default mode (RSEM)
        multiqc_html = f"RSEM_count_multiqc{assay_suffix}.html"
        multiqc_data = f"RSEM_count_multiqc{assay_suffix}_data.zip"
    
    # Join the files with commas
    combined_files = f"{multiqc_data},{multiqc_html}"
    
    # Add the column to the dataframe with the same value for all rows
    df = update_column(df, column_name, combined_files)
    
    return df

def add_raw_multiqc_reports_column(df, glds_prefix, assay_suffix):
    """Add the raw sequence data MultiQC reports column to the dataframe."""
    column_name = "Parameter Value[MultiQC File Names]"
    
    # Create the multiqc report filenames - both data zip and html
    multiqc_html = f"raw_multiqc{assay_suffix}.html"
    multiqc_data = f"raw_multiqc{assay_suffix}_data.zip"
    
    # Join the files with commas
    combined_files = f"{multiqc_data},{multiqc_html}"
    
    # Add the column to the dataframe with the same value for all rows
    df = update_column(df, column_name, combined_files)
    
    return df

def add_raw_counts_tables_column(df, glds_prefix, assay_suffix, mode=""):
    """Add the Raw Counts Tables column to the dataframe."""
    column_name = "Parameter Value[Raw Counts Tables]"
    
    # Create the unnormalized counts filename based on mode
    if mode == "microbes":
        # Microbes mode (FeatureCounts)
        counts_file = f"{glds_prefix}FeatureCounts_Unnormalized_Counts{assay_suffix}.csv"
    else:
        # Default mode (STAR/RSEM)
        rsem_file = f"{glds_prefix}RSEM_Unnormalized_Counts{assay_suffix}.csv"
        star_file = f"{glds_prefix}STAR_Unnormalized_Counts{assay_suffix}.csv"
        counts_file = f"{rsem_file},{star_file}"
    
    # Add the column to the dataframe with the same value for all rows
    df = update_column(df, column_name, counts_file)
    
    return df

def add_raw_counts_tables_rrnarm_column(df, glds_prefix, assay_suffix, mode=""):
    """Add the Raw Counts Tables rRNArm column to the dataframe."""
    column_name = "Parameter Value[Raw Counts Tables rRNArm]"
    
    # Create the unnormalized counts filename based on mode
    if mode == "microbes":
        # Microbes mode (FeatureCounts)
        counts_file = f"{glds_prefix}FeatureCounts_Unnormalized_Counts_rRNArm{assay_suffix}.csv"
    else:
        # Default mode (STAR/RSEM)
        counts_file = f"{glds_prefix}RSEM_Unnormalized_Counts_rRNArm{assay_suffix}.csv"
    
    # Add the column to the dataframe with the same value for all rows
    df = update_column(df, column_name, counts_file)
    
    return df

def add_normalized_counts_data_column(df, glds_prefix, assay_suffix):
    """Add the Normalized Counts Data column to the dataframe."""
    column_name = "Parameter Value[Normalized Counts Data]"
    
    # Create the normalized counts filenames - same for all samples
    normalized_files = [
        f"{glds_prefix}Normalized_Counts{assay_suffix}.csv",
        f"{glds_prefix}VST_Counts{assay_suffix}.csv"
    ]
    
    # Join the files with commas
    combined_files = ",".join(normalized_files)
    
    # Add the column to the dataframe with the same value for all rows
    df = update_column(df, column_name, combined_files)
    
    return df

def add_normalized_counts_data_rrnarm_column(df, glds_prefix, assay_suffix):
    """Add the Normalized Counts Data rRNArm column to the dataframe."""
    column_name = "Parameter Value[Normalized Counts Data rRNArm]"
    
    # Create the normalized counts filenames with rRNArm - same for all samples
    normalized_files = [
        f"{glds_prefix}Normalized_Counts_rRNArm{assay_suffix}.csv",
        f"{glds_prefix}VST_Counts_rRNArm{assay_suffix}.csv"
    ]
    
    # Join the files with commas
    combined_files = ",".join(normalized_files)
    
    # Add the column to the dataframe with the same value for all rows
    df = update_column(df, column_name, combined_files)
    
    return df

def add_differential_expression_column(df, glds_prefix, assay_suffix):
    """Add the Differential Expression Analysis Data column to the dataframe."""
    column_name = "Parameter Value[Differential Expression Analysis Data]"
    
    # Create the differential expression filenames - same for all samples
    de_files = [
        f"{glds_prefix}SampleTable{assay_suffix}.csv",
        f"{glds_prefix}contrasts{assay_suffix}.csv",
        f"{glds_prefix}differential_expression{assay_suffix}.csv"
    ]
    
    # Join the files with commas
    combined_files = ",".join(de_files)
    
    # Add the column to the dataframe with the same value for all rows
    df = update_column(df, column_name, combined_files)
    
    return df

def add_differential_expression_rrnarm_column(df, glds_prefix, assay_suffix):
    """Add the Differential Expression Analysis Data rRNArm column to the dataframe."""
    column_name = "Parameter Value[Differential Expression Analysis Data rRNArm]"
    
    # Create the differential expression filenames with rRNArm - same for all samples
    de_files = [
        f"{glds_prefix}SampleTable_rRNArm{assay_suffix}.csv",
        f"{glds_prefix}contrasts_rRNArm{assay_suffix}.csv",
        f"{glds_prefix}differential_expression_rRNArm{assay_suffix}.csv"
    ]
    
    # Join the files with commas
    combined_files = ",".join(de_files)
    
    # Add the column to the dataframe with the same value for all rows
    df = update_column(df, column_name, combined_files)
    
    return df

def add_aligned_sequence_data_column(df, glds_prefix, runsheet_df=None, mode=""):
    """Add the aligned sequence data column to the dataframe."""
    # Title Case column name
    column_name = "Parameter Value[Aligned Sequence Data]"
    
    # Find the sample name column in the assay table
    sample_col = next((col for col in df.columns if 'Sample Name' in col), None)
    
    if not sample_col:
        print("Warning: Could not find Sample Name column in assay table")
        # If no sample column, just use placeholder values
        if mode == "microbes":
            # Microbes mode (Bowtie2)
            values = [f"{glds_prefix}sample{i+1}_sorted.bam,{glds_prefix}sample{i+1}_sorted.bam.bai" for i in range(len(df))]
        else:
            # Default mode (STAR)
            values = [
                f"{glds_prefix}sample{i+1}_Aligned.sortedByCoord_sorted.out.bam,"
                f"{glds_prefix}sample{i+1}_Aligned.sortedByCoord_sorted.out.bam.bai,"
                f"{glds_prefix}sample{i+1}_Aligned.toTranscriptome.out.bam,"
                f"{glds_prefix}sample{i+1}_SJ.out.tab" 
                for i in range(len(df))
            ]
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
            
            if mode == "microbes":
                # Microbes mode (Bowtie2)
                files = [
                    f"{glds_prefix}{sample}_sorted.bam",
                    f"{glds_prefix}{sample}_sorted.bam.bai"
                ]
            else:
                # Default mode (STAR)
                files = [
                    f"{glds_prefix}{sample}_Aligned.sortedByCoord_sorted.out.bam",
                    f"{glds_prefix}{sample}_Aligned.sortedByCoord_sorted.out.bam.bai",
                    f"{glds_prefix}{sample}_Aligned.toTranscriptome.out.bam",
                    f"{glds_prefix}{sample}_SJ.out.tab"
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

def add_alignment_logs_column(df, glds_prefix, runsheet_df=None, mode=""):
    """Add the Alignment Logs column to the dataframe."""
    column_name = "Parameter Value[Aligned Sequence Data/Alignment Logs]"
    
    # Find the sample name column in the assay table
    sample_col = next((col for col in df.columns if 'Sample Name' in col), None)
    
    if not sample_col:
        print("Warning: Could not find Sample Name column in assay table")
        # If no sample column, just use placeholder values
        if mode == "microbes":
            # Microbes mode (Bowtie2)
            values = [f"{glds_prefix}sample{i+1}.bowtie2.log" for i in range(len(df))]
        else:
            # Default mode (STAR)
            values = [f"{glds_prefix}sample{i+1}_Log.final.out" for i in range(len(df))]
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
            
            if mode == "microbes":
                # Microbes mode (Bowtie2)
                values.append(f"{glds_prefix}{sample}.bowtie2.log")
            else:
                # Default mode (STAR)
                values.append(f"{glds_prefix}{sample}_Log.final.out")
    
    # Add the column to the dataframe
    if column_name not in df.columns:
        print(f"Adding column: {column_name}")
        df[column_name] = values
    else:
        print(f"Updating column: {column_name}")
        df[column_name] = values
    
    return df

def add_ercc_analyses_column(df, glds_prefix, assay_suffix):
    """Add the ERCC Analyses column to the dataframe."""
    column_name = "Parameter Value[ERCC Analyses]"
    
    # Create the ercc analyses filename - same for all samples
    ercc_analyses_file = f"{glds_prefix}ERCC_analysis{assay_suffix}.html"
    
    # Add the column to the dataframe with the same value for all rows
    df = update_column(df, column_name, ercc_analyses_file)
    
    return df

def clean_comma_space(df):
    """Remove spaces after commas in all string columns of the dataframe."""
    # Loop through all columns in the dataframe
    for col in df.columns:
        # Only process string (object) columns
        if df[col].dtype == 'object':
            # Replace comma-space with just comma
            df[col] = df[col].str.replace(", ", ",", regex=False)
            
    print("Removed spaces after commas in all string columns")
    return df

def clean_column_names(df):
    """Clean column names by removing any .# suffixes pandas adds to duplicates.
    
    Args:
        df: The DataFrame to clean column names
        
    Returns:
        The DataFrame with cleaned column names
    """
    # Create a mapping of old_name -> new_name (without .# suffix)
    name_mapping = {}
    for col in df.columns:
        # Use regex to match column names with .digits suffix
        if re.search(r'\.\d+$', col):
            # Remove the .# suffix
            base_name = re.sub(r'\.\d+$', '', col)
            name_mapping[col] = base_name
    
    # Rename columns using the mapping if any found
    if name_mapping:
        print(f"Cleaning {len(name_mapping)} column names by removing .# suffixes:")
        for old_name, new_name in name_mapping.items():
            print(f"  - {old_name} -> {new_name}")
        df = df.rename(columns=name_mapping)
    
    return df

def main():
    args = parse_args()
    
    # Find and parse the runsheet
    runsheet_df = find_runsheet(args.outdir, args.glds_accession)
    
    # Find RNA-Seq assay file and get its contents
    assay_df, assay_filename = extract_and_find_assay(args.outdir, args.glds_accession)
    print(f"Original assay table has {len(assay_df)} rows and {len(assay_df.columns)} columns")
    
    # Process and save assay file
    try:
        # Create GLDS prefix for filenames
        glds_prefix = f"{args.glds_accession.upper()}_rna_seq_"
        
        # Check if ERCC spike-ins are used
        ercc_used = has_ercc_spikes(runsheet_df)
        
        # Following the original star workflow order:
        
        # 0. Add read counts from MultiQC report (new)
        assay_df = add_read_counts(assay_df, args.outdir, args.glds_accession, args.assay_suffix, runsheet_df)
        
        # 0.5. Add Raw MultiQC reports column (new)
        assay_df = add_raw_multiqc_reports_column(assay_df, glds_prefix, args.assay_suffix)
        
        # 1. Trimmed Sequence Data column
        assay_df = add_trimmed_data_column(assay_df, glds_prefix, runsheet_df=runsheet_df)
        
        # 2. Trimmed Sequence Data/MultiQC Reports column
        assay_df = add_trimmed_multiqc_reports_column(assay_df, glds_prefix, args.assay_suffix)
        
        # 3. Trimmed Sequence Data/Trimming Reports column - commented out as trimming reports are now in MultiQC
        # assay_df = add_trimming_reports_column(assay_df, glds_prefix, runsheet_df=runsheet_df)
        
        # 4. Aligned Sequence Data column
        assay_df = add_aligned_sequence_data_column(assay_df, glds_prefix, runsheet_df=runsheet_df, mode=args.mode)
        
        # 5. Unmapped Reads column (moved here as requested)
        assay_df = add_unmapped_reads_column(assay_df, glds_prefix, runsheet_df=runsheet_df, mode=args.mode)
        
        # 6. Aligned Sequence Data/Alignment Logs column
        assay_df = add_alignment_logs_column(assay_df, glds_prefix, runsheet_df=runsheet_df, mode=args.mode)
        
        # 7. Aligned Sequence Data/MultiQC Reports column
        assay_df = add_align_multiqc_reports_column(assay_df, glds_prefix, args.assay_suffix)
        
        # 8. RSeQC/MultiQC Reports column
        assay_df = add_rseqc_multiqc_reports_column(assay_df, glds_prefix, args.assay_suffix)
        
        # 9. Raw Counts Data column
        assay_df = add_raw_counts_data_column(assay_df, glds_prefix, args.assay_suffix, args.mode)
        
        # 10. Raw Counts Data/MultiQC Reports column
        assay_df = add_raw_counts_multiqc_column(assay_df, glds_prefix, args.assay_suffix, args.mode)
        
        # 11. Raw Counts Tables column
        assay_df = add_raw_counts_tables_column(assay_df, glds_prefix, args.assay_suffix, args.mode)
        
        # 12. Normalized Counts Data column
        assay_df = add_normalized_counts_data_column(assay_df, glds_prefix, args.assay_suffix)
        
        # 13. Differential Expression Analysis Data column
        assay_df = add_differential_expression_column(assay_df, glds_prefix, args.assay_suffix)
        
        # All rRNArm columns
        
        # 14. Raw Counts Tables rRNArm column
        assay_df = add_raw_counts_tables_rrnarm_column(assay_df, glds_prefix, args.assay_suffix, args.mode)
        
        # 15. Normalized Counts Data rRNArm column
        assay_df = add_normalized_counts_data_rrnarm_column(assay_df, glds_prefix, args.assay_suffix)
        
        # 16. Differential Expression Analysis Data rRNArm column
        assay_df = add_differential_expression_rrnarm_column(assay_df, glds_prefix, args.assay_suffix)
        
        # ERCC column (conditionally added only if ERCC spike-ins are used)
        if ercc_used:
            print("Adding ERCC Analyses column")
            # 17. ERCC Analyses column
            assay_df = add_ercc_analyses_column(assay_df, glds_prefix, args.assay_suffix)
        
        # Clean comma-space in all string columns
        assay_df = clean_comma_space(assay_df)
        
        # Clean column names by removing any .# suffixes pandas adds
        assay_df = clean_column_names(assay_df)
        
        # Use the filename we found in extract_and_find_assay
        orig_filename = assay_filename
        
        # Only save the original filename version
        assay_df.to_csv(orig_filename, sep='\t', index=False)
        print(f"Assay table saved as: {orig_filename}")
        
        # Print summary of changes
        print("\n=== SUMMARY OF CHANGES ===")
        print(f"Processed assay table for {args.glds_accession} with {len(column_changes)} column operations:")
        for change in column_changes:
            print(f"  - {change}")
        print(f"Final assay table has {len(assay_df)} rows and {len(assay_df.columns)} columns")
        
    except Exception as e:
        print(f"Error processing assay file: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()