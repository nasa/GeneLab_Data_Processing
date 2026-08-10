#!/usr/bin/env python

import sys
import argparse
import zipfile
import pandas as pd
import json
import re


def parse_args():
    parser = argparse.ArgumentParser(
        prog='update_assay_table',
        description='Update Microarray Affymetrix assay table from ISA.zip with processed data file information.')
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--runsheet', required=True, 
                          help='Runsheet')
    required.add_argument('--glds_accession', required=True, 
                          help='GLDS accession number (e.g. GLDS-123)')
    required.add_argument('--isa_zip', action='store', default='',
                          help='Appropriate ISA file for the dataset')
    return parser.parse_args()


tty_colors = {
    'green': '\033[0;32m%s\033[0m',
    'yellow': '\033[0;33m%s\033[0m',
    'red': '\033[0;31m%s\033[0m'
}


def color_text(text, color='green'):
    """
    Colors text for output in terminal

    Args:
        text (str): input text
        color (str): a valid tty color ('red', 'yellow', or 'green')

    Returns:
        str: colored text

    """
    if sys.stdout.isatty():
        return tty_colors[color] % text
    else:
        return text


def report_failure_and_exit(message, color="red"):
    """
    Reports a failure and exits with status '1'.

    Args:
        message (str): Error message to report.
        color (str): Color in which to render the error message, default = 'red'
    """
    print("")
    print(color_text(f"Error: {message}", color))
    print("\nAssay table update failed.\n")

    sys.exit(1)


def report_warning(message, color="yellow"):
    """
    Reports are warning message.

    Args:
        message (str): Error message to report
        color (str): Color in which to render the error message, default = 'yellow'
    """
    print("")
    print(color_text(f"Warning: {message}", color))


def load_runsheet(runsheet_file):
    """
    Load the runsheet as a pandas.DataFrame.

    Args:
        runsheet_file (PathLike[str]): a file containing the assay runsheet used to generate the processed data

    Returns:
        pandas.DataFrame: sample information from the runsheet
    """
    try:
        runsheet_df = pd.read_csv(runsheet_file)
        print(f"Runsheet has {len(runsheet_df)} rows and {len(runsheet_df.columns)} columns")
        return runsheet_df
    except Exception as e:
        report_warning(f"Cannot read runsheet, proceeding without it: {e}")
        return None


def get_runsheet_sample_name_map(runsheet_df, assay_sample_names):
    """
    Generates a mapping of sample names in the assay table to the samplenames in the runsheet

    Args:
        runsheet_df (pandas.DataFrame): runsheet sample information
        assay_sample_names (list): sample names from the assay table

    Returns:
        dict: sample name mapping
    """
    sample_name_map = {}
    if 'Sample Name' in runsheet_df.columns:
        # Check for 'Original Sample Name' column to map between assay table and runsheet
        if 'Original Sample Name' in runsheet_df.columns:
            for _, row in runsheet_df.iterrows():
                orig_name = row['Original Sample Name']
                rs_name = row['Sample Name']
                if orig_name in assay_sample_names:
                    sample_name_map[orig_name] = rs_name
    return sample_name_map



def get_assay_table_from_isa(isa_file):
    """
    tries to find an assay table in an ISA zip file that matches the type expected for the provided assay

    Args:
        isa_file (PathLike[str]): path to ISA zip file

    Returns:
        pandas.DataFrame: assay table from extracted from ISA zip
    """
    valid_measurement = "transcription profiling"
    valid_technology = "DNA microarray"
    valid_platform = "Affymetrix"

    zip_file = zipfile.ZipFile(isa_file)
    isa_files = zip_file.namelist()

    # Parse investigation file to build STUDY ASSAYS table
    study_assays_table = {}
    study_assays_section = False

    investigation_file = next((f for f in isa_files if f.startswith('i_') and f.endswith('.txt')), None)
    if not investigation_file:
        report_failure_and_exit(f"Investigation file not found in ISA zip: {isa_file}")

    with zip_file.open(investigation_file, 'r') as f:
        for line in f:
            line = line.decode('utf-8').strip()

            # Track STUDY ASSAYS section
            if line == 'STUDY ASSAYS':
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
    required_keys = ['Study Assay Measurement Type', 'Study Assay Technology Type', 'Study Assay Technology Platform', 'Study Assay File Name']
    if not all(key in study_assays_table for key in required_keys):
        report_failure_and_exit("Missing required keys in STUDY ASSAYS section")

    # Get the values from the table
    measurement_types = study_assays_table['Study Assay Measurement Type']
    technology_types = study_assays_table['Study Assay Technology Type']
    technology_platforms = study_assays_table['Study Assay Technology Platform']
    file_names = study_assays_table['Study Assay File Name']

    # Ensure all lists have equal length
    if not (len(measurement_types) == len(technology_types) == len(technology_platforms) == len(file_names)):
        report_failure_and_exit("Measurement types, technology types, technology platforms, and file names have different lengths")

    # Find matching assay file
    matched_file = ""
    for i in range(len(measurement_types)):
        if (measurement_types[i].lower() == valid_measurement.lower() and
                technology_types[i].lower() == valid_technology.lower() and
                technology_platforms[i].lower() == valid_platform.lower()):
            matched_file = file_names[i]
            break

    if not matched_file:
        report_failure_and_exit(f"No assay file matched for {valid_measurement} assay. "
                                f"Measurement types: {measurement_types}, Technology types: {technology_types}, Technology platforms: {technology_platforms}")
    elif matched_file not in isa_files:
        # Load the matched assay file
        report_failure_and_exit(f"Matched assay file doesn't exist in ISA zip: {matched_file}")
    else:
        return pd.read_csv(zip_file.open(matched_file), sep='\t'), matched_file

    return pd.DataFrame(), matched_file


def add_parameter_column(df, column_name, value, glds_prefix=None):
    """
    Add a parameter column to the dataframe if it doesn't exist already. Use the same value for all rows in the table.
    
    Args:
        df (pandas.DataFrame): current assay table
        glds_prefix (str): a prefix to add to the start of each value (only for values that are filenames)
        column_name (str): parameter column name to add (e.g., "Parameter Value[Entry]")
        value (str): value to set for all rows in the table

    Returns:
        pandas.DataFrame: update assay table
    """
    # Apply prefix to value if provided
    if glds_prefix and isinstance(value, str):
        # Check if value already has the prefix
        if not value.startswith(glds_prefix):
            prefixed_value = f"{glds_prefix}{value}"
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


def add_protocol_ref_column(df):
    """Add Protocol REF GeneLab microarray data processing protocol column"""
    value = "GeneLab microarray data processing protocol"
    
    # Insert after "Array Data File" or "Derived Array Data File" if present, otherwise at the end
    insert_position = None

    # Pass 1: prefer the first "Derived Array Data File" column
    for i, col in enumerate(df.columns):
        if "Derived Array Data File" in col:
            column = "Derived Array Data File"
            insert_position = i + 1
            break

    # Pass 2: fall back to the first "Array Data File" column
    if insert_position is None:
        for i, col in enumerate(df.columns):
            if "Array Data File" in col:
                column = "Array Data File"
                insert_position = i + 1
                break

    # Pass 3: fall back to end of dataframe
    if insert_position is None:
        column = "end of dataframe"
        insert_position = len(df.columns)
    
    # Remove any existing Protocol REF columns with the exact data processing protocol value
    drop_cols = []
    drop_indices = []
    for i, col in enumerate(df.columns):
        if "protocol ref" in col.lower():
            col_data = df.iloc[:, i].astype(str).str.strip()
            if (col_data.str.lower() == value.lower()).any():
                drop_cols.append(col)
                drop_indices.append(i)
    
    if drop_cols:
        print(f"Removing existing Protocol REF columns with data processing protocol: {', '.join(drop_cols)}")
        df = df.drop(columns=drop_cols)
        # Adjust insert position if needed
        removed_before = sum(1 for idx in drop_indices if idx < insert_position)
        insert_position -= removed_before
        if insert_position < 0:
            insert_position = 0
    
    # Insert with temp name, rename "Protocol REF"
    temp_name = "Protocol REF_DP"
    while temp_name in df.columns:
        temp_name = f"{temp_name}_temp"
    df.insert(insert_position, temp_name, value)
    cols = list(df.columns)
    cols[insert_position] = "Protocol REF"
    df.columns = cols
    
    print(f"Added: Protocol REF (value: {value}) after {column}")
    return df


def add_raw_intensities_table_column(df, glds_prefix):
    """
    Add or update the raw intensities table column to the dataframe.

    Args:
        df (pandas.DataFrame): current assay table
        glds_prefix (str): a prefix to add to the start of each filename

    Returns:
        pandas.DataFrame: updated assay table
    """
    column_name = "Parameter Value[Raw Intensities Table]"
    # Create the raw intensities filename - same for all samples
    raw_intensities = (f"{glds_prefix}_array_raw_intensities_probes_GLmicroarray.csv")

    # Look for an existing column matching this name, ignoring case,
    # so we update/rename it instead of creating a duplicate column
    existing_col = next((col for col in df.columns if col.lower() == column_name.lower()), None)

    # If a differently-cased version exists, rename it to the canonical name
    if existing_col and existing_col != column_name:
        print(f"Renaming column '{existing_col}' to '{column_name}'")
        df = df.rename(columns={existing_col: column_name})

    # Set the value for all rows (creates the column if it didn't exist)
    print(f"{'Updating' if existing_col else 'Adding'} column: {column_name}")
    df[column_name] = raw_intensities

    return df


def add_normalized_expression_table_column(df, glds_prefix):
    """
    Add or update the normalized expression table column to the dataframe.

    Args:
        df (pandas.DataFrame): current assay table
        glds_prefix (str): a prefix to add to the start of each filename

    Returns:
        pandas.DataFrame: updated assay table
    """
    column_name = "Parameter Value[Normalized Expression Table]"
    # Create the normalized expression filenames - same for all samples
    normalized_files = [
        f"{glds_prefix}_array_normalized_expression_probset_GLmicroarray.csv",
        f"{glds_prefix}_array_normalized_intensities_probe_GLmicroarray.csv"
    ]

    combined_files = ','.join(normalized_files)

    # Look for an existing column matching this name, ignoring case,
    # so we update/rename it instead of creating a duplicate column
    existing_col = next((col for col in df.columns if col.lower() == column_name.lower()), None)

    # If a differently-cased version exists, rename it to the canonical name
    if existing_col and existing_col != column_name:
        print(f"Renaming column '{existing_col}' to '{column_name}'")
        df = df.rename(columns={existing_col: column_name})

    # Set the value for all rows (creates the column if it didn't exist)
    print(f"{'Updating' if existing_col else 'Adding'} column: {column_name}")
    df[column_name] = combined_files

    return df


def add_differential_expression_analysis_data_column(df, glds_prefix):
    """
    Add or update the Differential Expression Analysis Data column to the dataframe.

    Args:
        df (pandas.DataFrame): current assay table
        glds_prefix (str): a prefix to add to the start of each filename

    Returns:
        pandas.DataFrame: updated assay table
    """
    column_name = "Parameter Value[Differential Expression Analysis Data]"
    # Create the differential expression filenames - same for all samples
    de_files = [
        f"{glds_prefix}_array_SampleTable_GLmicroarray.csv",
        f"{glds_prefix}_array_contrasts_GLmicroarray.csv",
        f"{glds_prefix}_array_differential_expression_GLmicroarray.csv"
    ]

    # Join the files with commas
    combined_files = ','.join(de_files)

    # Look for an existing column matching this name, ignoring case,
    # so we update/rename it instead of creating a duplicate column
    existing_col = next((col for col in df.columns if col.lower() == column_name.lower()), None)

    # If a differently-cased version exists, rename it to the canonical name
    if existing_col and existing_col != column_name:
        print(f"Renaming column '{existing_col}' to '{column_name}'")
        df = df.rename(columns={existing_col: column_name})

    # Set the value for all rows (creates the column if it didn't exist)
    print(f"{'Updating' if existing_col else 'Adding'} column: {column_name}")
    df[column_name] = combined_files

    return df


def clean_comma_space(df):
    """Remove spaces after commas in all string columns of the dataframe."""
    for col in df.columns:
        try:
            # If there are duplicate columns, df[col] is a DataFrame, not Series
            if hasattr(df[col], 'dtype') and df[col].dtype == 'object':
                df[col] = df[col].str.replace(", ", ",", regex=False)
        except Exception:
            # If df[col] is a DataFrame (duplicate columns), apply to each
            for c in range(df.columns.get_loc(col), len(df.columns)):
                if df.columns[c] == col:
                    s = df.iloc[:, c]
                    if s.dtype == 'object':
                        df.iloc[:, c] = s.str.replace(", ", ",", regex=False)
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
    runsheet_df = load_runsheet(args.runsheet)

    # Find Microarray sequencing assay file and get its contents - use assay table if directly provided, else extract it from ISA.zip
    print(f"Extracting assay table from {args.isa_zip}")
    assay_df, assay_filename = get_assay_table_from_isa(args.isa_zip)
    print(f"Original assay table has {len(assay_df)} rows and {len(assay_df.columns)} columns")


    # Create a mapping from assay table sample names to runsheet sample names, exit if no sample column found
    sample_col = next((col for col in assay_df.columns if 'Sample Name' in col), None)
    if sample_col is None:
        report_failure_and_exit(f"Could not find 'Sample Name' column in assay table '{assay_filename}'")

    assay_sample_names = assay_df[sample_col].tolist()
    if runsheet_df is not None and 'Sample Name' in runsheet_df.columns:
        sample_name_map = get_runsheet_sample_name_map(runsheet_df, assay_sample_names)
    else:
        sample_name_map = {}

    # Process and save assay file
    try:

        
        assay_df = add_protocol_ref_column(assay_df)

        # Raw Intensities Table column
        assay_df = add_raw_intensities_table_column(assay_df, args.glds_accession)

        # Normalized Expression Table column
        assay_df = add_normalized_expression_table_column(assay_df, args.glds_accession)        

        # Differential Expression Analysis Data column
        assay_df = add_differential_expression_analysis_data_column(assay_df, args.glds_accession)

        # Clean comma-space in all string columns
        assay_df = clean_comma_space(assay_df)

        # Clean column names by removing any .# suffixes pandas adds
        assay_df = clean_column_names(assay_df)

        # Use the filename we found in extract_and_find_assay
        orig_filename = assay_filename

        # Create both original and modified output files
        # Original file (preserving the original name)
        assay_df.to_csv(orig_filename, sep='\t', index=False)
        print(f"Original assay table saved as: {orig_filename}")

        # Modified file with GLDS prefix
        if not orig_filename.startswith(args.glds_accession):
            mod_filename = f"{args.glds_accession}{orig_filename}"
        else:
            mod_filename = orig_filename

        assay_df.to_csv(mod_filename, sep='\t', index=False)
        print(f"Modified assay table saved as: {mod_filename}")

    except Exception as e:
        report_failure_and_exit(str(e))


if __name__ == "__main__":
    main()