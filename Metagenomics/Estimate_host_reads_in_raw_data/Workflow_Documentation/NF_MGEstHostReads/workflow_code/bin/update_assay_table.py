#!/usr/bin/env python

import sys
import argparse
import zipfile
import pandas as pd
import re
import os


def parse_args():
    parser = argparse.ArgumentParser(
        prog='update_assay_table',
        description='Update Metagenomics assay table from ISA.zip with an added "Parameter Value: Host Contamination" column.')
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--runsheet', required=True, 
                          help='Runsheet')
    isa_file_input = required.add_mutually_exclusive_group(required=True)
    isa_file_input.add_argument('--isa_zip', action='store', default='',
                                help='Appropriate ISA file for the dataset (a zip archive, providing this instead of '
                                     'an assay table directly will attempt to extract the correct assay table given '
                                     'the provided assay and technology types.)')
    isa_file_input.add_argument('--assay_table', action='store', default='',
                                help='Assay table for the dataset provided directly instead of extracting from an ISA '
                                     'zip file. If both are passed, assay_table will be used')
    required.add_argument('--summary_file', required=True, 
                          help='Summary file to extract host reads percentages')
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

    valid_measurement = 'Metagenomic sequencing'
    valid_technology = 'Whole-Genome Shotgun Sequencing'

    zip_file = zipfile.ZipFile(isa_file)
    isa_files = zip_file.namelist()

    # Parse investigation file to build STUDY ASSAYS table
    study_assays_table = {}
    study_assays_section = False

    with zip_file.open('i_Investigation.txt', 'r') as f:
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
    required_keys = ['Study Assay Measurement Type', 'Study Assay Technology Type', 'Study Assay File Name']
    if not all(key in study_assays_table for key in required_keys):
        report_failure_and_exit("Missing required keys in STUDY ASSAYS section")

    # Get the values from the table
    measurement_types = study_assays_table['Study Assay Measurement Type']
    technology_types = study_assays_table['Study Assay Technology Type']
    file_names = study_assays_table['Study Assay File Name']

    # Ensure all lists have equal length
    if not (len(measurement_types) == len(technology_types) == len(file_names)):
        report_failure_and_exit("Measurement types, technology types, and file names have different lengths")

    # Find matching assay file
    matched_file = ""
    for i in range(len(measurement_types)):
        if (measurement_types[i].lower() == valid_measurement.lower() and
                technology_types[i].lower() == valid_technology.lower()):
            matched_file = file_names[i]
            break

    if not matched_file:
        report_failure_and_exit(f"No assay file matched for {valid_measurement}. "
                                f"Measurement types: {measurement_types}, Technology types: {technology_types}")
    elif matched_file not in isa_files:
        # Load the matched assay file
        report_failure_and_exit(f"Matched assay file doesn't exist in ISA zip: {matched_file}")
    else:
        return pd.read_csv(zip_file.open(matched_file), sep='\t'), matched_file

    return pd.DataFrame(), matched_file


def load_summary_file(summary_file):
    """
    Load the host read removal summary file.

    Args:
        summary_file (PathLike[str]): path to summary TSV file

    Returns:
        pandas.DataFrame: summary data with sample IDs and host contamination percentages
    """
    try:
        summary_df = pd.read_csv(summary_file, sep='\t')
        print(f"Summary file has {len(summary_df)} rows and {len(summary_df.columns)} columns")
        return summary_df
    except Exception as e:
        report_failure_and_exit(f"Cannot read summary file: {e}")


def get_host_contamination_percentage(summary_df, sample_name):
    """
    Extract host contamination percentage for a given sample from the summary dataframe.

    Args:
        summary_df (pandas.DataFrame): summary data loaded from host removal summary file
        sample_name (str): sample name to look up

    Returns:
        float or str: contamination percentage, or 'N/A' if not found
    """
    # Find the Percent column (named Percent_<host>_reads)
    percent_col = next((col for col in summary_df.columns if col.startswith("Percent_")), None)
    if percent_col is None:
        report_warning(f"Could not find Percent column in summary file")
        return "N/A"

    match = summary_df[summary_df["Sample_ID"] == sample_name]
    if match.empty:
        report_warning(f"Sample '{sample_name}' not found in summary file")
        return "N/A"

    return match[percent_col].values[0]


def insert_or_update_unit_column(df, after_column, value="percent"):
    """
    Insert a Unit column immediately after a specified column, or update it if already present.

    Args:
        df (pandas.DataFrame): current assay table
        after_column (str): column name after which to insert the Unit column
        value (str): value to set for all rows (default: 'percent')

    Returns:
        pandas.DataFrame: updated assay table
    """
    col_idx  = df.columns.get_loc(after_column)
    next_col = df.columns[col_idx + 1] if col_idx + 1 < len(df.columns) else None

    if next_col and (next_col == "Unit" or re.match(r"^Unit\.\d+$", next_col)):
        print(f"Updating existing '{next_col}' column after '{after_column}'")
        df[next_col] = value
    else:
        df.insert(col_idx + 1, "Unit", value, allow_duplicates=True)
        print(f"Added 'Unit' column after '{after_column}'")

    return df


def add_host_contamination_column(df, sample_name_map, assay_sample_names, summary_df):
    """
    Add or update the host contamination column in the assay table.

    Args:
        df (pandas.DataFrame): current assay table
        sample_name_map (dict): mapping of assay table sample names to runsheet sample names
        assay_sample_names (list): list of sample names in the assay table
        summary_df (pandas.DataFrame): summary data from host removal summary file

    Returns:
        pandas.DataFrame: updated assay table
    """
    column_name = "Parameter Value[Host Contamination]"

    values = []
    for assay_sample in assay_sample_names:
        # Use mapped name if available, otherwise use the assay table name
        sample = sample_name_map.get(assay_sample, assay_sample)
        # Normalize name for summary file lookup (replace spaces with underscores)
        sample_normalized = sample.replace(" ", "_")
        percentage = get_host_contamination_percentage(summary_df, sample_normalized)
        values.append(str(percentage))

    if column_name not in df.columns:
        print(f"Adding column: {column_name}")
    else:
        print(f"Updating column: {column_name}")
    
    df[column_name] = values
    df = insert_or_update_unit_column(df, after_column=column_name)

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

    # Load summary file
    summary_df = load_summary_file(args.summary_file)

    # Find Metagenomic sequencing assay file and get its contents - use assay table if directly provided, else extract it from ISA.zip
    if args.assay_table:
        print(f"Reading user-provided assay table from file: {args.assay_table}")
        assay_filename = args.assay_table
        assay_df = pd.read_csv(open(assay_filename), sep='\t')
        print(f"Original assay table has {len(assay_df)} rows and {len(assay_df.columns)} columns")
    else:
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
        # Add host contamination
        assay_df = add_host_contamination_column(assay_df, sample_name_map=sample_name_map, assay_sample_names=assay_sample_names, summary_df=summary_df)
        
        # Clean comma-space in all string columns
        assay_df = clean_comma_space(assay_df)

        # Clean column names by removing any .# suffixes pandas adds
        assay_df = clean_column_names(assay_df)

        # Use the filename found in get_assay_table_from_isa() or passed from params.assay_table
        orig_filename = os.path.basename(assay_filename)

        # Save under original filename
        assay_df.to_csv(orig_filename, sep='\t', index=False)
        print(f"Original assay table saved as: {orig_filename}")

    except Exception as e:
        report_failure_and_exit(str(e))


if __name__ == "__main__":
    main()