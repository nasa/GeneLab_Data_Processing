#!/usr/bin/env python

import argparse
import subprocess
import os
import sys
import tempfile
import re 
import shutil
import pandas as pd
import requests 


####################
## 1.  For OSD ARG #
####################
# 1. Process the OSD arg to proper format
# 2. Download the ISA file
# 3. Convert to runsheet(s)
# 4. Select which runsheet to use

########################
## 1. For runsheet arg #
########################
# 1. Select which runsheet to use

##########################
## 2. Neutral flow after #
##########################
# 1. Validate schema of runsheet
# 2. Check if read_paths are URLs, prompt for download


# Process OSD arg: if numeric, append OSD-, if OSD-# or GLDS-#, leave it
def process_osd_argument(osd_arg):
    # Check if the argument is just numeric
    if osd_arg.isdigit():
        return f"OSD-{osd_arg}"
    # Check if it's already in the correct format (OSD-numeric or GLDS-numeric)
    elif re.match(r'^(OSD|GLDS)-\d+$', osd_arg):
        return osd_arg
    else:
        print("Invalid format for --OSD argument. Use 'numeric', 'OSD-numeric', or 'GLDS-numeric'.")
        sys.exit(1)

# Check provided OSD/GLDS is not on the list of those that can't be autoprocessed
def check_provided_osd_or_glds(osd_arg):
    # dictionaries of OSD/GLDS accessions and reason for not running, key = ID: value = reason
    # there are 3 because ID can be provided prefixed with "OSD-", "GLDS-", or nothing - not the most efficient here, but ¯\_(ツ)_/¯
    not_autoprocessable_OSD_dict = {
        "OSD-65": "This dataset has multiple different primers mixed in different orientations in each individual sample, and the workflow is unable to handle it in an automated fashion.",
        "OSD-66": "This dataset is not a standard amplicon dataset. It is comprised of hundreds of different primers targeting different regions of specific organisms, and the workflow is unable to handle it.",
        "OSD-82": "This dataset is still multiplexed, and we don't yet have the mapping information to split the samples apart appropriately."
    }

    not_autoprocessable_GLDS_dict = {
        "GLDS-65": "This dataset has multiple different primers mixed in different orientations in each individual sample, and the workflow is unable to handle it in an automated fashion.",
        "GLDS-66": "This dataset is not a standard amplicon dataset. It is comprised of hundreds of different primers targeting different regions of specific organisms, and the workflow is unable to handle it.",
        "GLDS-82": "This dataset is still multiplexed, and we don't yet have the mapping information to split the samples apart appropriately."
    }

    not_autoprocessable_dict = {
        "65": "This dataset has multiple different primers mixed in different orientations in each individual sample, and the workflow is unable to handle it in an automated fashion.",
        "66": "This dataset is not a standard amplicon dataset. It is comprised of hundreds of different primers targeting different regions of specific organisms, and the workflow is unable to handle it.",
        "82": "This dataset is still multiplexed, and we don't yet have the mapping information to split the samples apart appropriately."
    }

    # Checking based on OSD IDs
    if osd_arg in not_autoprocessable_OSD_dict:
        print(f"\nThe specified dataset {osd_arg} is unable to be processed with this workflow.")
        print(f"    Reason: {not_autoprocessable_OSD_dict[osd_arg]}\n")
        sys.exit(1)

    # checking based on GLDS IDs
    if osd_arg in not_autoprocessable_GLDS_dict:
        print(f"\n The specified dataset {osd_arg} is unable to be processed with this workflow.")
        print(f"    Reason: {not_autoprocessable_GLDS_dict[osd_arg]}\n")
        sys.exit(1)

    # checking based on plain IDs
    if osd_arg in not_autoprocessable_dict:
        print(f"\n The specified dataset {osd_arg} is unable to be processed with this workflow.")
        print(f"    Reason: {not_autoprocessable_dict[osd_arg]}\n")
        sys.exit(1)

# Run dpt-get-isa-archive in a temp folder, move it back to cd, return the filename
def download_isa_archive(accession_number):
    with tempfile.TemporaryDirectory() as temp_dir:
        try:
            # Run the command in the temporary directory
            subprocess.run(
                ["dpt-get-isa-archive", "--accession", str(accession_number)],
                check=True,
                text=True,
                cwd=temp_dir
            )

            # Find the downloaded zip file in the temp directory
            downloaded_files = [f for f in os.listdir(temp_dir) if f.endswith('.zip')]
            if not downloaded_files:
                print("No ISA archive file was downloaded.", file=sys.stderr)
                return None

            # Assuming there's only one file, get its name
            downloaded_file = downloaded_files[0]

            # Move the file back to the current directory
            shutil.move(os.path.join(temp_dir, downloaded_file), downloaded_file)

            full_path = os.path.abspath(downloaded_file)
            return full_path

        except subprocess.CalledProcessError as e:
            print("An error occurred while downloading ISA archive.", file=sys.stderr)
            sys.exit(1)

# Run dpt-isa-to-runsheet in a temp folder, move runsheet(s) back to cd, return list of runsheet(s)
def convert_isa_to_runsheet(accession_number, isa_zip):
    with tempfile.TemporaryDirectory() as temp_dir:
        # Copy the ISA archive to the temporary directory
        temp_isa_zip_path = shutil.copy(isa_zip, temp_dir)

        try:
            # Run the dpt-isa-to-runsheet command in the temporary directory
            subprocess.run(
                ["dpt-isa-to-runsheet", "--accession", accession_number, "--config-type", "amplicon", "--config-version", "Latest", "--isa-archive", os.path.basename(temp_isa_zip_path)],
                check=True,
                cwd=temp_dir,
                stdout=sys.stdout,
                stderr=sys.stderr
            )

            # Get the list of created files in the temp directory
            created_files = [f for f in os.listdir(temp_dir) if os.path.isfile(os.path.join(temp_dir, f)) and f != os.path.basename(temp_isa_zip_path)]

            # Move the created files back to the current directory
            moved_files = []
            for file in created_files:
                shutil.move(os.path.join(temp_dir, file), file)
                moved_files.append(file)

            return moved_files

        except subprocess.CalledProcessError as e:
            print("An error occurred while converting ISA archive to runsheet.", file=sys.stderr)
            sys.exit(1)


def handle_runsheet_selection(runsheet_files, target=None, specified_runsheet=None):
    selected_runsheet = None

    # Use the specified runsheet if provided
    if specified_runsheet and specified_runsheet in runsheet_files:
        selected_runsheet = specified_runsheet
        print(f"Using specified runsheet: {selected_runsheet}")
        return selected_runsheet

    if len(runsheet_files) == 1:
        if target:
            runsheet = runsheet_files[0]
            try:
                runsheet_df = pd.read_csv(runsheet)
                target_region = runsheet_df['Parameter Value[Library Selection]'].unique()[0]
                if target.lower() == target_region.lower():
                    selected_runsheet = runsheet
            except Exception as e:
                print(f"Error reading {runsheet}: {e}")
            print(f"Using runsheet: {selected_runsheet}")

    elif len(runsheet_files) > 1:
        if target:
            matching_runsheets = []
            for runsheet in runsheet_files:
                try:
                    runsheet_df = pd.read_csv(runsheet)
                    target_region = runsheet_df['Parameter Value[Library Selection]'].unique()[0]
                    if target.lower() == target_region.lower():
                        matching_runsheets.append(runsheet)
                except Exception as e:
                    print(f"Error reading {runsheet}: {e}")

            if len(matching_runsheets) == 1:
                # One matching runsheet found
                selected_runsheet = matching_runsheets[0]
                print(f"Using runsheet: {selected_runsheet}")

            elif len(matching_runsheets) > 1:
                # Multiple matching runsheets found
                print("The study contains multiple assays with the same target. Please specify one of the following runsheet names as a parameter for the --specify-runsheet argument:")
                for rs in matching_runsheets:
                    print(rs)
                return None

            else:
                # No matching runsheets found
                print("No runsheet matches the specified genomic target. Please check the target or specify a runsheet using --specify-runsheet.")
                return None

        else:
            # No target specified and multiple runsheets are available
            print("Multiple runsheets found but no genomic target specified. Cannot proceed. Use -t {16S, 18S, ITS} or --target {16S, 18S, ITS} to specify which assay/dataset to use.")
            return None

    # Remove unselected runsheet files if a runsheet was selected
    if selected_runsheet:
        unselected_runsheets = [file for file in runsheet_files if file != selected_runsheet]
        for file in unselected_runsheets:
            try:
                os.remove(file)
            except Exception as e:
                pass

    return selected_runsheet

def check_runsheet_read_paths(runsheet_df):
    # Check if a string is a URL / genelab URL
    def is_url(s):
        return "http://" in s or "https://" in s or "genelab-data.ndc.nasa.gov" in s


    # Check if 'read2_path' column exists 
    paired_end = runsheet_df['paired_end'].eq(True).all()

    # Check the first row to determine if the paths are URLs or local paths
    first_row = runsheet_df.iloc[0]

    uses_url = is_url(first_row['read1_path'])
    if uses_url:
        print("Runsheet references URLs.")
    else:
        print("Runsheet references local read files.")

    return uses_url

def sample_IDs_from_local(runsheet_df, output_file='unique-sample-IDs.txt'):
    # Check if the DataFrame is paired-end
    paired_end = runsheet_df['paired_end'].eq(True).all()

    with open(output_file, 'w') as file:
        for index, row in runsheet_df.iterrows():
            # Extract base names minus the suffixes
            base_read1 = os.path.basename(row['read1_path']).replace(row['raw_R1_suffix'], '')

            if paired_end:
                base_read2 = os.path.basename(row['read2_path']).replace(row['raw_R2_suffix'], '')
                # Check if base names match for paired-end data, necessary for snakemake arg expansion
                if base_read1 != base_read2:
                    print(f"Mismatch in sample IDs in row {index}: {base_read1} vs {base_read2}")
                    sys.exit(1)
            
            # Write the base name to the file
            file.write(f"{base_read1}\n")
    
    print(f"Unique sample IDs written to {output_file}")

def handle_url_downloads(runsheet_df, output_file='unique-sample-IDs.txt'):
    print("Downloading read files...")
    # Check if the DataFrame is paired-end
    paired_end = runsheet_df['paired_end'].eq(True).all()
    # Write 'Sample Name' into unique-sample-IDs.txt
    with open(output_file, 'w') as file:
        for sample_name in runsheet_df['Sample Name']:
            file.write(sample_name + '\n')

    # Create ./raw_reads/ directory if it does not exist
    raw_reads_dir = os.path.abspath('./raw_reads/')
    if not os.path.exists(raw_reads_dir):
        os.makedirs(raw_reads_dir)

    # Initialize count for skipped downloads
    skipped_downloads_count = 0
    # Iterate over each row and download files if they don't exist
    for _, row in runsheet_df.iterrows():
        sample_id = row['Sample Name']
        read1_path = os.path.join(raw_reads_dir, sample_id + row['raw_R1_suffix'])
        read2_path = os.path.join(raw_reads_dir, sample_id + row['raw_R2_suffix']) if paired_end else None

        # Download Read 1 if it doesn't exist
        if not os.path.exists(read1_path):
            download_url_to_file(row['read1_path'], read1_path)
        else:
            skipped_downloads_count += 1

        # Download Read 2 if it doesn't exist and if paired_end
        if paired_end and read2_path and not os.path.exists(read2_path):
            download_url_to_file(row['read2_path'], read2_path)
        elif paired_end and read2_path:
            skipped_downloads_count += 1

    # Print the number of skipped downloads
    if skipped_downloads_count > 0:
        print(f"{skipped_downloads_count} read file(s) were already present and were not downloaded.")

def download_url_to_file(url, file_path, max_retries=3, timeout_seconds=120):
    retries = 0
    success = False

    while retries < max_retries and not success:
        try:
            response = requests.get(url, stream=True, timeout=timeout_seconds)
            response.raise_for_status()  # Raises an HTTPError for bad status codes

            with open(file_path, 'wb') as file:
                shutil.copyfileobj(response.raw, file)
            success = True

        except (requests.exceptions.HTTPError, requests.exceptions.ConnectionError, requests.exceptions.Timeout) as e:
            retries += 1
            print(f"Attempt {retries}: Error occurred: {e}")

        except requests.exceptions.RequestException as e:
            print(f"An unexpected error occurred: {e}")
            break

    if not success:
        print("Failed to download the read files.")


def write_params(runsheet_df, uses_urls):
    
    # Extract necessary variables from runsheet_df
    data_type = "PE" if runsheet_df['paired_end'].eq(True).all() else "SE"
    raw_R1_suffix = runsheet_df['raw_R1_suffix'].unique()[0]
    raw_R2_suffix = runsheet_df['raw_R2_suffix'].unique()[0] if data_type == "PE" else ""
    f_primer = runsheet_df['F_Primer'].unique()[0]
    r_primer = runsheet_df['R_Primer'].unique()[0] if data_type == "PE" else ""
    target_region = runsheet_df['Parameter Value[Library Selection]'].unique()[0]

    # Determine raw_reads_directory
    if uses_urls:
        raw_reads_directory = os.path.abspath('./raw_reads/') + '/'
    else:
        read1_path_dir = os.path.dirname(runsheet_df['read1_path'].iloc[0])
        raw_reads_directory = os.path.abspath(read1_path_dir) + '/' if read1_path_dir else "./"

    with open("GLparams_file.csv", "w") as f:
        f.write("raw_reads_directory,raw_R1_suffix,raw_R2_suffix,f_primer,r_primer,target_region,data_type\n")
        if data_type == "PE":
            f.write(f"{raw_reads_directory},{raw_R1_suffix},{raw_R2_suffix},{f_primer},{r_primer},{target_region},{data_type}\n")
        else:
            f.write(f"{raw_reads_directory},{raw_R1_suffix},{f_primer},{r_primer},{target_region},{data_type}\n")

 

def write_input_file(runsheet_df):
    """ Write input file for the workflow..."""

    print("writing out GLfile.csv...")
    # Check if the DataFrame is paired-end
    paired_end = runsheet_df['paired_end'].eq(True).all()

    # Create ./raw_reads/ directory if it does not exist
    raw_reads_dir = os.path.abspath('./raw_reads/')
    if not os.path.exists(raw_reads_dir):
        os.makedirs(raw_reads_dir)

    # Create input file
    with open("GLfile.csv", 'w') as file:
        
        if paired_end:
            file.write(f"sample_id,forward,reverse,paired\n")
            # Iterate over each row and download files if they don't exist
            for _, row in runsheet_df.iterrows():
                sample_id = row['Sample Name']
                read1_path = os.path.join(raw_reads_dir, sample_id + row['raw_R1_suffix'])
                read2_path = os.path.join(raw_reads_dir, sample_id + row['raw_R2_suffix'])
                file.write(f"{sample_id},{read1_path},{read2_path},true\n")
        else:
            file.write(f"sample_id,forward,paired\n")
            for _, row in runsheet_df.iterrows():
                sample_id = row['Sample Name']
                read1_path = os.path.join(raw_reads_dir, sample_id + row['raw_R1_suffix'])
                file.write(f"{sample_id},{read1_path},false\n")


# Check for single primer set, also check for invalid characters in primers used, exit if either
def validate_primer_sequences(runsheet_df):
    errors = []

    # Check that there is only 1 entry in each primer column
    if len(runsheet_df['F_Primer'].unique()) > 1:
        errors.append(f"Multiple primer sequences present in F_Primer: {runsheet_df['F_Primer'].unique()}.")

    if len(runsheet_df['R_Primer'].unique()) > 1:
         errors.append(f"Multiple primer sequences present in R_primer: {runsheet_df['R_Primer'].unique()}.")
    

    # Check for non-letter characters in primer sequences
    def has_non_letter_characters(primer):
        # Pattern to find any character that is not a letter
        non_letter_pattern = re.compile(r'[^A-Za-z]')
        return non_letter_pattern.search(primer)

    # Check each unique primer in the F_Primer and R_Primer columns
    for f_primer in runsheet_df['F_Primer'].unique():
        if has_non_letter_characters(f_primer):
            errors.append(f"Non-letter characters detected in F_Primer: '{f_primer}'")

    for r_primer in runsheet_df['R_Primer'].unique():
        if has_non_letter_characters(r_primer):
            errors.append(f"Non-letter characters detected in R_Primer: '{r_primer}'")

    if errors:
        print("Error: Invalid primer sequence(s) detected in the runsheet.")
        for error in errors:
            print(f"  - {error}")
        print("Correct the primer sequences in the runsheet and rerun the workflow from the runsheet using the --runsheetPath argument.")
        sys.exit(1)


def main():
    # Argument parser setup with short argument names and an automatic help option
    parser = argparse.ArgumentParser(
        description='Create Runsheet from Genelab ID.',
        add_help=True,
        usage='%(prog)s [options]'  # Custom usage message
    )
    
    parser.add_argument('-o', '--OSD',
                        metavar='osd_number',
                        help='A GeneLab OSD dataset accession number to pull its read files and associated metadata. Acceptable formats: ###, OSD-###, GLDS-###',
                        type=str)
    
    parser.add_argument('-t', '--target',
                        choices=['16S', '18S', 'ITS'],
                        help='Specify the amplicon target for the assay. Options: 16S, 18S, ITS. This is used to select the appropriate dataset from an OSD study when multiple options are available.',
                        type=str)
    
    parser.add_argument('-r', '--runsheetPath',
                        metavar='/path/to/runsheet.csv',
                        help='Set up the Snakemake workflow using a specified runsheet file.',
                        type=str)

    
    parser.add_argument('--specify-runsheet',
                        help='Specifies the runsheet for an OSD dataset by name. Only used if there are multiple datasets with the same target in the study.',
                        metavar='runsheet_name',
                        type=str)
    
      
    # Check if no arguments were provided
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    try:
        args = parser.parse_args()
    except SystemExit:
        parser.print_help()
        sys.exit(1)

    target = args.target
    isa_zip = ""

    # If OSD is used, pull ISA metadata for the study, create and select the runsheet
    if args.OSD:
        accession_number = process_osd_argument(args.OSD)

        # checking OSD/GLDS ID is not on the list of those the workflow definitely can't handle
        check_provided_osd_or_glds(args.OSD)

        isa_zip = download_isa_archive(accession_number)
        if isa_zip:
            runsheet_files = convert_isa_to_runsheet(accession_number, isa_zip)
            if runsheet_files:
                runsheet_file = handle_runsheet_selection(runsheet_files, target, args.specify_runsheet)
                if runsheet_file is None:
                    sys.exit()
            else:
                print("No runsheet files were created.")
        else:
            print("No ISA archive was downloaded. Cannot proceed to runsheet conversion.", file=sys.stderr)
            sys.exit(1)
    
    # If a runsheet is specified, use that runsheet
    elif args.runsheetPath:
        runsheet_file = args.runsheetPath

    # Load the runsheet if a file is specified
    # Create unique-sample-IDs.txt based on filenames or 'Sample Name' if URLs
    # Download files if necessary
    if args.OSD or args.runsheetPath:
        if runsheet_file:
            #runsheet_df = validate_runsheet_schema(runsheet_file)
            runsheet_df = pd.read_csv(runsheet_file)
            if runsheet_df is not None:
                uses_urls = check_runsheet_read_paths(runsheet_df)

                # Check for primer file / invalid primers
                validate_primer_sequences(runsheet_df)

                # Create the 'unique-sample-IDs.txt' file and download read files if necessary
                if uses_urls:
                    handle_url_downloads(runsheet_df, output_file='unique-sample-IDs.txt')
                else:
                    sample_IDs_from_local(runsheet_df, output_file='unique-sample-IDs.txt')

                # Create the config.yaml file
                write_params(runsheet_df=runsheet_df, uses_urls=uses_urls)
                # Create input file required by the workflow
                write_input_file(runsheet_df=runsheet_df)
            else:
                print("Failed to validate the runsheet file.", file=sys.stderr)
                sys.exit(1)
        else:
            print("No runsheet file specified.", file=sys.stderr)
            sys.exit(1)
    



if __name__ == "__main__":
    main()
