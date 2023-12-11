import argparse
import subprocess
import os
import sys
import tempfile
import re 
import shutil
import pandas as pd
#import pandera as pa
import requests 
import yaml
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
# 3. Create config.yaml and unique-sample-IDs.txt
# 4. If --run is used: run the workflow

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

    # Change specified_runsheet to a basename in case a path is used as an arg for run_workflow.py
    if specified_runsheet:
        specified_runsheet_basename = os.path.basename(specified_runsheet)
    else:
        specified_runsheet_basename = None

    # Use the specified runsheet if provided
    if specified_runsheet and specified_runsheet in runsheet_files:
        selected_runsheet = specified_runsheet
        print(f"Using specified runsheet: {selected_runsheet}")
        return selected_runsheet

    if len(runsheet_files) == 1:
        # Automatically use the single runsheet file
        selected_runsheet = runsheet_files[0]
        print(f"Using runsheet: {selected_runsheet}")

    elif len(runsheet_files) > 1:
        if target:
            matching_runsheets = []
            for runsheet in runsheet_files:
                try:
                    runsheet_df = pd.read_csv(runsheet)
                    target_region = runsheet_df['Parameter Value[Library Selection]'].unique()[0]
                    if target.lower() in target_region.lower():
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


def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 
                  'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 
                  'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B', 
                  'D': 'H', 'H': 'D', 'N': 'N'}
    return ''.join(complement.get(base, base) for base in reversed(seq))

def create_config_yaml(isa_zip,
                       runsheet_file,
                       runsheet_df,
                       uses_urls,
                       output_dir,
                       min_trimmed_length,
                       trim_primers,
                       primers_linked,
                       anchor_primers,
                       discard_untrimmed,
                       left_trunc,
                       right_trunc,
                       left_maxEE,
                       right_maxEE,
                       concatenate_reads_only,
                       output_prefix):
    
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


    # Other default values
    output_dir = os.path.abspath(output_dir) + '/'
    primer_anchor = "^" if anchor_primers is True else ""

    f_linked_primer = f"{primer_anchor}{f_primer}...{reverse_complement(r_primer)}"
    r_linked_primer = f"{primer_anchor}{r_primer}...{reverse_complement(f_primer)}"

    # Make output_dir if it doesn't exist 
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    info_out_dir = os.path.join(output_dir, output_prefix + "Processing_Info") + os.sep
    fastqc_out_dir = os.path.join(output_dir, "FastQC_Outputs") + os.sep
    trimmed_reads_dir = os.path.join(output_dir, "Trimmed_Sequence_Data") + os.sep
    filtered_reads_dir = os.path.join(output_dir, "Filtered_Sequence_Data") + os.sep
    final_outputs_dir = os.path.join(output_dir, "Final_Outputs") + os.sep
    plots_dir = final_outputs_dir + "Plots" + os.sep

    # Write to config.yaml
    with open('config.yaml', 'w') as file:
        file.write("############################################################################################\n")
        file.write("## Configuration file for GeneLab Illumina amplicon processing workflow                   ##\n")
        file.write("## Developed by Michael D. Lee (Mike.Lee@nasa.gov)                                        ##\n")
        file.write("############################################################################################\n\n")

        file.write("############################################################\n")
        file.write("##################### VARIABLES TO SET #####################\n")
        file.write("############################################################\n\n")

        file.write("###########################################################################\n")
        file.write("##### These need to match what is specific to our system and our data #####\n")
        file.write("###########################################################################\n\n")

        file.write("## Path to ISA archive, only needed for saving a copy as metadata:\n")
        file.write(f"isa_archive:\n   \"{isa_zip}\"\n\n")

        file.write("## Path to runsheet:\n")
        file.write(f"runsheet:\n    \"{os.path.abspath(runsheet_file)}\"\n\n")

        file.write("## Set to \"PE\" for paired-end, \"SE\" for single-end.\n")
        file.write(f"data_type:\n    \"{data_type}\"\n\n")

        file.write("## single-column file with unique sample identifiers:\n")
        file.write("sample_info_file:\n    \"unique-sample-IDs.txt\"\n\n")

        file.write("## input reads directory (can be relative to workflow directory, or needs to be full path):\n")
        file.write(f"raw_reads_dir:\n    \"{raw_reads_directory}\"\n\n")

        file.write("## raw read suffixes:\n")
        file.write("  # e.g. for paired-end data, Sample-1_R1_raw.fastq.gz would be _R1_raw.fastq.gz for 'raw_R1_suffix' below\n")
        file.write("  # e.g. if single-end, Sample-1.fastq.gz would be .fastq.gz for 'raw_R1_suffix' below, and 'raw_R2_suffix' won't be used\n")
        file.write(f"raw_R1_suffix:\n    \"{raw_R1_suffix}\"\n")
        file.write(f"raw_R2_suffix:\n    \"{raw_R2_suffix}\"\n\n")

        file.write("## if we are trimming primers or not (\"TRUE\", or \"FALSE\")\n")
        file.write(f"trim_primers:\n    \"{trim_primers}\"\n\n")

        file.write("## primer sequences if we are trimming them (include anchoring symbols, e.g. '^', as needed, see: https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types)\n")
        file.write(f"F_primer:\n    \"{primer_anchor}{f_primer}\"\n")
        file.write(f"R_primer:\n    \"{primer_anchor}{r_primer}\"\n\n")

        # For linked primers
        file.write("## should cutadapt treat these as linked primers? (https://cutadapt.readthedocs.io/en/stable/recipes.html#trimming-amplicon-primers-from-paired-end-reads)\n")
        file.write(f"primers_linked:\n    \"{primers_linked}\"\n\n")
        file.write("## if primers are linked, we need to provide them as below, where the second half, following three periods, is the other primer reverse-complemented\n")
        file.write(f"  # (can reverse complement while retaining ambiguous bases at this site: http://arep.med.harvard.edu/labgc/adnan/projects/Utilities/revcomp.html)\n")
        file.write(f"  # include anchoring symbols, e.g. '^', as needed, see: https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types\n")
        file.write(f"F_linked_primer:\n    \"{f_linked_primer}\"\n")
        file.write(f"R_linked_primer:\n    \"{r_linked_primer}\"\n\n")

        file.write("## discard untrimmed, sets the \"--discard-untrimmed\" option if TRUE\n")
        file.write(f"discard_untrimmed:\n    \"{discard_untrimmed}\"\n\n")

        file.write("## target region (16S, 18S, or ITS is acceptable)\n")
        file.write("  # this determines which reference database is used for taxonomic classification\n")
        file.write("  # all are pulled from the pre-packaged DECIPHER downloads page here: http://www2.decipher.codes/Downloads.html\n")
        file.write("  # 16S uses SILVA\n")
        file.write("  # ITS uses UNITE\n")
        file.write("  # 18S uses PR2\n")
        file.write(f"target_region:\n    \"{target_region}\"\n\n")

        file.write("## concatenate only with dada2 instead of merging paired reads if TRUE\n")
        file.write("  # this is typically used with primers like 515-926, that captured 18S fragments that are typically too long to merge\n")
        file.write("  # note that 16S and 18S should have been separated already prior to running this workflow\n")
        file.write("  # this should likely be left as FALSE for any option other than \"18S\" above\n\n")

        file.write(f"concatenate_reads_only:\n    \"{concatenate_reads_only}\"\n\n")
        file.write(f"## values to be passed to dada2's filterAndTrim() function:\n")
        file.write(f"left_trunc:\n    {left_trunc}\n")
        file.write(f"right_trunc:\n    {right_trunc}\n")
        file.write(f"left_maxEE:\n    {left_maxEE}\n")
        file.write(f"right_maxEE:\n    {right_maxEE}\n\n")

        file.write("## minimum length threshold for cutadapt\n")
        file.write(f"min_cutadapt_len:\n    {min_trimmed_length}\n\n")

        file.write("######################################################################\n")
        file.write("##### The rest only need to be altered if we want to change them #####\n")
        file.write("######################################################################\n\n")

        file.write("## filename suffixes\n")
        file.write("primer_trimmed_R1_suffix:\n    \"_R1_trimmed.fastq.gz\"\n")
        file.write("primer_trimmed_R2_suffix:\n    \"_R2_trimmed.fastq.gz\"\n\n")

        file.write("filtered_R1_suffix:\n    \"_R1_filtered.fastq.gz\"\n")
        file.write("filtered_R2_suffix:\n    \"_R2_filtered.fastq.gz\"\n\n")

        file.write("## output prefix (if needed to distinguish from multiple primer sets, leave as empty string if not, include connecting symbol if adding, e.g. \"ITS-\")\n")
        file.write(f"output_prefix:\n    \"{output_prefix}\"\n\n")

        file.write("## output directories (all relative to processing directory, they will be created if needed)\n")
        file.write(f"info_out_dir:\n    \"{info_out_dir}\"\n")
        file.write(f"fastqc_out_dir:\n    \"{fastqc_out_dir}\"\n")
        file.write(f"trimmed_reads_dir:\n    \"{trimmed_reads_dir}\"\n")
        file.write(f"filtered_reads_dir:\n    \"{filtered_reads_dir}\"\n")
        file.write(f"final_outputs_dir:\n    \"{final_outputs_dir}\"\n")
        file.write(f"plots_dir:\n    \"{plots_dir}\"\n\n")

        # For general info and example usage command
        file.write("############################################################\n")
        file.write("###################### GENERAL INFO ########################\n")
        file.write("############################################################\n")
        file.write("# Workflow is currently equipped to work with paired-end data only, and reads are expected to be gzipped\n\n")
        file.write("## example usage command ##\n")
        file.write("# snakemake --use-conda --conda-prefix ${CONDA_PREFIX}/envs -j 2 -p\n\n")
        file.write("# `--use-conda` – this specifies to use the conda environments included in the workflow\n")
        file.write("# `--conda-prefix` – this allows us to point to where the needed conda environments should be stored...\n")
        file.write("# `-j` – this lets us set how many jobs Snakemake should run concurrently...\n")
        file.write("# `-p` – specifies to print out each command being run to the screen\n\n")
        file.write("# See `snakemake -h` for more options and details.\n")
    print("config.yaml was successfully created.")

# Example usage
# create_config_yaml(runsheet_df, uses_urls)


def main():
    # Argument parser setup with short argument names and an automatic help option
    parser = argparse.ArgumentParser(
        description='Run workflow for GeneLab data processing.',
        add_help=True,
        usage='%(prog)s [options]'  # Custom usage message
    )
    
    parser.add_argument('-o', '--OSD',
                        metavar='osd_number',
                        help='Set up the Snakemake workflow for a GeneLab OSD dataset and pull necessary read files and metadata. Acceptable formats: ###, OSD-###, GLDS-###',
                        type=str)
    
    parser.add_argument('-t', '--target',
                        choices=['16S', '18S', 'ITS'],
                        help='Specify the genomic target for the assay. Options: 16S, 18S, ITS. This is used to select the appropriate dataset from an OSD study when multiple options are available.',
                        type=str)
    
    parser.add_argument('-r', '--runsheetPath',
                        metavar='/path/to/runsheet.csv',
                        help='Set up the Snakemake workflow using a specified runsheet file.',
                        type=str)

    parser.add_argument('-x', '--run',
                        metavar='command',
                        nargs='?',
                        const="snakemake --use-conda --conda-prefix ${CONDA_PREFIX}/envs -j 2 -p",
                        type=str,
                        help='Specifies the command used to execute the snakemake workflow; Default: "snakemake --use-conda --conda-prefix ${CONDA_PREFIX}/envs -j 2 -p"')

    parser.add_argument('-d', '--outputDir',
                        metavar='/path/to/outputDir/',
                        default='./workflow_output/',  # Default value
                        help='Specifies the output directory for the output files generated by the workflow. Default: ./workflow_output/',
                        type=str)
    
    parser.add_argument('--specify-runsheet',
                        help='Specifies the runsheet for an OSD dataset by name. Only used if there are multiple datasets with the same target in the study.',
                        metavar='runsheet_name',
                        type=str)
    
    parser.add_argument('--trim-primers',
                        choices=['TRUE', 'FALSE'],
                        default='TRUE',
                        help='Specifies to trim primers (TRUE) or not (FALSE). Default: TRUE',
                        type=str)

    parser.add_argument('-m', '--min_trimmed_length',
                        metavar='length',
                        default=130,  # Default value
                        help='Specifies the MINIMUM length of trimmed reads. For paired-end data: if one read gets filtered, both reads are discarded. Default: 130',
                        type=int)
    
    parser.add_argument('--primers-linked',
                        choices=['TRUE', 'FALSE'],
                        default='TRUE',
                        help='If set to TRUE, instructs cutadapt to treat the primers as linked. Default: TRUE',
                        type=str)

    parser.add_argument('--anchor-primers',
                        choices=['TRUE', 'FALSE'],
                        default='TRUE',
                        help='Indicates if primers should be anchored (TRUE) or not (FALSE). Default: TRUE',
                        type=str)

    parser.add_argument('--discard-untrimmed',
                        choices=['TRUE', 'FALSE'],
                        default='TRUE',
                        help='If set to TRUE, instructs cutadapt to remove reads if the primers were not found in the expected location; if FALSE, these reads are kept. Default: TRUE',
                        type=str)

    parser.add_argument('--left-trunc',
                        default=0,
                        help='Specifies the length of the forwards reads, bases beyond this length will be truncated and reads shorter than this length are discarded. Default: 0 (no truncation)',
                        metavar='length',
                        type=int)

    parser.add_argument('--right-trunc',
                        default=0,
                        help='Specifies the length of the reverse reads, bases beyond this length will be truncated and reads shorter than this length are discarded. Default: 0 (no truncation)',
                        metavar='length',
                        type=int)

    parser.add_argument('--left-maxEE',
                        default=1,
                        help='Specifies the maximum expected error (maxEE) allowed for each forward read, reads with higher than maxEE will be discarded. Default: 1',
                        metavar='max_error',
                        type=int)

    parser.add_argument('--right-maxEE',
                        default=1,
                        help='Specifies the maximum expected error (maxEE) allowed for each forward read, reads with higher than maxEE will be discarded. Default: 1',
                        metavar='max_error',
                        type=int)

    parser.add_argument('--concatenate_reads_only',
                        choices=['TRUE', 'FALSE'],
                        default='FALSE',
                        help='If set to TRUE, specifies to concatenate forward and reverse reads only with dada2 instead of merging paired reads. Default: FALSE',
                        type=str)

    parser.add_argument('--output-prefix',
                        default='',
                        help='Specifies the prefix to use on all output files to distinguish multiple primer sets, leave as an empty string if only one primer set is being processed. Default: ""',
                        metavar='prefix',
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
    
    output_dir = args.outputDir
    min_trimmed_length = args.min_trimmed_length
    target = args.target
    isa_zip = ""

    # If OSD is used, pull ISA metadata for the study, create and select the runsheet
    if args.OSD:
        accession_number = process_osd_argument(args.OSD)
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

                # Create the 'unique-sample-IDs.txt' file and download read files if necessary
                if uses_urls:
                    handle_url_downloads(runsheet_df, output_file='unique-sample-IDs.txt')
                else:
                    sample_IDs_from_local(runsheet_df, output_file='unique-sample-IDs.txt')

                # Create the config.yaml file
                create_config_yaml(isa_zip=isa_zip,  
                                    runsheet_file=runsheet_file, 
                                    runsheet_df=runsheet_df, 
                                    uses_urls=uses_urls,
                                    output_dir=output_dir,
                                    min_trimmed_length=args.min_trimmed_length,
                                    trim_primers=args.trim_primers,
                                    primers_linked=args.primers_linked,
                                    anchor_primers=args.anchor_primers,
                                    discard_untrimmed=args.discard_untrimmed,
                                    left_trunc=args.left_trunc,
                                    right_trunc=args.right_trunc,
                                    left_maxEE=args.left_maxEE,
                                    right_maxEE=args.right_maxEE,
                                    concatenate_reads_only=args.concatenate_reads_only,
                                    output_prefix=args.output_prefix
                )
                
                print("Snakemake workflow setup is complete.")
            else:
                print("Failed to validate the runsheet file.", file=sys.stderr)
                sys.exit(1)
        else:
            print("No runsheet file specified.", file=sys.stderr)
            sys.exit(1)
    
    # Run the snakemake workflow if --run is used
    if args.run:
        snakemake_command = args.run if args.run is not None else "snakemake --use-conda --conda-prefix ${CONDA_PREFIX}/envs -j 2 -p"
        print(f"Running Snakemake command: {snakemake_command}")
        subprocess.run(snakemake_command, shell=True, check=True)
    
        # # Remove sample ID file
        # with open('config.yaml', 'r') as file:
        #     config_data = yaml.safe_load(file)
        #     sample_info_file = config_data.get('sample_info_file', '')  # Default to empty string if not found

        # if sample_info_file and os.path.exists(sample_info_file):
        #     os.remove(sample_info_file)
        
        # if isa_zip:
        #     try:
        #         os.remove(isa_zip)
        #     except FileNotFoundError:
        #         pass  # Ignore file not found error silently
        #     except Exception:
        #         pass 
        # # Remove all files if OSD run
        # if args.OSD:
        #     os.remove(runsheet_file)  # Assuming runsheet_file is a variable holding the file name
        #     os.remove("config.yaml")  # Ensure this is the correct file name

        # if args.runsheetPath:
        #     os.remove("config.yaml")  # Ensure this is the correct file name



if __name__ == "__main__":
    main()