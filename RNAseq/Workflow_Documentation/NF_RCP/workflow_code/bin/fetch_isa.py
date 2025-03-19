#!/usr/bin/env python3

"""
Downloads the ISA archive for a given dataset using the provided OSD ID.
"""

import argparse
import os
import sys
import requests

def main():
    parser = argparse.ArgumentParser(
        description="Download ISA archive for a given dataset using OSD ID."
    )
    parser.add_argument("--osd", required=True, help="OSD ID (e.g., OSD-576)")
    parser.add_argument("--outdir", required=True, help="Output directory to save the ISA archive")
    args = parser.parse_args()

    outdir = args.outdir

    if not args.osd.startswith('OSD-'):
        sys.exit(f"OSD accession ({args.osd}) was not provided in the correct format. It must start with 'OSD-'")

    # Build the JSON URL to get file information for a file with ISA in the name and a .zip extension
    json_url = f"https://visualization.osdr.nasa.gov/biodata/api/v2/dataset/{args.osd}/files/*ISA*.zip"

    # Fetch the JSON data
    response = requests.get(json_url)
    if response.status_code != 200:
        sys.exit(f"Error: Failed to retrieve file information. HTTP status code: {response.status_code}")

    # Parse the JSON data
    data = response.json()
    study_key = f"{args.osd}"
    study_data = data.get(study_key)
    if not study_data:
        sys.exit(f"Error: Study with OSD ID '{args.osd}' not found in the response.")

    # Get the list of study files
    files = study_data.get('files')
    if not files:
        sys.exit("Error: No files found in the response.")

    # Find the ISA archive file
    try:
        isa_file_name = next(filename for filename in files 
                         if 'ISA' in filename and filename.endswith('.zip'))
    except StopIteration:
        sys.exit("Error: ISA archive not found in the file list.")

    # Construct the full download URL
    download_url = files.get(isa_file_name).get('URL')
    if not download_url:
        sys.exit(f"Error: Download URL for ISA archive {isa_file_name} not found.")

    # Make the output directory and file path
    os.makedirs(outdir, exist_ok=True)
    output_file = os.path.join(outdir, isa_file_name)

    # Download the ISA archive
    with requests.get(download_url, stream=True) as r:
        if r.status_code == 200:
            with open(output_file, 'wb') as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
            print(f"ISA archive downloaded and saved to {output_file}")
        else:
            sys.exit(f"Error: Failed to download ISA archive. HTTP status code: {r.status_code}")

if __name__ == "__main__":
    main()