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

    osd_id = args.osd.split('-')[-1] if '-' in args.osd else args.osd
    outdir = args.outdir

    # Build the JSON URL to get file information
    json_url = f"https://osdr.nasa.gov/osdr/data/osd/files/{osd_id}"

    # Fetch the JSON data
    response = requests.get(json_url)
    if response.status_code != 200:
        sys.exit(f"Error: Failed to retrieve file information. HTTP status code: {response.status_code}")

    # Parse the JSON data
    data = response.json()
    study_key = f"OSD-{osd_id}"
    study_data = data.get("studies", {}).get(study_key)
    if not study_data:
        sys.exit(f"Error: Study with OSD ID {osd_id} not found in the response.")

    # Get the list of study files
    file_info_list = study_data.get("study_files", [])
    if not file_info_list:
        sys.exit("Error: No files found in the response.")

    # Find the ISA archive file
    isa_file_info = next(
        (file_info for file_info in file_info_list
         if 'ISA' in file_info.get('file_name', '') and file_info['file_name'].endswith('.zip')),
        None
    )
    if not isa_file_info:
        sys.exit("Error: ISA archive not found in the file list.")

    # Construct the full download URL
    remote_url = isa_file_info.get('remote_url')
    if not remote_url:
        sys.exit("Error: Download URL for ISA archive not found.")
    download_url = f"https://osdr.nasa.gov{remote_url}"

    # Make the output directory and file path
    os.makedirs(outdir, exist_ok=True)
    output_file = os.path.join(outdir, isa_file_info['file_name'])

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