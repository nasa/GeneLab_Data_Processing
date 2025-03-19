#!/usr/bin/env python

import requests
import argparse
import re
import sys
import json

def get_osd_and_glds(accession, api_url):
    # Fetch data from the API
    try:
        response = requests.get(api_url)
        response.raise_for_status()
        data = response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error fetching data from API: {e}", file=sys.stderr)
        sys.exit(1)
    except json.JSONDecodeError:
        print("Error decoding JSON response from API", file=sys.stderr)
        sys.exit(1)

    osd_accession = None
    glds_accessions = []

    # Check if the accession is OSD or GLDS
    if accession.startswith('OSD-'):
        osd_accession = accession
        
        # Find this OSD in the data
        # Search in wildcard endpoint results
        for osd_id, osd_data in data.items():
            if osd_id == accession:
                metadata = osd_data.get("metadata", {})
                identifiers = metadata.get("identifiers", "")
                glds_accessions = re.findall(r'GLDS-\d+', identifiers)
                break
    
    elif accession.startswith('GLDS-'):
        glds_accessions = [accession]
        
        # Find the OSD associated with this GLDS
        for osd_id, osd_data in data.items():
            metadata = osd_data.get("metadata", {})
            identifiers = metadata.get("identifiers", "")
            if accession in identifiers:
                osd_accession = metadata.get("accession")
                break
    else:
        print("Invalid accession format. Please use 'OSD-###' or 'GLDS-###'.", file=sys.stderr)
        sys.exit(1)

    if not osd_accession or not glds_accessions:
        print(f"No data found for {accession}", file=sys.stderr)
        sys.exit(1)

    return osd_accession, glds_accessions

def main():
    parser = argparse.ArgumentParser(description="Retrieve OSD and GLDS accessions.")
    parser.add_argument('--accession', required=True, help="Accession in the format 'OSD-###' or 'GLDS-###'")
    parser.add_argument('--api_url', required=True, help="OSDR API URL")
    args = parser.parse_args()

    osd_accession, glds_accessions = get_osd_and_glds(args.accession, args.api_url)

    # Output the results in a way that Nextflow can capture
    print(f"{osd_accession}")
    print(f"{','.join(glds_accessions)}")

if __name__ == "__main__":
    main()
