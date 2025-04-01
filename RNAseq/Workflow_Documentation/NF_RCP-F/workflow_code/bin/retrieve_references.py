#! /usr/bin/env python
""" Retrieves reference fasta and gtf annotations from ensembl.
Script should be compatible with releases from 73 and onward.
Older releases have different FTP directory formats and/or are missing checksums.
"""
from ftplib import FTP
from ftplib import error_perm
import subprocess
import argparse
from pathlib import Path

import pandas as pd
from tqdm import tqdm

def _parse_args():
  """ Parse command line args.
  """
  parser = argparse.ArgumentParser()
  parser.add_argument('--ensembl_version', metavar='103', required=True,
                      help='Ensembl release number')
  parser.add_argument('--organism', metavar='animal_animalus', required=True,
                      help='Scientific name for organism.')
  parser.add_argument('--target', metavar='toplevel', required=True, choices=["toplevel","primary_assembly"],
                      help='Fetch the toplevel or primary_assembly.')
  args = parser.parse_args()
  return args
args = _parse_args()


# parsed from CLI args
ENSEMBL_VERSION = int(args.ensembl_version)
ORGANISM = args.organism.capitalize()

# sets string modifiers needed to access the correct ftp server
if ORGANISM.lower() in ["arabidopsis_thaliana"]:
    non_vertebrate_mod1 = "genomes"
    non_vertebrate_mod2 = "plants/"
else:
    non_vertebrate_mod1 = ''
    non_vertebrate_mod2 = ''
# constants
SERVER = f"ftp.ensembl{non_vertebrate_mod1}.org"
RELEASE_FOLDER = f"pub/{non_vertebrate_mod2}release-{ENSEMBL_VERSION}"
FASTA_FOLDER = f"fasta/{ORGANISM.lower()}/dna"
TARGET_FASTA_SUFFIX = f"dna.{args.target}.fa.gz"
GTF_FOLDER = f"gtf/{ORGANISM.lower()}"
TARGET_GTF_SUFFIX = f"{ENSEMBL_VERSION}.gtf.gz"
 # although named txt, this is tab separated
TEMP_FILE = "tmp_ftp_data.txt"

def download_from_ftp(target_file_suffix, target_file_folder, verbose = False):
    """ Downloads both fasta and gtf annotations as well as verifies by checksum
    """
    ################################################################################
    # Login To FTP Server And Navigate To Fasta Folder
    ################################################################################
    ftp = FTP(SERVER)
    print(f"FTP: Logging into {SERVER}") if verbose else None
    ftp.login()

    # try,except block to capture any bad ftp naviations.
    try:
        print(f"FTP: Navigating into folder {RELEASE_FOLDER}")
        ftp.cwd(RELEASE_FOLDER)

        # cd into the fasta folder
        print(f"FTP: Searching in {target_file_folder}") if verbose else None
        ftp.cwd(target_file_folder)
    except error_perm as e:
        raise ValueError("FTPLIB ERROR RAISED: {e}.  This usually indicates the organism '{ORGANISM}' could not be located on Ensembl using this script.  This script may need updating if the organism is on ensembl as this indicates the ftp path isn't formed correctly.")

    ################################################################################
    # Find Fasta File
    ################################################################################
    with open(TEMP_FILE, "w") as f:
        def _write_file_list(line):
            f.write(f"{line}\n")
        ftp.retrlines("NLST", _write_file_list)
        #z = f.readlines()
    with open(TEMP_FILE, "r") as f:
        target_file = [filename for filename
                      in f.read().split("\n") if
                      target_file_suffix in filename]
        assert len(target_file) == 1, \
               f"Only one should match but {len(target_file)} matched"

        # extract from list
        target_file = target_file[0]

    ################################################################################
    # Download Fasta File
    ################################################################################
    file_size = ftp.size(target_file)
    file_size_str = f"{file_size/float(1<<30):.2f} GB"

    with open(target_file, 'wb') as f:
        with tqdm(total=file_size,
                  unit='B', unit_scale=True, unit_divisor=1024,
                  disable= not verbose) as pbar:
            def cb(data):
                pbar.update(len(data))
                f.write(data)

            print(f"FTP: Found target file: {target_file}") if verbose else None
            print("Downloading") if verbose else None
            ftp.retrbinary(f"RETR {target_file}", cb)
            print(f"Successfully downloaded: {target_file}") if verbose else None


    ################################################################################
    # Verify by CHECKSUM
    ################################################################################

    ### retrieve ensembl recorded checksum ###
    checksum_file = "CHECKSUMS"
    checksum = dict()
    def _checksum_for(line):
        if target_file in line:
            print(f"this one: {line}") if verbose else None
            # lines look like 42368 823064 Mus_musculus.GRCm38.dna.toplevel.fa.gz
            # we only want to first two numbers (w/o lead/trail whitespace)
            # should look like: [42368,823064]
            checksum[target_file] = [int(val) for val in line.split()[0:2]]
    ftp.retrlines(f"RETR {checksum_file}", _checksum_for)
    ensembl_checksum = checksum[target_file]
    print(f"Ensembl check sum for {target_file}: {ensembl_checksum}")

    ### compute checksum for downloaded file ###
    computed_checksum = subprocess.run(["sum", target_file], capture_output=True).stdout
    # convert to string from bytes
    # strip whitespace
    computed_checksum = [int(val) for val in computed_checksum.decode("utf-8").strip().split()[0:2]]

    ### Check if they match ###
    print(f"Successfully downloaded recorded checksums: {ensembl_checksum}") if verbose else None
    print(f"Successfully computed checksums: {computed_checksum}") if verbose else None
    checksum_match = ensembl_checksum == computed_checksum
    if checksum_match:
        print(f"Checksum matches: {checksum_match}") if verbose else None
    else:
        raise ValueError(f"Checksum for downloaded file ({target_file}) and the one recorded by ensembl FTP server do not match!")
    ftp.close()

if __name__ == "__main__":
    try:
        download_from_ftp(TARGET_FASTA_SUFFIX, FASTA_FOLDER, verbose = True)
    except AssertionError:
        dummy_path = f"D.N.E.{TARGET_FASTA_SUFFIX}".replace(".gz","") # remove gzipped extension, as would occur if the file existed
        print(f"Creating dummy file: '{dummy_path}'.  This will be downstream checked by Nextflow, specifically, the D.N.E. prefix")
        Path(dummy_path).touch()
    download_from_ftp(TARGET_GTF_SUFFIX, GTF_FOLDER, verbose = True)
