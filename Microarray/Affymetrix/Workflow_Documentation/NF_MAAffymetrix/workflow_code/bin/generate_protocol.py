#!/usr/bin/env python
"""
This script generates a protocol text file for GeneLab Affymetrix Microarray data processing.
It reads software versions from a YAML file and incorporates other parameters.
"""

import argparse
import os
import sys
import re
from datetime import datetime
from pathlib import Path
import pandas as pd
import yaml

def parse_args():
    """Sets up and parses the input parameters using argparse

    Returns:
        NameSpace: object holding input parameters 
    """
    parser = argparse.ArgumentParser(description='Generate protocol file for GeneLab Affymetrix Microarray pipeline')
    parser.add_argument('--outdir', required=True, type=Path,
                        help='Output directory for the protocol file')
    parser.add_argument('--software_table', required=True, type=Path,
                        help='Path to YAML file containing software versions')
    parser.add_argument('--assay_suffix', default='',
                        help='Suffix for the GeneLab assay type')
    parser.add_argument('--workflow_version', default='unknown',
                        help='Version of the NF_MAAffymetrix workflow manifest')
    parser.add_argument('--organism', required=False,
                        help='Organism name in the format "homo_sapiens"')
    parser.add_argument('--reference_source', required=False,
                        help='Source of the reference annotation')
    parser.add_argument('--reference_version', required=False,
                        help='Version of the reference annotation')
    parser.add_argument('--biomart_attribute', required=False,
                        help='Attribute for biomart query (if applicable)')
    parser.add_argument('--bioconductor_annotations', required=False, type=str,
                        help='bioconductor annotation database name')
    parser.add_argument('--annotations_db_info', required=False, type=Path,
                        help='Annotation DB info file from GeneLab Reference Annotation database used')
    parser.add_argument('--custom_annot_design', required=False, type=Path,
                        help='Path to the custom probe annotation design info file')
    parser.add_argument('--skip-DGE', type=bool, required=False,
                        help="Was DGE performed.")
    return parser.parse_args()


def read_annotation_versions(annot_info_file: Path, bioc_annot: str) -> dict[str, str]:
    """Reads the version information from the annotation DB info file into a dictionary mapping software name to version

    Args:
        annot_info_file (Path): the annotation info file produced by the GeneLab Reference Annotation pipeline

    Returns:
        dict: a dictionary mapping software name to version
    """
    versions = dict()
    try:
        with open(annot_info_file, 'r') as f:
            bioc_annot_found = False
            software_name = ''
            version_found = False
            doc_name_found = False
            for line in f:
                if bioc_annot_found:
                    versions[bioc_annot] = line.strip()
                    bioc_annot_found = False
                if version_found:
                    versions[software_name] = line.strip()
                    software_name = ''
                    version_found = False
                if doc_name_found:
                    versions['annot_doc'] = line.strip()
                    doc_name_found = False
                if line.startswith('Used'):
                    if re.search(bioc_annot, line):
                        bioc_annot_found = True
                    else:
                        version_found = True
                        software_name = line.split(' ')[1]
                if line.startswith('Based on:'):
                    doc_name_found = True
    except Exception as e:
        sys.stderr.write(f"Error reading annotation info file: {e}\n")
        sys.exit(1)
    
    return versions


def read_software_versions(yaml_file: Path) -> dict:
    """Reads the software versions YAML file into a dictionary mapping software name to version

    Args:
        yaml_file (Path): a YAML formatted file containing software version info produced by a Nextflow workflow

    Returns:
        dict: a dictionary mapping software name to version
    """
    try:
        with open(yaml_file, 'r') as f:
            return pd.DataFrame(yaml.safe_load(f)).set_index('name').to_dict()['version']
    except Exception as e:
        sys.stderr.write(f"Error reading software versions file: {e}\n")
        sys.exit(1)


def create_header(assay_suffix, workflow_version):
    """Create protocol header

    Args:
        assay_suffix (str): GeneLab assay suffix
        workflow_version (str): Workflow version

    Returns:
        str: Protocol header section
    """
    # Get current date
    current_date = datetime.now().strftime("%Y-%m-%d")

    # Create header
    header  = f"# GeneLab Microarray Pipeline Protocol{assay_suffix}\n"
    header += f"# Date: {current_date}\n\n"

    # Add appropriate protocol reference based on mode
    header += "Data were processed as described in GL-DPPD-7114-A "
    header += "(https://github.com/nasa/GeneLab_Data_Processing/blob/master/Microarray/Affymetrix/Pipeline_GL-DPPD-7114_Versions/GL-DPPD-7114-A.md), "
    
    header += f"using NF_MAAffymetrix version {workflow_version} "
    header += f"(https://github.com/nasa/GeneLab_Data_Processing/tree/NF_MAAffymetrix_{workflow_version}/Microarray/Affymetrix/Workflow_Documentation/NF_MAAffymetrix)."

    return header


def create_annot_and_dge_section(organism, biomart_attribute, ensembl_version, 
                                bioconductor_annotations, r_version, limma_version, annotation_versions,
                                custom_annots=pd.DataFrame(), skip_DGE=False):
    """Generate the probe annotation and DGE protocol sections

    Args:
        organism (str): full organism name lowercase with underscores instead of spaces (e.g., "homo_sapiens")
        biomart_attribute (str): Attribute for biomart query (name of the Agilent array as present in biomart or the custom design info file, if applicable)
        ensembl_version (str): Ensembl reference version
        bioconductor_annotations (str): Bioconductor annotation package name
        r_version (str): R software version
        limma_version (str): limma software version
        custom_annots (pd.DataFrame, optional): path to custom probe annotation design info file used during processing. Defaults to an empty path
        skip_DGE (bool, optional): Whether or not DGE was skipped during processing. Defaults to False.

    Returns:
        str: Protocol string for probe annotation and DGE
    """

    custom_annot_source_name = None
    custom_annot_download_link = None
    custom_annot_download_date = None
    custom_annot_create_date = None
    custom_annot_filename = None

    if not custom_annots.empty and biomart_attribute in custom_annots.index.values:
        custom_annot_source_name = custom_annots['annot_type'][biomart_attribute]
        custom_annot_filename = custom_annots['annot_filename'][biomart_attribute]
        if 'download_link' in custom_annots.columns:
            custom_annot_download_link = custom_annots['download_link'][biomart_attribute]
            custom_annot_download_date = custom_annots['download_date'][biomart_attribute]
        if 'create_date' in custom_annots.columns:
            custom_annot_create_date = custom_annots['create_date'][biomart_attribute]

    # Define versions for annotation package generation
    annot_doc = annotation_versions['annot_doc'] if 'annot_doc' in annotation_versions else 'GL-DPPD-7110-A'
    stringdb_version = annotation_versions['STRINGdb'] if 'STRINGdb' in annotation_versions else "2.16.4"
    pantherdb_version = annotation_versions['PANTHER.db'] if 'PANTHER.db' in annotation_versions else "1.0.12"
    
    # Define organism to annotation package mapping using scientific names
    organism_annotation_package = annotation_versions[bioconductor_annotations] if bioconductor_annotations in annotation_versions else "3.19.1"

    # Check if DGE was performed
    de_step = ""
    if not skip_DGE:
        de_step = f"Differential expression analysis was performed in R (version {r_version}) using limma (version {limma_version}); "
        de_step += "all groups were compared pairwise for each probeset to generate a moderated t-statistic and associated p- and adjusted p-value."
   
    organism_list=("homo_sapiens", "mus_musculus", "rattus_norvegicus", "drosophila_melanogaster", "caenorhabditis_elegans", "danio_rerio", "saccharomyces_cerevisiae")
    gene_mapping_step = ""

    # Case 1: a custom annotation source was used for this array design
    # no Ensembl FTP lookup happens in this case, which never calls
    # the Ensembl FTP helpers when annot_type is 'custom' or '3prime-IVT'.
    if biomart_attribute in custom_annots.index.values and custom_annot_source_name is not None:
        if '3prime-IVT' in custom_annot_source_name:
            gene_mapping_step = f"Gene annotations "
            annot_source = "3'-IVT"
        else:
            gene_mapping_step = "Annotations "
            annot_source = custom_annot_source_name.replace('_', ' ').title()
        gene_mapping_step += f"were retrieved for each probeset from {custom_annot_filename}, source: {annot_source}"
        if custom_annot_download_link is not None:
            created_text = ""
            if custom_annot_create_date is not None:
                created_text = f"created {custom_annot_create_date}, "
            gene_mapping_step += f" ({custom_annot_download_link}, {created_text}accessed {custom_annot_download_date})."
        else:
            gene_mapping_step += "."

    # Case 2: Ensembl FTP mart dump was used (no custom annotation)
    else:
        if organism == "arabidopsis_thaliana":
            database_name = "Plants Ensembl database"
            database_url = "plants.ensembl.org"
        elif organism in organism_list:
            database_name = "Ensembl database"
            database_url = "ensembl.org"
        else: # what should we do if the organism is not in the list??
            database_name = "TBD"
            database_url = "TBD"
        gene_mapping_step = f"Ensembl gene ID mappings were retrieved for each probeset using the {database_name} ftp server ({database_url}, release {ensembl_version})."

    # Gene annotations (STRINGdb/PANTHER/bioconductor merge) only runs when the
    # Ensembl FTP path succeeded (use_custom_annot == False in the QMD); every
    # custom-annotation branch (3prime-IVT, custom, NO_CUSTOM_ANNOT) skips this.
    annot_step = ""
    if custom_annot_source_name is None:
        annot_step = f"Gene annotations were assigned using the custom annotation tables generated in-house as detailed in {annot_doc} "
        annot_step += f"(https://github.com/nasa/GeneLab_Data_Processing/blob/master/GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/{annot_doc}/{annot_doc}.md), "
        annot_step += f"with STRINGdb (version {stringdb_version}), PANTHER.db (version {pantherdb_version}), and {bioconductor_annotations} (version {organism_annotation_package})."
    
    return " ".join((gene_mapping_step, de_step, annot_step))


def generate_protocol_content(args:argparse.Namespace, software_versions:dict, annotation_versions:dict, assay_suffix:str = "") -> str:
    """Generates a protocol string based on the input parameters and software versions

    Args:
        args (argparse): script input parameters
        software_versions (dict): a dictionary mapping software name to version
        assay_suffix (str): the Genelab assay suffix

    Returns:
        str: protocol text
    """

    header = create_header(assay_suffix, args.workflow_version)
        
    # Start building the description as a single paragraph
    # Add processing description with software versions
    oligo_version = software_versions.get('oligo', 'unknown')

    custom_annots = pd.DataFrame()
    if args.custom_annot_design.exists() and args.custom_annot_design.is_file():
        custom_annots = pd.read_csv(open(args.custom_annot_design, "r")).set_index('array_design')

    description = f"In short, a runsheet containing raw data file location and processing metadata from the study's *ISA.zip file was generated using dp_tools (version {software_versions.get('dp_tools', 'unknown')}). "
    description += f"The raw array data files were loaded into R (version {software_versions.get('R', 'unknown')}) using oligo (version {oligo_version}). "
    description += f"Raw data quality assurance density, pseudo image, MA, and boxplots were generated using oligo (version {oligo_version}). "
    description += f"The raw intensity data was background corrected and normalized across arrays via the oligo (version {oligo_version}) quantile method. "
    description += f"Normalized data quality assurance density, pseudo image, MA plots, and boxplots were generated using oligo (version {oligo_version}). "
    description += f"Normalized probe level data was summarized to the probeset level using the oligo (version {oligo_version}) RMA method."
    
    R_version = software_versions.get('R', 'unknown')
    limma_version = software_versions.get('limma', 'unknown')

    annot_and_dge_section = create_annot_and_dge_section(args.organism, args.biomart_attribute,
                                                         args.reference_version, args.bioconductor_annotations,
                                                         R_version, limma_version, 
                                                         annotation_versions, custom_annots, args.skip_DGE)

    # Create configuration section
    config = "\n\n\n## Configuration\n\n"

    # Add reference information if provided
    if hasattr(args, 'organism') and args.organism:
        config += f"- Organism: {args.organism}\n"
    if hasattr(args, 'reference_source') and args.reference_source:
        config += f"- Reference source: {args.reference_source}\n"
    if hasattr(args, 'reference_version') and args.reference_version:
        config += f"- Reference version: {args.reference_version}\n"
    
    # Add custom reference information if provided
    if hasattr(args, 'biomart_attribute') and args.biomart_attribute in custom_annots.index.values:
        config += f"- Array design: {args.biomart_attribute}\n"
        config += f"- Custom annotation source: {custom_annots['annot_type'][args.biomart_attribute]}\n"
        config += f"- Custom annotation file: {custom_annots['annot_filename'][args.biomart_attribute]}\n"
        if custom_annots['download_link'][args.biomart_attribute] is not None:
            config += f"- Custom annotation download link: {custom_annots['download_link'][args.biomart_attribute]}\n"
            config += f"- Custom annotation download date: {custom_annots['download_date'][args.biomart_attribute]}\n"
    
    # Create software versions section
    sw_section = "\n\n## All Software Versions\n\n"
    for software, version in software_versions.items():
        sw_section += f"- {software}: {version}\n"
    
    # Combine all sections
    return header + " " + description + " " + annot_and_dge_section + config + sw_section


def main():
    args = parse_args()
    
    # Read software versions from YAML file
    software_versions = read_software_versions(args.software_table)

    annotation_versions = dict()
    if args.annotations_db_info is not None and args.bioconductor_annotations is not None:
        annotation_versions = read_annotation_versions(args.annotations_db_info, args.bioconductor_annotations)

    assay_suffix = args.assay_suffix
    if not assay_suffix.startswith('_'):
        assay_suffix = f"_{assay_suffix}"
    
    # Generate protocol content
    protocol_content = generate_protocol_content(args, software_versions, annotation_versions, assay_suffix)
    
    
    # Write to output file
    output_file = os.path.join(args.outdir, f"protocol{assay_suffix}.txt")
    try:
        with open(output_file, 'w') as f:
            f.write(protocol_content)
        print(f"Protocol file generated successfully: {output_file}")
    except Exception as e:
        sys.stderr.write(f"Error writing protocol file: {e}\n")
        sys.exit(1)

if __name__ == "__main__":
    main()
