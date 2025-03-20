#!/usr/bin/env python

import os
import sys
import pandas as pd
import argparse
import glob
import csv
import datetime
import string
import enum
from pathlib import Path
import itertools
import numpy as np
import re

#############################################################################
# Differential Gene Expression (DGE) Validation Checks
#############################################################################
# Checks from dp_tools; ERCC results checks replace with ERCC gene presence checks
# - expect presence in unnormalized counts and absence in normalized counts
#
# File Existence Checks:
# - check_deseq2_normcounts_existence: DESeq2 normalized counts files (.tsv)
# - check_deseq2_dge_existence: DESeq2 DGE output files (.tsv)
# - check_ercc_presence: ERCC genes in unnormalized counts, absent in normalized counts
#
# Metadata Validation:
# - check_sample_table_against_runsheet: All runsheet samples in sample table
# - check_sample_table_for_correct_group_assignments: Sample group assignments match runsheet
# - check_contrasts_table_headers: Contrast headers match expected comparisons (e.g., A vs B)
# - check_contrasts_table_rows: Contrast rows contain correct group names
#
# DGE Table Content Validation:
# - check_dge_table_annotation_columns_exist: Required gene ID and annotation columns
# - check_dge_table_sample_columns_exist: All sample count columns present
# - check_dge_table_sample_columns_constraints: Sample counts ≥ 0
# - check_dge_table_group_columns_exist: Group mean/stdev columns for each condition
# - check_dge_table_group_columns_constraints: Group stats match manual calculation
# - check_dge_table_comparison_statistical_columns_exist: Log2fc/Stat/P.value/Adj.p.value columns
# - check_dge_table_group_statistical_columns_constraints: No nulls in Log2fc/Stat, no negatives in p-values
# - check_dge_table_fixed_statistical_columns_exist: All.mean/All.stdev/LRT.p.value columns
# - check_dge_table_fixed_statistical_columns_constraints: No nulls/negatives in means/stdevs, no negative p-values
# - check_dge_table_log2fc_within_reason: Log2fc values match direct calculation, signs correct
#
#############################################################################

def parse_runsheet(runsheet_path):
    """
    Parse the runsheet to extract sample information and dataset metadata.
    
    Args:
        runsheet_path: Path to the runsheet CSV file
        
    Returns:
        tuple: (sample_names, paired_end_values, metadata_dict)
            - sample_names: List of sample names
            - paired_end_values: Dictionary of paired-end values by sample
            - metadata_dict: Dictionary of metadata values
    """
    if not os.path.exists(runsheet_path):
        print(f"Error: Runsheet file not found: {runsheet_path}")
        sys.exit(1)
    
    try:
        # Try to read the runsheet using pandas
        df = pd.read_csv(runsheet_path)
        
        # Check for required columns
        required_columns = ['Sample Name', 'paired_end']
        missing_columns = [col for col in required_columns if col not in df.columns]
        
        if missing_columns:
            print(f"Error: Runsheet missing required columns: {', '.join(missing_columns)}")
            sys.exit(1)
        
        # Extract sample names and paired_end values
        sample_names = df['Sample Name'].tolist()
        paired_end_values = dict(zip(df['Sample Name'], df['paired_end']))
        
        # Check for tissue/group column (future support)
        tissue_groups = {}
        if 'tissue' in df.columns:
            tissue_groups = dict(zip(df['Sample Name'], df['tissue']))
        
        # Additional metadata that might be useful
        metadata = {}
        for col in df.columns:
            if col not in ['Sample Name']:
                # For columns that have the same value for all samples, store a single value
                unique_values = df[col].unique()
                if len(unique_values) == 1:
                    metadata[col] = unique_values[0]
        
        return sample_names, paired_end_values, tissue_groups, metadata
    
    except Exception as e:
        print(f"Error parsing runsheet: {str(e)}")
        sys.exit(1)

def check_directory_structure(outdir):
    """
    Check if the expected directory structure exists.
    
    Args:
        outdir: Output directory path
        
    Returns:
        bool: True if structure is valid, False otherwise
    """
    # Check for main directories
    required_dirs = [
        os.path.join(outdir, "04-DESeq2_NormCounts"),
        os.path.join(outdir, "05-DESeq2_DGE")
    ]
    
    # Optional rRNArm directories
    optional_dirs = [
        os.path.join(outdir, "04-DESeq2_NormCounts_rRNArm"),
        os.path.join(outdir, "05-DESeq2_DGE_rRNArm")
    ]
    
    missing_dirs = [d for d in required_dirs if not os.path.isdir(d)]
    
    if missing_dirs:
        print(f"Warning: Missing required directories: {', '.join(missing_dirs)}")
        return False
    
    return True

def initialize_vv_log(outdir):
    """Initialize or append to the VV_log.csv file."""
    vv_log_path = os.path.join(outdir, "VV_log.csv")
    
    # Create new file with header - force creation of header
    with open(vv_log_path, 'w') as f:
        f.write("component,sample_id,check_name,status,message,details\n")
    
    return vv_log_path

def log_check_result(log_path, component, sample_id, check_name, status, message="", details=""):
    """Log check result to the VV_log.csv file."""
    def escape_field(field, is_details=False):
        # Convert to string if not already
        field_str = str(field)
        
        # For details field, replace commas with semicolons to avoid CSV quoting
        if is_details and ',' in field_str:
            field_str = field_str.replace(', ', '; ')
        
        # Always apply consistent format - minimal quoting
        if any(c in field_str for c in ',"\n'):
            # If the field contains commas, quotes, or newlines, wrap in quotes and escape internal quotes
            return '"' + field_str.replace('"', '""') + '"'
        
        # Don't quote simple fields without special characters
        return field_str
    
    with open(log_path, 'a') as f:
        f.write(f"{escape_field(component)},{escape_field(sample_id)},{escape_field(check_name)},"
                f"{escape_field(status)},{escape_field(message)},{escape_field(details, True)}\n")

def check_deseq2_normcounts_existence(outdir, log_path, assay_suffix="_GLbulkRNAseq", mode="default"):
    """Check if DESeq2 normalized counts files exist.
    
    Args:
        outdir: Path to the NormCounts directory (not the parent outdir)
        log_path: Path to the log file
        assay_suffix: Suffix to append to filenames
        mode: Processing mode, either 'default' or 'microbes'
        
    Returns:
        bool: True if all required files exist, False otherwise
    """
    component = "dge_NormCounts"
    
    # Strip any stratification suffix from component name
    if "_" in os.path.basename(outdir):
        stratum = os.path.basename(outdir).split("_")[-1]
        if "rRNArm" in stratum:
            component = f"{component}_{stratum}"
        else:
            component = f"{component}_{stratum}"
    
    # Set the unnormalized counts filename based on mode
    if mode == "microbes":
        unnorm_file = f"FeatureCounts_Unnormalized_Counts{assay_suffix}.csv"
    else:
        unnorm_file = f"RSEM_Unnormalized_Counts{assay_suffix}.csv"
    
    # List of expected files
    expected_files = [
        unnorm_file,
        f"Normalized_Counts{assay_suffix}.csv",
        f"VST_Normalized_Counts{assay_suffix}.csv"
    ]
    
    missing_files = []
    
    for file in expected_files:
        file_path = os.path.join(outdir, file)
        if not os.path.exists(file_path):
            missing_files.append(file_path)
    
    if missing_files:
        print(f"WARNING: Missing DESeq2 normalized counts files: {', '.join(missing_files)}")
        log_check_result(log_path, component, "all", "check_deseq2_normcounts_existence", "RED", 
                         f"Missing DESeq2 normalized counts files", f"Missing files: {missing_files}")
        return False
    else:
        print(f"All DESeq2 normalized counts files found")
        log_check_result(log_path, component, "all", "check_deseq2_normcounts_existence", "GREEN", 
                         f"All DESeq2 normalized counts files found", "")
        return True

def check_deseq2_dge_existence(outdir, log_path, assay_suffix="_GLbulkRNAseq"):
    """Check if DESeq2 DGE files exist.
    
    Args:
        outdir: Path to the DGE directory (not the parent outdir)
        log_path: Path to the log file
        assay_suffix: Suffix to append to filenames
        
    Returns:
        bool: True if all required files exist, False otherwise
    """
    component = "dge"
    
    # Strip any stratification suffix from component name
    if "_" in os.path.basename(outdir):
        stratum = os.path.basename(outdir).split("_")[-1]
        if "rRNArm" in stratum:
            component = f"{component}_{stratum}"
        else:
            component = f"{component}_{stratum}"
    
    expected_files = [
        f"contrasts{assay_suffix}.csv",
        f"differential_expression{assay_suffix}.csv",
        f"SampleTable{assay_suffix}.csv"
    ]
    
    missing_files = []
    
    for file in expected_files:
        file_path = os.path.join(outdir, file)
        if not os.path.exists(file_path):
            missing_files.append(file_path)
    
    if missing_files:
        print(f"WARNING: Missing DESeq2 DGE files: {', '.join(missing_files)}")
        log_check_result(log_path, component, "all", "check_deseq2_dge_existence", "RED", 
                         f"Missing DESeq2 DGE files", f"Missing files: {missing_files}")
        return False
    else:
        print(f"All DESeq2 DGE files found")
        log_check_result(log_path, component, "all", "check_deseq2_dge_existence", "GREEN", 
                         "All DESeq2 DGE files found", "")
        return True

def r_style_make_names(s: str) -> str:
    """Recreates R's make.names function for individual strings.
    This function is often used to create syntactically valid names in R which are then saved in R outputs.
    Source: https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/make.names

    Args:
        s (str): A string to convert

    Returns:
        str: A string converted in the same way as R's make.names function
    """
    EXTRA_WHITELIST_CHARACTERS = "_ΩπϴλθijkuΑΒΓΔΕΖΗΘΙΚΛΜΝΞΟΠΡΣΤΥΦΧΨΩαβγδεζηθικλμνξοπρστυφχψω_µ" # Note: there are two "μμ" like characters one is greek letter mu, the other is the micro sign
    VALID_CHARACTERS = string.ascii_letters + string.digits + "." + EXTRA_WHITELIST_CHARACTERS
    REPLACEMENT_CHAR = "."
    new_string_chars = list()
    for char in s:
        if char in VALID_CHARACTERS:
            new_string_chars.append(char)
        else:
            new_string_chars.append(REPLACEMENT_CHAR)
    return "".join(new_string_chars)


class GroupFormatting(enum.Enum):
    r_make_names = enum.auto()
    ampersand_join = enum.auto()


def check_sample_table_against_runsheet(outdir, runsheet_path, log_path, assay_suffix="_GLbulkRNAseq", 
                                     stratum_factor="", stratum_value=""):
    """Check if the sample table includes all samples as denoted in the runsheet.
    
    Args:
        outdir: Path to the DGE directory (not the parent outdir)
        runsheet_path: Path to the runsheet CSV file
        log_path: Path to the log file
        assay_suffix: Suffix to append to the file name
        stratum_factor: Factor used for stratification (if any)
        stratum_value: Value of the stratum to check (if any)
    
    Returns:
        bool: True if check passes, False otherwise
    """
    component = "dge"
    check_name = "check_sample_table_against_runsheet"
    
    # Adjust the suffix and component name if stratified
    if stratum_value:
        check_suffix = f"_{stratum_value}"
        component_name = f"{component}{check_suffix}"
    else:
        check_suffix = ""
        component_name = component
    
    sample_table_path = os.path.join(outdir, f"SampleTable{assay_suffix}.csv")
    
    # Check if sample table exists
    if not os.path.exists(sample_table_path):
        print(f"WARNING: Sample table not found: {sample_table_path}")
        log_check_result(log_path, component_name, "all", check_name, "RED", 
                        f"Sample table not found", f"Expected at: {sample_table_path}")
        return False
    
    try:
        # Data specific preprocess
        df_rs = pd.read_csv(runsheet_path)
        
        # Filter runsheet by stratum if needed
        if stratum_factor and stratum_value:
            factor_col = f"Factor Value[{stratum_factor}]"
            if factor_col in df_rs.columns:
                # Remove "_rRNArm" suffix if present for filtering
                filter_value = stratum_value.replace("_rRNArm", "")
                
                # Filter samples that match this stratum value
                df_rs = df_rs[df_rs[factor_col].apply(
                    lambda x: r_style_make_names(str(x)) == filter_value
                )]
                
                print(f"Filtered runsheet to {len(df_rs)} samples for {stratum_factor}={filter_value}")
                
                if len(df_rs) == 0:
                    print(f"WARNING: No samples found in runsheet for {stratum_factor}={filter_value}")
                    log_check_result(log_path, component_name, "all", check_name, "RED", 
                                   f"No samples found in runsheet for stratum", 
                                   f"Stratum: {stratum_factor}={filter_value}")
                    return False
        
        df_rs = df_rs.set_index("Sample Name").sort_index()
        df_sample = pd.read_csv(sample_table_path, index_col=0).sort_index()

        # First check if there are any technical replicates in the runsheet
        tech_rep_pattern = "_techrep\\d+$"
        has_tech_reps = any(re.search(tech_rep_pattern, sample) for sample in df_rs.index)
        
        if has_tech_reps:
            # Function to get base sample name (remove _techrepX if present)
            def get_base_name(sample_name):
                return re.sub(tech_rep_pattern, "", sample_name)
            
            # Apply the function to get base sample names
            runsheet_base_names = set(get_base_name(sample) for sample in df_rs.index)
            sample_table_base_names = set(get_base_name(sample) for sample in df_sample.index)
            
            extra_samples = {
                "unique_to_runsheet": runsheet_base_names - sample_table_base_names,
                "unique_to_sampleTable": sample_table_base_names - runsheet_base_names,
            }
            
            print(f"Tech replicates detected in runsheet - checking base sample names")
        else:
            # No tech reps detected - use exact sample name matching
            extra_samples = {
                "unique_to_runsheet": set(df_rs.index) - set(df_sample.index),
                "unique_to_sampleTable": set(df_sample.index) - set(df_rs.index),
            }
            
            print(f"No tech replicates detected - checking exact sample names")

        # Check logic - all samples must be included
        all_samples_required = True
        if any([
            (extra_samples["unique_to_runsheet"] and all_samples_required),
            (extra_samples["unique_to_sampleTable"]),
        ]):
            print(f"WARNING: Sample mismatch between runsheet and sample table")
            for entry, values in extra_samples.items():
                if values:
                    print(f"  - {entry}: {values}")
            
            log_check_result(log_path, component_name, "all", check_name, "RED", 
                            f"Samples mismatched", f"Mismatched samples: {[f'{entry}:{v}' for entry, v in extra_samples.items() if v]}")
            return False
        else:
            samples_list = sorted(list(df_sample.index))
            samples_count = len(samples_list)
            print(f"All {samples_count} samples accounted for based on runsheet")
            
            # Include stratum info in the message if applicable
            stratum_info = f" for {stratum_factor}={stratum_value}" if stratum_factor and stratum_value else ""
            tech_rep_info = " (using base sample names)" if has_tech_reps else ""
            
            log_check_result(log_path, component_name, "all", check_name, "GREEN", 
                            f"All samples accounted for based on runsheet{stratum_info}{tech_rep_info}", 
                            f"Total samples: {samples_count}. Sample list: {'; '.join(samples_list)}")
            return True
    
    except Exception as e:
        print(f"ERROR checking sample table against runsheet: {str(e)}")
        log_check_result(log_path, component_name, "all", check_name, "RED", 
                        f"Error checking sample table", f"Error: {str(e)}")
        return False

def check_sample_table_for_correct_group_assignments(outdir, runsheet_path, log_path, assay_suffix="_GLbulkRNAseq",
                                                  stratum_factor="", stratum_value=""):
    """Check if the sample table is assigned to the correct experimental group.
    An experimental group is defined by the Factor Value columns found in the runsheet.
    
    Args:
        outdir: Path to the DGE directory (not the parent outdir)
        runsheet_path: Path to the runsheet CSV file
        log_path: Path to the log file
        assay_suffix: Suffix to append to the file name
        stratum_factor: Factor used for stratification (if any)
        stratum_value: Value of the stratum to check (if any)
    
    Returns:
        bool: True if check passes, False otherwise
    """
    component = "dge"
    check_name = "check_sample_table_for_correct_group_assignments"
    
    # Adjust the suffix and component name if stratified
    if stratum_value:
        check_suffix = f"_{stratum_value}"
        component_name = f"{component}{check_suffix}"
    else:
        check_suffix = ""
        component_name = component
    
    sample_table_path = os.path.join(outdir, f"SampleTable{assay_suffix}.csv")
    
    # Check if sample table exists
    if not os.path.exists(sample_table_path):
        print(f"WARNING: Sample table not found: {sample_table_path}")
        log_check_result(log_path, component_name, "all", check_name, "RED", 
                        f"Sample table not found", f"Expected at: {sample_table_path}")
        return False
    
    try:
        # Data specific preprocess
        df_sample = pd.read_csv(sample_table_path, index_col=0).sort_index()
        
        # Get factor values from runsheet
        df_rs = pd.read_csv(runsheet_path, dtype=str)
        
        # Filter runsheet by stratum if needed
        if stratum_factor and stratum_value:
            factor_col = f"Factor Value[{stratum_factor}]"
            if factor_col in df_rs.columns:
                # Remove "_rRNArm" suffix if present for filtering
                filter_value = stratum_value.replace("_rRNArm", "")
                
                # Filter samples that match this stratum value
                df_rs = df_rs[df_rs[factor_col].apply(
                    lambda x: r_style_make_names(str(x)) == filter_value
                )]
                
                print(f"Filtered runsheet to {len(df_rs)} samples for {stratum_factor}={filter_value}")
                
                if len(df_rs) == 0:
                    print(f"WARNING: No samples found in runsheet for {stratum_factor}={filter_value}")
                    log_check_result(log_path, component_name, "all", check_name, "RED", 
                                   f"No samples found in runsheet for stratum", 
                                   f"Stratum: {stratum_factor}={filter_value}")
                    return False
        
        df_rs = df_rs.set_index("Sample Name", drop=True)
        
        # Filter only Factor Value columns
        factor_cols = [col for col in df_rs.columns if col.startswith("Factor Value[")]
        
        if not factor_cols:
            print("WARNING: No Factor Value columns found in runsheet")
            log_check_result(log_path, component_name, "all", check_name, "RED", 
                            "No Factor Value columns found in runsheet", "")
            return False
        
        # Filter out the stratification factor if we're using it
        if stratum_factor and stratum_value:
            strat_factor_col = f"Factor Value[{stratum_factor}]"
            if strat_factor_col in factor_cols:
                # Remove the stratification factor from consideration
                factor_cols = [col for col in factor_cols if col != strat_factor_col]
                print(f"Excluding stratification factor '{stratum_factor}' from condition calculation")
        
        # If we've excluded all factors, add a warning
        if not factor_cols:
            print("WARNING: No remaining Factor Value columns found after excluding stratification factor")
            log_check_result(log_path, component_name, "all", check_name, "YELLOW", 
                           "No remaining factors for condition calculation after excluding stratification factor", 
                           f"Stratification factor '{stratum_factor}' was the only Factor Value column")
            return True
        
        df_rs = df_rs[factor_cols].loc[df_sample.index].sort_index()
        
        # Create expected conditions based on runsheet
        expected_conditions_based_on_runsheet = df_rs.apply(
            lambda x: "...".join(x), axis="columns"
        ).apply(r_style_make_names)
        
        # Check if conditions match
        mismatched_rows = expected_conditions_based_on_runsheet != df_sample["condition"]
        
        if not any(mismatched_rows):
            # Group samples by condition for reporting
            condition_to_samples = {}
            for sample, row in df_sample.iterrows():
                condition = row['condition']
                if condition not in condition_to_samples:
                    condition_to_samples[condition] = []
                condition_to_samples[condition].append(sample)
            
            # Format the details about group assignments
            group_details = []
            for condition, samples in sorted(condition_to_samples.items()):
                group_details.append(f"{condition}: {'; '.join(sorted(samples))}")
            
            # Include stratum info in the message if applicable
            stratum_info = f" for {stratum_factor}={stratum_value}" if stratum_factor and stratum_value else ""
            
            print(f"Conditions are formatted and assigned correctly for all {len(df_sample)} samples{stratum_info}")
            log_check_result(log_path, component_name, "all", check_name, "GREEN", 
                            f"Conditions are formatted and assigned correctly{stratum_info}", 
                            f"Sample to group assignments: {'; '.join(group_details)}")
            return True
        else:
            print("WARNING: Mismatch in expected conditions based on runsheet")
            mismatch_description = (
                df_sample[mismatched_rows]["condition"]
                + " <--SAMPLETABLE : RUNSHEET--> "
                + expected_conditions_based_on_runsheet[mismatched_rows]
            ).to_dict()
            
            for sample, mismatch in mismatch_description.items():
                print(f"  - {sample}: {mismatch}")
            
            log_check_result(log_path, component_name, "all", check_name, "RED", 
                            "Mismatch in expected conditions", f"Mismatched rows: {mismatch_description}")
            return False
    
    except Exception as e:
        print(f"ERROR checking sample table group assignments: {str(e)}")
        log_check_result(log_path, component_name, "all", check_name, "RED", 
                        f"Error checking group assignments", f"Error: {str(e)}")
        return False

def detect_stratification_factors(outdir, runsheet_path):
    """Auto-detect if the analysis is stratified based on directory structure.
    
    Args:
        outdir: Output directory path
        runsheet_path: Path to the runsheet CSV file
        
    Returns:
        tuple: (factor_name, factor_values) or (None, []) if no stratification detected
    """
    df = pd.read_csv(runsheet_path)
    
    # Get all Factor Value columns
    factor_cols = [col for col in df.columns if col.startswith("Factor Value[")]
    potential_factors = {}
    
    for col in factor_cols:
        # Extract factor name
        factor_name = col.replace("Factor Value[", "").replace("]", "")
        unique_values = df[col].unique()
        
        # Check if directories exist for each unique value
        for val in unique_values:
            safe_val = r_style_make_names(str(val))
            norm_dir = os.path.join(outdir, f"04-DESeq2_NormCounts_{safe_val}")
            dge_dir = os.path.join(outdir, f"05-DESeq2_DGE_{safe_val}")
            
            if os.path.exists(norm_dir) or os.path.exists(dge_dir):
                if factor_name not in potential_factors:
                    potential_factors[factor_name] = []
                potential_factors[factor_name].append(safe_val)
    
    # Return the factor with the most matching directories
    if potential_factors:
        factor, values = max(potential_factors.items(), key=lambda x: len(x[1]))
        return factor, values
    return None, []

def get_factor_stratified_paths(outdir, runsheet_path, target_factor=""):
    """Identify stratification factors and return the corresponding directory paths.
    
    Args:
        outdir: Output directory path
        runsheet_path: Path to the runsheet CSV file
        target_factor: Factor name to stratify by (if empty, auto-detect or use standard paths)
    
    Returns:
        dict: Dictionary mapping factor values to their directory paths.
              If no stratification is needed, returns a single entry with empty key.
    """
    # If no target factor specified, try to auto-detect
    if not target_factor:
        factor, values = detect_stratification_factors(outdir, runsheet_path)
        if factor:
            print(f"Auto-detected stratification by factor: '{factor}' with values: {values}")
            target_factor = factor
    
    # If still no factor or couldn't detect, return standard paths
    if not target_factor:
        return {"": {"norm_counts": os.path.join(outdir, "04-DESeq2_NormCounts"),
                     "dge": os.path.join(outdir, "05-DESeq2_DGE")}}
    
    # Read runsheet to identify values for target factor
    df = pd.read_csv(runsheet_path)
    factor_col = f"Factor Value[{target_factor}]"
    
    if factor_col not in df.columns:
        print(f"Warning: Stratification factor '{target_factor}' not found in runsheet")
        return {"": {"norm_counts": os.path.join(outdir, "04-DESeq2_NormCounts"),
                     "dge": os.path.join(outdir, "05-DESeq2_DGE")}}
    
    # Get unique values for the factor
    factor_values = df[factor_col].unique()
    result = {}
    
    # Create path mappings for each factor value
    for value in factor_values:
        safe_value = r_style_make_names(str(value))
        result[safe_value] = {
            "norm_counts": os.path.join(outdir, f"04-DESeq2_NormCounts_{safe_value}"),
            "dge": os.path.join(outdir, f"05-DESeq2_DGE_{safe_value}")
        }
        
        # Also check for rRNA removed paths if they exist
        rrna_norm_path = os.path.join(outdir, f"04-DESeq2_NormCounts_{safe_value}_rRNArm")
        rrna_dge_path = os.path.join(outdir, f"05-DESeq2_DGE_{safe_value}_rRNArm")
        
        if os.path.exists(rrna_norm_path) or os.path.exists(rrna_dge_path):
            result[f"{safe_value}_rRNArm"] = {
                "norm_counts": rrna_norm_path,
                "dge": rrna_dge_path
            }
    
    # Always include standard paths as fallback
    std_norm_path = os.path.join(outdir, "04-DESeq2_NormCounts")
    std_dge_path = os.path.join(outdir, "05-DESeq2_DGE")
    
    if os.path.exists(std_norm_path) or os.path.exists(std_dge_path):
        result[""] = {
            "norm_counts": std_norm_path,
            "dge": std_dge_path
        }
        
    return result

def print_summary(check_results, vv_log_path, overall_status="GREEN"):
    """Print a summary of the verification and validation checks.
    
    Args:
        check_results: Dictionary of check results
        vv_log_path: Path to the log file
        overall_status: Overall status (GREEN, YELLOW, RED)
        
    Returns:
        bool: True if there are no RED status checks, False otherwise
    """
    print("\n" + "="*80)
    print("VERIFICATION AND VALIDATION SUMMARY")
    print("="*80)
    
    # Read and process the VV log
    try:
        with open(vv_log_path, 'r') as f:
            reader = csv.DictReader(f)
            log_entries = list(reader)
            
        # Count statuses by component and check
        status_counts = {
            "GREEN": 0,
            "YELLOW": 0,
            "RED": 0
        }
        
        component_status = {}
        check_status = {}
        
        for entry in log_entries:
            component = entry['component']
            check = entry['check_name']
            status = entry['status']
            
            if status in status_counts:
                status_counts[status] += 1
            
            # Track status by component
            if component not in component_status:
                component_status[component] = {
                    "GREEN": 0,
                    "YELLOW": 0,
                    "RED": 0
                }
            component_status[component][status] += 1
            
            # Track status by check
            if check not in check_status:
                check_status[check] = {
                    "GREEN": 0,
                    "YELLOW": 0,
                    "RED": 0
                }
            check_status[check][status] += 1
        
        # Print status counts
        print(f"\nOverall Status: {overall_status}")
        print(f"Total Checks: {sum(status_counts.values())}")
        print(f"GREEN: {status_counts['GREEN']}")
        print(f"YELLOW: {status_counts['YELLOW']}")
        print(f"RED: {status_counts['RED']}")
        
        # Print component status
        print("\nStatus by Component:")
        for component, counts in component_status.items():
            total = sum(counts.values())
            status = "GREEN" if counts["RED"] == 0 else "RED"
            print(f"  {component}: {status} (GREEN: {counts['GREEN']}, YELLOW: {counts['YELLOW']}, RED: {counts['RED']})")
        
        # Print check status
        print("\nStatus by Check:")
        for check, counts in check_status.items():
            total = sum(counts.values())
            status = "GREEN" if counts["RED"] == 0 else "RED"
            print(f"  {check}: {status} (GREEN: {counts['GREEN']}, YELLOW: {counts['YELLOW']}, RED: {counts['RED']})")
        
        # Print specific check results
        print("\nSpecific Check Results:")
        
        check_display_names = {
            "file_existence_check": "DESeq2 NormCounts Existence",
            "dge_existence_check": "DESeq2 DGE Results Existence",
            "sample_table_check": "Sample Table Against Runsheet",
            "group_assignments_check": "Sample Table Group Assignments",
            "annotation_columns_check": "DGE Table Annotation Columns",
            "sample_columns_check": "DGE Table Sample Columns",
            "sample_columns_constraints_check": "DGE Table Sample Columns Constraints",
            "contrasts_headers_check": "Contrasts Table Headers",
            "contrasts_rows_check": "Contrasts Table Rows",
        }
        
        for check, result in check_results.items():
            if check != "overall":  # Skip the overall result
                display_name = check_display_names.get(check, check)
                result_str = "PASS" if result else "FAIL"
                status = "GREEN" if result else "RED"
                print(f"  {display_name}: {result_str} ({status})")
        
    except Exception as e:
        print(f"Error processing VV log: {str(e)}")
    
    print("\n" + "="*80)
    print(f"V&V log written to: {vv_log_path}")
    print("="*80 + "\n")
    
    # Return True if no RED statuses
    return status_counts.get('RED', 0) == 0

def check_contrasts_table_headers(outdir, runsheet_path, log_path, assay_suffix="_GLbulkRNAseq",
                               stratum_factor="", stratum_value=""):
    """Check if the contrasts table headers include expected comparisons from the runsheet.
    
    Args:
        outdir: Path to the DGE directory (not the parent outdir)
        runsheet_path: Path to the runsheet CSV file
        log_path: Path to the log file
        assay_suffix: Suffix to append to the file name
        stratum_factor: Factor used for stratification (if any)
        stratum_value: Value of the stratum to check (if any)
        
    Returns:
        bool: True if check passes, False otherwise
    """
    component = "dge"
    check_name = "check_contrasts_table_headers"
    
    # Adjust the suffix and component name if stratified
    if stratum_value:
        check_suffix = f"_{stratum_value}"
        component_name = f"{component}{check_suffix}"
    else:
        check_suffix = ""
        component_name = component
    
    # Look for contrasts table with the correct filename pattern
    contrasts_table_path = os.path.join(outdir, f"contrasts{check_suffix}{assay_suffix}.csv")
    
    # Check if contrasts table exists
    if not os.path.exists(contrasts_table_path):
        message = f"Contrasts table not found"
        log_check_result(log_path, component_name, "all", check_name, "RED", 
                        message, f"Expected at: {contrasts_table_path}")
        return False
    
    try:
        # Get expected groups from runsheet
        df_rs = pd.read_csv(runsheet_path)
        df_rs_factor_cols = df_rs[[col for col in df_rs.columns if col.startswith("Factor Value[")]]
        
        # Filter runsheet by stratum if needed
        if stratum_factor and stratum_value:
            factor_col = f"Factor Value[{stratum_factor}]"
            if factor_col not in df_rs.columns:
                log_check_result(log_path, component_name, "all", check_name, "RED", 
                                f"Stratification factor column not found in runsheet", 
                                f"Expected column: {factor_col}")
                return False
            df_rs = df_rs[df_rs[factor_col] == stratum_value]
            df_rs_factor_cols = df_rs_factor_cols[df_rs[factor_col] == stratum_value]
        
        # Check if we have any samples after filtering
        if df_rs.empty:
            log_check_result(log_path, component_name, "all", check_name, "RED", 
                            "No samples found in runsheet after filtering", 
                            f"Stratum factor: {stratum_factor}, Stratum value: {stratum_value}")
            return False
        
        # Get unique groups using the same formatting as dp_tools (ampersand_join)
        # First, get all unique combinations of factor values
        unique_factor_combinations = df_rs_factor_cols.drop_duplicates()
        
        # Format each combination like dp_tools does
        formatted_groups = []
        for _, row in unique_factor_combinations.iterrows():
            group = f"({' & '.join(row)})"
            formatted_groups.append(group)
        
        # Generate expected comparisons (all pairwise permutations)
        import itertools
        expected_comparisons = [
            f"{g1}v{g2}"
            for g1, g2 in itertools.permutations(formatted_groups, 2)
        ]
        
        # Read the contrasts table
        df_contrasts = pd.read_csv(contrasts_table_path, index_col=0)
        actual_comparisons = list(df_contrasts.columns)
        
        # Check if expected comparisons match actual comparisons
        differences = set(expected_comparisons).symmetric_difference(set(df_contrasts.columns))
        
        if not differences:
            status = "GREEN"
            message = "Contrasts table headers match expected comparisons"
            
            details = f"Found {len(expected_comparisons)} expected comparisons. Expected comparisons: {'; '.join(expected_comparisons)}. Actual comparisons: {'; '.join(actual_comparisons)}. All expected comparisons were found in the contrasts table."
            
            log_check_result(log_path, component_name, "all", check_name, status, message, details)
            return True
        else:
            print(f"WARNING: Contrasts table headers do not match expected comparisons")
            print(f"  - Expected: {expected_comparisons}")
            print(f"  - Actual: {actual_comparisons}")
            print(f"  - Missing: {differences}")
            print(f"  - Extra: {differences}")
            
            details = f"Differences found between expected and actual comparisons. Expected comparisons: {'; '.join(expected_comparisons)}. Actual comparisons: {'; '.join(actual_comparisons)}. Missing comparisons: {'; '.join(differences) if differences else 'None'}. Extra comparisons: {'; '.join(differences) if differences else 'None'}."
            
            log_check_result(log_path, component_name, "all", check_name, "RED", 
                             "Contrasts table headers do not match expected comparisons", details)
            return False
        
    except Exception as e:
        log_check_result(log_path, component_name, "all", check_name, "RED", 
                        "Error checking contrasts table headers", str(e))
        return False

def check_contrasts_table_rows(outdir, log_path, assay_suffix="_GLbulkRNAseq",
                           stratum_factor="", stratum_value=""):
    """Check if the contrasts table rows match expected formatting.
    
    Args:
        outdir: Path to the DGE directory (not the parent outdir)
        log_path: Path to the log file
        assay_suffix: Suffix to append to the file name
        stratum_factor: Factor used for stratification (if any)
        stratum_value: Value of the stratum to check (if any)
        
    Returns:
        bool: True if check passes, False otherwise
    """
    component = "dge"
    check_name = "check_contrasts_table_rows"
    
    # Adjust the suffix and component name if stratified
    if stratum_value:
        check_suffix = f"_{stratum_value}"
        component_name = f"{component}{check_suffix}"
    else:
        check_suffix = ""
        component_name = component
    
    # Look for contrasts table with the correct filename pattern
    contrasts_table_path = os.path.join(outdir, f"contrasts{check_suffix}{assay_suffix}.csv")
    
    # Check if contrasts table exists
    if not os.path.exists(contrasts_table_path):
        message = f"Contrasts table not found"
        log_check_result(log_path, component_name, "all", check_name, "RED", 
                        message, f"Expected at: {contrasts_table_path}")
        return False
    
    try:
        # Read the contrasts table
        df_contrasts = pd.read_csv(contrasts_table_path, index_col=0)
        
        def _get_groups_from_comparisons(s: str) -> set:
            """Extracts group names from a comparison string like 'G1vG2'.
            
            Args:
                s: Comparison string (e.g., '(Group1 & Factor2)v(Group2 & Factor2)')
                
            Returns:
                set: Set of group names
            """
            # Check if the format is '(G1)v(G2)' or just 'G1vG2'
            if '(' in s and ')v(' in s:
                g1, g2 = s.split(")v(")
                g1 = g1[1:]  # Remove leading parenthesis
                g2 = g2[:-1]  # Remove trailing parenthesis
            else:
                parts = s.split('v')
                # If there are multiple 'v's, we need to handle carefully
                if len(parts) != 2:
                    # This is likely malformed, but try to extract meaningful parts
                    g1, g2 = parts[0], parts[-1]
                else:
                    g1, g2 = parts
            
            # Handle special characters used in R formatting
            g1 = r_style_make_names(g1.replace(" & ", "..."))
            g2 = r_style_make_names(g2.replace(" & ", "..."))
            
            return {g1, g2}
        
        bad_columns = {}
        column_details = []
        
        for col_name, col_series in df_contrasts.items():
            try:
                # Get expected groups from the column name
                expected_values = _get_groups_from_comparisons(col_name)
                col_values = set(col_series.astype(str))
                
                # Record details for this column
                column_details.append(f"Column '{col_name}': Expected values: {', '.join(expected_values)}, Actual values: {', '.join(col_values)}")
                
                # Check if the column values match the expected values
                if expected_values != col_values:
                    # Sometimes R makes additional transformations to the names
                    # Let's check if the differences are just formatting issues
                    normalized_expected = {val.lower().replace('.', '') for val in expected_values}
                    normalized_actual = {val.lower().replace('.', '') for val in col_values}
                    
                    if normalized_expected != normalized_actual:
                        bad_columns[col_name] = {
                            "expected": expected_values,
                            "actual": col_values
                        }
            except Exception as e:
                # If parsing fails for a column, log it as a bad column
                column_details.append(f"Column '{col_name}': Error parsing expected values: {str(e)}")
                bad_columns[col_name] = {
                    "expected": f"Could not determine expected values: {str(e)}",
                    "actual": set(col_series)
                }
        
        if not bad_columns:
            print(f"Contrasts table rows match expected formatting")
            details = f"All {len(df_contrasts.columns)} comparisons have correct formatting. " + "; ".join(column_details)
            log_check_result(log_path, component_name, "all", check_name, "GREEN", 
                           "Contrasts table rows match expected formatting", details)
            return True
        else:
            print(f"WARNING: Contrasts table rows have formatting issues")
            
            # Prepare detailed error message
            error_details = []
            for col in bad_columns:
                info = bad_columns[col]
                error_details.append(
                    f"Column '{col}': Expected values: {'; '.join(str(x) for x in info['expected'])}, "
                    f"Actual values: {'; '.join(str(x) for x in info['actual'])}"
                )
            
            details = f"{len(bad_columns)} of {len(df_contrasts.columns)} columns have formatting issues: " + "; ".join(error_details) + "; All column details: " + "; ".join(column_details)
        
        log_check_result(log_path, component_name, "all", check_name, "RED", 
                        "Contrasts table rows do not match expected formatting", details)
        
        return False
        
    except Exception as e:
        log_check_result(log_path, component_name, "all", check_name, "RED", 
                        "Error checking contrasts table rows", str(e))
        return False

def check_dge_table_annotation_columns_exist(outdir, runsheet_path, log_path, assay_suffix="_GLbulkRNAseq",
                                        stratum_factor="", stratum_value=""):
    """Check if the DGE table includes annotation columns beyond just the gene ID.
    
    This verifies that the annotation module of the pipeline ran successfully by adding
    annotation columns (like SYMBOL, GENENAME, etc.) to the DGE table.
    
    Args:
        outdir: Path to the DGE directory (not the parent outdir)
        runsheet_path: Path to the runsheet CSV file
        log_path: Path to the log file
        assay_suffix: Suffix to append to the file name
        stratum_factor: Factor used for stratification (if any)
        stratum_value: Value of the stratum to check (if any)
        
    Returns:
        bool: True if check passes, False otherwise
    """
    component = "dge"
    check_name = "check_dge_table_annotation_columns_exist"
    
    # Adjust the suffix and component name if stratified
    if stratum_value:
        check_suffix = f"_{stratum_value}"
        component_name = f"{component}{check_suffix}"
    else:
        check_suffix = ""
        component_name = component
    
    # First get sample names from runsheet to identify which columns are samples
    try:
        df_rs = pd.read_csv(runsheet_path)
        
        # Filter runsheet by stratum if needed
        if stratum_factor and stratum_value:
            factor_col = f"Factor Value[{stratum_factor}]"
            if factor_col not in df_rs.columns:
                log_check_result(log_path, component_name, "all", check_name, "RED", 
                                f"Stratification factor column not found in runsheet", 
                                f"Expected column: {factor_col}")
                return False
            df_rs = df_rs[df_rs[factor_col] == stratum_value]
        
        # Check if we have any samples after filtering
        if df_rs.empty:
            log_check_result(log_path, component_name, "all", check_name, "RED", 
                            "No samples found in runsheet after filtering", 
                            f"Stratum factor: {stratum_factor}, Stratum value: {stratum_value}")
            return False
        
        # Get sample names
        sample_names = set(df_rs['Sample Name'])
        
    except Exception as e:
        log_check_result(log_path, component_name, "all", check_name, "RED", 
                        "Error parsing runsheet", str(e))
        return False
    
    # Get DGE table paths with the correct pattern
    dge_table_path = os.path.join(outdir, f"differential_expression{check_suffix}{assay_suffix}.csv")
    
    # Check if the DGE table exists
    if not os.path.exists(dge_table_path):
        message = f"DGE table not found"
        log_check_result(log_path, component_name, "all", check_name, "RED", 
                        message, f"Expected at: {dge_table_path}")
        return False
    
    try:
        # Read the DGE table
        df_dge = pd.read_csv(dge_table_path)
        
        # Identify columns that are not sample names
        all_columns = list(df_dge.columns)
        
        # The first column is typically the gene ID column
        gene_id_column = all_columns[0]
        
        # Identify sample columns by finding which column names match sample names from runsheet
        sample_columns = [col for col in all_columns if col in sample_names]
        
        # If we can't find any sample columns based on runsheet, fallback to heuristic
        if not sample_columns:
            # Heuristic: Look for columns that follow a common pattern for sample names
            sample_columns = [col for col in all_columns if "Rep" in col or "_R" in col]
        
        if not sample_columns:
            status = "RED"
            message = "Could not identify sample columns in DGE table"
            details = f"Available columns: {', '.join(all_columns)}"
            log_check_result(log_path, component_name, "all", check_name, status, message, details)
            return False
        
        # Find the index of the first sample column
        first_sample_idx = min([all_columns.index(col) for col in sample_columns])
        
        # Get annotation columns (between gene ID and first sample)
        annotation_columns = all_columns[1:first_sample_idx]
        
        if not annotation_columns:
            status = "RED"
            message = "No annotation columns found in DGE table"
            details = f"Gene ID column: {gene_id_column}; No annotation columns were found between gene ID and sample columns."
            log_check_result(log_path, component_name, "all", check_name, status, message, details)
            return False
        else:
            # Found mandatory annotation columns
            print(f"Found {len(annotation_columns)} annotation columns in DGE table")
            log_check_result(log_path, component_name, "all", check_name, "GREEN", 
                           f"Found {len(annotation_columns)} annotation columns in DGE table", 
                           f"Gene ID column: {gene_id_column}; Annotation columns: {'; '.join(annotation_columns)}")
            return True
        
    except Exception as e:
        log_check_result(log_path, component_name, "all", check_name, "RED", 
                        "Error checking annotation columns", str(e))
        return False

def check_dge_table_sample_columns_exist(outdir, runsheet_path, log_path, assay_suffix="_GLbulkRNAseq",
                                          stratum_factor="", stratum_value=""):
    """Check if all sample columns from the runsheet exist in the DGE table.
    
    Args:
        outdir: Path to the DGE directory (not the parent outdir)
        runsheet_path: Path to the runsheet CSV file
        log_path: Path to the log file
        assay_suffix: Suffix to append to the file name
        stratum_factor: Factor used for stratification (if any)
        stratum_value: Value of the stratum to check (if any)
        
    Returns:
        bool: True if check passes, False otherwise
    """
    component = "dge"
    check_name = "check_dge_table_sample_columns_exist"
    
    # Adjust the suffix and component name if stratified
    if stratum_value:
        check_suffix = f"_{stratum_value}"
        component_name = f"{component}{check_suffix}"
    else:
        check_suffix = ""
        component_name = component
    
    # First get sample names from runsheet
    try:
        df_rs = pd.read_csv(runsheet_path)
        
        # Filter runsheet by stratum if needed
        if stratum_factor and stratum_value:
            factor_col = f"Factor Value[{stratum_factor}]"
            if factor_col not in df_rs.columns:
                log_check_result(log_path, component_name, "all", check_name, "RED", 
                                f"Stratification factor column not found in runsheet", 
                                f"Expected column: {factor_col}")
                return False
            df_rs = df_rs[df_rs[factor_col] == stratum_value]
        
        # Check if we have any samples after filtering
        if df_rs.empty:
            log_check_result(log_path, component_name, "all", check_name, "RED", 
                            "No samples found in runsheet after filtering", 
                            f"Stratum factor: {stratum_factor}, Stratum value: {stratum_value}")
            return False
        
        # Get sample names
        sample_names = list(df_rs['Sample Name'])
        
        # Check if there are any technical replicates in the sample names
        tech_rep_pattern = "_techrep\\d+$"
        has_tech_reps = any(re.search(tech_rep_pattern, sample) for sample in sample_names)
        
        # Function to get base sample name (remove _techrepX if present)
        def get_base_name(sample_name):
            return re.sub(tech_rep_pattern, "", sample_name)
        
        # If tech reps are present, convert to base sample names
        if has_tech_reps:
            expected_samples = set(get_base_name(sample) for sample in sample_names)
        else:
            expected_samples = set(sample_names)
        
    except Exception as e:
        log_check_result(log_path, component_name, "all", check_name, "RED", 
                        "Error parsing runsheet", str(e))
        return False
    
    # Get DGE table path with the correct pattern
    dge_table_path = os.path.join(outdir, f"differential_expression{check_suffix}{assay_suffix}.csv")
    
    # Check if the DGE table exists
    if not os.path.exists(dge_table_path):
        message = f"DGE table not found"
        log_check_result(log_path, component_name, "all", check_name, "RED", 
                        message, f"Expected at: {dge_table_path}")
        return False
    
    try:
        # Read the DGE table
        df_dge = pd.read_csv(dge_table_path)
        dge_columns = list(df_dge.columns)
        
        # Find sample columns in the DGE table
        if has_tech_reps:
            # For each expected base sample name, check if it appears in DGE columns
            found_base_samples = set()
            for expected_base in expected_samples:
                for col in dge_columns:
                    # If the column matches exactly or its base name matches
                    if col == expected_base or get_base_name(col) == expected_base:
                        found_base_samples.add(expected_base)
                        break
            
            existing_samples = found_base_samples
        else:
            # Direct matching for non-tech rep datasets
            existing_samples = set(col for col in dge_columns if col in expected_samples)
        
        # Calculate missing samples
        missing_samples = expected_samples - existing_samples
        
        if not missing_samples:
            status = "GREEN"
            message = "All sample columns present in DGE table"
            details = f"Found all {len(expected_samples)} expected sample names in the DGE table."
            log_check_result(log_path, component_name, "all", check_name, status, message, details)
            return True
        else:
            status = "RED"
            message = "Some sample columns missing from DGE table"
            if has_tech_reps:
                details = f"{len(missing_samples)} of {len(expected_samples)} expected sample base names are missing from the DGE table. Missing samples: {'; '.join(missing_samples)}"
            else:
                details = f"{len(missing_samples)} of {len(expected_samples)} expected sample columns are missing from the DGE table. Missing columns: {'; '.join(missing_samples)}"
            log_check_result(log_path, component_name, "all", check_name, status, message, details)
            return False
        
    except Exception as e:
        log_check_result(log_path, component_name, "all", check_name, "RED", 
                        "Error checking sample columns", str(e))
        return False

def check_dge_table_sample_columns_constraints(outdir, runsheet_path, log_path, assay_suffix="_GLbulkRNAseq",
                                          stratum_factor="", stratum_value=""):
    """Check if the sample columns in the DGE table meet value constraints.
    
    Specifically, this checks that all count values in sample columns are >= 0.
    
    Args:
        outdir: Path to the DGE directory (not the parent outdir)
        runsheet_path: Path to the runsheet CSV file
        log_path: Path to the log file
        assay_suffix: Suffix to append to the file name
        stratum_factor: Factor used for stratification (if any)
        stratum_value: Value of the stratum to check (if any)
        
    Returns:
        bool: True if check passes, False otherwise
    """
    component = "dge"
    check_name = "check_dge_table_sample_columns_constraints"
    
    # Define the minimum allowed count value
    MINIMUM_COUNT = 0
    
    # Adjust the suffix and component name if stratified
    if stratum_value:
        check_suffix = f"_{stratum_value}"
        component_name = f"{component}{check_suffix}"
    else:
        check_suffix = ""
        component_name = component
    
    # First get sample names from runsheet
    try:
        df_rs = pd.read_csv(runsheet_path)
        
        # Filter runsheet by stratum if needed
        if stratum_factor and stratum_value:
            factor_col = f"Factor Value[{stratum_factor}]"
            if factor_col not in df_rs.columns:
                log_check_result(log_path, component_name, "all", check_name, "RED", 
                                f"Stratification factor column not found in runsheet", 
                                f"Expected column: {factor_col}")
                return False
            df_rs = df_rs[df_rs[factor_col] == stratum_value]
        
        # Check if we have any samples after filtering
        if df_rs.empty:
            log_check_result(log_path, component_name, "all", check_name, "RED", 
                            "No samples found in runsheet after filtering", 
                            f"Stratum factor: {stratum_factor}, Stratum value: {stratum_value}")
            return False
        
        # Get sample names
        sample_names = list(df_rs['Sample Name'])
        
        # Check if there are any technical replicates in the sample names
        tech_rep_pattern = "_techrep\\d+$"
        has_tech_reps = any(re.search(tech_rep_pattern, sample) for sample in sample_names)
        
    except Exception as e:
        log_check_result(log_path, component_name, "all", check_name, "RED", 
                        "Error parsing runsheet", str(e))
        return False
    
    # Get DGE table path with the correct pattern
    dge_table_path = os.path.join(outdir, f"differential_expression{check_suffix}{assay_suffix}.csv")
    
    # Check if the DGE table exists
    if not os.path.exists(dge_table_path):
        message = f"DGE table not found"
        log_check_result(log_path, component_name, "all", check_name, "RED", 
                        message, f"Expected at: {dge_table_path}")
        return False
    
    try:
        # Read the DGE table
        df_dge = pd.read_csv(dge_table_path)
        
        # Identify sample columns in the DGE table
        if has_tech_reps:
            # When there are tech reps, identify sample columns by looking for columns 
            # that either directly match a sample name or match a base sample name
            def get_base_name(sample_name):
                return re.sub(tech_rep_pattern, "", sample_name)
                
            # Extract base sample names from the runsheet
            base_sample_names = set(get_base_name(sample) for sample in sample_names)
            
            # Find columns in DGE table that correspond to samples
            sample_columns = []
            for col in df_dge.columns:
                # Check if this column is a sample column (matches sample name or base name)
                if col in sample_names or get_base_name(col) in base_sample_names:
                    sample_columns.append(col)
        else:
            # For datasets without tech reps, use direct matching
            sample_columns = [col for col in df_dge.columns if col in sample_names]
        
        if not sample_columns:
            status = "RED"
            message = "No sample columns found in DGE table"
            details = f"Could not identify sample columns among: {', '.join(df_dge.columns[:10])}..."
            log_check_result(log_path, component_name, "all", check_name, status, message, details)
            return False
        
        # Filter to only include sample columns
        df_dge_samples = df_dge[sample_columns]
        
        # Check if all values in each column meet the minimum count constraint
        columns_with_issues = []
        min_values = {}
        
        for col in sample_columns:
            if not (df_dge_samples[col] >= MINIMUM_COUNT).all():
                columns_with_issues.append(col)
                min_val = float(df_dge_samples[col].min())  # Convert numpy type to Python float
                min_values[col] = min_val
        
        constraint_description = f"All counts are greater than or equal to {MINIMUM_COUNT}"
        
        if not columns_with_issues:
            status = "GREEN"
            message = "All sample columns meet the minimum count constraint"
            details = f"All {len(sample_columns)} sample columns satisfy: {constraint_description}"
            log_check_result(log_path, component_name, "all", check_name, status, message, details)
            return True
        else:
            status = "RED"
            message = "Some sample columns have values below the minimum count constraint"
            details = f"{len(columns_with_issues)} of {len(sample_columns)} sample columns with counts below {MINIMUM_COUNT}: {'; '.join(columns_with_issues)}"
            log_check_result(log_path, component_name, "all", check_name, status, message, details)
            return False
        
    except Exception as e:
        log_check_result(log_path, component_name, "all", check_name, "RED", 
                        "Error checking sample columns", str(e))
        return False

def check_dge_table_group_columns_exist(outdir, runsheet_path, log_path, assay_suffix="_GLbulkRNAseq",
                                    stratum_factor="", stratum_value=""):
    """Check if all group summary statistic columns exist in the DGE table.
    
    This includes columns with prefixes 'Group.Mean_' and 'Group.Stdev_' for each experimental group.
    
    Args:
        outdir: Path to the DGE directory (not the parent outdir)
        runsheet_path: Path to the runsheet CSV file
        log_path: Path to the log file
        assay_suffix: Suffix to append to the file name
        stratum_factor: Factor used for stratification (if any)
        stratum_value: Value of the stratum to check (if any)
        
    Returns:
        bool: True if check passes, False otherwise
    """
    component = "dge"
    check_name = "check_dge_table_group_columns_exist"
    
    # Define the group column prefixes
    GROUP_PREFIXES = ["Group.Mean_", "Group.Stdev_"]
    
    # Adjust the suffix and component name if stratified
    if stratum_value:
        check_suffix = f"_{stratum_value}"
        component_name = f"{component}{check_suffix}"
    else:
        check_suffix = ""
        component_name = component
    
    # Get DGE table path with the correct pattern
    dge_table_path = os.path.join(outdir, f"differential_expression{check_suffix}{assay_suffix}.csv")
    
    # Check if the DGE table exists
    if not os.path.exists(dge_table_path):
        message = f"DGE table not found"
        log_check_result(log_path, component_name, "all", check_name, "RED", 
                      message, f"Expected at: {dge_table_path}")
        return False
    
    try:
        # First get expected groups from runsheet
        df_rs = pd.read_csv(runsheet_path)
        
        # Filter runsheet by stratum if needed
        if stratum_factor and stratum_value:
            factor_col = f"Factor Value[{stratum_factor}]"
            if factor_col not in df_rs.columns:
                log_check_result(log_path, component_name, "all", check_name, "RED", 
                                f"Stratification factor column not found in runsheet", 
                                f"Expected column: {factor_col}")
                return False
            df_rs = df_rs[df_rs[factor_col] == stratum_value]
        
        # Check if we have any samples after filtering
        if df_rs.empty:
            log_check_result(log_path, component_name, "all", check_name, "RED", 
                            "No samples found in runsheet after filtering", 
                            f"Stratum factor: {stratum_factor}, Stratum value: {stratum_value}")
            return False
        
        # Extract factor value columns from runsheet
        factor_cols = [col for col in df_rs.columns if col.startswith("Factor Value[")]
        
        # Group samples by their factor combinations
        groups = {}
        for _, row in df_rs.iterrows():
            factors = [row[col] for col in factor_cols]
            
            # Format in two ways: one with dots (for R-style) and one with parentheses and ampersands
            r_style_group = "...".join(factors)
            paren_style_group = f"({' & '.join(factors)})"
            
            if r_style_group not in groups:
                groups[r_style_group] = paren_style_group
        
        # Create a sorted list of group display names for reporting
        group_names = sorted([paren_style for paren_style in groups.values()])
        
        # Read DGE table and check for expected columns
        df_dge = pd.read_csv(dge_table_path)
        dge_columns = set(df_dge.columns)
        
        # Find missing columns - for each group, consider it missing only if both formats are missing
        missing_columns = []
        found_count = 0
        
        for prefix in GROUP_PREFIXES:
            for r_style_group, paren_style_group in groups.items():
                paren_col = f"{prefix}{paren_style_group}"
                dot_col = f"{prefix}{r_style_group}"
                
                if paren_col in dge_columns or dot_col in dge_columns:
                    found_count += 1
                else:
                    # Prefer to report the parentheses style in the error message
                    missing_columns.append(paren_col)
        
        if not missing_columns:
            status = "GREEN"
            message = "All group summary statistic columns present in DGE table"
            log_check_result(log_path, component_name, "all", check_name, "GREEN", 
                           "All group summary statistic columns present in DGE table",
                           f"Found {'; '.join(GROUP_PREFIXES)} columns for {len(group_names)} groups: {'; '.join(group_names)}")
            return True
        else:
            status = "RED"
            message = "Missing group summary statistic columns in DGE table"
            
            # Simplified error message
            details = f"Missing columns: {', '.join(missing_columns[:5])}"
            if len(missing_columns) > 5:
                details += f" and {len(missing_columns) - 5} more"
            
            log_check_result(log_path, component_name, "all", check_name, status, message, details)
            return False
        
    except Exception as e:
        log_check_result(log_path, component_name, "all", check_name, "RED", 
                        "Error checking group columns", str(e))
        return False

def check_dge_table_group_columns_constraints(outdir, runsheet_path, log_path, assay_suffix="_GLbulkRNAseq",
                                     stratum_factor="", stratum_value=""):
    """Check if the group statistics columns in the DGE table meet value constraints.
    
    Group statistics columns include Group.Mean_* and Group.Stdev_* for each group.
    
    Args:
        outdir: Path to the DGE directory (not the parent outdir)
        runsheet_path: Path to the runsheet CSV file
        log_path: Path to the log file
        assay_suffix: Suffix to append to the file name
        stratum_factor: Factor used for stratification (if any)
        stratum_value: Value of the stratum to check (if any)
        
    Returns:
        bool: True if check passes, False otherwise
    """
    component = "dge"
    check_name = "check_dge_table_group_columns_constraints"
    
    # Adjust the suffix and component name if stratified
    if stratum_value:
        check_suffix = f"_{stratum_value}"
        component_name = f"{component}{check_suffix}"
    else:
        check_suffix = ""
        component_name = component
    
    # First get sample names from runsheet
    try:
        df_rs = pd.read_csv(runsheet_path)
        
        # Filter runsheet by stratum if needed
        if stratum_factor and stratum_value:
            factor_col = f"Factor Value[{stratum_factor}]"
            if factor_col not in df_rs.columns:
                log_check_result(log_path, component_name, "all", check_name, "RED", 
                                f"Stratification factor column not found in runsheet", 
                                f"Expected column: {factor_col}")
                return False
            df_rs = df_rs[df_rs[factor_col] == stratum_value]
        
        # Check if we have any samples after filtering
        if df_rs.empty:
            log_check_result(log_path, component_name, "all", check_name, "RED", 
                            "No samples found in runsheet after filtering", 
                            f"Stratum factor: {stratum_factor}, Stratum value: {stratum_value}")
            return False
        
        # Get sample names
        sample_names = list(df_rs['Sample Name'])
        
        # Check if there are any technical replicates in the sample names
        tech_rep_pattern = "_techrep\\d+$"
        has_tech_reps = any(re.search(tech_rep_pattern, sample) for sample in sample_names)
        
    except Exception as e:
        log_check_result(log_path, component_name, "all", check_name, "RED", 
                        "Error parsing runsheet", str(e))
        return False
    
    # Get group information
    try:
        # Read sample table to get group assignments
        sample_table_path = os.path.join(outdir, f"SampleTable{check_suffix}{assay_suffix}.csv") 
        
        if not os.path.exists(sample_table_path):
            log_check_result(log_path, component_name, "all", check_name, "RED", 
                           "Sample table not found", 
                           f"Expected at: {sample_table_path}")
            return False
            
        df_sample_table = pd.read_csv(sample_table_path)
        
        # Check if 'condition' column exists
        if "condition" not in df_sample_table.columns:
            log_check_result(log_path, component_name, "all", check_name, "RED", 
                           "Condition column missing in sample table", 
                           f"Expected columns: condition")
            return False
        
        # Get group names
        unique_conditions = df_sample_table["condition"].unique()
        
    except Exception as e:
        log_check_result(log_path, component_name, "all", check_name, "RED", 
                       "Error processing sample table", str(e))
        return False
    
    # Get DGE table path with the correct pattern
    dge_table_path = os.path.join(outdir, f"differential_expression{check_suffix}{assay_suffix}.csv")
    
    # Check if the DGE table exists
    if not os.path.exists(dge_table_path):
        message = f"DGE table not found"
        log_check_result(log_path, component_name, "all", check_name, "RED", 
                       message, f"Expected at: {dge_table_path}")
        return False
    
    try:
        # Read the DGE table
        df_dge = pd.read_csv(dge_table_path)
        
        # Function to get base sample name (remove _techrepX if present)
        def get_base_name(sample_name):
            return re.sub(tech_rep_pattern, "", sample_name)
        
        # Check if all samples from the runsheet are in the DGE table
        missing_samples = []
        for sample in sample_names:
            if has_tech_reps:
                base_name = get_base_name(sample)
                # Check if either the sample or its base name is in the DGE columns
                if base_name not in df_dge.columns and sample not in df_dge.columns:
                    missing_samples.append(sample)
            else:
                if sample not in df_dge.columns:
                    missing_samples.append(sample)
        
        if missing_samples:
            log_check_result(log_path, component_name, "all", check_name, "RED", 
                           "Missing samples in DGE table", 
                           f"Cannot check group statistics; the following samples are missing: {'; '.join(missing_samples)}")
            return False
        
        # Create a set of expected columns for group statistics
        expected_group_prefixes = [
            "Group.Mean_", "Group.Stdev_"
        ]
        
        groups = []
        for condition in unique_conditions:
            # Convert condition to R-friendly group name parentheses format
            group_name = f"({condition.replace('.', ' ')})"
            groups.append(group_name)
        
        expected_columns = []
        for prefix in expected_group_prefixes:
            for group in groups:
                expected_columns.append(f"{prefix}{group}")
        
        # Check if all expected columns are present
        missing_columns = [col for col in expected_columns if col not in df_dge.columns]
        
        if missing_columns:
            log_check_result(log_path, component_name, "all", check_name, "RED", 
                           "Group summary statistics columns missing", 
                           f"Missing columns: {', '.join(missing_columns)}")
            return False
        
        # Check if the mean columns have no null values and are non-negative
        mean_columns = [col for col in expected_columns if col.startswith("Group.Mean_")]
        for col in mean_columns:
            # Check for nulls
            if df_dge[col].isnull().any():
                status = "RED"
                message = f"Group mean column contains null values"
                details = f"Column {col} has {df_dge[col].isnull().sum()} null values"
                log_check_result(log_path, component_name, "all", check_name, status, message, details)
                return False
            
            # Check for negative values
            if (df_dge[col] < 0).any():
                status = "RED"
                message = f"Group mean column contains negative values"
                details = f"Column {col} has {(df_dge[col] < 0).sum()} negative values"
                log_check_result(log_path, component_name, "all", check_name, status, message, details)
                return False
        
        # Check if the stdev columns have no null values and are non-negative
        stdev_columns = [col for col in expected_columns if col.startswith("Group.Stdev_")]
        for col in stdev_columns:
            # Check for nulls
            if df_dge[col].isnull().any():
                status = "RED"
                message = f"Group standard deviation column contains null values"
                details = f"Column {col} has {df_dge[col].isnull().sum()} null values"
                log_check_result(log_path, component_name, "all", check_name, status, message, details)
                return False
            
            # Check for negative values
            if (df_dge[col] < 0).any():
                status = "RED"
                message = f"Group standard deviation column contains negative values"
                details = f"Column {col} has {(df_dge[col] < 0).sum()} negative values"
                log_check_result(log_path, component_name, "all", check_name, status, message, details)
                return False
        
        status = "GREEN"
        message = "All group summary statistic columns meet constraints"
        details = f"Group mean and standard deviation columns for {len(groups)} groups have no null or negative values"
        log_check_result(log_path, component_name, "all", check_name, status, message, details)
        return True
        
    except Exception as e:
        log_check_result(log_path, component_name, "all", check_name, "RED", 
                       "Error checking group columns", str(e))
        return False

def check_dge_table_comparison_statistical_columns_exist(outdir, runsheet_path, log_path, assay_suffix="_GLbulkRNAseq",
                                                stratum_factor="", stratum_value=""):
    """Check if all comparison statistical columns exist in the DGE table.
    
    This includes columns with prefixes 'Log2fc_', 'Stat_', 'P.value_', and 'Adj.p.value_' 
    for each pairwise group comparison.
    
    Args:
        outdir: Path to the DGE directory (not the parent outdir)
        runsheet_path: Path to the runsheet CSV file
        log_path: Path to the log file
        assay_suffix: Suffix to append to the file name
        stratum_factor: Factor used for stratification (if any)
        stratum_value: Value of the stratum to check (if any)
        
    Returns:
        bool: True if check passes, False otherwise
    """
    component = "dge"
    check_name = "check_dge_table_comparison_statistical_columns_exist"
    
    # Define the comparison statistical column prefixes
    COMPARISON_PREFIXES = ["Log2fc_", "Stat_", "P.value_", "Adj.p.value_"]
    
    # Adjust the suffix and component name if stratified
    if stratum_value:
        check_suffix = f"_{stratum_value}"
        component_name = f"{component}{check_suffix}"
    else:
        check_suffix = ""
        component_name = component
    
    # Get DGE table path with the correct pattern
    dge_table_path = os.path.join(outdir, f"differential_expression{check_suffix}{assay_suffix}.csv")
    
    # Check if the DGE table exists
    if not os.path.exists(dge_table_path):
        message = f"DGE table not found"
        log_check_result(log_path, component_name, "all", check_name, "RED", 
                       message, f"Expected at: {dge_table_path}")
        return False
    
    try:
        # First get expected groups from runsheet
        df_rs = pd.read_csv(runsheet_path)
        
        # Filter runsheet by stratum if needed
        if stratum_factor and stratum_value:
            factor_col = f"Factor Value[{stratum_factor}]"
            if factor_col not in df_rs.columns:
                log_check_result(log_path, component_name, "all", check_name, "RED", 
                                f"Stratification factor column not found in runsheet", 
                                f"Expected column: {factor_col}")
                return False
            df_rs = df_rs[df_rs[factor_col] == stratum_value]
        
        # Check if we have any samples after filtering
        if df_rs.empty:
            log_check_result(log_path, component_name, "all", check_name, "RED", 
                            "No samples found in runsheet after filtering", 
                            f"Stratum factor: {stratum_factor}, Stratum value: {stratum_value}")
            return False
        
        # Extract factor value columns from runsheet
        factor_cols = [col for col in df_rs.columns if col.startswith("Factor Value[")]
        
        # Group samples by their factor combinations
        groups = {}
        for _, row in df_rs.iterrows():
            factors = [row[col] for col in factor_cols]
            
            # Format in two ways: one with dots (for R-style) and one with parentheses and ampersands
            r_style_group = "...".join(factors)
            paren_style_group = f"({' & '.join(factors)})"
            
            if r_style_group not in groups:
                groups[r_style_group] = paren_style_group
        
        # Create a sorted list of group display names for reporting
        group_names = sorted([paren_style for paren_style in groups.values()])
        
        # Generate all pairwise comparisons between groups
        expected_comparisons = []
        for g1, g2 in itertools.permutations(group_names, 2):
            comparison = f"{g1}v{g2}"
            expected_comparisons.append(comparison)
        
        # Generate expected column names
        expected_columns = []
        for prefix in COMPARISON_PREFIXES:
            for comparison in expected_comparisons:
                expected_columns.append(f"{prefix}{comparison}")
        
        # Read DGE table and check for expected columns
        df_dge = pd.read_csv(dge_table_path)
        dge_columns = set(df_dge.columns)
        
        # Find missing columns
        missing_columns = [col for col in expected_columns if col not in dge_columns]
        
        if not missing_columns:
            status = "GREEN"
            message = "All comparison statistical columns present in DGE table"
            log_check_result(log_path, component_name, "all", check_name, status, message, 
                            f"Found all statistical columns with prefixes: {'; '.join(COMPARISON_PREFIXES)} for {len(expected_comparisons)} comparisons: {'; '.join(expected_comparisons)}")
            return True
        else:
            status = "RED"
            message = "Missing comparison statistical columns in DGE table"
            details = f"Missing these statistical columns: {', '.join(missing_columns[:5])}"
            if len(missing_columns) > 5:
                details += f" and {len(missing_columns) - 5} more"
            log_check_result(log_path, component_name, "all", check_name, status, message, details)
            return False
        
    except Exception as e:
        log_check_result(log_path, component_name, "all", check_name, "RED", 
                        "Error checking comparison statistical columns", str(e))
        return False

def nonNull(series):
    """Check if a pandas Series has no null values.
    
    Args:
        series: Pandas Series to check
        
    Returns:
        bool: True if no null values, False otherwise
    """
    return not series.isnull().any()

def nonNegative(series):
    """Check if a pandas Series has only non-negative values (ignoring nulls).
    
    Args:
        series: Pandas Series to check
        
    Returns:
        bool: True if all non-null values are non-negative, False otherwise
    """
    return ((series >= 0) | (series.isnull())).all()

def onlyAllowedValues(series, allowed_values):
    """Check if a pandas Series contains only allowed values (ignoring nulls).
    
    Args:
        series: Pandas Series to check
        allowed_values: List of allowed values
        
    Returns:
        bool: True if all non-null values are in allowed_values, False otherwise
    """
    return ((series.isin(allowed_values)) | (series.isnull())).all()

def utils_common_constraints_on_dataframe(df, constraints):
    """Apply common constraints to a dataframe.
    
    Args:
        df: Pandas DataFrame
        constraints: List of tuples (column_set, constraint_dict)
            where column_set is a set of column names and constraint_dict defines the constraints
            
    Returns:
        dict: Dictionary of failed constraints with lists of column names
    """
    issues = {
        "Failed non null constraint": [],
        "Failed non negative constraint": [],
        "Failed allowed values constraint": []
    }
    
    for (col_set, col_constraints) in constraints:
        # Make a copy to avoid modifying the original
        constraints_copy = col_constraints.copy()
        
        # Only include columns that exist in the dataframe
        existing_cols = [col for col in col_set if col in df.columns]
        if not existing_cols:
            continue
            
        # Check each column against constraints
        for colname in existing_cols:
            colseries = df[colname]
            
            # Check non-null constraint
            if constraints_copy.get("nonNull", False) and not nonNull(colseries):
                issues["Failed non null constraint"].append(colname)
                
            # Check non-negative constraint
            if constraints_copy.get("nonNegative", False) and not nonNegative(colseries):
                issues["Failed non negative constraint"].append(colname)
                
            # Check allowed values constraint
            allowed_values = constraints_copy.get("allowedValues", None)
            if allowed_values and not onlyAllowedValues(colseries, allowed_values):
                issues["Failed allowed values constraint"].append(colname)
                
    return issues

def check_dge_table_group_statistical_columns_constraints(outdir, runsheet_path, log_path, assay_suffix="_GLbulkRNAseq",
                                                 stratum_factor="", stratum_value=""):
    """Check if the comparison statistical columns in the DGE table meet specific constraints.
    
    This validates that Log2fc, Stat, P.value, and Adj.p.value columns meet requirements:
    - Log2fc and Stat should not have any null values
    - P.value and Adj.p.value should not have negative values (but can have nulls)
    
    Args:
        outdir: Path to the DGE directory (not the parent outdir)
        runsheet_path: Path to the runsheet CSV file
        log_path: Path to the log file
        assay_suffix: Suffix to append to the file name
        stratum_factor: Factor used for stratification (if any)
        stratum_value: Value of the stratum to check (if any)
        
    Returns:
        bool: True if check passes, False otherwise
    """
    component = "dge"
    check_name = "check_dge_table_group_statistical_columns_constraints"
    
    # Adjust the suffix and component name if stratified
    if stratum_value:
        check_suffix = f"_{stratum_value}"
        component_name = f"{component}{check_suffix}"
    else:
        check_suffix = ""
        component_name = component
    
    # Get DGE table path with the correct pattern
    dge_table_path = os.path.join(outdir, f"differential_expression{check_suffix}{assay_suffix}.csv")
    
    # Check if the DGE table exists
    if not os.path.exists(dge_table_path):
        message = f"DGE table not found"
        log_check_result(log_path, component_name, "all", check_name, "RED", 
                       message, f"Expected at: {dge_table_path}")
        return False
    
    try:
        # First get expected groups from runsheet
        df_rs = pd.read_csv(runsheet_path)
        
        # Filter runsheet by stratum if needed
        if stratum_factor and stratum_value:
            factor_col = f"Factor Value[{stratum_factor}]"
            if factor_col not in df_rs.columns:
                log_check_result(log_path, component_name, "all", check_name, "RED", 
                                f"Stratification factor column not found in runsheet", 
                                f"Expected column: {factor_col}")
                return False
            df_rs = df_rs[df_rs[factor_col] == stratum_value]
        
        # Check if we have any samples after filtering
        if df_rs.empty:
            log_check_result(log_path, component_name, "all", check_name, "RED", 
                            "No samples found in runsheet after filtering", 
                            f"Stratum factor: {stratum_factor}, Stratum value: {stratum_value}")
            return False
        
        # Extract factor value columns from runsheet
        factor_cols = [col for col in df_rs.columns if col.startswith("Factor Value[")]
        
        # Group samples by their factor combinations
        groups = {}
        for _, row in df_rs.iterrows():
            factors = [row[col] for col in factor_cols]
            
            # Format in two ways: one with dots (for R-style) and one with parentheses and ampersands
            r_style_group = "...".join(factors)
            paren_style_group = f"({' & '.join(factors)})"
            
            if r_style_group not in groups:
                groups[r_style_group] = paren_style_group
        
        # Generate all pairwise comparisons between groups
        group_names = list(groups.values())
        expected_comparisons = []
        for g1, g2 in itertools.permutations(group_names, 2):
            comparison = f"{g1}v{g2}"
            expected_comparisons.append(comparison)
        
        # Define constraints for different types of columns
        constraints = [
            (set([f"Log2fc_{comp}" for comp in expected_comparisons]), {"nonNull": True}),
            (set([f"Stat_{comp}" for comp in expected_comparisons]), {"nonNull": True}),
            # P-values can be NA in DESeq2 for various reasons (low counts, outliers, etc.)
            (set([f"P.value_{comp}" for comp in expected_comparisons]), {"nonNegative": True}),
            (set([f"Adj.p.value_{comp}" for comp in expected_comparisons]), {"nonNegative": True}),
        ]
        
        # Read DGE table
        df_dge = pd.read_csv(dge_table_path)
        
        # Apply constraints
        issues = utils_common_constraints_on_dataframe(df_dge, constraints)
        
        # Check for any issues
        has_issues = any(len(issue_cols) > 0 for issue_cols in issues.values())
        
        if not has_issues:
            status = "GREEN"
            message = "All comparison statistical columns meet constraints"
            details = "Log2fc and Stat columns have no null values; P.value and Adj.p.value columns have no negative values"
            log_check_result(log_path, component_name, "all", check_name, status, message, details)
            return True
        else:
            status = "RED"
            message = "Comparison statistical columns fail constraints"
            
            # Format the details with issues
            details_parts = []
            for issue_type, columns in issues.items():
                if columns:
                    details_parts.append(f"{issue_type}: {', '.join(columns[:5])}")
                    if len(columns) > 5:
                        details_parts[-1] += f" and {len(columns) - 5} more"
            
            details = "; ".join(details_parts)
            log_check_result(log_path, component_name, "all", check_name, status, message, details)
            return False
            
    except Exception as e:
        log_check_result(log_path, component_name, "all", check_name, "RED", 
                        "Error checking comparison statistical columns constraints", str(e))
        return False

def check_dge_table_fixed_statistical_columns_exist(outdir, log_path, assay_suffix="_GLbulkRNAseq",
                                           stratum_factor="", stratum_value=""):
    """Check if the fixed statistical columns exist in the DGE table.
    
    These columns include dataset-level statistics like All.mean, All.stdev, and LRT.p.value.
    
    Args:
        outdir: Path to the DGE directory (not the parent outdir)
        log_path: Path to the log file
        assay_suffix: Suffix to append to the file name
        stratum_factor: Factor used for stratification (if any)
        stratum_value: Value of the stratum to check (if any)
        
    Returns:
        bool: True if check passes, False otherwise
    """
    component = "dge"
    check_name = "check_dge_table_fixed_statistical_columns_exist"
    
    # Define the expected fixed statistical columns
    expected_columns = {"All.mean", "All.stdev", "LRT.p.value"}
    
    # Adjust the suffix and component name if stratified
    if stratum_value:
        check_suffix = f"_{stratum_value}"
        component_name = f"{component}{check_suffix}"
    else:
        check_suffix = ""
        component_name = component
    
    # Get DGE table path with the correct pattern
    dge_table_path = os.path.join(outdir, f"differential_expression{check_suffix}{assay_suffix}.csv")
    
    # Check if the DGE table exists
    if not os.path.exists(dge_table_path):
        message = f"DGE table not found"
        log_check_result(log_path, component_name, "all", check_name, "RED", 
                       message, f"Expected at: {dge_table_path}")
        return False
    
    try:
        # Read DGE table and check for expected columns
        df_dge = pd.read_csv(dge_table_path)
        df_dge_columns = set(df_dge.columns)
        
        # Find missing columns
        missing_columns = expected_columns - df_dge_columns
        
        if not missing_columns:
            status = "GREEN"
            message = "All dataset summary statistic columns present in DGE table"
            log_check_result(log_path, component_name, "all", check_name, status, message, 
                            f"Found all dataset summary columns: {'; '.join(sorted(expected_columns))}")
            return True
        else:
            status = "RED"
            message = "Missing dataset summary statistic columns in DGE table"
            details = f"Missing these dataset summary columns: {', '.join(sorted(missing_columns))}"
            log_check_result(log_path, component_name, "all", check_name, status, message, details)
            return False
        
    except Exception as e:
        log_check_result(log_path, component_name, "all", check_name, "RED", 
                        "Error checking fixed statistical columns", str(e))
        return False

def check_dge_table_fixed_statistical_columns_constraints(outdir, log_path, assay_suffix="_GLbulkRNAseq",
                                                 stratum_factor="", stratum_value=""):
    """Check if the fixed statistical columns in the DGE table meet specific constraints.
    
    This validates that:
    - All.mean and All.stdev columns have no null values and no negative values
    - LRT.p.value column has no negative values (can have nulls)
    
    Args:
        outdir: Path to the DGE directory (not the parent outdir)
        log_path: Path to the log file
        assay_suffix: Suffix to append to the file name
        stratum_factor: Factor used for stratification (if any)
        stratum_value: Value of the stratum to check (if any)
        
    Returns:
        bool: True if check passes, False otherwise
    """
    component = "dge"
    check_name = "check_dge_table_fixed_statistical_columns_constraints"
    
    # Adjust the suffix and component name if stratified
    if stratum_value:
        check_suffix = f"_{stratum_value}"
        component_name = f"{component}{check_suffix}"
    else:
        check_suffix = ""
        component_name = component
    
    # Get DGE table path with the correct pattern
    dge_table_path = os.path.join(outdir, f"differential_expression{check_suffix}{assay_suffix}.csv")
    
    # Check if the DGE table exists
    if not os.path.exists(dge_table_path):
        message = f"DGE table not found"
        log_check_result(log_path, component_name, "all", check_name, "RED", 
                       message, f"Expected at: {dge_table_path}")
        return False
    
    try:
        # Define constraints for different types of columns
        constraints = [
            ({"All.mean", "All.stdev"}, {"nonNull": True, "nonNegative": True}),
            ({"LRT.p.value"}, {"nonNegative": True}),
        ]
        
        # Read DGE table
        df_dge = pd.read_csv(dge_table_path)
        
        # Apply constraints
        issues = utils_common_constraints_on_dataframe(df_dge, constraints)
        
        # Check for any issues
        has_issues = any(len(issue_cols) > 0 for issue_cols in issues.values())
        
        if not has_issues:
            status = "GREEN"
            message = "All dataset summary statistic columns meet constraints"
            details = "All.mean and All.stdev columns have no null or negative values; LRT.p.value column has no negative values"
            log_check_result(log_path, component_name, "all", check_name, status, message, details)
            return True
        else:
            status = "RED"
            message = "Dataset summary statistic columns fail constraints"
            
            # Format the details with issues
            details_parts = []
            for issue_type, columns in issues.items():
                if columns:
                    details_parts.append(f"{issue_type}: {', '.join(columns[:5])}")
                    if len(columns) > 5:
                        details_parts[-1] += f" and {len(columns) - 5} more"
            
            details = "; ".join(details_parts)
            log_check_result(log_path, component_name, "all", check_name, status, message, details)
            return False
            
    except Exception as e:
        log_check_result(log_path, component_name, "all", check_name, "RED", 
                        "Error checking fixed statistical columns constraints", str(e))
        return False

def check_dge_table_log2fc_within_reason(outdir, runsheet_path, log_path, assay_suffix="_GLbulkRNAseq",
                                       stratum_factor="", stratum_value=""):
    """Check if the log2fc values in the DGE table are within reasonable bounds compared to the group means.
    
    This function performs two key checks:
    1. Validates that the log2fc values computed by DESeq2 are reasonably close to directly computing log2fc
       from the group means (within a specified tolerance)
    2. Checks that the sign (positive/negative) of log2fc values matches the expected sign based on
       group mean comparisons for genes with substantial expression differences
    
    Args:
        outdir: Path to the DGE directory (not the parent outdir)
        runsheet_path: Path to the runsheet CSV file
        log_path: Path to the log file
        assay_suffix: Suffix to append to the file name
        stratum_factor: Factor used for stratification (if any)
        stratum_value: Value of the stratum to check (if any)
        
    Returns:
        bool: True if check passes, False otherwise
    """
    component = "dge"
    check_name = "check_dge_table_log2fc_within_reason"
    
    # Define thresholds for log2fc validation
    LOG2FC_CROSS_METHOD_PERCENT_DIFFERENCE_THRESHOLD = 10  # Percent
    LOG2FC_CROSS_METHOD_TOLERANCE_PERCENT = 50  # Percent
    THRESHOLD_PERCENT_MEANS_DIFFERENCE = 50  # percent
    SMALL_COUNTS_THRESHOLD = 20  # Minimum sum of Group.Mean values to consider
    
    # Adjust the suffix and component name if stratified
    if stratum_value:
        check_suffix = f"_{stratum_value}"
        component_name = f"{component}{check_suffix}"
    else:
        check_suffix = ""
        component_name = component
    
    # Get DGE table path with the correct pattern
    dge_table_path = os.path.join(outdir, f"differential_expression{check_suffix}{assay_suffix}.csv")
    
    # Check if the DGE table exists
    if not os.path.exists(dge_table_path):
        message = f"DGE table not found"
        log_check_result(log_path, component_name, "all", check_name, "RED", 
                      message, f"Expected at: {dge_table_path}")
        return False
    
    try:
        # Read DGE table to get actual column names
        df_dge = pd.read_csv(dge_table_path)
        original_count = len(df_dge)
        
        # Filter out low-count genes
        group_mean_cols = [col for col in df_dge.columns if col.startswith("Group.Mean_")]
        if group_mean_cols:
            # Sum of all Group.Mean columns
            df_dge["Group.Mean_SUM"] = df_dge[group_mean_cols].sum(axis=1)
            # Filter out genes with low counts
            df_dge_filtered = df_dge[df_dge["Group.Mean_SUM"] > SMALL_COUNTS_THRESHOLD]
            filtered_count = original_count - len(df_dge_filtered)
            filtering_info = f"Filtered out {filtered_count}/{original_count} ({filtered_count/original_count*100:.1f}%) genes with sum of Group.Mean values below {SMALL_COUNTS_THRESHOLD}."
            print(filtering_info)
            df_dge = df_dge_filtered
        else:
            filtering_info = "No Group.Mean columns found for filtering."
        
        # Find Log2fc columns that exist in the data
        log2fc_columns = [col for col in df_dge.columns if col.startswith("Log2fc_")]
        
        if not log2fc_columns:
            log_check_result(log_path, component_name, "all", check_name, "RED",
                           "No Log2fc columns found in DGE table", 
                           f"Available columns: {', '.join(df_dge.columns[:10])}...")
            return False
        
        # Extract comparisons from column names
        comparisons = [col[len("Log2fc_"):] for col in log2fc_columns]
        
        # Track issues
        err_msg_yellow = ""
        all_suspect_signs = {}
        last_percent_within_tolerance = 0
        
        # Process each comparison
        for comparison in comparisons:
            query_column = f"Log2fc_{comparison}"
            
            # Extract group names based on the comparison format
            try:
                group1_name = comparison.split(")v(")[0] + ")"
                group2_name = "(" + comparison.split(")v(")[1]
            except:
                try:
                    # Try simple "v" separator without parentheses
                    group1_name = comparison.split("v")[0]
                    group2_name = comparison.split("v")[1]
                except:
                    # Skip malformed comparison
                    continue
            
            # Get group mean columns
            group1_mean_col = f"Group.Mean_{group1_name}"
            group2_mean_col = f"Group.Mean_{group2_name}"
            
            # Skip if mean columns don't exist
            if group1_mean_col not in df_dge.columns or group2_mean_col not in df_dge.columns:
                continue
            
            # Compute log2fc directly from group means
            safe_denom = df_dge[group2_mean_col].replace(0, np.nan)  # Avoid division by zero
            ratio = df_dge[group1_mean_col] / safe_denom
            computed_log2fc = np.log2(ratio)
            
            # Calculate percent difference between methods
            safe_query_values = df_dge[query_column].replace(0, np.nan)  # Avoid division by zero
            abs_percent_difference = abs(((computed_log2fc - df_dge[query_column]) / safe_query_values) * 100)
            
            # Check what percentage is within tolerance (ignoring NaN values)
            valid_mask = ~abs_percent_difference.isna()
            if valid_mask.sum() > 0:
                within_tolerance_mask = abs_percent_difference < LOG2FC_CROSS_METHOD_PERCENT_DIFFERENCE_THRESHOLD
                percent_within_tolerance = (within_tolerance_mask[valid_mask].sum() / valid_mask.sum()) * 100
                last_percent_within_tolerance = percent_within_tolerance
                
                # Flag if not enough within tolerance
                if percent_within_tolerance < LOG2FC_CROSS_METHOD_TOLERANCE_PERCENT:
                    err_msg_yellow += (
                        f"For comparison '{comparison}': {percent_within_tolerance:.2f}% of genes have absolute percent differences "
                        f"(between log2fc direct computation and DESeq2's approach) "
                        f"less than {LOG2FC_CROSS_METHOD_PERCENT_DIFFERENCE_THRESHOLD}% which is below the minimum required "
                        f"({LOG2FC_CROSS_METHOD_TOLERANCE_PERCENT}%). "
                        f"This may indicate misassigned or misaligned columns. "
                    )
            
            # SIGN CHECK: Filter genes with substantial mean differences
            abs_mean_diffs = abs((df_dge[group1_mean_col] - df_dge[group2_mean_col]) / safe_denom) * 100
            mask = (abs_mean_diffs > THRESHOLD_PERCENT_MEANS_DIFFERENCE) & (~abs_mean_diffs.isna())
            
            if mask.sum() > 0:  # Only proceed if we have genes meeting criteria
                # Calculate expected sign
                positive_sign_expected = (df_dge[group1_mean_col] - df_dge[group2_mean_col])[mask] > 0
                
                # Check if actual log2fc matches expected sign
                actual_sign_positive = df_dge[query_column][mask] > 0
                matches_expected = (actual_sign_positive & positive_sign_expected) | (~actual_sign_positive & ~positive_sign_expected)
                
                # Find genes with wrong signs
                wrong_sign_mask = ~matches_expected.values
                if wrong_sign_mask.any():
                    wrong_sign_genes = df_dge[mask][wrong_sign_mask]
                    
                    # Get gene identifiers - use first column as gene ID
                    gene_id_col = df_dge.columns[0]
                    
                    # Get examples (max 5)
                    examples = []
                    for _, row in wrong_sign_genes.head(5).iterrows():
                        gene_id = row[gene_id_col]
                        examples.append(f"{gene_id}")
                    
                    # Store information about wrong sign genes
                    all_suspect_signs[comparison] = {
                        "count": wrong_sign_mask.sum(),
                        "total": mask.sum(),
                        "examples": examples
                    }
        
        # Determine status based on findings
        if all_suspect_signs:
            status = "YELLOW"  # Changed from RED to YELLOW
            message = "Log2fc signs do not match expected direction based on group means"
            
            # Collect all unique genes with wrong signs across all contrasts
            all_wrong_sign_genes = set()
            total_genes_checked = 0
            total_wrong_signs = 0
            wrong_sign_gene_ids = []
            
            # Track contrast-specific information for context
            contrasts_info = []
            
            # Process each comparison and collect consolidated info
            for comparison, info in all_suspect_signs.items():
                all_wrong_sign_genes.update(info["examples"])
                total_genes_checked += info["total"]
                total_wrong_signs += info["count"]
                
                # Get all gene IDs from this comparison (not just examples)
                # Need to fetch full list by re-querying the data
                wrong_sign_mask = ~matches_expected.values if 'matches_expected' in locals() else None
                if wrong_sign_mask is not None and any(wrong_sign_mask):
                    gene_id_col = df_dge.columns[0]
                    wrong_sign_gene_ids.extend(df_dge.loc[mask][wrong_sign_mask][gene_id_col].tolist())
            
            # Check how many of these problematic genes had low counts
            # Create a complete DGE table with the Group.Mean_SUM column
            if "Group.Mean_SUM" not in df_dge.columns:
                df_dge_with_sum = pd.read_csv(dge_table_path)
                df_dge_with_sum["Group.Mean_SUM"] = df_dge_with_sum[group_mean_cols].sum(axis=1)
            else:
                df_dge_with_sum = df_dge
                
            # Count low-expression genes among the problematic ones
            gene_id_col = df_dge.columns[0]
            low_count_wrong_sign_genes = 0
            
            # If we have gene IDs, check them
            if wrong_sign_gene_ids:
                low_count_mask = df_dge_with_sum[df_dge_with_sum[gene_id_col].isin(wrong_sign_gene_ids)]["Group.Mean_SUM"] <= SMALL_COUNTS_THRESHOLD
                low_count_wrong_sign_genes = low_count_mask.sum()
            
            # Limit to 10 gene examples maximum if there are more
            gene_examples = list(all_wrong_sign_genes)[:10]
            if len(all_wrong_sign_genes) > 10:
                gene_examples.append("...")
            
            # Calculate overall percentage
            overall_percent = (total_wrong_signs / total_genes_checked) * 100 if total_genes_checked > 0 else 0
            
            # Create consolidated details message
            details = f"Found {total_wrong_signs} genes with incorrect log2fc signs across all comparisons ({total_wrong_signs}/{total_genes_checked} genes). Genes: {'; '.join(gene_examples)}. {low_count_wrong_sign_genes}/{total_wrong_signs} genes with incorrect signs had Group.Mean values below {SMALL_COUNTS_THRESHOLD}."
            
            log_check_result(log_path, component_name, "all", check_name, status, message, details)
            return False
            
        elif err_msg_yellow:
            status = "YELLOW"
            message = "Log2fc values show significant differences from direct computation"
            details = f"Log2fc values differ significantly from direct calculation."
            log_check_result(log_path, component_name, "all", check_name, status, message, details)
            return True  # Yellow status still passes
            
        else:
            status = "GREEN"
            message = "All log2fc values are within reasonable bounds"
            details = f"Log2fc values match direct calculation and all signs match expected direction."
            log_check_result(log_path, component_name, "all", check_name, status, message, details)
            return True
            
    except Exception as e:
        import traceback
        error_details = traceback.format_exc()
        log_check_result(log_path, component_name, "all", check_name, "RED", 
                       "Error checking log2fc reasonableness", f"{str(e)}\n{error_details}")
        return False

def check_ercc_presence(outdir, runsheet_path, log_path, assay_suffix="_GLbulkRNAseq", mode="default"):
    """Check for the presence/absence of ERCC spike-ins in count matrices.
    
    If the runsheet indicates has_ercc=True, this function verifies:
    1. ERCC spike-ins are present in the unnormalized counts
    2. Whether they are present/removed in normalized counts
    
    Args:
        outdir: Path to the output directory
        runsheet_path: Path to the runsheet CSV
        log_path: Path to the log file
        assay_suffix: Suffix for the assay files
        mode: Processing mode, either 'default' or 'microbes'
        
    Returns:
        bool: True if check passes, False otherwise
    """
    component = "dge_NormCounts"
    check_name = "check_ercc_presence"
    
    try:
        # Read runsheet to check for has_ercc flag
        runsheet_df = pd.read_csv(runsheet_path)
        
        # Find has_ercc column - could be capitalized in different ways
        ercc_col = None
        for col in runsheet_df.columns:
            if col.lower() == 'has_ercc':
                ercc_col = col
                break
        
        # If no has_ercc column or it's false, this check is not applicable
        if not ercc_col or not runsheet_df[ercc_col].iloc[0]:
            log_check_result(log_path, component, "all", check_name, "GREEN", 
                           "ERCC spike-in check skipped", "No has_ercc=True found in runsheet")
            return True
        
        # Set the unnormalized counts file based on mode
        if mode == "microbes":
            unnorm_path = os.path.join(outdir, f"FeatureCounts_Unnormalized_Counts{assay_suffix}.csv")
        else:
            unnorm_path = os.path.join(outdir, f"RSEM_Unnormalized_Counts{assay_suffix}.csv")
            
        norm_path = os.path.join(outdir, f"Normalized_Counts{assay_suffix}.csv")
        vst_norm_path = os.path.join(outdir, f"VST_Normalized_Counts{assay_suffix}.csv")
        
        # Check if files exist
        files_missing = []
        if not os.path.exists(unnorm_path):
            files_missing.append(f"Unnormalized counts file: {unnorm_path}")
        if not os.path.exists(norm_path):
            files_missing.append(f"Normalized counts file: {norm_path}")
        
        if files_missing:
            log_check_result(log_path, component, "all", check_name, "RED", 
                           "Required count files not found", "; ".join(files_missing))
            return False
        
        # Count ERCC in unnormalized counts
        unnorm_df = pd.read_csv(unnorm_path)
        gene_id_col = unnorm_df.columns[0]  # First column should be gene IDs
        ercc_unnorm = sum(unnorm_df[gene_id_col].str.startswith('ERCC-'))
        
        # Count ERCC in normalized counts
        norm_df = pd.read_csv(norm_path)
        gene_id_col_norm = norm_df.columns[0]
        ercc_norm = sum(norm_df[gene_id_col_norm].str.startswith('ERCC-'))
        
        # Count ERCC in VST normalized counts (if exists)
        ercc_vst = 0
        if os.path.exists(vst_norm_path):
            vst_df = pd.read_csv(vst_norm_path)
            gene_id_col_vst = vst_df.columns[0]
            ercc_vst = sum(vst_df[gene_id_col_vst].str.startswith('ERCC-'))
        
        # Check results and build message
        if ercc_unnorm == 0:
            status = "RED"
            message = "ERCC spike-ins missing from unnormalized counts despite has_ercc=True"
            details = f"Runsheet indicates has_ercc={runsheet_df[ercc_col].iloc[0]}, but no ERCC IDs found in unnormalized counts"
            log_check_result(log_path, component, "all", check_name, status, message, details)
            return False
        
        # Determine if ERCCs were processed as expected
        status = "GREEN"
        message = "ERCC spike-ins found in appropriate files"
        
        details_parts = [
            f"Found {ercc_unnorm} ERCC spike-ins in unnormalized counts",
            f"Found {ercc_norm} ERCC spike-ins in normalized counts"
        ]
        
        if os.path.exists(vst_norm_path):
            details_parts.append(f"Found {ercc_vst} ERCC spike-ins in VST normalized counts")
        
        # Check if ERCCs were removed during normalization
        if ercc_unnorm > 0 and ercc_norm == 0:
            details_parts.append("ERCC spike-ins were correctly removed during normalization")
        elif ercc_unnorm > 0 and ercc_norm > 0:
            details_parts.append("ERCC spike-ins were retained during normalization")
            
        details = "; ".join(details_parts)
        log_check_result(log_path, component, "all", check_name, status, message, details)
        return True
        
    except Exception as e:
        import traceback
        error_details = traceback.format_exc()
        log_check_result(log_path, component, "all", check_name, "RED", 
                       f"Error checking ERCC spike-ins", f"{str(e)}\n{error_details}")
        return False
                       
    
def main():
    parser = argparse.ArgumentParser(description="Verify and validate DESeq2 normalization and DGE outputs")
    parser.add_argument("--outdir", required=True, help="Path to the output directory")
    parser.add_argument("--runsheet", required=True, help="Path to the runsheet CSV file")
    parser.add_argument("--assay_suffix", default="_GLbulkRNAseq", help="Assay suffix")
    parser.add_argument("--stratify_by", default="", help="Factor to stratify analysis by (e.g., 'Plant Part')")
    parser.add_argument("--mode", default="default", choices=["default", "microbes"], 
                       help="Processing mode: 'default' for RSEM or 'microbes' for FeatureCounts")
    args = parser.parse_args()
    
    # Initialize the VV log
    vv_log_path = initialize_vv_log(args.outdir)
    
    # Get stratified paths if needed
    stratified_paths = get_factor_stratified_paths(args.outdir, args.runsheet, args.stratify_by)
    
    # Track overall check results
    all_checks_passed = True
    check_results = {}
    
    # Run checks for each stratum
    for stratum, paths in stratified_paths.items():
        print(f"\n{'=' * 50}")
        if stratum:
            print(f"Running checks for stratum: {stratum}")
            final_suffix = f"{args.assay_suffix}_{stratum}"
            
            # Extract factor name and value from stratum
            factor_name = args.stratify_by
            factor_value = stratum
        else:
            print("Running checks for standard analysis (no stratification)")
            final_suffix = args.assay_suffix
            factor_name = ""
            factor_value = ""
        
        # Adjust paths based on stratum
        norm_counts_dir = paths["norm_counts"]
        dge_dir = paths["dge"]
        
        print(f"Checking norm counts directory: {norm_counts_dir}")
        print(f"Checking DGE directory: {dge_dir}")
        
        # Initialize check results for this stratum
        norm_counts_passed = False
        dge_passed = False
        sample_table_passed = False
        group_assignments_passed = False
        contrasts_headers_passed = False
        contrasts_rows_passed = False
        annotation_columns_passed = False
        sample_columns_passed = False
        sample_columns_constraints_passed = False
        
        # Check DESeq2 normalized counts files
        print("\nChecking DESeq2 normalized counts files...")
        try:
            norm_counts_passed = check_deseq2_normcounts_existence(norm_counts_dir, vv_log_path, final_suffix, args.mode)
        except Exception as e:
            print(f"Error during normalized counts check: {str(e)}")
            log_check_result(vv_log_path, "dge_NormCounts", "all", "check_deseq2_normcounts_existence", 
                           "RED", f"Exception during check", str(e))
        
        # Check ERCC spike-ins if applicable
        print("\nChecking ERCC spike-ins...")
        try:
            ercc_passed = check_ercc_presence(norm_counts_dir, args.runsheet, vv_log_path, final_suffix, args.mode)
        except Exception as e:
            print(f"Error during ERCC spike-in check: {str(e)}")
            log_check_result(vv_log_path, "dge_NormCounts", "all", "check_ercc_presence", 
                           "RED", f"Exception during check", str(e))
            ercc_passed = False
        
        # Check DESeq2 DGE files
        print("\nChecking DESeq2 DGE files...")
        try:
            dge_passed = check_deseq2_dge_existence(dge_dir, vv_log_path, final_suffix)
        except Exception as e:
            print(f"Error during DGE files check: {str(e)}")
            log_check_result(vv_log_path, "dge_DGE", "all", "check_deseq2_dge_existence", 
                           "RED", f"Exception during check", str(e))
        
        # Check sample table against runsheet
        print("\nChecking sample table against runsheet...")
        try:
            sample_table_passed = check_sample_table_against_runsheet(
                dge_dir, args.runsheet, vv_log_path, final_suffix, factor_name, factor_value
            )
        except Exception as e:
            print(f"Error during sample table check: {str(e)}")
            log_check_result(vv_log_path, "dge", "all", "check_sample_table_against_runsheet", 
                           "RED", f"Exception during check", str(e))
        
        # Check sample table for correct group assignments
        print("\nChecking sample table for correct group assignments...")
        try:
            group_assignments_passed = check_sample_table_for_correct_group_assignments(
                dge_dir, args.runsheet, vv_log_path, final_suffix, factor_name, factor_value
            )
        except Exception as e:
            print(f"Error during group assignments check: {str(e)}")
            log_check_result(vv_log_path, "dge", "all", "check_sample_table_for_correct_group_assignments", 
                           "RED", f"Exception during check", str(e))
        
        # Check contrasts table headers
        print("\nChecking contrasts table headers...")
        try:
            contrasts_headers_passed = check_contrasts_table_headers(
                dge_dir, args.runsheet, vv_log_path, final_suffix, factor_name, factor_value
            )
        except Exception as e:
            print(f"Error during contrasts table headers check: {str(e)}")
            log_check_result(vv_log_path, "dge", "all", "check_contrasts_table_headers", 
                           "RED", f"Exception during check", str(e))
        
        # Check contrasts table rows
        print("\nChecking contrasts table rows...")
        try:
            contrasts_rows_passed = check_contrasts_table_rows(
                dge_dir, vv_log_path, final_suffix, factor_name, factor_value
            )
        except Exception as e:
            print(f"Error during contrasts table rows check: {str(e)}")
            log_check_result(vv_log_path, "dge", "all", "check_contrasts_table_rows", 
                           "RED", f"Exception during check", str(e))
        
        # Check DGE table annotation columns
        print("\nChecking DGE table annotation columns...")
        try:
            annotation_columns_passed = check_dge_table_annotation_columns_exist(
                dge_dir, args.runsheet, vv_log_path, final_suffix, factor_name, factor_value
            )
        except Exception as e:
            print(f"Error during DGE table annotation columns check: {str(e)}")
            log_check_result(vv_log_path, "dge", "all", "check_dge_table_annotation_columns_exist", 
                           "RED", f"Exception during check", str(e))
        
        # Check DGE table sample columns existence
        print("\nChecking DGE table sample columns existence...")
        try:
            sample_columns_passed = check_dge_table_sample_columns_exist(
                dge_dir, args.runsheet, vv_log_path, final_suffix, factor_name, factor_value
            )
        except Exception as e:
            print(f"Error during DGE table sample columns check: {str(e)}")
            log_check_result(vv_log_path, "dge", "all", "check_dge_table_sample_columns_exist", 
                           "RED", f"Exception during check", str(e))
        
        # Check DGE table sample columns constraints
        print("\nChecking DGE table sample columns constraints...")
        try:
            sample_columns_constraints_passed = check_dge_table_sample_columns_constraints(
                dge_dir, args.runsheet, vv_log_path, final_suffix, factor_name, factor_value
            )
        except Exception as e:
            print(f"Error during DGE table sample columns constraints check: {str(e)}")
            log_check_result(vv_log_path, "dge", "all", "check_dge_table_sample_columns_constraints", 
                           "RED", f"Exception during check", str(e))
        
        # Check DGE table group columns existence
        print("\nChecking DGE table group columns existence...")
        try:
            group_columns_passed = check_dge_table_group_columns_exist(
                dge_dir, args.runsheet, vv_log_path, final_suffix, factor_name, factor_value
            )
        except Exception as e:
            print(f"Error during DGE table group columns check: {str(e)}")
            log_check_result(vv_log_path, "dge", "all", "check_dge_table_group_columns_exist", 
                           "RED", f"Exception during check", str(e))
        
        # Check DGE table group columns constraints
        print("\nChecking DGE table group columns constraints...")
        try:
            group_columns_constraints_passed = check_dge_table_group_columns_constraints(
                dge_dir, args.runsheet, vv_log_path, final_suffix, factor_name, factor_value
            )
        except Exception as e:
            print(f"Error during DGE table group columns constraints check: {str(e)}")
            log_check_result(vv_log_path, "dge", "all", "check_dge_table_group_columns_constraints", 
                           "RED", f"Exception during check", str(e))
        
        # Check DGE table comparison statistical columns
        print("\nChecking DGE table comparison statistical columns...")
        try:
            comparison_columns_passed = check_dge_table_comparison_statistical_columns_exist(
                dge_dir, args.runsheet, vv_log_path, final_suffix, factor_name, factor_value
            )
        except Exception as e:
            print(f"Error during DGE table comparison statistical columns check: {str(e)}")
            log_check_result(vv_log_path, "dge", "all", "check_dge_table_comparison_statistical_columns_exist", 
                           "RED", f"Exception during check", str(e))
        
        # Check DGE table group statistical columns constraints
        print("\nChecking DGE table group statistical columns constraints...")
        try:
            group_statistical_constraints_passed = check_dge_table_group_statistical_columns_constraints(
                dge_dir, args.runsheet, vv_log_path, final_suffix, factor_name, factor_value
            )
        except Exception as e:
            print(f"Error during DGE table group statistical columns constraints check: {str(e)}")
            log_check_result(vv_log_path, "dge", "all", "check_dge_table_group_statistical_columns_constraints", 
                           "RED", f"Exception during check", str(e))
        
        # Check DGE table fixed statistical columns exist
        print("\nChecking DGE table fixed statistical columns...")
        try:
            fixed_statistical_columns_passed = check_dge_table_fixed_statistical_columns_exist(
                dge_dir, vv_log_path, final_suffix, factor_name, factor_value
            )
        except Exception as e:
            print(f"Error during DGE table fixed statistical columns check: {str(e)}")
            log_check_result(vv_log_path, "dge", "all", "check_dge_table_fixed_statistical_columns_exist", 
                           "RED", f"Exception during check", str(e))
        
        # Check DGE table fixed statistical columns constraints
        print("\nChecking DGE table fixed statistical columns constraints...")
        try:
            fixed_statistical_constraints_passed = check_dge_table_fixed_statistical_columns_constraints(
                dge_dir, vv_log_path, final_suffix, factor_name, factor_value
            )
        except Exception as e:
            print(f"Error during DGE table fixed statistical columns constraints check: {str(e)}")
            log_check_result(vv_log_path, "dge", "all", "check_dge_table_fixed_statistical_columns_constraints", 
                           "RED", f"Exception during check", str(e))
        
        # Check DGE table log2fc within reason
        print("\nChecking DGE table log2fc within reason...")
        try:
            log2fc_within_reason_passed = check_dge_table_log2fc_within_reason(
                dge_dir, args.runsheet, vv_log_path, final_suffix, factor_name, factor_value
            )
        except Exception as e:
            print(f"Error during DGE table log2fc within reason check: {str(e)}")
            log_check_result(vv_log_path, "dge", "all", "check_dge_table_log2fc_within_reason", 
                           "RED", f"Exception during check", str(e))
        
        # Update overall status
        stratum_passed = all([
            norm_counts_passed, 
            ercc_passed,
            dge_passed, 
            sample_table_passed, 
            group_assignments_passed,
            contrasts_headers_passed,
            contrasts_rows_passed,
            annotation_columns_passed,
            sample_columns_passed,
            sample_columns_constraints_passed,
            group_columns_passed,
            group_columns_constraints_passed,
            comparison_columns_passed,
            group_statistical_constraints_passed,
            fixed_statistical_columns_passed,
            fixed_statistical_constraints_passed,
            log2fc_within_reason_passed
        ])
        
        # Store results for each check
        component_id = f"dge_{stratum}" if stratum else "dge"
        check_results[component_id] = {
            "normcounts_check": norm_counts_passed,
            "ercc_check": ercc_passed,
            "dge_check": dge_passed,
            "sample_table_check": sample_table_passed,
            "group_assignments_check": group_assignments_passed,
            "contrasts_headers_check": contrasts_headers_passed,
            "contrasts_rows_check": contrasts_rows_passed,
            "annotation_columns_check": annotation_columns_passed,
            "sample_columns_check": sample_columns_passed,
            "sample_columns_constraints_check": sample_columns_constraints_passed,
            "group_columns_check": group_columns_passed,
            "group_columns_constraints_check": group_columns_constraints_passed,
            "comparison_columns_check": comparison_columns_passed,
            "group_statistical_constraints_check": group_statistical_constraints_passed,
            "fixed_statistical_columns_check": fixed_statistical_columns_passed,
            "fixed_statistical_constraints_check": fixed_statistical_constraints_passed,
            "log2fc_within_reason_check": log2fc_within_reason_passed,
            "overall": stratum_passed
        }
        
        all_checks_passed = all_checks_passed and stratum_passed
        
        print(f"\nStatus for {stratum if stratum else 'standard analysis'}: {'PASSED' if stratum_passed else 'FAILED'}")
    
    # Determine overall status
    overall_status = "GREEN" if all_checks_passed else "RED"
    
    # Print the summary
    print_summary(check_results, vv_log_path, overall_status)
    
    print(f"\n{'=' * 50}")
    print(f"Overall DGE validation status: {'PASSED' if all_checks_passed else 'FAILED'}")
    
    # Always return success (0) regardless of check results
    return 0

if __name__ == "__main__":
    sys.exit(main()) 