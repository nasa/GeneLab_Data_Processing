""" Schemas for validation 
Uses Schema to allow usage of validation functions
"""
import pandera as pa

check_single_value = pa.Check(
    lambda x: len(x.unique()) == 1,
    title="Check that all values in the column are identical",
    description="Useful for columns that contain dataset level metadata like organism and paired_end.",
    error="Dataset level columns do NOT contain one unique value"
    )

runsheet = pa.DataFrameSchema(
        columns={
            "Original Sample Name": pa.Column(str),
            "Study Assay Measurement": pa.Column(str),
            "Study Assay Technology Type": pa.Column(str),
            "Study Assay Technology Platform": pa.Column(str),
            "Source Name": pa.Column(str),
            "Label": pa.Column(str),
            "Array Data File Name": pa.Column(str),
            "Array Data File Path": pa.Column(str),
            "Original Sample Name": pa.Column(str),
            "organism": pa.Column(str, check_single_value),
        },
        # define checks at the DataFrameSchema-level
        # checks=check_read2_path_populated_if_paired_end
    )
