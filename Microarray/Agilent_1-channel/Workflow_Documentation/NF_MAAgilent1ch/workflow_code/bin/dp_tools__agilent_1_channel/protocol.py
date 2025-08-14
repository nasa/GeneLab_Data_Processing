from pathlib import Path
import re
from typing import Union
import yaml
import logging
import textwrap

from dp_tools.core.entity_model import Dataset

log = logging.getLogger(__name__)

from dp_tools.core.check_model import ValidationProtocol, FlagCode
# ORDER is intentional to allow shared name functions to overload in favor of microarray
from dp_tools import bulkRNASeq

from . import checks

CONFIG = {
    "Metadata-check_metadata_attributes_exist": {
        "expected_attrs": ["biomart_attribute", "organism"]
    },
}


def validate(
    dataset: Dataset,
    config_path: Path = None,
    run_args: dict = None,
    report_args: dict = None,
    protocol_args: dict = None,
    defer_run: bool = False,
) -> Union[ValidationProtocol, ValidationProtocol.Report]:

    if config_path is not None:
        with open(config_path, "r") as f:
            config = yaml.safe_load(f)
    else:
        config = CONFIG

    if run_args is None:
        run_args = dict()

    if report_args is None:
        report_args = dict()

    if protocol_args is None:
        protocol_args = dict()
    # init validation protocol
    vp = ValidationProtocol(**protocol_args)
    # fmt: on
    with vp.component_start(
        name=dataset.name,
        description="Validate microarray processed data",
    ):
        with vp.component_start(
            name="Metadata", description="Metadata file validation"
        ):
            with vp.payload(payloads=[{"dataset": dataset}]):
                vp.add(
                    checks.check_metadata_attributes_exist,
                    config=config["Metadata-check_metadata_attributes_exist"],
                    full_description=textwrap.dedent(f"""
                        - Check: Runsheet includes required metadata columns, {config["Metadata-check_metadata_attributes_exist"]["expected_attrs"]}
                            - Reason:
                                - 'Organism' column is used to map the appropriate annotation table
                                - 'biomart_attribute' column is used to map the appropriate biomart attribute for initial Ensembl based gene mapping
                            - Potential Source of Problems:
                                - Runsheet does not include these columns. If automatically generated using 'dp_tools', report this to the maintainer of 'dp_tools'.  If manually generated, ensure {config["Metadata-check_metadata_attributes_exist"]["expected_attrs"]} columns are populated.
                            - Flag Condition:
                                - Green: Columns {config["Metadata-check_metadata_attributes_exist"]["expected_attrs"]} exist
                                - Halt: Columns {config["Metadata-check_metadata_attributes_exist"]["expected_attrs"]} do not exist
                    """)
                )
                vp.add_manual(
                    description = "Manually validate runsheet against ISA assay table",
                    start_instructions = "Open runsheet (in Metadata folder). Open OSD webpage, navigate to microarray assay table",
                    pass_fail_questions = [
                        "Does the runsheet open?",
                        "Is the number of samples in parity?",
                        "Do the file extensions in the 'Array Data File Name' column end in '.txt' or '.txt.gz'?",
                        "Do all factor values exist in both tables. Does the runsheet have units included as appropriate?",
                    ]
                )
        with vp.component_start(
            name="Raw Intensities",
            description="",
        ):
            with vp.payload(
                payloads=[
                    {
                        "raw_intensities": lambda: dataset.data_assets["raw intensities table"].path,
                        "samples": dataset.samples
                    }
                ]
            ):
                vp.add(
                    checks.check_raw_intensities_table,
                    full_description=textwrap.dedent(f"""
                        - Check: Ensure raw intensities table has all samples and values within [0,+inf)
                            - Reason:
                                - Part of processing output
                            - Potential Source of Problems:
                                - Bug in processing script or malformed raw data files
                            - Flag Condition:
                                - Green: All conditions met
                                - Halt: At least one condition failed
                    """)
                    )
            with vp.payload(
                payloads=[
                    {
                        "organism": lambda: dataset.metadata["organism"],
                        "samples": lambda: set(dataset.samples),
                        "dge_table": lambda: dataset.data_assets[
                            "raw intensities table"
                        ].path,
                        "runsheet": lambda: dataset.data_assets["runsheet"].path,
                    }
                ]
            ):
                vp.add(
                    bulkRNASeq.checks.check_dge_table_annotation_columns_exist,
                    full_description=textwrap.dedent(f"""
                            - Check: Ensure ['SYMBOL','GENENAME','REFSEQ','ENTREZID','STRING_id','GOSLIMS_IDS','ENSEMBL'] columns exist (note: 'ENSEMBL' is replaced by 'TAIR' for Arabidospis)
                                - Reason:
                                    - These columns should be populated during annotation for all single gene mapping probes
                                - Potential Source of Problems:
                                    - Bug in processing script during annotation step
                                - Flag Condition:
                                    - Green: All columns present
                                    - Halt: At least one column is missing
                        """)
                    )
        with vp.component_start(
            name="Normalized expression",
            description="",
        ):
            with vp.payload(
                payloads=[
                    {
                        "normalized_expression": lambda: dataset.data_assets["normalized expression table"].path,
                        "samples": dataset.samples
                    }
                ]
            ):
                vp.add(
                    checks.check_normalized_expression_table,
                    full_description=textwrap.dedent(f"""
                        - Check: Ensure normalized expression table has all samples and values within (-inf,+inf)
                            - Reason:
                                - Part of processing output. Note: Values are log2 transformed
                            - Potential Source of Problems:
                                - Bug in processing script
                            - Flag Condition:
                                - Green: All conditions met
                                - Halt: At least one condition failed
                    """)
                    )
            with vp.payload(
                payloads=[
                    {
                        "organism": lambda: dataset.metadata["organism"],
                        "samples": lambda: set(dataset.samples),
                        "dge_table": lambda: dataset.data_assets[
                            "normalized expression table"
                        ].path,
                        "runsheet": lambda: dataset.data_assets["runsheet"].path,
                    }
                ]
            ):
                vp.add(
                    bulkRNASeq.checks.check_dge_table_annotation_columns_exist,
                    full_description=textwrap.dedent(f"""
                            - Check: Ensure ['SYMBOL','GENENAME','REFSEQ','ENTREZID','STRING_id','GOSLIMS_IDS','ENSEMBL'] columns exist (note: 'ENSEMBL' is replaced by 'TAIR' for Arabidospis)
                                - Reason:
                                    - These columns should be populated during annotation for all single gene mapping probes
                                - Potential Source of Problems:
                                    - Bug in processing script during annotation step
                                - Flag Condition:
                                    - Green: All columns present
                                    - Halt: At least one column is missing
                        """)
                    )
        with vp.component_start(
            name="DE Metadata",
            description="",
        ):

            with vp.component_start(
                name="Sample Table",
                description="",
            ):
                with vp.payload(
                    payloads=[
                        {
                            "runsheet": lambda: dataset.data_assets["runsheet"].path,
                            "sampleTable": lambda: dataset.data_assets[
                                "sample table"
                            ].path,
                        }
                    ]
                ):
                    vp.add(
                        checks.check_sample_table_against_runsheet,
                        config={"all_samples_required": True},
                        full_description=textwrap.dedent(f"""
                            - Check: Ensure all samples denoted in the runsheet are present
                                - Reason:
                                    - Sample Table should be inclusive of all samples processed
                                - Potential Source of Problems:
                                    - Bug in processing script that results in missing sample table columns
                                - Flag Condition:
                                    - Green: All samples present
                                    - Halt: At least one sample indicated in runsheet is missing
                        """)
                    )
                    vp.add(
                        checks.check_sample_table_for_correct_group_assignments,
                        full_description=textwrap.dedent(f"""
                            - Check: Ensure sample to group mapping consistent with runsheet
                                - Reason:
                                    - Group mapping indicated by sample table. Mis-mapping will result in incorrect DE output
                                - Potential Source of Problems:
                                    - Bug in processing script
                                - Flag Condition:
                                    - Green: Consistent sample to group mapping based on runsheet
                                    - Halt: At least one inconsistency in mapping
                        """)
                        )

            with vp.component_start(
                name="Contrasts Table",
                description="",
            ):
                with vp.payload(
                    payloads=[
                        {
                            "runsheet": lambda: dataset.data_assets["runsheet"].path,
                            "contrasts_table": lambda: dataset.data_assets[
                                "DE contrasts table"
                            ].path,
                        }
                    ]
                ):
                    vp.add(
                        checks.check_contrasts_table_headers,
                        full_description=textwrap.dedent(f"""
                            - Check: Ensure contrast table header correctly formatted using runsheet as reference for groups
                                - Reason:
                                    - Incorrect contrasts will result in incorrect DE output
                                - Potential Source of Problems:
                                    - Bug in processing script
                                - Flag Condition:
                                    - Green: Consistent contrast header based on runsheet
                                    - Halt: At least one inconsistency in header
                        """)
                        )
                    vp.add(
                        bulkRNASeq.checks.check_contrasts_table_rows,
                        full_description=textwrap.dedent(f"""
                            - Check: Ensure contrast table rows correctly formatted using runsheet as reference for groups
                                - Reason:
                                    - Incorrect rows will result in incorrect DE output as groups will become mis-mapped in contrasts
                                - Potential Source of Problems:
                                    - Bug in processing script
                                - Flag Condition:
                                    - Green: Consistent contrast rows (i.e. groups) based on runsheet 
                                    - Halt: At least one inconsistency in rows
                        """)
                        )
        with vp.component_start(
            name="DE Output",
            description="",
        ):
            with vp.payload(
                payloads=[
                    {
                        "organism": lambda: dataset.metadata["organism"],
                        "samples": lambda: set(dataset.samples),
                        "dge_table": lambda: dataset.data_assets[
                            "DE annotated table"
                        ].path,
                        "runsheet": lambda: dataset.data_assets["runsheet"].path,
                    }
                ]
            ):
                vp.add(
                    bulkRNASeq.checks.check_dge_table_annotation_columns_exist,
                    full_description=textwrap.dedent(f"""
                            - Check: Ensure ['SYMBOL','GENENAME','REFSEQ','ENTREZID','STRING_id','GOSLIMS_IDS','ENSEMBL'] columns exist (note: 'ENSEMBL' is replaced by 'TAIR' for Arabidospis)
                                - Reason:
                                    - These columns should be populated during annotation for all single gene mapping probes
                                - Potential Source of Problems:
                                    - Bug in processing script during annotation step
                                - Flag Condition:
                                    - Green: All columns present
                                    - Halt: At least one column is missing
                        """)
                    )
                vp.add(
                    bulkRNASeq.checks.check_dge_table_sample_columns_exist,
                    full_description=textwrap.dedent(f"""
                            - Check: Ensure all samples denoted in the runsheet are present
                                - Reason:
                                    - DE Table should be inclusive of all samples processed, these columns hold the normalized log2 expression for each probe
                                - Potential Source of Problems:
                                    - Bug in processing script that results in missing sample columns
                                - Flag Condition:
                                    - Green: All samples present
                                    - Halt: Expected set of samples is the set of samples in the table (either missing or extra elements)
                        """)
                    )
                vp.add(
                    checks.check_dge_table_sample_columns_constraints,
                    full_description=textwrap.dedent(f"""
                            - Check: Ensure all sample columns contain float values, consistent with probe log2 transformed normalized expression values
                                - Reason:
                                    - All these values should be populated with a float value
                                - Potential Source of Problems:
                                    - Bug in processing script that results in either null values or improper value type
                                - Flag Condition:
                                    - Green: All values are float values
                                    - Halt: At least value is null or a non-float type
                        """)
                    )
                vp.add(
                    bulkRNASeq.checks.check_dge_table_group_columns_exist,
                    full_description=textwrap.dedent(f"""
                            - Check: Ensure all groups are present using runsheet as reference for groups
                                - Reason:
                                    - All groups MUST be processed during DE
                                - Potential Source of Problems:
                                    - Bug in processing script
                                - Flag Condition:
                                    - Green: All groups represented
                                    - Halt: Expected set of groups not found (either missing or extra elements)
                        """)
                    )
                vp.add(
                    bulkRNASeq.checks.check_dge_table_group_columns_constraints,
                    full_description=textwrap.dedent(f"""
                            - Check: Ensure ["Group.Stdev_", "Group.Mean_"] are properly computed
                                - Reason:
                                    - These columns are part of proper DE output
                                - Potential Source of Problems:
                                    - Bug in processing script, possibly a calculation formula/method bug
                                - Flag Condition:
                                    - Green: Within 0.001, all values are consistent with recomputation
                                    - Halt: At least one value is greater than 0.001 different than the recomputed value
                        """)
                    )
                vp.add(
                    checks.check_dge_table_comparison_statistical_columns_exist,
                    full_description=textwrap.dedent(f"""
                            - Check: Ensure all comparison related columns exist in both directions: ["Log2fc_(GROUP1vGROUP2)", "T.stat_(GROUP1vGROUP2)", "P.value_(GROUP1vGROUP2)", "Adj.p.value_(GROUP1vGROUP2)"]  
                                - Reason:
                                    - These columns are part of proper DE output
                                - Potential Source of Problems:
                                    - Bug in processing script, possibly a mis-mapping of sample to groups
                                - Flag Condition:
                                    - Green: All columns exist
                                    - Halt: At least one column is missing
                        """)
                    )
                vp.add(
                    checks.check_dge_table_group_statistical_columns_constraints,
                    full_description=textwrap.dedent(f"""
                            - Check: Ensure all comparison related columns met the following constraints: ["Log2fc: not Null", "T.stat: not Null", "P.value: Not negative, Not Null", "Adj.p.value: Not negative, Not Null"]  
                                - Reason:
                                    - These columns are part of proper DE output
                                - Potential Source of Problems:
                                    - Bug in processing script, possibly a miscalculation or column mapping issue
                                - Flag Condition:
                                    - Green: All columns met constraints
                                    - Halt: At least one constraint failed
                        """)
                    )
                vp.add(
                    checks.check_dge_table_fixed_statistical_columns_exist,
                    full_description=textwrap.dedent(f"""
                            - Check: Ensure all the following dataset wide columns exist ["All.mean", "All.stdev", "F.p.value"]  
                                - Reason:
                                    - These columns are part of proper DE output
                                - Potential Source of Problems:
                                    - Bug in processing script
                                - Flag Condition:
                                    - Green: All columns exist
                                    - Halt: At least one column is missing
                        """)
                    )
                vp.add(
                    checks.check_dge_table_fixed_statistical_columns_constraints,
                    full_description=textwrap.dedent(f"""
                            - Check: Ensure all the following dataset wide columns met the following constraints ["All.mean": nonNull & nonNegative, "All.stdev": nonNull & nonNegative, "F.p.value": nonNegative]  
                                - Reason:
                                    - These columns are part of proper DE output
                                - Potential Source of Problems:
                                    - Bug in processing script, possibly a miscalculation or column mapping issue
                                - Flag Condition:
                                    - Green: All columns met constraints
                                    - Halt: At least one constraint failed
                        """)
                    )
                vp.add(
                    checks.check_dge_table_log2fc_within_reason,
                    full_description=textwrap.dedent(f"""
                            - Check: Ensure log2fc for all comparisons is reasonable based on log2fc normalized expression values
                                - Reason:
                                    - These columns are part of proper DE output
                                - Potential Source of Problems:
                                    - Bug in processing script
                                - Flag Condition:
                                    - Green: All log2fc values are consistent with recomputation within 0.1% tolerance
                                    - Halt: At least 1 log2fc value inconsistent with recomputation (greater than tolerance)
                        """)
                    )

            with vp.component_start(
                name="Viz Tables",
                description="Extended from the dge tables",
            ):
                with vp.payload(
                    payloads=[
                        {
                            "samples": lambda: set(dataset.samples),
                            "pca_table": lambda: dataset.data_assets[
                                "viz PCA table"
                            ].path,
                        }
                    ]
                ):
                    vp.add(
                        bulkRNASeq.checks.check_viz_pca_table_index_and_columns_exist,
                        full_description=textwrap.dedent(f"""
                                - Check: Ensure all samples (row-indices) present and columns PC1, PC2 and PC3 are present
                                    - Reason:
                                        - PCA table should include all samples and PC1, PC2, PC3 (for 3D PCA viz)
                                    - Potential Source of Problems:
                                        - Bug in processing script
                                    - Flag Condition:
                                        - Green: All samples and all columns present
                                        - Halt: At least one sample or column is missing
                            """)
                        )
        with vp.component_start(
            name="Processing Report",
            description="",
        ):
            vp.add_manual(
                description = "Loading report",
                start_instructions = "Load html report into web browser",
                pass_fail_questions = [
                    "Does the report render without issue?",
                ]
            )
            vp.add_manual(
                description = "Section 2 Load Metadata and Raw Data",
                start_instructions = "Navigate to section 2",
                pass_fail_questions = [
                    "Does the content of the table match the runsheet?",
                    "Do the number of entries in the embedded runsheet match the number of samples/runsheet-rows?",
                ]
            )
            vp.add_manual(
                description = "Section 3 QA For Raw Data - Density Plots",
                start_instructions = "Navigate to section 3 - Density Plots",
                pass_fail_questions = [
                    "Is every sample included in the legend?",
                    "Are axes and axe titles clear?",
                ],
                pass_flag_questions = [
                    "Are all lines have similar one or two peak shape, each line not likely overlapping?",
                ]
            )
            vp.add_manual(
                description = "Section 3 QA For Raw Data - Pseudo Image Plots",
                start_instructions = "Navigate to section 3 - Pseudo Image Plots",
                pass_fail_questions = [
                    "Does every sample have its own plot?",
                    "Is each plot white-blue gradient?",
                ],
                pass_flag_questions = [
                    "Are there clear physical streaks, bright (i.e. white) circles or other features?",
                ]
            )
            vp.add_manual(
                description = "Section 3 QA For Raw Data - MA Plots",
                start_instructions = "Navigate to section 3 - MA Plots",
                pass_fail_questions = [
                    "Does every sample have its own plot?",
                ],
                pass_flag_questions = [
                    "Are positive control points distributed from left end to right end of MA cloud",
                    "Are negative control points concentrated at left end of MA cloud?",
                ]
            )
            vp.add_manual(
                description = "Section 3 QA For Raw Data - Foreground Background Plots",
                start_instructions = "Navigate to section 3 - Foreground Background Plots",
                pass_fail_questions = [
                    "Does every sample have its own plot?",
                ],
                pass_flag_questions = [
                    "Are the vast majority (i.e. 99.9%) of points are above the blue diagonal line",
                ]
            )
            vp.add_manual(
                description = "Section 3 QA For Raw Data - Boxplots",
                start_instructions = "Navigate to section 3 - Boxplots",
                pass_fail_questions = [
                    "Does every sample have its own plot?",
                ],
                pass_flag_questions = [
                    "Do the boxplots have overlapping distributions?",
                ]
            )
            ###################################################
            ### Normalized data plots
            ###################################################
            vp.add_manual(
                description = "Section 6 Normalized Data Quality Assessment - Density Plots",
                start_instructions = "Navigate to section 6 - Density Plots",
                pass_fail_questions = [
                    "Is every sample included in the legend?",
                    "Are axes and axe titles clear?",
                ],
                pass_flag_questions = [
                    "Are all lines nearly fully overlapping?",
                ]
            )
            vp.add_manual(
                description = "Section 6 Normalized Data Quality Assessment - Pseudo Image Plots",
                start_instructions = "Navigate to section 6 - Pseudo Image Plots",
                pass_fail_questions = [
                    "Does every sample have its own plot?",
                ],
                pass_flag_questions = [
                    "Is each plot white-blue gradient?",
                    "Are there clear physical streaks, bright (i.e. white) circles or other features?",
                ]
            )
            vp.add_manual(
                description = "Section 6 Normalized Data Quality Assessment - MA Plots",
                start_instructions = "Navigate to section 6 - MA Plots",
                pass_fail_questions = [
                    "Does every sample have its own plot?",
                ],
                pass_flag_questions = [
                    "Are positive control points distributed from left end to right end of MA 'cloud'",
                    "Are negative control points concentrated at left end of MA 'cloud'?",
                ]
            )
            vp.add_manual(
                description = "Section 6 Normalized Data Quality Assessment - Boxplots",
                start_instructions = "Navigate to section 6 - Boxplots",
                pass_fail_questions = [
                    "Does every sample have its own plot?",
                ],
                pass_flag_questions = [
                    "Do the boxplots have largely the same distributions? (i.e. are the boxes all nearly horizontally aligned)",
                ]
            )
    # return protocol object without running or generating a report
    if defer_run:
        return vp

    vp.run(**run_args)

    # return report
    return vp.report(**report_args, combine_with_flags=dataset.loaded_assets_dicts)
