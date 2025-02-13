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
