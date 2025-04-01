from pathlib import Path
import re
from typing import Union
import yaml
import logging

from dp_tools.core.entity_model import Dataset

log = logging.getLogger(__name__)

from dp_tools.core.check_model import ValidationProtocol

from .checks import *

CONFIG = {
    "Metadata-check_metadata_attributes_exist": {
        "expected_attrs": ["paired_end", "organism"]
    },
    "Raw Reads-check_for_outliers": {
        "mqc_module": "FastQC",
        "mqc_plot": "general_stats",
        "mqc_keys": [
            "percent_gc",
            "avg_sequence_length",
            "total_sequences",
            "percent_duplicates",
        ],
        "thresholds": [
            {"code": "YELLOW", "stdev_threshold": 2, "middle_fcn": "median"},
            {"code": "RED", "stdev_threshold": 4, "middle_fcn": "median"},
        ],
    },
    "Trim Reads-check_for_outliers": {
        "mqc_module": "FastQC",
        "mqc_plot": "general_stats",
        "mqc_keys": [
            "percent_gc",
            "avg_sequence_length",
            "total_sequences",
            "percent_duplicates",
        ],
        "thresholds": [
            {"code": "YELLOW", "stdev_threshold": 2, "middle_fcn": "median"},
            {"code": "RED", "stdev_threshold": 4, "middle_fcn": "median"},
        ],
    },
    "Raw Reads By Sample-check_fastqgz_file_contents": {
        "count_lines_to_check": 200000000
    },
    "Trim Reads By Sample-check_fastqgz_file_contents": {
        "count_lines_to_check": 200000000
    },
}

# Manual kept in sync for now
COMPONENTS_LIST = [
    "Metadata",  # for raw reads V&V
    "Raw Reads",  # for raw reads V&V
    "Raw Reads By Sample",  # for raw reads V&V
    "Trim Reads",  # for trim reads V&V
    "Trimmed Reads By Sample",  # for trim reads V&V
]


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

    # Modify protocol_args to convert run_components to skip_components based on COMPONENTS_LIST
    if (
        "run_components" in protocol_args
        and protocol_args.get("run_components") is not None
    ):
        protocol_args["skip_components"] = [
            c for c in COMPONENTS_LIST if c not in protocol_args["run_components"]
        ]
        # Check if any run components are not in COMPONENTS_LIST
        if set(protocol_args["run_components"]) - set(COMPONENTS_LIST):
            raise ValueError(
                f"run_components contains components not in COMPONENTS_LIST. Unique to run_components: {set(protocol_args['run_components']) - set(COMPONENTS_LIST)}. All Components: {COMPONENTS_LIST}"
            )
        del protocol_args["run_components"]

    # init validation protocol
    vp = ValidationProtocol(**protocol_args)
    # fmt: on
    with vp.component_start(
        name=dataset.name,
        description="Validate processing from trim reads through differential gene expression output",
    ):

        with vp.component_start(
            name="Metadata", description="Metadata file validation"
        ):
            with vp.payload(payloads=[{"dataset": dataset}]):
                vp.add(
                    check_metadata_attributes_exist,
                    config=config["Metadata-check_metadata_attributes_exist"],
                )

        with vp.component_start(
            name="Raw Reads", description="Raw Reads Outliers Detection"
        ):
            with vp.payload(
                payloads=[
                    {
                        "dataset": dataset,
                        "data_asset_keys": ["raw reads fastQC ZIP"],
                    }
                ]
                if not dataset.metadata["paired_end"]
                else [
                    {
                        "dataset": dataset,
                        "data_asset_keys": [
                            "raw forward reads fastQC ZIP",
                        ],
                    },
                    {
                        "dataset": dataset,
                        "data_asset_keys": [
                            "raw reverse reads fastQC ZIP",
                        ],
                    },
                ]
            ):
                vp.add(
                    check_for_outliers, config=config["Raw Reads-check_for_outliers"]
                )

            with vp.payload(
                payloads=[
                    {
                        "samples": list(dataset.samples),
                        "multiqc_report_path": lambda: dataset.data_assets[
                            "raw MultiQC directory"
                        ].path,
                        "name_reformat_func": lambda: lambda s: re.sub(
                            "_raw|_R1_raw|_R2_raw$", "", s
                        ),
                    },
                ]
            ):
                vp.add(
                    check_sample_in_multiqc_report,
                    description="Check all samples are present in raw reads multiQC report",
                )

        with vp.component_start(
            name="Trim Reads", description="Trimmed Reads Outliers Detection"
        ):
            with vp.payload(
                payloads=[
                    {
                        "dataset": dataset,
                        "data_asset_keys": ["trimmed reads fastQC ZIP"],
                    }
                ]
                if not dataset.metadata["paired_end"]
                else [
                    {
                        "dataset": dataset,
                        "data_asset_keys": [
                            "trimmed forward reads fastQC ZIP",
                        ],
                    },
                    {
                        "dataset": dataset,
                        "data_asset_keys": [
                            "trimmed reverse reads fastQC ZIP",
                        ],
                    },
                ]
            ):
                vp.add(
                    check_for_outliers, config=config["Trim Reads-check_for_outliers"]
                )
            with vp.payload(
                payloads=[
                    {
                        "samples": list(dataset.samples),
                        "multiqc_report_path": lambda: dataset.data_assets[
                            "trimmed fastQC MultiQC directory"
                        ].path,
                        "name_reformat_func": lambda: lambda s: re.sub(
                            "_R1|_R2$", "", s
                        ),
                    },
                    {
                        "samples": list(dataset.samples),
                        "multiqc_report_path": lambda: dataset.data_assets[
                            "trimming MultiQC directory"
                        ].path,
                        "name_reformat_func": lambda: lambda s: re.sub(
                            "_raw|_R1_raw|_R2_raw$", "", s
                        ),
                    },
                ]
            ):
                vp.add(
                    check_sample_in_multiqc_report,
                    description="Check that all samples are present in the trimmed FastQC and trimming report multiQC reports",
                )
        with vp.component_start(
            name="STAR Alignments",
            description="Dataset wide checks including outliers detection",
        ):
            with vp.payload(
                payloads=[
                    {
                        "dataset": dataset,
                        "data_asset_keys": ["aligned log Final"],
                    }
                ]
            ):
                vp.add(
                    check_for_outliers,
                    config=config["STAR Alignments-check_for_outliers"],
                )
            with vp.payload(
                payloads=[
                    {
                        "samples": list(dataset.samples),
                        "multiqc_report_path": lambda: dataset.data_assets[
                            "aligned MultiQC directory"
                        ].path,
                    },
                ]
            ):
                vp.add(
                    check_sample_in_multiqc_report,
                    description="Check all samples are present in STAR multiQC report",
                )
            
        for sample in dataset.samples.values():
            with vp.component_start(
                name=sample.name, description="Samples level checks"
            ):
                with vp.component_start(
                    name="Raw Reads By Sample", description="Raw reads"
                ):
                    with vp.payload(
                        payloads=(
                            [
                                {
                                    "file": lambda sample=sample: sample.data_assets[
                                        "raw forward reads fastq GZ"
                                    ].path
                                },
                                {
                                    "file": lambda sample=sample: sample.data_assets[
                                        "raw reverse reads fastq GZ"
                                    ].path
                                },
                            ]
                            if dataset.metadata["paired_end"]
                            else [
                                {
                                    "file": lambda sample=sample: sample.data_assets[
                                        "raw reads fastq GZ"
                                    ].path
                                },
                            ]
                        )
                    ):
                        vp.add(
                            check_fastqgz_file_contents,
                            config=config[
                                "Raw Reads By Sample-check_fastqgz_file_contents"
                            ],
                        )
                        vp.add(
                            check_gzip_file_integrity,
                        )
                    with vp.payload(
                        payloads=[
                            {
                                "sample": sample,
                                "reads_key_1": "raw forward reads fastQC ZIP",
                                "reads_key_2": "raw reverse reads fastQC ZIP",
                            },
                        ],
                    ):
                        vp.add(
                            check_forward_and_reverse_reads_counts_match,
                            skip=(not dataset.metadata["paired_end"]),
                        )
                with vp.component_start(
                    name="Trimmed Reads By Sample", description="Trimmed reads"
                ):
                    with vp.payload(
                        payloads=(
                            [
                                {
                                    "file": lambda sample=sample: sample.data_assets[
                                        "trimmed forward reads fastq GZ"
                                    ].path
                                },
                                {
                                    "file": lambda sample=sample: sample.data_assets[
                                        "trimmed reverse reads fastq GZ"
                                    ].path
                                },
                            ]
                            if dataset.metadata["paired_end"]
                            else [
                                {
                                    "file": lambda sample=sample: sample.data_assets[
                                        "trimmed reads fastq GZ"
                                    ].path
                                }
                            ]
                        )
                    ):
                        vp.add(check_file_exists, description="Check reads files exist")
                        vp.add(
                            check_fastqgz_file_contents,
                            config=config[
                                "Trim Reads By Sample-check_fastqgz_file_contents"
                            ],
                        )

                    with vp.payload(
                        payloads=[
                            {
                                "sample": sample,
                                "reads_key_1": "trimmed forward reads fastQC ZIP",
                                "reads_key_2": "trimmed reverse reads fastQC ZIP",
                            },
                        ],
                    ):
                        vp.add(
                            check_forward_and_reverse_reads_counts_match,
                            skip=(not dataset.metadata["paired_end"]),
                        )

    # return protocol object without running or generating a report
    if defer_run:
        return vp

    vp.run(**run_args)

    # return report
    return vp.report(**report_args, combine_with_flags=dataset.loaded_assets_dicts)