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
        "expected_attrs": ["paired_end", "has_ERCC", "organism"]
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
    "STAR Alignments By Sample-check_thresholds-Mapped": {
        "mqc_key": "STAR",
        "stat_string": "uniquely_mapped_percent + multimapped_percent",
        "thresholds": [
            {"code": "YELLOW", "type": "lower", "value": 70},
            {"code": "RED", "type": "lower", "value": 50},
        ],
    },
    "STAR Alignments By Sample-check_thresholds-MultiMapped": {
        "mqc_key": "STAR",
        "stat_string": "multimapped_toomany_percent + multimapped_percent",
        "thresholds": [
            {"code": "YELLOW", "type": "lower", "value": 30},
            {"code": "RED", "type": "lower", "value": 15},
        ],
    },
    "STAR Alignments-check_for_outliers": {
        "mqc_module": "STAR",
        "mqc_plot": "general_stats",
        "mqc_keys": [
            "uniquely_mapped_percent",
            "avg_mapped_read_length",
            "mismatch_rate",
            "deletion_rate",
            "deletion_length",
            "insertion_rate",
            "insertion_length",
            "multimapped_percent",
            "multimapped_toomany_percent",
            "unmapped_mismatches_percent",
            "unmapped_tooshort_percent",
            "unmapped_other_percent",
        ],
        "thresholds": [
            {"code": "YELLOW", "stdev_threshold": 2, "middle_fcn": "median"},
            {"code": "RED", "stdev_threshold": 4, "middle_fcn": "median"},
        ],
    },
    "RSeQC-check_for_outliers-geneBody_coverage": {
        "mqc_module": "RSeQC",
        "mqc_plot": "Gene Body Coverage",
        "mqc_keys": ["_ALL"],
        "thresholds": [
            {"code": "YELLOW", "stdev_threshold": 2, "middle_fcn": "median"},
            {"code": "RED", "stdev_threshold": 4, "middle_fcn": "median"},
        ],
    },
    "RSeQC-check_for_outliers-infer_experiment": {
        "mqc_module": "RSeQC",
        "mqc_plot": "Infer experiment",
        "mqc_keys": ["_ALL"],
        "thresholds": [
            {"code": "YELLOW", "stdev_threshold": 2, "middle_fcn": "median"},
            {"code": "RED", "stdev_threshold": 4, "middle_fcn": "median"},
        ],
    },
    "RSeQC-check_for_outliers-inner_distance": {
        "mqc_module": "RSeQC",
        "mqc_plot": "Inner Distance",
        "mqc_keys": ["_ALL"],
        "thresholds": [
            {"code": "YELLOW", "stdev_threshold": 2, "middle_fcn": "median"},
            {"code": "RED", "stdev_threshold": 4, "middle_fcn": "median"},
        ],
    },
    "RSeQC-check_for_outliers-read_distribution": {
        "mqc_module": "RSeQC",
        "mqc_plot": "Read Distribution",
        "mqc_keys": ["_ALL"],
        "thresholds": [
            {"code": "YELLOW", "stdev_threshold": 2, "middle_fcn": "median"},
            {"code": "RED", "stdev_threshold": 4, "middle_fcn": "median"},
        ],
    },
    "RSeQC-check_strandedness_assessable_from_infer_experiment": {
        "stranded_assessment_range": {"max": 100, "min": 75},
        "unstranded_assessment_range": {"min": 40, "max": 60},
        "valid_dominant_strandedness_assessments": [
            "Sense (% Tags)",
            "Antisense (% Tags)",
        ],
    },
    "RSEM Counts-check_for_outliers": {
        "mqc_module": "Rsem",
        "mqc_plot": "general_stats",
        "mqc_keys": [
            "Unalignable",
            "Alignable",
            "Filtered",
            "Total",
            "alignable_percent",
            "Unique",
            "Multi",
            "Uncertain",
        ],
        "thresholds": [
            {"code": "YELLOW", "stdev_threshold": 2, "middle_fcn": "median"},
            {"code": "RED", "stdev_threshold": 4, "middle_fcn": "median"},
        ],
    },
}

# Manual kept in sync for now
COMPONENTS_LIST = [
    "Metadata",  # for raw reads V&V
    "Raw Reads",  # for raw reads V&V
    "Raw Reads By Sample",  # for raw reads V&V
    "Trim Reads",  # for trim reads V&V
    "Trimmed Reads By Sample",  # for trim reads V&V
    "STAR Alignments",  # for star alignment V&V
    "STAR Alignments By Sample",  # for star alignment V&V
    "RSeQC By Sample",  # for RSeQC V&V
    "RSeQC",  # for RSeQC V&V
    "RSEM Counts",  # for after RSEM V&V
    "Unnormalized Gene Counts",  # for after RSEM V&V
    "DGE Metadata",  # for post DGE
    "DGE Metadata ERCC",  # for post DGE
    "DGE Output",  # for post DGE
    "DGE Output ERCC",  # for post DGE
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

        with vp.component_start(
            name="RSeQC",
            description="RSeQC submodule outliers checking and other submodule specific dataset wide checks",
        ):
            with vp.payload(
                payloads=[
                    {
                        "dataset": dataset,
                        "data_asset_keys": ["genebody coverage out"],
                    }
                ]
            ):
                vp.add(
                    check_for_outliers,
                    description="Check for outliers in geneBody Coverage",
                    config=config["RSeQC-check_for_outliers-geneBody_coverage"],
                )
            with vp.payload(
                payloads=[
                    {
                        "dataset": dataset,
                        "data_asset_keys": ["infer experiment out"],
                    }
                ]
            ):
                vp.add(
                    check_for_outliers,
                    description="Check for outliers in infer experiment",
                    config=config["RSeQC-check_for_outliers-infer_experiment"],
                )
            with vp.payload(
                payloads=[
                    {
                        "dataset": dataset,
                        "data_asset_keys": ["inner distance out"],
                    }
                ]
            ):
                vp.add(
                    check_for_outliers,
                    description="Check for outliers in inner distance",
                    config=config["RSeQC-check_for_outliers-inner_distance"],
                    skip=(not dataset.metadata["paired_end"]),
                )
            with vp.payload(
                payloads=[
                    {
                        "dataset": dataset,
                        "data_asset_keys": ["read distribution out"],
                    }
                ]
            ):
                vp.add(
                    check_for_outliers,
                    description="Check for outliers in read distribution",
                    config=config["RSeQC-check_for_outliers-read_distribution"],
                )

            with vp.payload(payloads=[{"dataset": dataset}]):
                vp.add(
                    check_strandedness_assessable_from_infer_experiment,
                    config=config[
                        "RSeQC-check_strandedness_assessable_from_infer_experiment"
                    ],
                )
            with vp.payload(
                payloads=[
                    {
                        "samples": list(dataset.samples),
                        "multiqc_report_path": lambda: dataset.data_assets[
                            "genebody coverage MultiQC directory"
                        ].path,
                    },
                    {
                        "samples": list(dataset.samples),
                        "multiqc_report_path": lambda: dataset.data_assets[
                            "infer experiment MultiQC directory"
                        ].path,
                        "name_reformat_func": lambda: lambda s: re.sub(
                            "_infer_expt$", "", s
                        ),
                    },
                    {
                        "samples": list(dataset.samples),
                        "multiqc_report_path": lambda: dataset.data_assets[
                            "read distribution MultiQC directory"
                        ].path,
                        "name_reformat_func": lambda: lambda s: re.sub(
                            "_read_dist$", "", s
                        ),
                    },
                ]
            ):
                vp.add(
                    check_sample_in_multiqc_report,
                    description="Check all samples are present in RSeQC multiQC reports",
                )
            with vp.payload(
                payloads=[
                    {
                        "samples": list(dataset.samples),
                        "multiqc_report_path": lambda: dataset.data_assets[
                            "inner distance MultiQC directory"
                        ].path,
                    },
                ]
            ):
                vp.add(
                    check_sample_in_multiqc_report,
                    description="Check all samples are present in RSeQC inner distance multiQC report (paired end only)",
                    skip=(not dataset.metadata["paired_end"]),
                )
        with vp.component_start(
            name="RSEM Counts",
            description="Dataset wide checks including outliers detection",
        ):
            with vp.payload(
                payloads=[
                    {
                        "dataset": dataset,
                        "data_asset_keys": ["sample counts stats directory"],
                    }
                ]
            ):
                vp.add(
                    check_for_outliers, config=config["RSEM Counts-check_for_outliers"]
                )
            with vp.payload(
                payloads=[
                    {
                        "samples": list(dataset.samples),
                        "multiqc_report_path": lambda: dataset.data_assets[
                            "RSEM counts MultiQC directory"
                        ].path,
                    },
                ]
            ):
                vp.add(
                    check_sample_in_multiqc_report,
                    description="Check all samples are present in RSEM multiQC report",
                )
        with vp.component_start(
            name="Unnormalized Gene Counts",
            description="Validate normalization related output",
        ):

            with vp.payload(
                payloads=[
                    {
                        "unnormalizedCountTable": lambda: dataset.data_assets[
                            "star unnormalized counts table"
                        ].path,
                        "samplewise_tables": lambda: {
                            s.name: s.data_assets["sample reads per gene table"].path
                            for s in dataset.samples.values()
                        },
                    },
                ]
            ):
                vp.add(
                    check_aggregate_star_unnormalized_counts_table_values_against_samplewise_tables
                )
            with vp.payload(
                payloads=[
                    {
                        "unnormalizedCountTable": lambda: dataset.data_assets[
                            "rsem unnormalized counts table"
                        ].path,
                        "samplewise_tables": lambda: {
                            s.name: s.data_assets["sample gene counts table"].path
                            for s in dataset.samples.values()
                        },
                    },
                ]
            ):
                vp.add(
                    check_aggregate_rsem_unnormalized_counts_table_values_against_samplewise_tables
                )
                vp.add(
                    check_ERCC_subgroup_representation,
                    skip=(not dataset.metadata["has_ERCC"]),
                )

        with vp.component_start(
            name="DGE Metadata",
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
                        check_sample_table_against_runsheet,
                        config={"all_samples_required": True},
                    )
                    vp.add(check_sample_table_for_correct_group_assignments)

            with vp.component_start(
                name="Contrasts Tables",
                description="",
            ):
                with vp.payload(
                    payloads=[
                        {
                            "runsheet": lambda: dataset.data_assets["runsheet"].path,
                            "contrasts_table": lambda: dataset.data_assets[
                                "DESeq2 contrasts table"
                            ].path,
                        }
                    ]
                ):
                    vp.add(check_contrasts_table_headers)
                    vp.add(check_contrasts_table_rows)

        with vp.component_start(
            name="DGE Metadata ERCC",
            description="",
            skip=(not dataset.metadata["has_ERCC"]),
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
                                "ERCC sample table"
                            ].path,
                        }
                    ]
                ):
                    vp.add(
                        check_sample_table_against_runsheet,
                        config={"all_samples_required": False},
                    )
                    vp.add(check_sample_table_for_correct_group_assignments)

            with vp.component_start(
                name="Contrasts Tables",
                description="",
            ):
                with vp.payload(
                    payloads=[
                        {
                            "runsheet": lambda: dataset.data_assets["runsheet"].path,
                            "contrasts_table": lambda: dataset.data_assets[
                                "ERCC normalized DESeq2 contrasts table"
                            ].path,
                        }
                    ]
                ):
                    vp.add(check_contrasts_table_headers)
                    vp.add(check_contrasts_table_rows)

        with vp.component_start(
            name="DGE Output",
            description="",
        ):
            with vp.payload(
                payloads=[
                    {
                        "rsem_table_path": lambda: dataset.data_assets[
                            "rsem unnormalized counts table"
                        ].path,
                        "deseq2_table_path": lambda: dataset.data_assets[
                            "DESeq2 unnormalized counts table"
                        ].path,
                    }
                ]
            ):
                vp.add(
                    check_rsem_counts_and_unnormalized_tables_parity,
                    skip=(
                        "rsem unnormalized counts table" not in dataset.data_assets
                        or "DESeq2 unnormalized counts table" not in dataset.data_assets
                    ),
                )

            with vp.payload(
                payloads=[
                    {
                        "organism": lambda: dataset.metadata["organism"],
                        "samples": lambda: set(dataset.samples),
                        "dge_table": lambda: dataset.data_assets[
                            "DESeq2 annotated DGE table"
                        ].path,
                        "runsheet": lambda: dataset.data_assets["runsheet"].path,
                    }
                ]
            ):
                vp.add(check_dge_table_annotation_columns_exist)
                vp.add(check_dge_table_sample_columns_exist)
                vp.add(check_dge_table_sample_columns_constraints)
                vp.add(check_dge_table_group_columns_exist)
                vp.add(check_dge_table_group_columns_constraints)
                vp.add(check_dge_table_comparison_statistical_columns_exist)
                vp.add(check_dge_table_group_statistical_columns_constraints)
                vp.add(check_dge_table_fixed_statistical_columns_exist)
                vp.add(check_dge_table_fixed_statistical_columns_constraints)
                vp.add(check_dge_table_log2fc_within_reason)

            with vp.component_start(
                name="Viz Tables",
                description="Extended from the dge tables",
            ):
                with vp.payload(
                    payloads=[
                        {
                            "organism": lambda: dataset.metadata["organism"],
                            "samples": lambda: set(dataset.samples),
                            "dge_table": lambda: dataset.data_assets[
                                "DESeq2 annotated DGE extended for viz table"
                            ].path,
                            "runsheet": lambda: dataset.data_assets["runsheet"].path,
                        }
                    ]
                ):
                    vp.add(check_dge_table_annotation_columns_exist)
                    vp.add(check_dge_table_sample_columns_exist)
                    vp.add(check_dge_table_sample_columns_constraints)
                    vp.add(check_dge_table_group_columns_exist)
                    vp.add(check_dge_table_group_columns_constraints)
                    vp.add(check_dge_table_comparison_statistical_columns_exist)
                    vp.add(check_dge_table_group_statistical_columns_constraints)
                    vp.add(check_dge_table_fixed_statistical_columns_exist)
                    vp.add(check_dge_table_fixed_statistical_columns_constraints)
                    vp.add(check_dge_table_log2fc_within_reason)
                    vp.add(check_viz_table_columns_exist)
                    vp.add(check_viz_table_columns_constraints)

                with vp.payload(
                    payloads=[
                        {
                            "samples": lambda: set(dataset.samples),
                            "pca_table": lambda: dataset.data_assets[
                                "DESeq2 viz PCA table"
                            ].path,
                        }
                    ]
                ):
                    vp.add(check_viz_pca_table_index_and_columns_exist)

        with vp.component_start(
            name="DGE Output ERCC",
            description="",
            skip=(not dataset.metadata["has_ERCC"]),
        ):
            with vp.payload(
                payloads=[
                    {
                        "organism": lambda: dataset.metadata["organism"],
                        "samples": lambda: set(
                            pd.read_csv(
                                dataset.data_assets["ERCC sample table"].path,
                                index_col=0,
                            ).index
                        ),
                        "dge_table": lambda: dataset.data_assets[
                            "ERCC normalized DESeq2 annotated DGE table"
                        ].path,
                        "runsheet": lambda: dataset.data_assets["runsheet"].path,
                    }
                ]
            ):
                vp.add(check_dge_table_annotation_columns_exist)
                vp.add(check_dge_table_sample_columns_exist)
                vp.add(check_dge_table_sample_columns_constraints)
                vp.add(check_dge_table_group_columns_exist)
                vp.add(check_dge_table_group_columns_constraints)
                vp.add(check_dge_table_comparison_statistical_columns_exist)
                vp.add(check_dge_table_group_statistical_columns_constraints)
                vp.add(check_dge_table_fixed_statistical_columns_exist)
                vp.add(check_dge_table_fixed_statistical_columns_constraints)
                vp.add(check_dge_table_log2fc_within_reason)

            with vp.component_start(
                name="Viz Tables",
                description="Extended from the dge tables",
            ):
                with vp.payload(
                    payloads=[
                        {
                            "organism": lambda: dataset.metadata["organism"],
                            "samples": lambda: set(
                                pd.read_csv(
                                    dataset.data_assets["ERCC sample table"].path,
                                    index_col=0,
                                ).index
                            ),
                            "dge_table": lambda: dataset.data_assets[
                                "ERCC normalized DESeq2 annotated DGE extended for viz table"
                            ].path,
                            "runsheet": lambda: dataset.data_assets["runsheet"].path,
                        }
                    ]
                ):
                    vp.add(check_dge_table_annotation_columns_exist)
                    vp.add(check_dge_table_sample_columns_exist)
                    vp.add(check_dge_table_sample_columns_constraints)
                    vp.add(check_dge_table_group_columns_exist)
                    vp.add(check_dge_table_group_columns_constraints)
                    vp.add(check_dge_table_comparison_statistical_columns_exist)
                    vp.add(check_dge_table_group_statistical_columns_constraints)
                    vp.add(check_dge_table_fixed_statistical_columns_exist)
                    vp.add(check_dge_table_fixed_statistical_columns_constraints)
                    vp.add(check_dge_table_log2fc_within_reason)
                    vp.add(check_viz_table_columns_exist)
                    vp.add(check_viz_table_columns_constraints)

                with vp.payload(
                    payloads=[
                        {
                            "samples": lambda: set(
                                pd.read_csv(
                                    dataset.data_assets["ERCC sample table"].path,
                                    index_col=0,
                                ).index
                            ),
                            "pca_table": lambda: dataset.data_assets[
                                "ERCC normalized DESeq2 viz PCA table"
                            ].path,
                        }
                    ]
                ):
                    vp.add(check_viz_pca_table_index_and_columns_exist)

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
                                    "file": lambda: sample.data_assets[
                                        "raw forward reads fastq GZ"
                                    ].path
                                },
                                {
                                    "file": lambda: sample.data_assets[
                                        "raw reverse reads fastq GZ"
                                    ].path
                                },
                            ]
                            if dataset.metadata["paired_end"]
                            else [
                                {
                                    "file": lambda: sample.data_assets[
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
                                    "file": lambda: sample.data_assets[
                                        "trimmed forward reads fastq GZ"
                                    ].path
                                },
                                {
                                    "file": lambda: sample.data_assets[
                                        "trimmed reverse reads fastq GZ"
                                    ].path
                                },
                            ]
                            if dataset.metadata["paired_end"]
                            else [
                                {
                                    "file": lambda: sample.data_assets[
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

                with vp.component_start(
                    name="STAR Alignments By Sample",
                    description="STAR Alignment outputs",
                ):

                    with vp.payload(
                        payloads=[
                            {
                                "file": lambda: sample.data_assets[
                                    "aligned ToTranscriptome Bam"
                                ].path,
                            },
                            {
                                "file": lambda: sample.data_assets[
                                    "aligned SortedByCoord Bam"
                                ].path,
                            },
                        ]
                    ):
                        vp.add(
                            check_bam_file_integrity,
                            config={
                                "samtools_bin": "samtools"
                            },  # assumes accessible on path already
                        )

                    with vp.payload(
                        payloads=[
                            {
                                "multiqc_inputs": lambda: [
                                    sample.data_assets["aligned log Final"].path
                                ],
                            },
                        ]
                    ):
                        vp.add(
                            check_thresholds,
                            config=config[
                                "STAR Alignments By Sample-check_thresholds-Mapped"
                            ],
                            description="Check that mapping rates are reasonable, specifically most reads map to the target genome",
                        )
                        vp.add(
                            check_thresholds,
                            config=config[
                                "STAR Alignments By Sample-check_thresholds-MultiMapped"
                            ],
                            description="Check that mapping rates are reasonable, specifically that a considerable amount of reads multimap to the target genome",
                        )

                with vp.component_start(
                    name="RSeQC By Sample",
                    description="RNASeq QA outputs",
                ):
                    with vp.component_start(
                        name="geneBody_coverage",
                        description="Assess integrity of transcripts and library prep signatures",
                    ):
                        with vp.payload(
                            payloads=[
                                {
                                    "input_dir": lambda: sample.data_assets[
                                        "genebody coverage out"
                                    ].path
                                },
                            ]
                        ):
                            vp.add(check_genebody_coverage_output)
                    with vp.component_start(
                        name="inner_distance",
                        description="Reports on distance between mate reads based on gene annotations",
                        skip=(not dataset.metadata["paired_end"]),
                    ):
                        with vp.payload(
                            payloads=[
                                {
                                    "input_dir": lambda: sample.data_assets[
                                        "inner distance out"
                                    ].path
                                },
                            ]
                        ):
                            vp.add(check_inner_distance_output)
    # return protocol object without running or generating a report
    if defer_run:
        return vp

    vp.run(**run_args)

    # return report
    return vp.report(**report_args, combine_with_flags=dataset.loaded_assets_dicts)