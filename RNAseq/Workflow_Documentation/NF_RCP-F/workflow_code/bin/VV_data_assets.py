#! /usr/bin/env python
""" Validation/Verification for raw reads in RNASeq Concensus Pipeline
"""
import argparse
from pathlib import Path

from dp_tools.bulkRNASeq.vv_protocols import validate_bulkRNASeq
from dp_tools.core.check_model import ValidationProtocol
from dp_tools.core.loaders import load_data

ASSET_KEY_SETS = [
    "glds metadata",
    "demuliplexed paired end raw data",  # paired only
    "qc reports for paired end raw data",  # paired only
    "paired end trimmed reads",  # paired only
    "qc reports for paired end trimmed reads data",  # paired only
    "RSeQC output for paired end data",  # paired only
    "demuliplexed single end raw data",  # single only
    "qc reports for single end raw data",  # single only
    "single end trimmed reads",  # single only
    "qc reports for single end trimmed reads data",  # single only
    "RSeQC output for single end data",  # single only
    "STAR alignments",
    "RSEM counts",
    "DGE Output",
    "ERCC DGE Output",  # ERCC Only
    "RSEM Output",
]

SKIP_LIST_COMPONENTS = [
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

##############################################################
# Utility Functions To Handle Logging, Config and CLI Arguments
##############################################################
def _parse_args():
    """Parse command line args."""
    parser = argparse.ArgumentParser()

    parser.add_argument("--root-path", required=True, help="Root data path")
    parser.add_argument(
        "--runsheet-path", required=True, help="Path to dataset runsheet.csv file"
    )

    parser.add_argument("--accession", required=True, help="Accession number")
    parser.add_argument(
        "--max-flag-code",
        default=80,
        help="Throw an exception if any flag code exceeds this value",
    )
    parser.add_argument(
        "--data-asset-sets",
        nargs="+",
        help=f"Data assets keys to load. Must be a subset of this list: {ASSET_KEY_SETS}",
        required=True,
    )
    parser.add_argument(
        "--run-components",
        nargs="+",
        help=f"V&V Components to run. Must be a subset of this list: {SKIP_LIST_COMPONENTS}",
        required=True,
    )

    args = parser.parse_args()
    return args


def main():
    args = _parse_args()
    max_flag_code = int(args.max_flag_code)

    # invert requested components
    # this is because these are supplied as skips
    assert set(args.run_components).issubset(
        set(SKIP_LIST_COMPONENTS)
    ), "Must supply valid components (see help menu)"
    skip_components = list(set(SKIP_LIST_COMPONENTS) - set(args.run_components))

    datasystem = load_data(
        key_sets=args.data_asset_sets,
        config=("bulkRNASeq", "Latest"),
        root_path=Path(args.root_path),
        runsheet_path=Path(args.runsheet_path),
    )

    vp = validate_bulkRNASeq(
        datasystem.dataset,
        report_args={"include_skipped": True},
        protocol_args={"skip_components": skip_components},
        defer_run=True,
    )

    print(f"{'QUEUED CHECK COMPONENT TREE':*^60}")
    print(vp.queued_checks(include_individual_checks=False))
    print(f"{'QUEUED CHECK COMPONENT TREE WITH SKIPPED COMPONENTS':*^60}")
    print(
        vp.queued_checks(
            include_individual_checks=False, include_skipped_components=True
        )
    )
    print(f"{'QUEUED CHECK COMPONENT TREE WITH INVIDUAL CHECKS':*^60}")
    print(vp.queued_checks(include_individual_checks=True))

    vp.run()
    report = vp.report(
        include_skipped=False, combine_with_flags=datasystem.dataset.loaded_assets_dicts
    )
    # data assets loading information

    # output default dataframe
    samples = list(datasystem.dataset.samples)
    ValidationProtocol.append_sample_column(report["flag_table"], samples=samples)
    df = report["flag_table"]
    output_fn = f"VV_log.tsv"
    df.to_csv(output_fn, sep="\t")

    # halt on error
    flagged_messages = "\n".join(
        [msg for msg in df.loc[df["code_level"] >= max_flag_code]["message"]]
    )
    assert (
        df["code_level"].max() < max_flag_code
    ), f"Maximum flag code exceeded: {max_flag_code}. Printing flag messages that caused this halt: {flagged_messages}"


if __name__ == "__main__":
    main()
