#! /usr/bin/env python
""" Validation/Verification for microarray processed data
"""
import argparse
from pathlib import Path

from dp_tools.core.check_model import ValidationProtocol
from dp_tools.core.loaders import load_data
from dp_tools import plugin_api

ASSET_KEY_SETS = [
    "glds metadata",
    "processed",
]

SKIP_LIST_COMPONENTS = [
    "Metadata",  # for raw reads V&V
    "DGE Metadata",  # for post DGE
    "DGE Output",  # for post DGE
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
    parser.add_argument(
        "--plugin-dir",
        help=f"Plugin directory to load",
        required=True,
    )

    args = parser.parse_args()
    return args


def main():
    args = _parse_args()
    max_flag_code = int(args.max_flag_code)

    plugin = plugin_api.load_plugin(Path(args.plugin_dir))

    # invert requested components
    # this is because these are supplied as skips
    assert set(args.run_components).issubset(
        set(SKIP_LIST_COMPONENTS)
    ), "Must supply valid components (see help menu)"
    skip_components = list(set(SKIP_LIST_COMPONENTS) - set(args.run_components))

    datasystem = load_data(
        key_sets=args.data_asset_sets,
        config=plugin.config,
        root_path=Path(args.root_path),
        runsheet_path=Path(args.runsheet_path),
    )

    vp = plugin.protocol.validate(
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
    print(f"{'QUEUED CHECK COMPONENT TREE WITH INDIVIDUAL CHECKS':*^60}")
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
        ["--".join(list(row.values())) for row in df.loc[df["code_level"] >= max_flag_code][["description","message"]].to_dict(orient="records")]
    )
    assert (
        df["code_level"].max() < max_flag_code
    ), f"Maximum flag code exceeded: {max_flag_code}. Printing flag messages that caused this halt: {flagged_messages}"


if __name__ == "__main__":
    main()
