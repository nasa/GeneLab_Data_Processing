#! /usr/bin/env python
import argparse
from pathlib import Path

from dp_tools.core.loaders import load_data
from dp_tools.core.post_processing import generate_md5sum_table
from dp_tools.plugin_api import load_plugin


##############################################################
# Utility Functions To Handle Logging, Config and CLI Arguments
##############################################################
def _parse_args():
    """Parse command line args."""
    parser = argparse.ArgumentParser()

    parser.add_argument("--root-path", required=True, help="Root data path")

    parser.add_argument("--runsheet-path", required=True, help="Runsheet path")

    parser.add_argument("--plug-in-dir", required=True, help="Plugin path")

    args = parser.parse_args()
    return args


def main(root_dir: Path, runsheet_path: Path, plug_in_dir: Path):
    plugin = load_plugin(Path(plug_in_dir))
    
    ds = load_data(
        config=plugin.config,
        root_path=(root_dir),
        runsheet_path=runsheet_path,
    )

    df = generate_md5sum_table(
        ds.dataset,
        config=plugin.config,
        include_tags=True,
    )

    unique_tags = set(df["tags"].sum())
    for tag in unique_tags:
        df_subset = df.loc[df["tags"].apply(lambda l: tag in l)].drop(
            "tags", axis="columns"
        )
        df_subset.to_csv(f"{tag}_md5sum.tsv", sep="\t", index=False)

    # Log missing files
    print(df.columns)
    missing_files = df.loc[df['md5sum'] == "USER MUST ADD MANUALLY!"]["filename"].to_list()
    if missing_files:
        with open("Missing_md5sum_files.txt", "w") as f:
            for missing in missing_files:
                f.write(missing+"\n")

if __name__ == "__main__":
    import logging

    logging.basicConfig(level=logging.DEBUG)
    log = logging.getLogger(__name__)
    args = _parse_args()
    main(Path(args.root_path), runsheet_path=Path(args.runsheet_path), plug_in_dir=args.plug_in_dir)
