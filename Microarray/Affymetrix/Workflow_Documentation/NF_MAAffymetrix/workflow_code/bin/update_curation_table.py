#! /usr/bin/env python
""" Validation/Verification for raw reads in RNASeq Concensus Pipeline
"""
import argparse
from pathlib import Path

import pandas as pd

from dp_tools.core.post_processing import update_curation_tables
from dp_tools.core.loaders import load_data
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

    parser.add_argument("--isa-path", required=True, help="ISA Archive Path")

    args = parser.parse_args()
    return args


def main(root_dir: Path, runsheet_path: Path, plug_in_dir: Path, isa_path: Path):
    plugin = load_plugin(Path(plug_in_dir))
    
    ds = load_data(
        config=plugin.config,
        root_path=(root_dir),
        runsheet_path=runsheet_path,
    )

    # inject ISA data asset
    class ISA_ASSET:
        path = isa_path
        config = {
            "resource categories": {"publish to repo": False}
        }

    ds.dataset.data_assets["ISA Archive"] = ISA_ASSET()
    isa_path = ds.dataset.data_assets["ISA Archive"].path
    
    update_curation_tables(ds.dataset, config=plugin.config)


if __name__ == "__main__":
    import logging

    logging.basicConfig(level=logging.DEBUG)
    log = logging.getLogger(__name__)
    args = _parse_args()
    main(Path(args.root_path), 
         runsheet_path=Path(args.runsheet_path), 
         plug_in_dir=args.plug_in_dir,
         isa_path=args.isa_path
         )
