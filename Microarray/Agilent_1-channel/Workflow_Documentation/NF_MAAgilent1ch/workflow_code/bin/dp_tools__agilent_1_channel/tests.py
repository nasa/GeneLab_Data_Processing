import pytest
from pathlib import Path
import os

import pandas as pd
import hashlib


from dp_tools import plugin_api

from dp_tools.core.loaders import load_data
from dp_tools.scripts.convert import isa_to_runsheet

# set for testing
@pytest.fixture
def root_test_dir():
    """This should be development machine specific, path should be set by env variable for privacy"""
    return Path(os.environ["TEST_ASSETS_DIR"])

@pytest.fixture()
def plugins():
    return plugin_api.load_plugin(Path(__file__).parent)

@pytest.fixture
def glds367_test_dir(root_test_dir):
    return root_test_dir / "GLDS-367"

def pseudo_fingerprint(df):
    return (
        df["code"].apply(lambda v: v.value).sum()
        + df["code"].apply(lambda v: v.value).mean()
        + df.shape[0] * df.shape[1]
    )

@pytest.mark.parametrize(
    "data_asset_keys,components,expected_flag_table_shape,expected_outlier_table_shape,expected_flag_table_fingerprint,expected_halt_flag_count",
    [
        pytest.param(
            ("glds metadata", "processed"),
            ["Metadata", "DGE Metadata", "DGE Output"],
            (35, 7),
            (0, 0),
            1047.2857142857142,
            1, # Just missing 'Array Design REF' column expected for legacy runsheet
            id="Run all checks",
        ),
    ],
)
def test_protocol_GLDS367_Agile1Channel(
    glds367_test_dir,
    data_asset_keys,
    components,
    expected_flag_table_shape,
    expected_outlier_table_shape,
    expected_flag_table_fingerprint,
    expected_halt_flag_count,
    plugins
):

    

    datasystem = load_data(
        key_sets=data_asset_keys,
        config=plugins['agilent_1_channel'].config,
        root_path=(glds367_test_dir),
        runsheet_path=(
            glds367_test_dir / "Metadata/GLDS-367_microarray_v0_runsheet.csv"
        ),
    )

    report = plugins['agilent_1_channel'].protocol.validate(
        datasystem.dataset,
        report_args={"include_skipped": False},
        protocol_args={"run_components": components},
    )

    assert (
        sum(report["flag_table"]["code_level"] >= 80) == expected_halt_flag_count
    ), f"Found more than expected HALT+ flags: {report['flag_table'].loc[report['flag_table']['code_level'] >= 80].to_dict(orient = 'records')}"

    assert (
        report["flag_table"].shape,
        report["outliers"].shape,
        pseudo_fingerprint(report["flag_table"]),
    ) == (
        expected_flag_table_shape,
        expected_outlier_table_shape,
        expected_flag_table_fingerprint,
    )


def test_microarray_glds367_isa_to_runsheet(glds367_test_dir, plugins):
    """This tests validation as it would be run on dataset after demultiplexing"""
    df_runsheet = isa_to_runsheet(
        "GLDS-367", glds367_test_dir / "Metadata" / "GLDS-367_metadata_GLDS-367-ISA.zip", 
        config=plugins['agilent_1_channel'].config,
        schema=plugins['agilent_1_channel'].schemas.runsheet,
    )

    assert df_runsheet.shape == (16, 14)  # 2 factor values
    assert (
        hashlib.sha1(pd.util.hash_pandas_object(df_runsheet).values).hexdigest()
        == "08ab92a374fdf522b825483918bbf8608c0dbb53"
    ), "Hash did not match, the means the contents changed. Manually validation and reset of test hash is in order"

def test_microarray_glds367_isa_to_runsheet_with_inject(glds367_test_dir, plugins):
    """This tests validation as it would be run on dataset after demultiplexing"""
    df_runsheet = isa_to_runsheet(
        "GLDS-367", glds367_test_dir / "Metadata" / "GLDS-367_metadata_GLDS-367-ISA.zip", 
        config=plugins['agilent_1_channel'].config,
        schema=plugins['agilent_1_channel'].schemas.runsheet,
        inject={'biomart_attribute':'agilent_wholegenome_4x44k_v1'}
    )

    assert df_runsheet.shape == (16, 14)  # 2 factor values
    assert (
        hashlib.sha1(pd.util.hash_pandas_object(df_runsheet).values).hexdigest()
        == "4e2b1e747d4b9a65f2e25293cfee7cbf1a5871a9"
    ), "Hash did not match, the means the contents changed. Manually validation and reset of test hash is in order"



def test_key_based_loader_GLDS367_microarray(root_test_dir, plugins):

    ds = load_data(
        key_sets=("glds metadata", "processed"),
        config=plugins['agilent_1_channel'].config,
        root_path=(root_test_dir / "GLDS-367"),
        runsheet_path=(
            root_test_dir / "GLDS-367/Metadata/GLDS-367_microarray_v0_runsheet.csv"
        ),
    )