import string
from pathlib import Path
import logging
import enum
from typing import Union
import itertools

log = logging.getLogger(__name__)

import pandas as pd

from dp_tools.core.check_model import FlagCode, FlagEntry, FlagEntryWithOutliers

class GroupFormatting(enum.Enum):
    r_make_names = enum.auto()
    ampersand_join = enum.auto()

def utils_runsheet_to_expected_groups(
    runsheet: Path,
    formatting: GroupFormatting = GroupFormatting.ampersand_join,
    limit_to_samples: list = None,
    map_to_lists: bool = False,
) -> Union[dict[str, str], dict[str, list[str]]]:
    df_rs = (
        pd.read_csv(runsheet, index_col="Sample Name", dtype=str)
        .filter(regex="^Factor Value\[.*\]")
        .sort_index()
    )  # using only Factor Value columns

    if limit_to_samples:
        df_rs = df_rs.filter(items=limit_to_samples, axis="rows")

    match formatting:
        case GroupFormatting.r_make_names:
            expected_conditions_based_on_runsheet = (
                df_rs.apply(lambda x: "...".join(x), axis="columns")
                .apply(r_style_make_names)  # join factors with '...'
                .to_dict()
            )  # reformat entire group in the R style
        case GroupFormatting.ampersand_join:
            expected_conditions_based_on_runsheet = df_rs.apply(
                lambda x: f"({' & '.join(x)})", axis="columns"
            ).to_dict()
        case _:
            raise ValueError(
                f"Formatting method invalid, must be one of the following: {list(GroupFormatting)}"
            )

    # convert from {sample: group} dict
    #         to {group: [samples]} dict
    if map_to_lists:
        unique_groups = set(expected_conditions_based_on_runsheet.values())
        reformatted_dict: dict[str, list[str]] = dict()
        for query_group in unique_groups:
            reformatted_dict[query_group] = [
                sample
                for sample, group in expected_conditions_based_on_runsheet.items()
                if group == query_group
            ]
        expected_conditions_based_on_runsheet: dict[str, list[str]] = reformatted_dict

    return expected_conditions_based_on_runsheet

## Dataframe and Series specific helper functions
def nonNull(df: pd.DataFrame) -> bool:
    # negation since it checks if any are null
    return ~df.isnull().any(axis=None)


def nonNegative(df: pd.DataFrame) -> bool:
    """This ignores null values, use nonNull to validate that condition"""
    return ((df >= 0) | (df.isnull())).all(axis=None)


def onlyAllowedValues(df: pd.DataFrame, allowed_values: list) -> bool:
    """This ignores null values, use nonNull to validate that condition"""
    return ((df.isin(allowed_values)) | (df.isnull())).all(axis=None)


def utils_common_constraints_on_dataframe(
    df: pd.DataFrame, constraints: tuple[tuple[set, dict], ...]
) -> dict:

    issues: dict[str, list[str]] = {
        "Failed non null constraint": list(),
        "Failed non negative constraint": list(),
    }

    for (col_set, col_constraints) in constraints:
        # this will avoid overriding the original constraints dictionary
        # which is likely used in the check message
        col_constraints = col_constraints.copy()

        # limit to only columns of interest
        query_df = df[col_set]
        for (colname, colseries) in query_df.iteritems():
            # check non null constraint
            if col_constraints.pop("nonNull", False) and nonNull(colseries) == False:
                issues["Failed non null constraint"].append(colname)
            # check non negative constraint
            if (
                col_constraints.pop("nonNegative", False)
                and nonNegative(colseries) == False
            ):
                issues["Failed non negative constraint"].append(colname)
            # check allowed values constraint
            if allowedValues := col_constraints.pop("allowedValues", False):
                if onlyAllowedValues(colseries, allowedValues) == False:
                    issues["Failed non negative constraint"].append(colname)

            # raise exception if there are unhandled constraint keys
            if col_constraints:
                raise ValueError(f"Unhandled constraint types: {col_constraints}")

    return issues

def check_dge_table_group_statistical_columns_constraints(
    dge_table: Path, runsheet: Path, **_
) -> FlagEntry:
    expected_groups = utils_runsheet_to_expected_groups(runsheet, map_to_lists=True)
    expected_comparisons = [
        "v".join(paired_groups)
        for paired_groups in itertools.permutations(expected_groups, 2)
    ]

    resolved_constraints = (
        ({f"Log2fc_{comp}" for comp in expected_comparisons}, {"nonNull": True}),
        ({f"T.stat_{comp}" for comp in expected_comparisons}, {"nonNull": True}),
        # can be removed from analysis before p-value and adj-p-value assessed
        # ref: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#why-are-some-p-values-set-to-na
        (
            {f"P.value_{comp}" for comp in expected_comparisons},
            {"nonNegative": True, "nonNull": False},
        ),
        (
            {f"Adj.p.value_{comp}" for comp in expected_comparisons},
            {"nonNegative": True, "nonNull": False},
        ),
    )

    df_dge = pd.read_csv(dge_table)

    # issue trackers
    # here: {prefix+constraint: [failed_columns]}
    issues: dict[str, list[str]] = dict()

    issues = utils_common_constraints_on_dataframe(df_dge, resolved_constraints)

    # check logic
    if not any([issue_type for issue_type in issues.values()]):
        code = FlagCode.GREEN
        message = f"All values in columns met constraint: {resolved_constraints}"
    else:
        code = FlagCode.HALT
        message = (
            f"Issues found {issues} that" f"fail the contraint: {resolved_constraints}."
        )
    return {"code": code, "message": message}

def enclose_in_parens(s: str) -> str:
    """Recreates R's make.names function for individual strings.
    This function is often used to create syntactically valid names in R which are then saved in R outputs.
    Source: https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/make.names

    Args:
        s (str): A string to convert

    Returns:
        str: A string enclosed in parentheses
    """
    return f"({s})"

def r_style_make_names(s: str) -> str:
    """Recreates R's make.names function for individual strings.
    This function is often used to create syntactically valid names in R which are then saved in R outputs.
    Source: https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/make.names

    Args:
        s (str): A string to convert

    Returns:
        str: A string converted in the same way as R's make.names function
    """
    EXTRA_WHITELIST_CHARACTERS = "_ΩπϴλθijkuΑΒΓΔΕΖΗΘΙΚΛΜΝΞΟΠΡΣΤΥΦΧΨΩαβγδεζηθικλμνξοπρστυφχψω_µ" # Note: there are two "μμ" like characters one is greek letter mu, the other is the micro sign
    VALID_CHARACTERS = string.ascii_letters + string.digits + "." + EXTRA_WHITELIST_CHARACTERS
    REPLACEMENT_CHAR = "."
    new_string_chars = list()
    for char in s:
        if char in VALID_CHARACTERS:
            new_string_chars.append(char)
        else:
            new_string_chars.append(REPLACEMENT_CHAR)
    return "".join(new_string_chars)

def check_if_valid_extensions(file, valid_ext):
    """description
    params (description + data type)
    return (description + data type)"""
    pass
def check_factor():
    """description"""
    pass

def check_if_extensions_valid(file, valid_ext):
    """ Description of the function
    Description Line 2

    :param file: Input raw data file
    :type file: Path
    :param valid_ext: Extensions that are allow for the raw data files
    :type valid_ext: list[str]
    :return: A required fields-only flag entry dictionary
    :rtype: FlagEntry
    """
    pass # Does nothing (... is equivalent)

def check_factor_values_in_runsheet():
    """ Description of the function
    Description Line 2

    :param file: Input raw data file
    :type file: Path
    :param valid_ext: Extensions that are allow for the raw data files
    :type valid_ext: list[str]
    :return: A required fields-only flag entry dictionary
    :rtype: FlagEntry
    """
    pass # Does nothing (... is equivalent)

def check_sample_table_against_runsheet(
    runsheet: Path, sampleTable: Path, all_samples_required: bool
) -> FlagEntry:
    """Check the sample table includes all samples as denoted in the runsheet.

    Args:
        runsheet (Path): csv file used for processing, the index denotes all samples
        sampleTable (Path): csv file that pairs each sample with resolved experimental group (called condition within the table)
        all_samples_required (bool): denotes if all samples must be shared or if a subset of samples from the runsheet is okay.

    Returns:
        FlagEntry: A check result
    """
    # data specific preprocess
    df_rs = pd.read_csv(runsheet, index_col="Sample Name").sort_index()
    df_sample = pd.read_csv(sampleTable, index_col="sample").sort_index()

    extra_samples: dict[str, set[str]] = {
        "unique_to_runsheet": set(df_rs.index) - set(df_sample.index),
        "unique_to_sampleTable": set(df_sample.index) - set(df_rs.index),
    }

    # check logic
    if any(
        [
            (extra_samples["unique_to_runsheet"] and all_samples_required),
            (extra_samples["unique_to_sampleTable"]),
        ]
    ):
        code = FlagCode.HALT
        message = f"Samples mismatched: {[f'{entry}:{v}' for entry, v  in extra_samples.items() if v]}"
    else:
        code = FlagCode.GREEN
        message = f"All samples accounted for based on runsheet (All samples required?: {all_samples_required})"
    return {"code": code, "message": message}

def check_dge_table_fixed_statistical_columns_constraints(
    dge_table: Path, **_
) -> FlagEntry:
    # data specific preprocess
    fixed_stats_columns = (
        ({"All.mean", "All.stdev"}, {"nonNull": True, "nonNegative": True}),
        ({"F.p.value"}, {"nonNull": False, "nonNegative": True}),
    )

    df_dge = pd.read_csv(dge_table)

    # issue trackers
    # here: {prefix+constraint: [failed_columns]}
    issues: dict[str, list[str]] = dict()

    issues = utils_common_constraints_on_dataframe(df_dge, fixed_stats_columns)

    # check logic
    if not any([issue_type for issue_type in issues.values()]):
        code = FlagCode.GREEN
        message = f"All values in columns met constraint: {fixed_stats_columns}"
    else:
        code = FlagCode.HALT
        message = (
            f"Issues found {issues} that" f"fail the constraint: {fixed_stats_columns}."
        )
    return {"code": code, "message": message}

def check_dge_table_comparison_statistical_columns_exist(
    dge_table: Path, runsheet: Path, **_
) -> FlagEntry:
    # data specific preprocess
    COMPARISON_PREFIXES = ["Log2fc_", "T.stat_", "P.value_", "Adj.p.value_"]
    expected_groups = utils_runsheet_to_expected_groups(runsheet, map_to_lists=True)
    expected_comparisons = [
        "v".join(paired_groups)
        for paired_groups in itertools.permutations(expected_groups, 2)
    ]
    expected_columns = {
        "".join(comb)
        for comb in itertools.product(COMPARISON_PREFIXES, expected_comparisons)
    }
    df_dge_columns = set(pd.read_csv(dge_table).columns)
    missing_cols = expected_columns - df_dge_columns

    # check logic
    if not missing_cols:
        code = FlagCode.GREEN
        message = f"All comparison summary statistic columns (Prefixes: {COMPARISON_PREFIXES}) present. {sorted(list(expected_columns))}"
    else:
        code = FlagCode.HALT
        message = f"Missing these comparison summary statistic columns (Prefixes: {COMPARISON_PREFIXES}): {sorted(list(missing_cols))}"
    return {"code": code, "message": message}

def check_dge_table_fixed_statistical_columns_exist(dge_table: Path, **_) -> FlagEntry:
    # data specific preprocess
    fixed_stats_columns = {
        "All.mean": {"nonNull": True, "nonNegative": True},
        "All.stdev": {"nonNull": True, "nonNegative": True},
        "F.p.value": {"nonNull": False, "nonNegative": True},
    }
    expected_columns = set(fixed_stats_columns)
    df_dge_columns = set(pd.read_csv(dge_table).columns)
    missing_cols = expected_columns - df_dge_columns

    # check logic
    if not missing_cols:
        code = FlagCode.GREEN
        message = f"All dataset summary stat columns present. {sorted(list(expected_columns))}"
    else:
        code = FlagCode.HALT
        message = (
            f"Missing these dataset summary stat columns: {sorted(list(missing_cols))}"
        )
    return {"code": code, "message": message}

def check_sample_table_for_correct_group_assignments(
    runsheet: Path, sampleTable: Path
) -> FlagEntry:
    """Check the sample table is assigned to the correct experimental group.
    An experimental group is defined by the Factor Value columns found in the runsheet.

    Args:
        runsheet (Path): csv file used for processing, includes metadata used for experimental group designation
        sampleTable (Path): csv file that pairs each sample with resolved experimental group (called condition within the table)

    Returns:
        FlagEntry: A check result
    """
    df_sample = pd.read_csv(sampleTable, index_col=0).set_index("sample", ).sort_index()
    # data specific preprocess
    df_rs = (
        pd.read_csv(runsheet, index_col="Sample Name", dtype=str) # Ensure no factor value columns are misinterpreted as numeric
        .filter(regex="^Factor Value\[.*\]")
        .loc[df_sample.index]  # ensure only sampleTable groups are checked
        .sort_index()
    )  # using only Factor Value columns

    # TODO: refactor with utils_runsheet_to_expected_groups
    expected_conditions_based_on_runsheet = df_rs.apply(
        lambda x: " & ".join(x), axis="columns"
    ).apply(  # join factors with '...'
        enclose_in_parens
    )  # reformat entire group in the R style

    mismatched_rows = expected_conditions_based_on_runsheet != df_sample["group"]

    # check logic
    if not any(mismatched_rows):
        code = FlagCode.GREEN
        message = f"Conditions are formatted and assigned correctly based on runsheet for all {len(df_sample)} samples in sample table: {list(df_sample.index)}"
    else:
        code = FlagCode.HALT
        mismatch_description = (
            df_sample[mismatched_rows]["group"]
            + " <--SAMPLETABLE : RUNSHEET--> "
            + expected_conditions_based_on_runsheet[mismatched_rows]
        ).to_dict()
        message = f"Mismatch in expected conditions based on runsheet for these rows: {mismatch_description}"
    return {"code": code, "message": message}