import string
from pathlib import Path
import logging
import enum
from typing import Union
import itertools
from statistics import mean
import re

from loguru import logger as log
import pandas as pd
import pandera as pa

from dp_tools.core.check_model import FlagCode, FlagEntry, FlagEntryWithOutliers
from dp_tools.core.entity_model import Dataset

class GroupFormatting(enum.Enum):
    r_make_names = enum.auto()
    ampersand_join = enum.auto()
    ampersand_join_and_remove_non_ascii = enum.auto()

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
        case GroupFormatting.ampersand_join_and_remove_non_ascii:
            def remove_non_ascii(s) -> str:
                if isinstance(s, str):
                    new = re.sub(r'[^\x00-\x7F]+', '', s)
                    print(new, s)
                    return new
                else:
                    return s

            df_rs = df_rs.applymap(remove_non_ascii)

            expected_conditions_based_on_runsheet = df_rs.apply(
                lambda x: f"({' & '.join(x)})", axis="columns"
            ).to_dict()
        case _:
            raise ValueError(
                f"Formatting method invalid, must be one of the following: {list(GroupFormatting)}. Supplied: {formatting}"
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

def check_dge_table_group_columns_constraints(
    dge_table: Path, runsheet: Path, samples: set[str], **_
) -> FlagEntry:
    FLOAT_TOLERANCE = (
        0.001  # Percent allowed difference due to float precision differences
    )
    # data specific preprocess
    GROUP_PREFIXES = ["Group.Stdev_", "Group.Mean_"]
    expected_groups = utils_runsheet_to_expected_groups(runsheet, formatting=GroupFormatting.ampersand_join_and_remove_non_ascii, limit_to_samples=samples)
    query_columns = {
        "".join(comb)
        for comb in itertools.product(GROUP_PREFIXES, expected_groups.values())
    }

    expected_group_lists = utils_runsheet_to_expected_groups(
        runsheet, formatting=GroupFormatting.ampersand_join_and_remove_non_ascii, map_to_lists=True, limit_to_samples=samples
    )
    df_dge = pd.read_csv(dge_table)

    # issue trackers
    issues: dict[str, list[str]] = {
        f"mean computation deviates by more than {FLOAT_TOLERANCE} percent": [],
        f"standard deviation deviates by more than {FLOAT_TOLERANCE} percent": [],
    }

    group: str
    sample_set: list[str]
    for group, sample_set in expected_group_lists.items():
        abs_percent_differences = abs(
            (df_dge[f"Group.Mean_{group}"] - df_dge[sample_set].mean(axis="columns"))
            / df_dge[sample_set].mean(axis="columns")
            * 100
        )
        if any(abs_percent_differences > FLOAT_TOLERANCE):
            issues[
                f"mean computation deviates by more than {FLOAT_TOLERANCE} percent"
            ].append(group)

        abs_percent_differences = abs(
            (df_dge[f"Group.Stdev_{group}"] - df_dge[sample_set].std(axis="columns"))
            / df_dge[sample_set].mean(axis="columns")
            * 100
        )
        if any(abs_percent_differences > FLOAT_TOLERANCE):
            issues[
                f"standard deviation deviates by more than {FLOAT_TOLERANCE} percent"
            ].append(group)

    # check logic
    contraint_description = f"Group mean and standard deviations are correctly computed from samplewise normalized counts within a tolerance of {FLOAT_TOLERANCE} percent (to accomodate minor float related differences )"
    if not any([issue_type for issue_type in issues.values()]):
        code = FlagCode.GREEN
        message = f"All values in columns: {query_columns} met constraint: {contraint_description}"
    else:
        code = FlagCode.HALT
        message = (
            f"Issues found {issues} that"
            f"fail the contraint: {contraint_description}."
        )
    return {"code": code, "message": message}

def check_contrasts_table_headers(contrasts_table: Path, runsheet: Path) -> FlagEntry:
    # data specific preprocess
    expected_groups = utils_runsheet_to_expected_groups(runsheet, map_to_lists=True)
    expected_comparisons = [
        "v".join(paired_groups)
        for paired_groups in itertools.permutations(expected_groups, 2)
    ]
    df_contrasts = pd.read_csv(contrasts_table, index_col=0)

    # check logic
    differences = set(expected_comparisons).symmetric_difference(
        set(df_contrasts.columns)
    )
    if not differences:
        code = FlagCode.GREEN
        message = f"Contrasts header includes expected comparisons as determined runsheet Factor Value Columns: {sorted(set(expected_comparisons))}"
    else:
        code = FlagCode.HALT
        message = f"Contrasts header does not match expected comparisons as determined runsheet Factor Value Columns: {sorted(differences)}"
    return {"code": code, "message": message}

def check_metadata_attributes_exist(
    dataset: Dataset, expected_attrs: list[str]
) -> FlagEntry:
    missing_metadata_fields = list(set(expected_attrs) - set(dataset.metadata))

    # check if any missing_metadata_fields are present
    # check logic
    if not missing_metadata_fields:
        code = FlagCode.GREEN
        message = f"All expected metadata keys found: Expected {expected_attrs}, Found {sorted(set(dataset.metadata))}"
    else:
        code = FlagCode.HALT
        message = f"Missing dataset metadata (source from Runsheet): {sorted(missing_metadata_fields)}"
    return {"code": code, "message": message}

def check_raw_intensities_table(
    raw_intensities: Path, samples: list[str]
) -> FlagEntry:
    schema = pa.DataFrameSchema(
        columns = {sample: pa.Column(float, pa.Check.ge(0)) for sample in samples}
    )

    log.trace(schema)

    df = pd.read_csv(raw_intensities)

    try:
        schema.validate(df, lazy=True)
        error_message = None
    except pa.errors.SchemaErrors as err:
        log.trace(err)
        error_message = err.schema_errors
    if error_message is None:
        code = FlagCode.GREEN
        message = (
            f"Table conforms to schema: {repr(schema)}"
        )
    else:
        code = FlagCode.HALT
        message = (
            f"{error_message}"
        )
    return {"code": code, "message": message}

def check_normalized_expression_table(
    normalized_expression: Path, samples: list[str]
) -> FlagEntry:
    schema = pa.DataFrameSchema(
        columns = {sample: pa.Column(float) for sample in samples}
    )

    df = pd.read_csv(normalized_expression)

    try:
        schema.validate(df, lazy=True)
        error_message = None
    except pa.errors.SchemaErrors as err:
        log.trace(err)
        error_message = err.schema_errors
    if error_message is None:
        code = FlagCode.GREEN
        message = (
            f"Table conforms to schema: {repr(schema)}"
        )
    else:
        code = FlagCode.HALT
        message = (
            f"{error_message}"
        )
    return {"code": code, "message": message}

def check_viz_table_columns_constraints(
    dge_table: Path, runsheet: Path, **_
) -> FlagEntry:
    # data specific preprocess
    expected_groups = utils_runsheet_to_expected_groups(runsheet, map_to_lists=True)
    expected_comparisons = [
        "v".join(paired_groups)
        for paired_groups in itertools.permutations(expected_groups, 2)
    ]
    viz_pairwise_columns_constraints = (
        (
            {f"Log2_Adj.p.value_{comp}" for comp in expected_comparisons},
            {"nonNull": False},
        ),
        (
            {f"Sig.1_{comp}" for comp in expected_comparisons},
            {"allowedValues": [False, True], "nonNull": False},
        ),
        (
            {f"Sig.05_{comp}" for comp in expected_comparisons},
            {"allowedValues": [False, True], "nonNull": False},
        ),
        (
            {f"Log2_P.value_{comp}" for comp in expected_comparisons},
            {"nonNegative": False, "nonNull": False},
        ),
        (
            {f"Updown_{comp}" for comp in expected_comparisons},
            {"allowedValues": [1, -1], "nonNull": True},
        ),
    )

    df_viz = pd.read_csv(dge_table)

    # issue trackers
    # here: {prefix+constraint: [failed_columns]}
    issues: dict[str, list[str]] = dict()

    issues = utils_common_constraints_on_dataframe(
        df_viz, viz_pairwise_columns_constraints
    )

    # check logic
    if not any([issue_type for issue_type in issues.values()]):
        code = FlagCode.GREEN
        message = (
            f"All values in columns met constraint: {viz_pairwise_columns_constraints}"
        )
    else:
        code = FlagCode.HALT
        message = (
            f"Issues found {issues} that"
            f"fail the contraint: {viz_pairwise_columns_constraints}."
        )
    return {"code": code, "message": message}

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

def check_dge_table_log2fc_within_reason(
    dge_table: Path, runsheet: Path, **_
) -> FlagEntry:
    """ Note: This function assumes the normalized expression values are log2 transformed
    """
    LOG2FC_CROSS_METHOD_PERCENT_DIFFERENCE_THRESHOLD = 1  # Percent
    PERCENT_ROWS_WITHIN_TOLERANCE = 99.5  # Percent

    # data specific preprocess
    expected_groups = utils_runsheet_to_expected_groups(runsheet, formatting = GroupFormatting.ampersand_join_and_remove_non_ascii, map_to_lists=True)
    expected_comparisons = [
        "v".join(paired_groups)
        for paired_groups in itertools.permutations(expected_groups, 2)
    ]
    df_dge = pd.read_csv(dge_table)

    # Track error messages
    error_list: list[tuple[str,float]] = list()
    for comparison in expected_comparisons:
        query_column = f"Log2fc_{comparison}"
        group1_mean_col = (
            "Group.Mean_" + comparison.split(")v(")[0] + ")"
        )  # Uses parens and adds them back to prevent slicing on 'v' within factor names
        group2_mean_col = "Group.Mean_" + "(" + comparison.split(")v(")[1]
        print(df_dge[group1_mean_col].describe())
        print(df_dge[group2_mean_col].describe())
        computed_log2fc = df_dge[group1_mean_col] - df_dge[group2_mean_col]
        abs_percent_difference = abs(
            ((computed_log2fc - df_dge[query_column]) / df_dge[query_column]) * 100
        )
        percent_within_tolerance = (
            mean(
                abs_percent_difference
                < LOG2FC_CROSS_METHOD_PERCENT_DIFFERENCE_THRESHOLD
            )
            * 100
        )

        if percent_within_tolerance < PERCENT_ROWS_WITHIN_TOLERANCE: # add current query column to error list
            error_list.append((query_column,percent_within_tolerance,f"First index out of tolerance: {abs_percent_difference.idxmin()}"))

    # inplace sort error list for deterministic order
    error_list.sort()
    if error_list == list():
        code = FlagCode.GREEN
        message = f"All log2fc values recomputed and consistent (within {LOG2FC_CROSS_METHOD_PERCENT_DIFFERENCE_THRESHOLD})"
    else:
        code = FlagCode.HALT
        message = f"At least one log2fc values recomputed and is not consistent (within {LOG2FC_CROSS_METHOD_PERCENT_DIFFERENCE_THRESHOLD}). These columns were flagged: {error_list}"

    return {"code": code, "message": message}

def check_dge_table_sample_columns_constraints(
    dge_table: Path, samples: set[str], **_
) -> FlagEntry:
    MINIMUM_COUNT = 0
    # data specific preprocess
    df_dge = pd.read_csv(dge_table)[samples]

    schema = pa.DataFrameSchema({
        sample: pa.Column(float)
        for sample in samples
    })

    try:
        schema.validate(df_dge, lazy=True)
        error_cases = None
        error_data = None
    except pa.errors.SchemaErrors as err:
        error_cases = err.failure_cases
        error_data = err.data
    if error_cases == error_data == None:
        code = FlagCode.GREEN
        message = (
            f"All values in columns: {samples} met constraints: {repr(schema)}"
        )
    else:
        code = FlagCode.HALT
        message = (
            f"{error_cases}:::{error_data}"
        )
    return {"code": code, "message": message}

def check_dge_table_group_statistical_columns_constraints(
    dge_table: Path, runsheet: Path, **_
) -> FlagEntry:
    expected_groups = utils_runsheet_to_expected_groups(runsheet, formatting = GroupFormatting.ampersand_join_and_remove_non_ascii, map_to_lists=True)
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
    expected_groups = utils_runsheet_to_expected_groups(runsheet, formatting = GroupFormatting.ampersand_join_and_remove_non_ascii, map_to_lists=True)
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
    df_sample = pd.read_csv(sampleTable, index_col="sample").sort_index()
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