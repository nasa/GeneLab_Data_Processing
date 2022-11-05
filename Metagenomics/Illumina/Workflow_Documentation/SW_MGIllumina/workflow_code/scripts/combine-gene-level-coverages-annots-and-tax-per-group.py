#!/usr/bin/env python

"""
This is an ad hoc script for the corresponding workflow. It combines the sample gene-level coverage files, taxonomy, and KO annotations into one table for each group of samples.
It produces 2 output files, one normalized to coverage-per-million, one not normalized.

Modified from my `bit-GL-combine-KO-and-tax-tables`: https://github.com/AstrobioMike/bioinf_tools
"""

import os
import sys
import argparse
import textwrap
import pandas as pd
from math import isnan
from numpy import NaN

parser = argparse.ArgumentParser(description="This is an ad hoc script for the corresponding workflow. It combines the individual sample gene-level coverage files and KO annotations into one table. \
                                              It produces 2 output files, one normalized to coverage-per-million, one not normalized.")

required = parser.add_argument_group('required arguments')


required.add_argument("input_coverage_tables", metavar="input-coverage-tables", type=str, nargs="+", help="Input gene-level coverage tables (as written, expected to end with extension '.tsv'.")
required.add_argument("-a", "--KO-annotations-file", help="Input KO annotation table", action="store", dest="KO_tab")
required.add_argument("-t", "--taxonomy-file", help="Input taxonomy table", action="store", dest="tax_tab")
required.add_argument("-g", "--group-ID", help="Group ID", action="store", dest="group")
parser.add_argument("-o", "--output-dir", help="Output directory", action="store", dest="output_dir", default="./")

if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(0)

args = parser.parse_args()

################################################################################


def main():

    check_all_inputs_exist(args.input_coverage_tables, args.KO_tab, args.tax_tab)

    input_files, sample_names = setup_input_lists(args.input_coverage_tables, args.group)

    building_df = combine_KO_and_tax_tab(args.KO_tab, args.tax_tab)

    unnormd_combined_tab, normd_combined_tab = process_and_combine_each_coverage_table(input_files, sample_names)

    unnormd_final_tab, normd_final_tab = combine_all(building_df, unnormd_combined_tab, normd_combined_tab)

    # writing out
    unnormd_final_tab.to_csv(args.output_dir + args.group + "-gene-coverages-annotations-and-tax.tsv", index=False, sep="\t", na_rep = "NA")
    normd_final_tab.to_csv(args.output_dir + args.group + "-CPM-normalized-gene-coverages-annotations-and-tax.tsv", index=False, sep="\t", na_rep = "NA")


################################################################################


# setting some colors
tty_colors = {
    'green' : '\033[0;32m%s\033[0m',
    'yellow' : '\033[0;33m%s\033[0m',
    'red' : '\033[0;31m%s\033[0m'
}


### functions ###
def color_text(text, color='green'):
    if sys.stdout.isatty():
        return tty_colors[color] % text
    else:
        return text


def wprint(text):
    """ print wrapper """

    print(textwrap.fill(text, width=80, initial_indent="  ",
          subsequent_indent="  ", break_on_hyphens=False))


def check_all_inputs_exist(input_tables, KO_tab, tax_tab):

    for file in input_tables + [KO_tab] + [tax_tab]:
        if not os.path.exists(file):
            print("")
            wprint(color_text("It seems the specified input file '" + str(file) + "' can't be found.", "yellow"))
            print("\nExiting for now.\n")
            sys.exit(1)


def setup_input_lists(input_tables, group):
    """ setting up input lists for file locations and sample names """

    input_files = []
    sample_names = []

    for sample in input_tables:
        input_files.append(sample)
        sample_names.append(os.path.splitext(os.path.basename(sample))[0].replace("-" + group + "-gene-coverages", ""))

    return(input_files, sample_names)


def combine_KO_and_tax_tab(KO_tab, tax_tab):

    KO_df = pd.read_csv(KO_tab, sep="\t")
    tax_df = pd.read_csv(tax_tab, sep="\t", dtype = {"taxid": pd.Int64Dtype()}) # this is needed to keep the taxids without a decimal while handling those that are NA
    combined_df = pd.merge(KO_df, tax_df)
    return(combined_df)

def process_and_combine_each_coverage_table(input_files, sample_names):
    """ reads in each table, creates combined tables, one normalized to coverage-per-million, one not normalized  """

    normd_tabs = []
    unnormd_tabs = []

    # iterator to access the same input file and sample name
    for i in range(len(input_files)):

        unnormd_tab = pd.read_csv(input_files[i], sep="\t")

        # generating a normalized version
        normd_tab = unnormd_tab.copy()
        normd_tab.coverage = normd_tab.coverage / normd_tab.coverage.sum() * 1000000

        # changing coverage column headers to be sample name
        unnormd_tab.rename(columns = {"coverage":sample_names[i]}, inplace = True)
        normd_tab.rename(columns = {"coverage":sample_names[i]}, inplace = True)

        # adding to lists
        unnormd_tabs.append(unnormd_tab)
        normd_tabs.append(normd_tab)


    # combining tables
    unnormd_combined_tab = pd.concat(unnormd_tabs, axis = 1).T.drop_duplicates().T
    normd_combined_tab = pd.concat(normd_tabs, axis = 1).T.drop_duplicates().T

    return(unnormd_combined_tab, normd_combined_tab)

def combine_all(building_df, unnormd_combined_tab, normd_combined_tab):
    """ combines KO annotations and tax with coverage tables """

    final_unnormd_tab = pd.merge(building_df, unnormd_combined_tab)
    final_normd_tab = pd.merge(building_df, normd_combined_tab)

    return(final_unnormd_tab, final_normd_tab)

if __name__ == "__main__":
    main()
