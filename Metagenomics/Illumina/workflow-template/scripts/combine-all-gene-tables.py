#!/usr/bin/env python

"""
This is an ad hoc script for the corresponding workflow. It combines and summarizes all the group sample gene-level coverages by KO annotations (1 table) and taxonomy (another table),
including not annotated and not classified, and also produces CPM (coverage-per-million) normalized versions.

Modified from my `bit-GL-combine-KO-and-tax-tables`: https://github.com/AstrobioMike/bioinf_tools
"""

import os
import sys
import argparse
import textwrap
import pandas as pd
from math import isnan
from numpy import NaN

parser = argparse.ArgumentParser(description="This is an ad hoc script for the corresponding workflow. It combines and summarizes \
                                              all the group sample gene-level coverages by KO annotations (1 table) and taxonomy (another table),\
                                              including not annotated and not classified, and also produces CPM (coverage-per-million) normalized versions.")

required = parser.add_argument_group('required arguments')

required.add_argument("input_tables", metavar="input-tables", type=str, nargs="+", help="Input coverage, annotation, and tax tables (as written, expected to end with extension '.tsv'.")
parser.add_argument("-o", "--output-dir", help='Desired output directory (default: "./")', action="store", default="./", dest="output_dir")

if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(0)

args = parser.parse_args()

################################################################################


def main():

    check_all_inputs_exist(args.input_tables)

    input_files, group_names = setup_input_lists(args.input_tables)

    KO_dict, tax_dict = {}, {}

    na_taxids = []

    KO_collapsed_tabs, tax_collapsed_tabs, KO_dict, tax_dict = process_each_table(input_files, KO_dict, tax_dict, na_taxids)

    combined_KO_tab, combined_norm_KO_tab, combined_tax_tab, combined_norm_tax_tab = combine_tabs(KO_collapsed_tabs, tax_collapsed_tabs, KO_dict, tax_dict)

    # writing out tables
    combined_KO_tab.to_csv(args.output_dir + "All-combined-KO-function-coverages.tsv", index=False, sep="\t")
    combined_norm_KO_tab.to_csv(args.output_dir + "All-combined-KO-function-CPM-normalized-coverages.tsv", index=False, sep="\t")

    combined_tax_tab.to_csv(args.output_dir + "All-combined-taxonomy-coverages.tsv", index=False, sep="\t")
    combined_norm_tax_tab.to_csv(args.output_dir + "All-combined-taxonomy-CPM-normalized-coverages.tsv", index=False, sep="\t")


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


def check_all_inputs_exist(input_tables):

    for file in input_tables:
        if not os.path.exists(file):
            print("")
            wprint(color_text("It seems the specified input file '" + str(file) + "' can't be found.", "yellow"))
            print("\nExiting for now.\n")
            sys.exit(1)


def setup_input_lists(input_tables):
    """ setting up input lists for file locations and sample names """

    input_files = []
    group_names = []

    for group in input_tables:
        input_files.append(group)
        group_names.append(os.path.splitext(os.path.basename(group))[0].replace("-gene-coverages-annotations-and-tax", ""))

    return(input_files, group_names)


def add_to_KO_dict(table, KO_dict):
    """ function for building KO mapping dictionary """

    for index, row in table.iterrows():

        if str(row["KO_ID"]).startswith("K"):

            if str(row["KO_ID"]) not in KO_dict:

                KO_dict[row["KO_ID"]] = row["KO_function"]

    return(KO_dict)


def add_to_tax_dict(table, tax_dict, na_taxids):
    """ function for building tax mapping dictionary """

    for index, row in table.iterrows():

        # skipping if not classified
        if not pd.isna(row["taxid"]):

            if not row["taxid"] in tax_dict:
                tax_dict[row["taxid"]] = row[["domain", "phylum", "class", "order", "family", "genus", "species"]].tolist()

                # some taxids have all NA for these ranks (like 1 and 131567), keep track so can sum with not classified
                if len(set(row[["domain", "phylum", "class", "order", "family", "genus", "species"]].tolist())) == 1:
                    na_taxids.append(row["taxid"])


    return(tax_dict, na_taxids)


def get_na_taxids(tax_dict):
    """ some taxids have all NA for these ranks (like 1 and 131567), keeping track of those so can sum together with the "Not classified" row """

    na_taxids = []

    for key in tax_dict:
        if len(set(mock_dict[key])) == 1:
            na_taxids.append(key)

    return(na_taxids)


def process_each_table(input_files, KO_dict, tax_dict, na_taxids):
    """ reads in each table, normalizes coverage values, collapses based on KO annotations """

    KO_collapsed_tabs = []
    tax_collapsed_tabs = []

    for i in range(len(input_files)):

        tab = pd.read_csv(input_files[i], sep="\t", dtype = {'taxid': str})

        # getting sample column names
        sample_cols = tab.columns[11:].tolist()

        # building dictionaries that will hold all KO terms and taxa from all input files
        KO_dict = add_to_KO_dict(tab, KO_dict)
        tax_dict, na_taxids = add_to_tax_dict(tab, tax_dict, na_taxids)

        # making collapsed KO and tax tables
        KO_tab = tab[['KO_ID'] + sample_cols].copy()
            # collapsing based on KO terms
        KO_tab = KO_tab.groupby(by = ['KO_ID'], dropna = False).sum()

        tax_tab = tab[['taxid'] + sample_cols].copy()
          # setting any taxids that are all NA at these standard ranks to "NA" (some have an assigned taxid, but don't have a D/P/C/O/F/G/S taxid, like 1 and 131567)
        tax_tab.replace(na_taxids, NaN, inplace = True)
            # collapsing based on tax
        tax_tab = tax_tab.groupby(by = ['taxid'], dropna = False).sum()

        # appending to lists of tables
        KO_collapsed_tabs.append(KO_tab)
        tax_collapsed_tabs.append(tax_tab)

    return(KO_collapsed_tabs, tax_collapsed_tabs, KO_dict, tax_dict)


def add_KO_functions(tab, KO_dict):
    """ adds KO functions to combined table based on KO_ID and KO_dict object holding mappings """

    KO_functions = []

    for KO in tab.KO_ID:

        if KO in KO_dict:

            KO_functions.append(str(KO_dict[KO]))

        else:

            KO_functions.append("Not annotated")

    tab.insert(1, "KO_function", KO_functions)

    return(tab)

def add_tax_info(tab, tax_dict):
    """ adds lineage info back to combined table based on taxid and tax_dict object holding mappings """

    domain_list, phylum_list, class_list, order_list, family_list, genus_list, species_list = [], [], [], [], [], [], []

    for taxid in tab.taxid:

        if taxid in tax_dict:

            if isinstance(tax_dict[taxid][0], str):
                domain_list.append(tax_dict[taxid][0])
            else:
                domain_list.append("NA")

            if isinstance(tax_dict[taxid][1], str):
                phylum_list.append(tax_dict[taxid][1])
            else:
                phylum_list.append("NA")

            if isinstance(tax_dict[taxid][2], str):
                class_list.append(tax_dict[taxid][2])
            else:
                class_list.append("NA")

            if isinstance(tax_dict[taxid][3], str):
                order_list.append(tax_dict[taxid][3])
            else:
                order_list.append("NA")

            if isinstance(tax_dict[taxid][4], str):
                family_list.append(tax_dict[taxid][4])
            else:
                family_list.append("NA")

            if isinstance(tax_dict[taxid][5], str):
                genus_list.append(tax_dict[taxid][5])
            else:
                genus_list.append("NA")


            if isinstance(tax_dict[taxid][6], str):
                species_list.append(tax_dict[taxid][6])
            else:
                species_list.append("NA")

        else:
            domain_list.append("NA")
            phylum_list.append("NA")
            class_list.append("NA")
            order_list.append("NA")
            family_list.append("NA")
            genus_list.append("NA")
            species_list.append("NA")

    tab.insert(1, "domain", domain_list)
    tab.insert(2, "phylum", phylum_list)
    tab.insert(3, "class", class_list)
    tab.insert(4, "order", order_list)
    tab.insert(5, "family", family_list)
    tab.insert(6, "genus", genus_list)
    tab.insert(7, "species", species_list)

    return(tab)


def combine_tabs(KO_tab_list, tax_tab_list, KO_dict, tax_dict):
    """ combines all KO tables into one and all tax tables into one """

    # combining KO tabs
    KO_combined_tab = pd.concat(KO_tab_list, axis=1).drop_duplicates().fillna(0).sort_index()

    # moving index to be column and changing that NaN to be "Not annotated", and naming column back to KO_ID
    KO_combined_tab = KO_combined_tab.reset_index().fillna("Not annotated")
    KO_combined_tab.rename(columns = {"index":'KO_ID'}, inplace = True)

    # adding KO functions
    KO_combined_tab = add_KO_functions(KO_combined_tab, KO_dict)

    # combining tax tabs
    tax_combined_tab = pd.concat(tax_tab_list, axis=1).drop_duplicates().fillna(0).sort_index()

    # moving index to be column and naming column back to taxid
    tax_combined_tab = tax_combined_tab.reset_index()

    # changing the NaN to be "Not annotated" and naming column back to taxid
    tax_combined_tab['index'] = tax_combined_tab['index'].fillna("Not classified")
    tax_combined_tab.rename(columns = {"index":'taxid'}, inplace = True)

    # adding tax full lineage info
    tax_combined_tab = add_tax_info(tax_combined_tab, tax_dict)

    # making CPM-normalized versions of each
    KO_combined_norm_tab = KO_combined_tab.copy()
    tax_combined_norm_tab = tax_combined_tab.copy()

    # getting sample column names
    sample_cols = KO_combined_norm_tab.columns[2:].tolist()

    # making normalized versions
    for col in sample_cols:
        KO_combined_norm_tab[col] = KO_combined_norm_tab[col] / KO_combined_norm_tab[col].sum() * 1000000
        tax_combined_norm_tab[col] = tax_combined_norm_tab[col] / tax_combined_norm_tab[col].sum() * 1000000

    return(KO_combined_tab, KO_combined_norm_tab, tax_combined_tab, tax_combined_norm_tab)

if __name__ == "__main__":
    main()
