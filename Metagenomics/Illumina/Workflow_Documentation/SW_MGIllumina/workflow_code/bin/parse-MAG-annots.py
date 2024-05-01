#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description='This script does whatever it needs to do.')

required = parser.add_argument_group('required arguments')

required.add_argument("-i", "--input-tsv", help='no help for you, come back, 2 years!', action="store", required=True)
required.add_argument("-w", "--wanted-things", help="what'd i tell you?", action="store", required=True)
required.add_argument("-M", "--MAG-ID", action="store", required=True)

parser.add_argument("-o", "--output_tsv", help='Default: "out.tsv"', action="store", dest="output_tsv", default="out.tsv")

args = parser.parse_args()

targets_set = set(line.strip() for line in open(args.wanted_things))

out_tab = open(args.output_tsv, "a")

for line in open(args.input_tsv):
    line = line.strip().split("\t")
    if line[2] != "NA":

        # dropping last coding seq # field so matches contig ID
        if line[0].rsplit('_', 1)[0] in targets_set:

            out_tab.write(str(args.MAG_ID) + "\t" + line[2] + "\n")

out_tab.close()
