#! /usr/bin/env python

################################################################################
# Adapted from create_table_v2.py
# Sent by Ana Uriarte Acuna, 2021_04_02
################################################################################

import re
import os, sys
import numpy as np
import codecs
import warnings
import argparse
import tempfile
import zipfile

##############################################################
# Utility Functions To Handle CLI Arguments
##############################################################
def _parse_args():
    """Parse command line args."""
    parser = argparse.ArgumentParser(description="Create metasheet for GeneLab dataset")
    parser.add_argument(
        "--accession",
        metavar="GLDS-194",
        required=True,
        help="GeneLab dataset accession number",
    )
    parser.add_argument(
        "--isa-zip",
        metavar="ISA.ZIP",
        required=True,
        help="Path to GeneLab ISA zip folder",
    )
    parser.add_argument(
        "--output-dir",
        metavar="out/Dir",
        required=True,
        help="Path create the metasheet table",
    )
    parser.add_argument(
        "--runsheet",
        metavar="AST*****.csv",
        required=True,
        help="Path runsheet to merge data from",
    )

    args = parser.parse_args()
    return args


##############################################################
# Additional Functions
##############################################################
def _unzip_ISA(isa_zip_path: str) -> str:
    """Unzips ISA and places into a tmp contents folder.
    Returns path to temporary directory holding ISA zip file contents.

    :param isa_zip_path: path to isa zip file
    """
    temp_dir = tempfile.mkdtemp()
    with zipfile.ZipFile(isa_zip_path, "r") as zip_ref:
        zip_ref.extractall(temp_dir)
    return temp_dir


# *******************************************************************************************************
# Function that opens sample and investigation files
# input : none
# output : s_file -> string with all data from sample file
#         i_file -> string with all data from investigation file
#         file_name_s -> string name of sample file
#         file_name_i -> string name of investigation file
#         metadata_directory -> absolute path of directory with metadata files
# ********************************************************************************************************
def open_files(success, isa_zip_path):
    success = -1
    encoding_flag = 0

    try:
        os.system("clear")
        # metadata_directory = input("Metadata directory: ")
        metadata_directory = _unzip_ISA(isa_zip_path)
        if metadata_directory.lower() == "exit":
            print("Goodbye")
            sys.exit()
        # For now thw output directory is the same as the input one
        outdir = metadata_directory
    except Exception as e:
        # print(e, "\nPress enter to try another directory...")
        # non-interactive version
        raise e
    else:
        file_name_s = ""
        file_name_i = ""
        sample_file = ""
        investigation_file = ""
        try:
            for filename in os.listdir(metadata_directory):
                if filename.endswith(".txt"):
                    if filename[0] == "s" and filename[1] == "_":
                        file_name_s = outdir + "/" + filename
                        try:
                            f1 = open(
                                os.path.join(metadata_directory, file_name_s), "r"
                            )
                            sample_file = f1.read()
                            f1.close()
                        except Exception as e:
                            success += 1
                            print(e)
                    if filename[0] == "i" and filename[1] == "_":
                        file_name_i = outdir + "/" + filename
                        try:
                            f1 = open(
                                os.path.join(metadata_directory, file_name_i), "r"
                            )
                            investigation_file = f1.read()
                            f1.close()
                        except Exception as e:
                            success += 1
                            error = str(e)
                            if "can't decode byte" in error:
                                print("ERROR: Need to change file encoding.")
                                encoding_flag = 1
            if file_name_s == "" and encoding_flag == 0:
                print(
                    "ERROR: No sample file found. Please make sure the metadata files are complete and named correctly."
                )
                sys.exit()
            else:
                if sample_file == "" and encoding_flag == 0:
                    print("ERROR: The study file is empty.")
                    sys.exit()
            if file_name_i == "" and encoding_flag == 0:
                print(
                    "ERROR: No investigation file found. Please make sure the metadata files are complete and named correctly."
                )
                sys.exit()
            else:
                if investigation_file == "" and encoding_flag == 0:
                    print("Error: The investigation file is empty.")
                    sys.exit()

        except Exception as e:
            print(e, "FILE PATH DOES NOT EXIST")
            success += 1
            # input()

    return [
        success,
        sample_file,
        investigation_file,
        file_name_s,
        file_name_i,
        metadata_directory,
    ]


# *******************************************************************************
# Funciton that removes extra tabs and line breaks from raw data file
# input : file
# output : clean file
# *******************************************************************************
def clean_file(file):
    str_file = ""
    new_file_data = ""
    new_line = ""

    quote = 0
    for character in file:
        if character == '"':
            quote += 1
        if quote % 2 == 1:
            if character == "\t" or character == "\n":
                new_file_data += " "
            else:
                new_file_data += character
        else:
            new_file_data += character

    return new_file_data


# *******************************************************************************
# Funciton that parses investigation file to get needed information
# For the metadata table we only need assay file name(s) and techonoly type(s)
# input : i_file object
# output : none
# *******************************************************************************
def i_file_parse(i_file):
    for line in i_file.raw_data.split("\n"):
        if "Study Assay File Name" in line:
            num_col = 0
            for column in line.split("\t"):
                if num_col > 0:
                    i_file.assay_file_names.append(column.strip('"'))
                num_col += 1
        if "Study Assay Technology Type" in line:
            num_col = 0
            read = 0
            for column in line.split("\t"):
                if num_col == 0:
                    if "Study Assay Technology Type" == column:
                        read = 1
                    else:
                        read = 0
                else:
                    if read == 1:
                        i_file.assay_tech_types.append(column.strip('"'))
                num_col += 1


# **************************************************************************************************************
# Function that parses sample file and pulls out factors, characteristics that change and parameters that change
# For plant datasets it also pulls out required characteristcs/parameters
# input : s_file object
# output : column_index -> list with index of charactersitcs, factors and parameters to be included in table
#         data_matrix -> matrix with columns of characteristcs, factors and parameters to be included in table
# **************************************************************************************************************
def s_file_parse(s_file):

    # Split the first line by tabs to get  headers and number of lines
    for i, line in enumerate(s_file.raw_data.split("\n")):
        if i == 0:
            s_file.headers = line.split("\t")
    s_file.num_lines = i

    # create temporary lists, only characteristcs and parameters that change will be added to final
    temp_characteristics_idx = []
    temp_parameters_idx = []
    for idx, title in enumerate(s_file.headers):
        # characterists that are always added
        if (
            "organism" in title.lower()
            or "strain" in title.lower()
            or "cell line" in title.lower()
            or "ecotype" in title.lower()
            or "cultivar" in title.lower()
        ):
            s_file.final_characteristics_idx.append(idx)
            # save organism idx to check if plant dataset
            if "organism" in title.lower():
                org_idx = idx
        if "characteristic" in title.lower():
            temp_characteristics_idx.append(idx)
        if "factor" in title.lower():
            s_file.factors_idx.append(idx)
        if "sample name" in title.lower():
            if "original submitted sample name" in title.lower():
                s_file.original_sample_name_idx.append(idx)
            else:
                s_file.sample_name_idx.append(idx)
        if "unit" in title.lower():
            s_file.units_idx.append(idx)
        if "parameter" in title.lower():
            temp_parameters_idx.append(idx)
        # parameters that are always added
        if "exposure duration" in title.lower():
            s_file.final_parameters_idx.append(idx)
        # characteristics that will only be included if organism is plant
        if (
            "hardware" in title.lower()
            or "material type" in title.lower()
            or "age" in title.lower()
            or "genotype" in title.lower()
            or "temperature" in title.lower()
            or "preservation method" in title.lower()
        ):
            s_file.plant_params_idx.append(idx)
        # had to add second condition to light to remove "spacefLIGHT"
        if "light" in title.lower():
            if "space" not in title.lower():
                s_file.plant_params_idx.append(idx)

    # add units for factors
    if len(s_file.factors_idx) > 0:
        for i in range(len(s_file.factors_idx)):
            s_file.final_factors_idx.append(s_file.factors_idx[i])
            for element in s_file.units_idx:
                if s_file.factors_idx[i] == element - 1:
                    s_file.final_factors_idx.append(element)

    # check what characterists change through samples
    for idx in temp_characteristics_idx:
        whole_column = []
        for line in s_file.raw_data.split("\n"):
            for col_idx, column in enumerate(line.split("\t")):
                if col_idx == idx:
                    whole_column.append(column)
        if len(set(whole_column)) > 2:
            # at least 2 different values
            s_file.characteristics_change_idx.append(idx)

    # add units to characterists
    if len(s_file.characteristics_change_idx) > 0:
        for i in range(len(s_file.characteristics_change_idx)):
            s_file.final_characteristics_idx.append(
                s_file.characteristics_change_idx[i]
            )
            for element in s_file.units_idx:
                if s_file.characteristics_change_idx[i] == element - 1:
                    s_file.final_characteristics_idx.append(element)

    # add units to parameters
    # need to do this because exposure duration is added before and it has units
    if len(s_file.final_parameters_idx) > 0:
        for i in range(len(s_file.final_parameters_idx)):
            for element in s_file.units_idx:
                if s_file.final_parameters_idx[i] == element - 1:
                    s_file.final_parameters_idx.append(element)

    # check what parameters change through samples
    for idx in temp_parameters_idx:
        whole_column = []
        for line in s_file.raw_data.split("\n"):
            for col_idx, column in enumerate(line.split("\t")):
                if col_idx == idx:
                    whole_column.append(column)
        if len(set(whole_column)) > 2:
            # at least 2 different values
            s_file.parameters_change_idx.append(idx)

    # add units to parameters that change
    if len(s_file.parameters_change_idx) > 0:
        for i in range(len(s_file.parameters_change_idx)):
            s_file.final_parameters_idx.append(s_file.parameters_change_idx[i])
            for element in s_file.units_idx:
                if s_file.parameters_change_idx[i] == element - 1:
                    s_file.final_parameters_idx.append(element)

    # add units to required paramters for plants
    if len(s_file.plant_params_idx) > 0:
        # print("Factors: ")
        for i in range(len(s_file.plant_params_idx)):
            # print(s_file.factors_idx[i], ":", s_file.factors[i])
            s_file.final_plant_params_idx.append(s_file.plant_params_idx[i])
            for element in s_file.units_idx:
                if s_file.plant_params_idx[i] == element - 1:
                    s_file.final_plant_params_idx.append(element)

    # get organism
    organism_column = []
    for i, line in enumerate(s_file.raw_data.split("\n")):
        for j, col in enumerate(line.split("\t")):
            if i > 0 and j == org_idx:
                col = col.strip('"')
                organism_column.append(col.lower())

    # for now only plant datasets have required params
    if (
        "arabidopsis thaliana" in organism_column
        or "brassica rapa" in organism_column
        or "ceratopteris richardii" in organism_column
        or "chlamydomonas reinhardtii" in organism_column
        or "eruca vesicaria" in organism_column
    ):
        s_file.organism_type = "plant"

    column_index = []
    if s_file.organism_type == "plant":
        column_index.extend(
            s_file.sample_name_idx
            + s_file.original_sample_name_idx
            + s_file.final_characteristics_idx
            + s_file.final_plant_params_idx
            + s_file.final_parameters_idx
            + s_file.final_factors_idx
        )
    else:
        column_index.extend(
            s_file.sample_name_idx
            + s_file.original_sample_name_idx
            + s_file.final_characteristics_idx
            + s_file.final_parameters_idx
            + s_file.final_factors_idx
        )

    # remove repetitions
    column_index = list(dict.fromkeys(column_index))

    # Get raw data into matrix
    data_matrix = []
    all_factors_list = []
    for i, idx in enumerate(column_index):
        whole_column = []
        for line in s_file.raw_data.split("\n"):
            all_factors = ""
            for col_idx, column in enumerate(line.split("\t")):
                if col_idx == idx:
                    whole_column.append(column)
                if col_idx in s_file.final_factors_idx:
                    all_factors += "_" + column.strip('"')
            if i == 0:
                all_factors_list.append(all_factors)
        whole_column = whole_column[1:]
        data_matrix.append(whole_column)

    all_factors_list = all_factors_list[1 : len(all_factors_list) - 1]

    # create list of repetition for each sample
    rep_list = []
    name_dict = {}
    count = 1
    for i, element in enumerate(all_factors_list):
        if element in name_dict:
            count = name_dict[element]
            name_dict[element] = count + 1
        else:
            name_dict[element] = 1
        rep_list.append(name_dict[element])

    s_file.rep_list = rep_list

    return column_index, data_matrix


# **************************************************************************************************************
# Function that opens assay file(s) and gets samples in the file and raw data file names
# input : a_file object
#        i_file object (that has assay file names)
#        metadata_directory (to open files)
# output : none
# **************************************************************************************************************
def a_file_parse(a_file, i_file, metadata_directory):
    for filename, device_type in i_file.assay_map.items():
        if device_type == "":
            warnings.warn("Found 'null' device type, likely due to empty field, skipping")
            continue

        if device_type == "Image analysis":
            warnings.warn("Image analysis data file parsing is not yet supported.")
            continue

        # replace spaces with _ to match file name
        device_type = device_type.replace(" ", "_")

        file_name_a = metadata_directory + "/" + filename
        if file_name_a == "":
            print(
                "ERROR: No assay file found for device type "
                + device_type
                + ". Please make sure the metadata files are complete and named correctly."
            )
            print("The device type must be included in the assay file name.")
            sys.exit()

        f1 = open(
            os.path.join(metadata_directory, file_name_a),
            "r",
            encoding="windows-1252",
        )
        assay_file = f1.read()
        f1.close()
        assay_file = clean_file(assay_file)

        if assay_file == "":
            print("ERROR: The assay file for " + device_type + "is empty.")
            sys.exit()

        # get sample names to check later which samples from sample file are in assay
        names = []
        file_names = []
        idx = -1
        for i, line in enumerate(assay_file.split("\n")):
            # flag to get only first appearance of raw data file
            flag = 0
            for j, column in enumerate(line.split("\t")):
                if i == 0:
                    if (
                        "raw data file" in column.lower()
                        or "array data file" in column.lower()
                        or "raw spectral data file" in column.lower()
                    ):
                        flag += 1
                        if flag == 1:
                            idx = j
                    if (
                        "sample name" in column.lower()
                        and "comment" not in column.lower()
                    ):
                        sample_name_idx = j
                else:
                    if j == idx:
                        file_names.append(column)
                    if j == sample_name_idx:
                        names.append(column)

        a_file.sample_names.append(names)
        a_file.data_file.append(file_names)


# **************************************************************************************************************
# Function that creates metadata table with tsv format
# input : i_file object
#        s_file object
#        a_file object
#        data_columns -> matrix with columns of values to be included
#        column_index -> list with indexes of columns from sample table to be included
#        glds_num
# output : data -> formatted table to be written in file
# **************************************************************************************************************
def create_table(i_file, s_file, a_file, data_columns, column_index, glds_num):

    data = []
    # create header to keep track of number of columns, idxs without line formatting
    header = []
    # list of units to concatenate to value instead of having separate columns
    unit_idx = []

    # First line
    line = ""

    # First column is sample name
    header.append("Sample Name")
    line += "Sample Name\t"

    # Next columns are file names, for headers add file type and assay type between ()
    for assay_type in i_file.assay_tech_types:
        # for nucleotide sequencing, bisulfite sequencing and rna sequencing Raw Data File is included
        if "sequencing" in assay_type.lower():
            header.append("Raw Data File (" + assay_type + ")")
            line += "Raw Data File (" + assay_type + ")\t"
        # for microarray array data file is included
        elif "microarray" in assay_type.lower():
            header.append("Array Data File (" + assay_type + ")")
            line += "Array Data File (" + assay_type + ")\t"
        # for mass spectrometry raw spectral data file is included
        elif "spectrometry" in assay_type.lower():
            header.append("Raw Spectral Data File (" + assay_type + ")")
            line += "Raw Spectral Data File (" + assay_type + ")\t"

    # GLDS #
    header.append("GLDS-#")
    line += "GLDS-#\t"

    # Device type
    header.append("Device Type")
    line += "Device Type\t"

    # Characteristics, factors and parameters
    for i in range(1, len(column_index)):
        # if column corresponds to a unit save it to be concatenated in same column
        if "Unit" in s_file.headers[:][column_index[i]]:
            unit_idx.append(i)
        else:
            header.append(s_file.headers[:][column_index[i]].strip('"'))
            line += s_file.headers[:][column_index[i]].strip('"') + "\t"

    # Repetition number
    header.append("Rep #")
    line += "Rep #\t"

    # Add required headers for plant datasets
    extra_columns = 0
    if s_file.organism_type == "plant":
        required_fields = [
            "Hardware",
            "Material Type",
            "Age",
            "Ecotype",
            "Genotype",
            "Time",
            " Light",
            "Temperature",
            "Treatment type",
            "Dose",
            "Preservation method",
        ]
        for field in required_fields:
            pattern = ".*" + field.lower() + ".*"
            r = re.compile(pattern)

            if not (any(r.match(line.lower()) for line in header)):
                header.append(field)
                line += field + "\t"
                extra_columns += 1

    # End of first line (column headers)
    line += "\n"
    data.append(line)

    # Add values for every column for each line
    for i in range(0, s_file.num_lines - 1):
        line = ""

        for j in range(len(header) + len(unit_idx)):
            # column 1 : sample name
            if j == 0:
                line += data_columns[j][i] + "\t"

            # column 2 -> n = number of device types : file names
            if j > 0 and j <= len(i_file.assay_tech_types):
                if i_file.assay_tech_types[j - 1] == 'Image analysis':
                    warnings.warn("Image analysis data file parsing is not yet supported.")
                    continue
                if i_file.assay_tech_types[j - 1] == '':
                    warnings.warn("Found 'null' device type, likely due to empty field, skipping")
                    continue

                assay_row = -1
                for row, element in enumerate(a_file.sample_names[j - 1]):
                    # check if sample is included in specific assay file
                    if data_columns[0][i] == element:
                        # save row# from assay file in case order is different
                        assay_row = row
                # if sample was found in assay table
                if assay_row >= 0:
                    line += a_file.data_file[j - 1][assay_row] + "\t"
                else:
                    line += "NA\t"

            # column n+1 : glds number
            if j == len(i_file.assay_tech_types) + 1:
                line += "GLDS-" + glds_num + "\t"

            # column n+2 : device type
            if j == len(i_file.assay_tech_types) + 2:
                if i_file.assay_tech_types[j - 1 - 2] == '':
                    warnings.warn("Found 'null' device type, likely due to empty field, skipping")
                    continue
                if i_file.assay_tech_types[j - 1 - 2] == 'Image analysis':
                    warnings.warn("Image analysis data file parsing is not yet supported.")
                    continue

                # flag to see if it's first value added
                flag = 0
                for k in range(len(i_file.assay_tech_types)):
                    if data_columns[0][i] in a_file.sample_names[k]:
                        if flag == 0:
                            line += i_file.assay_tech_types[k]
                            flag += 1
                        else:
                            # if it's not the first value, separate values with comma
                            line += ", " + i_file.assay_tech_types[k]
                line += "\t"

            # column n+3 -> m = len(headers) - extra_columns : characteristcs, parameters and factors
            if j > (len(i_file.assay_tech_types) + 2) and j < (
                len(header) + len(unit_idx) - extra_columns - 1
            ):
                if j - (len(i_file.assay_tech_types) + 2) not in unit_idx:
                    # if the next value in the list is a unit, concatenate values in same column
                    if j + 1 - (len(i_file.assay_tech_types) + 2) in unit_idx:
                        line += (
                            data_columns[j - (len(i_file.assay_tech_types) + 2)][
                                i
                            ].strip('"')
                            + " "
                            + data_columns[j + 1 - (len(i_file.assay_tech_types) + 2)][
                                i
                            ].strip('"')
                            + "\t"
                        )
                    # if there is no unit, just add value
                    else:
                        line += (
                            data_columns[j - (len(i_file.assay_tech_types) + 2)][
                                i
                            ].strip('"')
                            + "\t"
                        )

            # column n+3+m = repetition number
            if j == len(header) + len(unit_idx) - extra_columns - 1:
                line += "Rep " + str(s_file.rep_list[i]) + "\t"

            # column n+4+m -> len(headers) : required characterists for plant datasets
            if j > len(header) + len(unit_idx) - extra_columns - 1:
                # add empty values
                line += " \t"

        # end of line
        line += "\n"
        data.append(line)

    return data


# **************************************************************************************************************
# Function that writes metadata table into metadata directory
# input : table -> formatted table to be written in file
#        metadata_table -> directory where file will be written
#        glds_num -> glds number to be included in file name
# output : none
# **************************************************************************************************************
def write_table(table, metadata_directory, glds_num, output_dir):
    file_name = "GLDS-" + glds_num + "_metadata_table.txt"
    try:
        newf = open(os.path.join(output_dir, file_name), "w+")
        for i in range(len(table)):
            print(table[i])
            newf.write(table[i])
        newf.close()
        print(file_name, " file created.")
        # input()
    except Exception as e:
        print(e)
        print("ERROR: Error writing " + file_name + " file.")
        # input()


# **************************************************************************************************************
# Class for sample file
# attributes : info needed to create metadata table (descriptions below)
# methods : none
# **************************************************************************************************************
class sample_file_class:
    def __init__(self):
        # FILE INFO
        # name of sample file
        self.file_name = []
        # raw data read from file
        self.raw_data = []
        # number of lines in the file
        self.num_lines = 0
        # number of columns in the file
        self.num_cols = 0

        # SAMPLE INFO
        # headers of sample table
        self.headers = []
        # column index of sample names
        self.sample_name_idx = []
        # column index for original sample name
        self.original_sample_name_idx = []
        # type of organism ie plant
        self.organism_type = ""

        # list of column index for units
        self.units_idx = []
        # list of column index of characteristcs that change through samples
        self.characteristics_change_idx = []
        # list of column index for required characteristics, charateristics that change through samples and their units
        self.final_characteristics_idx = []

        # list of column index for factors
        self.factors_idx = []
        # list of column index for factors and their units
        self.final_factors_idx = []

        # list of column index for required parameters, parameters that change through samples and their units
        self.final_parameters_idx = []
        # list of column index for parameters that change through samples
        self.parameters_change_idx = []

        # list of repetition number for each sample
        self.rep_list = []

        # characteristics for plants only
        self.plant_params_idx = []
        self.final_plant_params_idx = []


# **************************************************************************************************************
# Class for sample file
# attributes : info needed to create metadata table (descriptions below)
# methods : none
# **************************************************************************************************************
class assay_file_class:
    def __init__(self):
        # file name and raw data not used since there can be more than one assay file
        self.file_name = []
        self.raw_data = []

        # number of assay files
        self.num_files = []

        # These lists follow the order of tech types from investigation file
        # list of lists containing raw data file names for each assay file
        self.data_file = []

        # list of lists containing sample names for each assay file
        self.sample_names = []


# **************************************************************************************************************
# Class for investigation file
# attributes : info needed to create metadata table (descriptions below)
# methods : none
# **************************************************************************************************************
class investigation_file_class:
    def __init__(self):
        # FILE INFO
        # file name
        self.file_name = []
        # raw data read from file
        self.raw_data = []

        # assay files and device type info
        # list of assay file names
        self.assay_file_names = []
        # list of device types (for each device type there must be an assay file)
        self.assay_tech_types = []

    @property
    def assay_map(self):
        return {assay_filename:assay_type for  assay_filename, assay_type in zip(self.assay_file_names, self.assay_tech_types)}


if __name__ == "__main__":
    args = _parse_args()
    # t1 = time.time()

    # try to open files until successful
    # success = 0
    # while success >= 0:
    # No longer interactive, removing this while loop
    success = 0
    [
        success,
        sample_file,
        investigation_file,
        file_name_s,
        file_name_i,
        metadata_directory,
    ] = open_files(success, args.isa_zip)

    # get glds# from curator (not in ISA files)
    glds_num = args.accession.replace("GLDS-", "")

    # create objects that contain all info of sample, assay and investigation files
    s_file = sample_file_class()
    a_file = assay_file_class()
    i_file = investigation_file_class()

    # parse investigation file to get assay file name and tech type
    i_file.file_name = file_name_i
    i_file.raw_data = clean_file(investigation_file)
    i_file_parse(i_file)

    # parse sample file to get columns to be included in metadata table
    s_file.raw_data = clean_file(sample_file)
    s_file.file_name = file_name_s
    [column_index, data_columns] = s_file_parse(s_file)

    # parse assay file(s) to get raw data file names
    a_file_parse(a_file, i_file, metadata_directory)

    # create table
    table = create_table(i_file, s_file, a_file, data_columns, column_index, glds_num)

    # write table
    write_table(table, metadata_directory, glds_num, output_dir=args.output_dir)

    #############################################################################
    # merge additional data already parsed in runsheet, added by Jonathan Oribello
    #############################################################################

    # load tabular files as dataframes
    import pandas as pd

    df_meta = pd.read_csv(f"{args.accession}_metadata_table.txt", sep="\t")
    df_runsheet = pd.read_csv(args.runsheet)

    # determine columns that need to be merged into meta table
    cols_to_merge_from_runsheet = [
        "Sample Name",
        "organism",
        "paired_end",
        "has_ERCC",
    ] + [col for col in df_runsheet.columns if col.startswith("Factor Value")]

    # preprocess before merging, keep only columns to merge, adjust sample_name to match df_meta as a merge key
    df_runsheet = df_runsheet[cols_to_merge_from_runsheet]

    # merge
    df_merged = df_runsheet.merge(
        df_meta, on="Sample Name", suffixes=[None, "_metatable"]
    )

    # drop redudant columns from metatable (probably factor values)
    df_merged = df_merged.drop(
        columns=[col for col in df_merged.columns if col.endswith("_metatable")]
    )

    # cleanup and "unnamed:" columns as imported from original meta_table
    df_merged = df_merged.drop(
        columns=[col for col in df_merged.columns if col.startswith("Unnamed:")]
    )

    # write back to file
    df_merged.to_csv(f"{args.accession}_metadata_table.txt", sep="\t", index=False)
