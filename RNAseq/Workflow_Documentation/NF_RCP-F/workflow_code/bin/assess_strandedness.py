#! /usr/bin/env python

# a script to assess the dataset-wise strand selection from the sample-wise output from RSeQC infer_experiment
from pathlib import Path
from typing import Tuple
from statistics import median
import warnings

# DONE: Update this to be in sync with check configuration, as of 04/25/2022
STRANDEDNESS_ASSIGNMENT_THRESHOLD_MAX = 1.00 # values between above this will be assigned strandedness
STRANDEDNESS_ASSIGNMENT_THRESHOLD_MIN = 0.75
AMBIGUOUS_ASSIGNMENT_THRESHOLD_MAX = 0.75 # values above this, but below strandedness assignment will raise an warning
AMBIGUOUS_ASSIGNMENT_THRESHOLD_MIN = 0.60
UNSTRANDEDNESS_ASSIGMENT_THRESHOLD_MAX = 0.60 # values between this minimum and the ambiguous threshold are assigned unstranded
UNSTRANDEDNESS_ASSIGMENT_THRESHOLD_MIN = 0.40 # values between this minimum and the ambiguous threshold are assigned unstranded


def main(root_dir: str):
    results = dict()
    sense_results = list()
    antisense_results = list()
    undetermined_results = list()
    for file in Path(root_dir).glob("*"):
        with open(file, "r") as f:
            contents = f.read()
        result = _get_stranded_tuple(contents)
        undetermined_results.append(result[0])
        antisense_results.append(result[1])
        sense_results.append(result[2])

        results[file.name] = result
    
    # determine average strandedness
    median_sense  =  median(sense_results)
    median_antisense = median(antisense_results)
    median_undetermined = median(undetermined_results)


    if median_sense > median_antisense:
        dominant, value = "sense", median_sense
    elif median_antisense >= median_sense:
        dominant, value = "antisense", median_antisense

    # assess strandedness
    if STRANDEDNESS_ASSIGNMENT_THRESHOLD_MAX > value > STRANDEDNESS_ASSIGNMENT_THRESHOLD_MIN:
        assignment = dominant
    elif AMBIGUOUS_ASSIGNMENT_THRESHOLD_MAX  > value  > AMBIGUOUS_ASSIGNMENT_THRESHOLD_MIN:
        warnings.warn(f"Strandedness assignment is not well defined for this dataset. median sense: {median_sense}, median antisense: {median_antisense}")
        assignment = dominant
    elif UNSTRANDEDNESS_ASSIGMENT_THRESHOLD_MAX  > value > UNSTRANDEDNESS_ASSIGMENT_THRESHOLD_MIN:
        assignment = "unstranded"

    with open("result.txt", "w") as f:
        f.write(f"{assignment}:{value}")
    


def _get_stranded_tuple(text: str) -> Tuple[float,float,float]:
    """ Parses stdout from infer_experiment """
    for line in text.split("\n"):
        if line.startswith("Fraction of reads failed to determine:"):
            undetermined = float(line.split()[-1])
        elif line.startswith('Fraction of reads explained by "1++,1--,2+-,2-+":') or line.startswith('Fraction of reads explained by "++,--":') :
            antisense = float(line.split()[-1])
        elif line.startswith('Fraction of reads explained by "1+-,1-+,2++,2--":') or line.startswith('Fraction of reads explained by "+-,-+":'):
            sense = float(line.split()[-1])
    return (undetermined, sense, antisense)

if __name__ == "__main__":
    main("infer_out")
