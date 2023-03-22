#! /usr/bin/env python
""" This script generates excel files indicating the files to publish
"""
import sys
from pathlib import Path

import pandas as pd
from openpyxl.utils import get_column_letter

from displayablepaths import DisplayablePath

def infer_sample(path: Path, sample_list: list) -> str:
    """ Infers the sample based on the file path
    """
    inferred_sample = None
    for sample in sample_list:
        if sample in str(path.path):
            if not inferred_sample:
                inferred_sample = sample
            else:
                raise ValueError(f"Inferred another sample: {sample} but already inferred {inferred_sample}")
    if inferred_sample:
        return inferred_sample
    else:
        #print(f"No sample inferred from {str(path.path)}")
        #print("Returning '?'")
        return 'All samples'


def get_samples(runsheet_path: Path) -> list:
    """ Returns a list of sample names from a given runsheet
    """
    df = pd.read_csv(runsheet_path)
    sample_col = [col for col in df.columns if "sample" in col.lower()][0]
    return list(df[sample_col])


def main(root_path, runsheet_path, template, outputDir: Path = None):
    if template == "RNASeq":
        samples = get_samples(runsheet_path)
        #print(samples)
        paths = DisplayablePath.make_tree(root_path)
        pandas_rows = list()
        for path in paths:
            if path.path.is_file():
                #print(path.path)
                sub_root_dir = path.path.parents[0]
                try:
                    #print(f"{infer_sample(path, samples)}, {path.path.parents[len(path.path.parents)-3].name}, {path.path.name}")
                    
                    parent_dir = path.path.relative_to(root_path).parents
                    print(path.path)
                    print(list(parent_dir))
                    print(list(parent_dir)[0])
                    pandas_rows.append({'Sample Name': infer_sample(path, samples),
                                        'ParentDir':str(list(parent_dir)[0]),
                                        'FileName':path.path.name})
                except IndexError: #files in the root dir GLDS-373/software_versions.txt
                    #print(f"{infer_sample(path, samples)}, {'ROOT'}, {path.path.name}")
                    pandas_rows.append({'Sample Name': infer_sample(path, samples),
                                        'ParentDir':'ROOT',
                                        'FileName':path.path.name})
            else:
                #print(f"Skipping directory: {path}")
                ...

        # create table
        df = pd.DataFrame(pandas_rows).set_index(keys="Sample Name")
        df.to_csv("test.out")

        def categorize_excel_file(fname: str) -> str:
            """ """
            if any(
                    (
                      (fname.endswith("_raw.fastq.gz")),
                      (fname == "raw_multiqc_report.zip"),
                    )
                  ):
                return "raw"
            elif any(
                      (
                        # trimmed files
                        (fname.endswith("_trimming_report.txt")),
                        (fname.endswith("_trimmed.fastq.gz")),
                        (fname == "trimmed_multiqc_report.zip"),
                        (fname == "trimming_multiqc_report.zip"),

                        # alignment files
                        (fname.endswith("_Aligned.sortedByCoord.out.bam")),
                        (fname.endswith("_Aligned.toTranscriptome.out.bam")),
                        (fname.endswith("_SJ.out.tab")),

                        (fname.endswith("_Log.final.out")),
                        (fname == "align_multiqc_report.zip"),

                        # rseqc files
                        (fname == "rseqc-genebody_coverage_multiqc_report.zip"),
                        (fname == "rseqc-infer_experiment_multiqc_report.zip"),
                        (fname == "rseqc-inner_distance_multiqc_report.zip"),
                        (fname == "rseqc-read_distribution_multiqc_report.zip"),
                        (fname.endswith(".geneBodyCoverage.txt")),
                        (fname.endswith(".inner_distance_freq.txt")),
                        (fname.endswith("_read_dist.out")),
                        (fname.endswith("_infer_expt.out")),

                        # raw counts files
                        (fname.endswith(".genes.results")),
                        (fname.endswith(".isoforms.results")),
                        (fname == "Unnormalized_Counts.csv"),

                        # Normalized Counts
                        (fname == "Normalized_Counts.csv"),
                        (fname == "ERCC_Normalized_Counts.csv"),
                        (fname == "count_multiqc_report.zip"),

                        # DGE Files
                        (fname == "contrasts.csv"),
                        (fname == "differential_expression.csv"),
                        (fname == "visualization_output_table.csv"),
                        (fname == "visualization_PCA_table.csv"),
                        (fname == "ERCCnorm_contrasts.csv"),
                        (fname == "ERCCnorm_differential_expression.csv"),
                        (fname == "visualization_output_table_ERCCnorm.csv"),
                        (fname == "visualization_PCA_table_ERCCnorm.csv"),
                        (fname == "SampleTable.csv"),
                      )
                    ):
                return "processed"
            return "unpublished"

        def categorize_excel_sheet(fname: str) -> str:
            """ """
            if any(
                    (
                      (fname.endswith("_raw.fastq.gz")),
                      (fname == "raw_multiqc_report.zip"),
                    )
                  ):
                return "raw"
            elif any(
                      (
                        (fname.endswith("_trimming_report.txt")),
                        (fname.endswith("_trimmed.fastq.gz")),
                        (fname == "trimmed_multiqc_report.zip"),
                        (fname == "trimming_multiqc_report.zip"),
                      )
                    ):
                return "Trimmed_Files"
            elif any(
                      (
                        (fname.endswith("_Aligned.sortedByCoord.out.bam")),
                        (fname.endswith("_Aligned.toTranscriptome.out.bam")),
                        (fname.endswith("_SJ.out.tab")),

                        (fname.endswith("_Log.final.out")),
                        (fname == "align_multiqc_report.zip"),
                      )
                    ):
                return "Alignment_Files"
            elif any(
                      (
                        (fname.endswith(".genes.results")),
                        (fname.endswith(".isoforms.results")),
                        (fname == "Unnormalized_Counts.csv"),
                      )
                    ):
                return "Raw_Counts_Files"
            elif any(
                      (
                        (fname == "Normalized_Counts.csv"),
                        (fname == "ERCC_Normalized_Counts.csv"),
                        (fname == "count_multiqc_report.zip"),
                      )
                    ):
                return "Normalized_Counts_Files"
            elif any(
                      (
                        # DGE Files
                        (fname == "contrasts.csv"),
                        (fname == "differential_expression.csv"),
                        (fname == "visualization_output_table.csv"),
                        (fname == "visualization_PCA_table.csv"),
                        (fname == "ERCCnorm_contrasts.csv"),
                        (fname == "ERCCnorm_differential_expression.csv"),
                        (fname == "visualization_output_table_ERCCnorm.csv"),
                        (fname == "visualization_PCA_table_ERCCnorm.csv"),
                        (fname == "SampleTable.csv"),
                      )
                    ):
                return "DGE_Files"
            elif any(
                      (
                        # rseqc files
                        (fname == "rseqc-genebody_coverage_multiqc_report.zip"),
                        (fname == "rseqc-infer_experiment_multiqc_report.zip"),
                        (fname == "rseqc-inner_distance_multiqc_report.zip"),
                        (fname == "rseqc-read_distribution_multiqc_report.zip"),
                        (fname.endswith(".geneBodyCoverage.txt")),
                        (fname.endswith(".inner_distance_freq.txt")),
                        (fname.endswith("_read_dist.out")),
                        (fname.endswith("_infer_expt.out")),
                      )
                    ):
                return "RSeQC_Analysis"
            return "Other"

        def categorize_excel_column(fname: str) -> str:
            """ """
            if fname.endswith("_raw.fastq.gz"):
                return "Raw Fastq Files"
            elif fname.endswith("_trimmed.fastq.gz"):
                return "Trimmed Fastq Files"
            elif fname.endswith("_multiqc_report.zip"):
                return "MultiQC Files"
            elif any(
                      (
                        (fname.endswith("_trimming_report.txt")),
                      )
                    ):
                return "Trimmed Reports"
            elif any(
                      (
                        (fname.endswith("_Aligned.sortedByCoord.out.bam")),
                        (fname.endswith("_Aligned.toTranscriptome.out.bam")),
                        (fname.endswith("_SJ.out.tab")),
                      )
                    ):
                return "Alignment Data"
            elif any(
                      (
                        # rseqc files
                        (fname.endswith(".geneBodyCoverage.txt")),
                      )
                    ):
                return "GeneBody Coverage"
            elif any(
                      (
                        # rseqc files
                        (fname.endswith(".inner_distance_freq.txt")),
                      )
                    ):
                return "Inner Distance"
            elif any(
                      (
                        # rseqc files
                        (fname.endswith("_read_dist.out")),
                      )
                    ):
                return "Read Distribution"
            elif any(
                      (
                        # rseqc files
                        (fname.endswith("_infer_expt.out")),
                      )
                    ):
                return "Infer Experiment"
            elif any(
                      (
                        (fname.endswith("_Log.final.out")),
                      )
                    ):
                return "Alignment Logs"
            elif any(
                      (
                        (fname.endswith(".genes.results")),
                        (fname.endswith(".isoforms.results")),
                        (fname == "Unnormalized_Counts.csv"),
                      )
                    ):
                return "Raw Counts Data"
            elif any(
                      (
                        (fname == "Normalized_Counts.csv"),
                        (fname == "ERCC_Normalized_Counts.csv"),
                      )
                    ):
                return "Normalized Counts Data"
            elif any(
                      (
                        # DGE Files
                        (fname == "contrasts.csv"),
                        (fname == "differential_expression.csv"),
                        (fname == "visualization_output_table.csv"),
                        (fname == "visualization_PCA_table.csv"),
                        (fname == "ERCCnorm_contrasts.csv"),
                        (fname == "ERCCnorm_differential_expression.csv"),
                        (fname == "visualization_output_table_ERCCnorm.csv"),
                        (fname == "visualization_PCA_table_ERCCnorm.csv"),
                        (fname == "SampleTable.csv"),
                      )
                    ):
                return "DGE Data"
            return "Other_Col"

        def expand_all_samples(df):
            df_all_samples = df.loc[df["Sample Name"] == "All samples"]
            df_sample_wise = df.loc[df["Sample Name"] != "All samples"]
            for sample in df_sample_wise["Sample Name"].unique():
                for fname in df_all_samples["FileName"].unique():
                    df_sample_wise = df_sample_wise.append({"Sample Name":sample, "FileName":fname}, ignore_index=True)
            return df_sample_wise


        ## categorize each file
        df_assigned = df.drop(labels="ParentDir",axis=1)
        df_assigned["Sample Name"] = df_assigned.index.to_list()
        df_assigned = expand_all_samples(df_assigned)
        df_assigned["Excel File"] = [categorize_excel_file(fname) for fname in df_assigned["FileName"]]
        df_assigned["Excel Sheet"] = [categorize_excel_sheet(fname) for fname in df_assigned["FileName"]]
        df_assigned["Excel Column"] = [categorize_excel_column(fname) for fname in df_assigned["FileName"]]

        # testing zip for blanks
        # write to excel file
        SPACES_BETWEEN_SAMPLES = 1
        COLUMN_LENGTH_BUFFER = 8
        for excel_file in df_assigned["Excel File"].unique():
            print(excel_file)
            output_file = f"{excel_file}_file_names.xlsx"
            options = {}
            options['strings_to_formulas'] = False
            options['strings_to_urls'] = False
            writer = pd.ExcelWriter(output_file, engine='openpyxl', options=options)
            df_file = df_assigned.loc[df_assigned["Excel File"] == excel_file]
            for SHEETNAME in df_file["Excel Sheet"].unique():
                print(f"Working on sheet name: {SHEETNAME}")
                df_file_sheet = df_file.loc[df_file["Excel Sheet"] == SHEETNAME]
                sample_dfs = list()
                dataset_df = None
                for i, (sample, gb) in enumerate(df_file_sheet.groupby("Sample Name")):
                    columns = dict()
                    for col in gb["Excel Column"].unique():
                        columns[col] = gb.loc[gb["Excel Column"] == col]["FileName"].to_list()
                    # extend columns to match length of longest column
                    max_col = max([len(col) for col in columns.values()])
                    for key, col in columns.items():
                        columns[key] = col + [''] * (max_col - len(col) + SPACES_BETWEEN_SAMPLES)
                    # create sample column
                    columns["Sample Name"] = [sample] +  [''] * (max_col - 1 + SPACES_BETWEEN_SAMPLES)
                    
                    # convert to dataframe
                    sample_dfs.append(pd.DataFrame(columns))
                # collapse into All samples if all sample_dfs are equivalent (excluding sample name)
                if (all(sample_df.drop("Sample Name", axis=1).equals(sample_dfs[0].drop("Sample Name", axis=1)) for sample_df in sample_dfs)):
                    all_samples_df = sample_dfs[0]
                    all_samples_df["Sample Name"] = ["All Samples"] + all_samples_df["Sample Name"].to_list()[1:]
                    sample_dfs = [all_samples_df]
                for sample_df in sample_dfs:
                    dataset_df = sample_df if not isinstance(dataset_df, pd.DataFrame)  else pd.concat([dataset_df, sample_df])
                
                dataset_df = dataset_df.reset_index(drop=True)
                proper_order = ["Sample Name"] + [col for col in dataset_df.columns if col != "Sample Name"]
                dataset_df = dataset_df.reindex(columns=proper_order)
                print(f"Writing sheet {SHEETNAME} in file {output_file}")
                print(dataset_df.head(20))
                dataset_df.to_excel(writer, sheet_name=SHEETNAME, index=False, header=True)
                # autoresize based on items
                worksheet = writer.sheets[SHEETNAME]
                for idx, col in enumerate(dataset_df.columns):
                    max_len = max((dataset_df[col].astype(str).map(len).max()), len(str(col))) + COLUMN_LENGTH_BUFFER # max of either longest column value or column header plus a buffer character
                    worksheet.column_dimensions[get_column_letter(idx+1)].width = max_len
                writer.save()
        return 

if __name__ == '__main__':
    main(root_path = Path(sys.argv[1]), runsheet_path = Path(sys.argv[2]), template = "RNASeq")
