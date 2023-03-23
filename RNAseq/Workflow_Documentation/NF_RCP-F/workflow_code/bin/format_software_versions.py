#! /usr/bin/env python
""" This script reformats a software_versions.txt file as output by the
RCP Nextflow Pipeline and outputs as an md table
"""
from __future__ import annotations
import sys
from pathlib import Path

import pandas as pd

PUBLISH_TABLE_ORDER = [
    "Nextflow",
    "FastQC",
    "MultiQC",
    "Cutadapt",
    "TrimGalore!",
    "gtfToGenePred",
    "genePredToBed",
    "STAR",
    "RSEM",
    "R",
    "Bioconductor",
    "DESeq2",
    "tximport",
    # "tidyverse",
    # "STRINGdb",
    # "PANTHER.db",
]

# Certain tools lack a way to print version information and thus cannot be capture dynamically as part of the workflow
# As a workaround, those tool versions should be manually specified here
MANUALLY_DEFINED = [
    {
        "Program": "gtfToGenePred",
        "Version": "377",
        "Relevant Links": "https://anaconda.org/bioconda/ucsc-gtftogenepred",
    },
    {
        "Program": "genePredToBed",
        "Version": "377",
        "Relevant Links": "https://anaconda.org/bioconda/ucsc-genepredtobed",
    },
]

def _parse_samtools_block(text) -> dict:
    """Parses an Samtools version output"""
    for line in text.splitlines():
        if line.startswith("samtools version:"):
            return {
                "Program": "Samtools",
                "Version": line.split(":")[1].strip(),
                "Relevant Links": "https://github.com/samtools/samtools",
            }
    raise ValueError


def _parse_Nextflow_block(text) -> dict:
    """Parses an Nextflow version output"""
    for line in text.splitlines():
        if line.startswith("Nextflow Version:"):
            return {
                "Program": "Nextflow",
                "Version": line.split(":")[1],
                "Relevant Links": "https://github.com/nextflow-io/nextflow",
            }
    raise ValueError


def _parse_rseqc_block(text) -> dict:
    """Parses an RSEQC version output.  Note, these this operates on individual tools version output."""
    results = list()
    for line in text.splitlines():
        if (len(line.split()) == 2) and (".py" in line):
            tool, version = line.split()
            return {
                "Program": tool,
                "Version": version,
                "Relevant Links": "https://sourceforge.net/projects/rseqc",
            }
    raise ValueError("NO VERSIONS SUCCESSFULLY PARSED")


def _parse_RSEM_block(text) -> dict:
    """Parses an RSEM version output"""
    for line in text.splitlines():
        if line.startswith("COUNT_RSEM_version: Current version: RSEM"):
            return {
                "Program": "RSEM",
                "Version": line.split(":")[-1].split()[-1].lstrip("v"),
                "Relevant Links": "https://github.com/deweylab/RSEM",
            }


def _parse_STAR_block(text) -> dict:
    """Parses an STAR version output"""
    for line in text.splitlines():
        if line.startswith("ALIGN_STAR_version:"):
            return {
                "Program": "STAR",
                "Version": line.split()[-1],
                "Relevant Links": "https://github.com/alexdobin/STAR",
            }


def _parse_FastQC_block(text) -> dict:
    """Parses an Fastqc version output"""
    for line in text.splitlines():
        if line.startswith("FastQC v"):
            return {
                "Program": "FastQC",
                "Version": line.split()[-1].lstrip("v"),
                "Relevant Links": "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/",
            }


def _parse_multiqc_block(text) -> dict:
    """Parses an multiqc version output"""
    for line in text.splitlines():
        if line.startswith("multiqc, version"):
            return {
                "Program": "MultiQC",
                "Version": line.split()[-1],
                "Relevant Links": "https://multiqc.info/",
            }


def _parse_Trimgalore_block(text) -> list[dict, dict]:
    """Parses out Trimgalore and Cutadapt versions"""
    lines = text.splitlines()
    for i, line in enumerate(lines):
        if "[powered by Cutadapt]" in line:
            # next line has version
            version = lines[i + 1].split()[-1]
            version_trimgalore = {
                "Program": "TrimGalore!",
                "Version": version,
                "Relevant Links": "https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/",
            }
        if line.startswith("cutadapt version:"):
            version = line.split(":")[-1]
            version_cutadapt = {
                "Program": "Cutadapt",
                "Version": version,
                "Relevant Links": "https://cutadapt.readthedocs.io/en/stable/",
            }
    return [version_trimgalore, version_cutadapt]


R_VERSION_TABLE_DICT = {
    "DESeq2": {
        "table_name": "DESeq2",
        "link": "https://bioconductor.org/packages/release/bioc/html/DESeq2.html",
    },
    "PANTHER.db": {
        "table_name": "PANTHER.db",
        "link": "https://bioconductor.org/packages/release/data/annotation/html/PANTHER.db.html",
    },
    "STRINGdb": {
        "table_name": "STRINGdb",
        "link": "https://www.bioconductor.org/packages/release/bioc/html/STRINGdb.html",
    },
    "tidyverse": {"table_name": "tidyverse", "link": "https://www.tidyverse.org"},
    "tximport": {
        "table_name": "tximport",
        "link": "https://bioconductor.org/packages/release/bioc/html/tximport.html",
    },
    "org.Mm.eg.db": {
        "table_name": "org.Mm.eg.db",
        "link": "https://bioconductor.org/packages/release/data/annotation/html/org.Mm.eg.db.html",
    },
    "org.Rn.eg.db": {
        "table_name": "org.Rn.eg.db",
        "link": "https://bioconductor.org/packages/release/data/annotation/html/org.Rn.eg.db.html",
    },
    "org.Dr.eg.db": {
        "table_name": "org.Dr.eg.db",
        "link": "https://bioconductor.org/packages/release/data/annotation/html/org.Dr.eg.db.html",
    },
    "org.Dm.eg.db": {
        "table_name": "org.Dm.eg.db",
        "link": "https://bioconductor.org/packages/release/data/annotation/html/org.Dm.eg.db.html",
    },
    "org.Ce.eg.db": {
        "table_name": "org.Ce.eg.db",
        "link": "https://bioconductor.org/packages/release/data/annotation/html/org.Ce.eg.db.html",
    },
    "org.Sc.sgd.db": {
        "table_name": "org.Sc.sgd.db",
        "link": "https://bioconductor.org/packages/release/data/annotation/html/org.Sc.sgd.db.html",
    },
    "org.At.tair.db": {
        "table_name": "org.At.tair.db",
        "link": "https://bioconductor.org/packages/release/data/annotation/html/org.At.tair.db.html",
    },
    "org.EcK12.eg.db": {
        "table_name": "org.EcK12.eg.db",
        "link": "https://bioconductor.org/packages/release/data/annotation/html/org.EcK12.eg.db.html",
    },
}


def _parse_R_block(
    text: str, filter_to_rename_dict: dict = R_VERSION_TABLE_DICT, filter_to_rename=True
) -> list:
    """Parses out R package versions"""

    def _is_package(token: str):
        """Returns true if this token represents a package"""
        if all(
            (
                len(token.split("_")) == 2,
                "=" not in token,
            )
        ):
            return True
        return False

    lines = text.splitlines()
    versions = list()
    for i, line in enumerate(lines):
        # parse
        #  [1] org.Mm.eg.db_3.12.0         PANTHER.db_1.0.10
        tokens = line.split()
        if not len(tokens) >= 1:
            continue
        if line.startswith("R version"):
            versions.append(
                {
                    "Program": "R",
                    "Version": line.split()[2],
                    "Relevant Links": "https://www.r-project.org",
                }
            )

        if line.startswith("BioC_version_associated_with_R_version"):
            # next line is Bioconductor version
            versions.append(
                {
                    "Program": "Bioconductor",
                    "Version": lines[i + 1].strip(),
                    "Relevant Links": "https://bioconductor.org",
                }
            )

        if tokens[0].startswith("[") and tokens[0].endswith("]"):
            for potential_package_token in tokens[1:]:
                if _is_package(potential_package_token):
                    package, version = potential_package_token.split("_")
                    if package not in filter_to_rename_dict.keys():
                        continue
                    else:
                        versions.append(
                            {
                                "Program": filter_to_rename_dict[package]["table_name"],
                                "Version": version,
                                "Relevant Links": filter_to_rename_dict[package][
                                    "link"
                                ],
                            }
                        )
    return versions


def main(software_versions_path: Path):
    with software_versions_path.open() as f:
        text_blocks = f.read().split("<><><>")
    results = list()
    for text_block in text_blocks:
        if "COUNT_RSEM_version: Current version: RSEM" in text_block:
            results.append(_parse_RSEM_block(text_block))
        elif "FastQC v" in text_block:
            results.append(_parse_FastQC_block(text_block))
        elif "multiqc, version" in text_block:
            results.append(_parse_multiqc_block(text_block))
        elif "ALIGN_STAR_version:" in text_block:
            results.append(_parse_STAR_block(text_block))
        elif "samtools version:" in text_block:
            results.append(_parse_samtools_block(text_block))
        elif "Nextflow Version:" in text_block:
            results.append(_parse_Nextflow_block(text_block))
        elif "R version" in text_block:
            results.extend(_parse_R_block(text_block))
        elif "Quality-/Adapter-/RRBS-/Speciality-Trimming" in text_block:
            results.extend(_parse_Trimgalore_block(text_block))
        elif "RSeQC" in text_block:
            results.append(_parse_rseqc_block(text_block))
        else:
            # raise NotImplementedError(f"Scripts does not know how to parse: {text_block}")
            print(f"WARNING: Script does not know how to parse: {text_block}")
            pass
    # print(results)
    # Add manually defined versions
    results.extend(MANUALLY_DEFINED)

    df = pd.DataFrame(results)
    df = df.set_index(keys="Program")
    # print(df.head())
    # find any software not in conserved list (e.g. organism specific annotation databases)
    extra_software = list(set(df.index).difference(set(PUBLISH_TABLE_ORDER)))
    # sort to ensure reproducible output
    extra_software.sort()
    # and add to end of the table
    PUBLISH_TABLE_ORDER.extend(extra_software)
    # print(PUBLISH_TABLE_ORDER)
    print(df)
    print(PUBLISH_TABLE_ORDER)
    # deduplicate any entries
    df = df[~df.index.duplicated(keep="first")]
    print(df)
    print(df.index)
    df = df.reindex(PUBLISH_TABLE_ORDER)
    print(df.index)
    print(df)
    output_file = software_versions_path.parent / Path("software_versions.md")
    df.to_markdown(output_file, index=True)
    print(f"Wrote {output_file}")


if __name__ == "__main__":
    main(software_versions_path=Path(sys.argv[1]))
