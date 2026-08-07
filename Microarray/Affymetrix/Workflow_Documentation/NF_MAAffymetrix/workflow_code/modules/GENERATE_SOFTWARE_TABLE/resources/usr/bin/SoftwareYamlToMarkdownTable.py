#! /usr/bin/env python
from pathlib import Path

import yaml
import click
import pandas as pd


AFFYMETRIX_SOFTWARE_DPPD = [
    "R",
    "Bioconductor",
    "DT",
    "dplyr",
    "tibble",
    "stringr",
    "purrr",
    "oligo",
    "limma",
    "glue",
    "matrixStats",
    "statmod",
    "dp_tools",
    "singularity",
    "Quarto",
    "nextflow"
]

AFFYMETRIX_SOFTWARE_DPPD = [s.lower() for s in AFFYMETRIX_SOFTWARE_DPPD]

ASSUMED_SOFTWARE = [{
    "name": "singularity",
    "version": 3.9,
    "homepage": "https://sylabs.io"
}]


## Used when the R library metadata doesn't encode any URLS
HOMEPAGE_MAP = {
    "statmod":"https://cran.r-project.org/web/packages/statmod/index.html",
    "oligo":"https://www.bioconductor.org/packages/3.22/bioc/html/oligo.html", # UPDATE ON bioconductor version update
}

## Used when certain packages are conditionally used, and therefore dropped from the software table
NOT_USED_SENTINEL = "(Not used for this dataset)"

@click.command()
@click.argument("input_yaml", type=click.Path(exists=True))
@click.argument("filename")
@click.argument("skip_de", type=click.BOOL)
def yamlToMarkdown(input_yaml: Path, filename: str, skip_de: bool):
    """ Using a software versions """
    with open(input_yaml, "r") as f:
        data = yaml.safe_load(f)

    data.extend(ASSUMED_SOFTWARE)
    df = pd.DataFrame(data)

    if skip_de:
        AFFYMETRIX_SOFTWARE_DPPD.remove('limma')
        AFFYMETRIX_SOFTWARE_DPPD.remove('statmod')
        AFFYMETRIX_SOFTWARE_DPPD.remove('matrixstats')

    # Drop software explicitly marked as unused for this dataset (e.g. purrr, only invoked on the 3prime-IVT custom-annotation branch) 
    # and remove them from the software table too, so the completeness assert below doesn't demand software that legitimately never ran
    unused_mask = df["version"].astype(str) == NOT_USED_SENTINEL
    unused_software = set(df.loc[unused_mask, "name"].str.lower())
    for name in unused_software:
        if name in AFFYMETRIX_SOFTWARE_DPPD:
            AFFYMETRIX_SOFTWARE_DPPD.remove(name)
    df = df.loc[~unused_mask]

    # Filter to direct software used (i.e. exclude dependencies of the software)
    df = df.loc[df["name"].str.lower().isin(AFFYMETRIX_SOFTWARE_DPPD)]

    assert len(AFFYMETRIX_SOFTWARE_DPPD) == len(df), f"Not all software accounted for! Missing: {set(AFFYMETRIX_SOFTWARE_DPPD) - set(df['name'].str.lower())}"

    print(df.apply(lambda row: print(row) , axis="columns"))
    df['homepage'] = df.apply(lambda row: HOMEPAGE_MAP[row['name']] if row['homepage'] == "NO URLS ENCODED" else row['homepage'], axis="columns")

    print(df[['name','version','homepage']])

    df = df.rename({"name":"Program","version":"Version","homepage":"Relevant Links"}, axis="columns")

    # Sort by program name for deterministic output
    df = df.sort_values("Program")

    with open("software_versions_GLmicroarray.md", "w") as f:
        f.write(df[["Program","Version","Relevant Links"]].to_markdown(index = False))
    
    
if __name__ == '__main__':
    yamlToMarkdown()
