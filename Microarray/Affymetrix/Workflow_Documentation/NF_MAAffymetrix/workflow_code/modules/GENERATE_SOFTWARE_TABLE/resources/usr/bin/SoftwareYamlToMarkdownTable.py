#! /usr/bin/env python
from pathlib import Path

import yaml
import click
import pandas as pd


AFFYMETRIX_SOFTWARE_DPPD = [
    "R",
    "DT",
    "dplyr",
    "tibble",
    "stringr",
    "R.utils",
    "oligo",
    "limma",
    "glue",
    "biomaRt",
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
    "biomaRt":"https://bioconductor.org/packages/3.14/bioc/html/biomaRt.html", # UPDATE ON biomaRt version update
    "oligo":"https://www.bioconductor.org/packages/3.14/bioc/html/oligo.html", # UPDATE ON biomaRt version update
}

@click.command()
@click.argument("input_yaml", type=click.Path(exists=True))
def yamlToMarkdown(input_yaml: Path):
    """ Using a software versions """
    with open(input_yaml, "r") as f:
        data = yaml.safe_load(f)

    data.extend(ASSUMED_SOFTWARE)
    df = pd.DataFrame(data)

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
