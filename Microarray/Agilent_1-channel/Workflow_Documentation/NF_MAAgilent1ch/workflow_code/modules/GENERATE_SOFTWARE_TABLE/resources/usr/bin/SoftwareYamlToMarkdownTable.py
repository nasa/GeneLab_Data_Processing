#! /usr/bin/env python
from pathlib import Path

import yaml
import click
import pandas as pd


AGILENT_SOFTWARE_DPPD = [
    "R",
    "DT",
    "dplyr",
    "stringr",
    "R.utils",
    "limma",
    "glue",
    "ggplot2",
    "biomaRt",
    "matrixStats",
    "statmod",
    "dp_tools",
    "singularity",
    "Quarto",
]

AGILENT_SOFTWARE_DPPD = [s.lower() for s in AGILENT_SOFTWARE_DPPD]

ASSUMED_SOFTWARE = [{
    "name": "singularity",
    "version": 3.9,
    "homepage": "https://sylabs.io"
}]


## Used when the R library metadata doesn't encode any URLS
HOMEPAGE_MAP = {
    "statmod":"https://cran.r-project.org/web/packages/statmod/index.html",
    "biomaRt":"https://bioconductor.org/packages/3.14/bioc/html/biomaRt.html", # UPDATE ON biomaRt version update
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
    df = df.loc[df["name"].str.lower().isin(AGILENT_SOFTWARE_DPPD)]

    assert len(AGILENT_SOFTWARE_DPPD) == len(df), f"Not all software accounted for! Missing: {set(AGILENT_SOFTWARE_DPPD) - set(df['name'].str.lower())}"

    print(df.apply(lambda row: print(row) , axis="columns"))
    df['homepage'] = df.apply(lambda row: HOMEPAGE_MAP[row['name']] if row['homepage'] == "NO URLS ENCODED" else row['homepage'], axis="columns")

    print(df[['name','version','homepage']])

    df = df.rename({"name":"Program","version":"Version","homepage":"Relevant Links"}, axis="columns")

    with open("software_versions.md", "w") as f:
        f.write(df[["Program","Version","Relevant Links"]].to_markdown(index = False))
    
    
if __name__ == '__main__':
    yamlToMarkdown()

