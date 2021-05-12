# GeneLab bioinformatics processing protocol for Illumina amplicon sequencing data

> **The document [`GL-DPPD-7104-A.md`](GL-DPPD-7104-A.md) holds an overview and some example commands of how GeneLab processes Illumina amplicon datasets. Exact processing commands for specific datasets that have been released is available in the [GLDS_Processing_Scripts](GLDS_Processing_Scripts) sub-directory and is also provided with their processed data in the [GeneLab Data Systems (GLDS) repository](https://genelab-data.ndc.nasa.gov/genelab/projects).**  

**Developed and maintained by:**  
Michael D. Lee (Mike.Lee@nasa.gov)

---

# General Info
The processing is implemented as a [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow and utilizes [conda](https://docs.conda.io/en/latest/) environments. If helpful, an introduction to conda with installation help and links to other resources can be found [here at Happy Belly Bioinformatics](https://astrobiomike.github.io/unix/conda-intro).

Once we have conda, we can install snakemake like so:

```bash
conda install -c conda-forge -c bioconda -c defaults snakemake
```

All files required for utilizing the workflow are in the [workflow-template](workflow-template) directory. Copying the address of that directory and pasting it into [GitZip here](http://kinolien.github.io/gitzip/) will conveniently allow us to download just that directory. We can then modify the variables in the [config.yaml](workflow-template/config.yaml) file as needed, like pointing to the directory where our starting reads are located, or specifying what our primers are, and then run the workflow. 

**Example command to run workflow**
```bash
snakemake --use-conda --conda-prefix ${CONDA_PREFIX}/envs -j 2 -p
```

* `--use-conda` – this specifies to use the conda environments included in the workflow
* `--conda-prefix` – this allows us to point to where the needed conda environments should be stored. Including this means if we use the workflow on a different dataset somewhere else in the future, it will re-use the same conda environments rather than make new ones. The value listed here, `${CONDA_PREFIX}/envs`, is the default location for conda environments (the variable `${CONDA_PREFIX}` will be expanded to the appropriate location on whichever system it is run on).
* `-j` – this lets us set how many jobs Snakemake should run concurrently (keep in mind that many of the thread and cpu parameters set in the config.yaml file will be multiplied by this)
* `-p` – specifies to print out each command being run to the screen

See `snakemake -h` for more options and details.
