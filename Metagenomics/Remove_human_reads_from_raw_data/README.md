# GeneLab removal of human reads from metagenomic datasets

> **It is NASA's policy that any human reads are to be removed from metagenomics datasets prior to being hosted in [GeneLab's data repository](https://genelab-data.ndc.nasa.gov/genelab/projects). As such, all metagenomics datasets are screened against a human reference-genome [kraken2](https://github.com/DerrickWood/kraken2/wiki) database. The document [`GL-DPPD-7105-A.md`](GL-DPPD-7105-A.md) holds an overview and some example commands of how GeneLab peforms this.**  

---

# General Info
The processing is implemented as a [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow and utilizes [conda](https://docs.conda.io/en/latest/) environments. If helpful, an introduction to conda with installation help and links to other resources can be found [here at Happy Belly Bioinformatics](https://astrobiomike.github.io/unix/conda-intro).

Once we have conda, we can install snakemake like so:

```bash
conda install -c conda-forge -c bioconda -c defaults snakemake
```

All files required for utilizing the workflow are in the [workflow-template](workflow-template) directory. Copying the address of that directory and pasting it into [GitZip here](http://kinolien.github.io/gitzip/) will conveniently allow us to download just that directory. We can then modify the variables in the [config.yaml](workflow-template/config.yaml) file as needed, like pointing to the directory where our starting reads are located, and then run the workflow. 

A quick example (though involves downloading the ~4GB reference database) can be run with the files included in the [workflow-template](workflow-template) directory after specifying a location for the reference database in the [config.yaml](workflow-template/config.yaml) file.

**Example command to run workflow**
```bash
snakemake --use-conda --conda-prefix ${CONDA_PREFIX}/envs -j 2 -p
```

* `--use-conda` – this specifies to use the conda environments included in the workflow
* `--conda-prefix` – this allows us to point to where the needed conda environments should be stored. Including this means if we use the workflow on a different dataset somewhere else in the future, it will re-use the same conda environments rather than make new ones. The value listed here, `${CONDA_PREFIX}/envs`, is the default location for conda environments (the variable `${CONDA_PREFIX}` will be expanded to the appropriate location on whichever system it is run on).
* `-j` – this lets us set how many jobs Snakemake should run concurrently (keep in mind that many of the thread and cpu parameters set in the config.yaml file will be multiplied by this)
* `-p` – specifies to print out each command being run to the screen

See `snakemake -h` for more options and details.


## Database Info
The database we use was built with kraken2 v2.1.1 as detailed below, and by default will be downloaded to run with the above workflow (it's ~4.3 GB uncompressed). 

---

### Kraken2 human database build

> The following was performed on 29-Nov-2020 with kraken v2.1.1.

**Downloading human reference (takes ~2 minutes as run here):**

```bash
kraken2-build --download-library human --db kraken2-human-db --threads 30 --no-masking
```

**Downloading NCBI taxonomy info needed (takes ~10 minutes):**

```bash
kraken2-build --download-taxonomy --db kraken2-human-db/
```

**Building database (takes ~20 minutes as run here):**

```bash
kraken2-build --build --db kraken2-human-db/ --threads 30
```

**Removing intermediate files:**

```bash
kraken2-build --clean --db kraken2-human-db/
```

---

### Download database as built on 29-Nov-2020
The reference database 3GB compressed and ~4.3GB uncompressed. If using the accompanying Snakemake workflow, it will download and use this reference database. It can also be downloaded and unpacked by itself with the following:

```bash
curl -L -o kraken2-human-db.tar.gz https://ndownloader.figshare.com/files/25627058

tar -xzvf kraken2-human-db.tar.gz
```

---
