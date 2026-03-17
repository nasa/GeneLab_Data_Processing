# Reference database info
The database we use was built with kraken2 v2.1.1 as detailed below, and by default will be downloaded automatically the first time the workflow is run (it's ~4.3 GB uncompressed).

---

## Kraken2 human database build

> The following was performed on 29-Nov-2020 with kraken v2.1.1, which at that time pulled RefSeq human genome reference GCF_000001405.39, GRCh38.p13.

**Download human reference (takes ~2 minutes as run here):**

```bash
kraken2-build --download-library human --db kraken2-human-db --threads 30 --no-masking
```

**Download NCBI taxonomy info needed (takes ~10 minutes):**

```bash
kraken2-build --download-taxonomy --db kraken2-human-db/
```

**Build the database (takes ~20 minutes as run here):**

```bash
kraken2-build --build --db kraken2-human-db/ --threads 30
```

**Remove intermediate files:**

```bash
kraken2-build --clean --db kraken2-human-db/
```

---

## Download database as built on 29-Nov-2020
The reference database is 3GB compressed and ~4.3GB uncompressed. If using the accompanying Snakemake workflow, it will download and use this reference database. It can also be downloaded and unpacked independently by running the following commands:

```bash
curl -L -o kraken2-human-db.tar.gz https://ndownloader.figshare.com/files/25627058

tar -xzvf kraken2-human-db.tar.gz
```

---
