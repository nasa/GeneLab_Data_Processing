# Reference database info
The database used will depend on the host. The ones that have been created thus far are detailed and available below.


## Mouse ([GRCm39 | GCF_000001635.27](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.27)) database build
This database was built with kraken2 v2.1.1 on 26-Jan-2022.

**Download NCBI taxonomy info needed (takes ~10 minutes):**

```bash
kraken2-build --download-taxonomy --db kraken2-mouse-db/
```

**Downloading mouse reference genome:**

```bash
curl -LO https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz

gunzip GCF_000001635.27_GRCm39_genomic.fna.gz
```


**Adding mouse fasta to database:**

```bash
kraken2-build --add-to-library GCF_000001635.27_GRCm39_genomic.fna.gz --no-masking --db kraken2-mouse-db/
```

**Build the database (takes ~20 minutes as run here):**

```bash
kraken2-build --build --db kraken2-mouse-db/ --threads 30 --no-masking
```

**Remove intermediate files:**

```bash
kraken2-build --clean --db kraken2-mouse-db/
```

### Download mouse kraken2 db


The reference database is ~2.6GB compressed and ~3.8GB uncompressed. It can be downloaded and unpacked with the following:

```bash
curl -L -o kraken2-mouse-db.tar.gz https://figshare.com/ndownloader/files/33900572

tar -xzvf kraken2-mouse-db.tar.gz
```

--- 

