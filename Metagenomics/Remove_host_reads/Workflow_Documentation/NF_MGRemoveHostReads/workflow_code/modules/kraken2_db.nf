process KRAKEN2_DB {
    tag "Creating host reads database in ${params.ref_dbs_dir}"
    storeDir "${params.ref_dbs_dir}"
        
    input:
    val(host_name)
    val(host_fasta)

    output:
    path("kraken2-${host_name}-db/"), emit: build
    path("versions.txt"), emit: version, optional: true

    script:
    if (params.db_url)
        """
        echo "Downloading and unpacking database from ${params.db_url}"
        wget -O kraken2-${host_name}-db.tar.gz --timeout=3600 --tries=0 --continue ${params.db_url}

        tar -zxvf kraken2-${host_name}-db.tar.gz

        # Cleaning up
        [ -f  kraken2-${host_name}-db.tar.gz ] && rm -rf kraken2-${host_name}-db.tar.gz            
        """
    else if (host_fasta)
        """
        echo "Attempting to build a custom ${host_name} reference database from ${host_fasta}"

        # install taxonomy
        k2 download-taxonomy --db kraken2-${host_name}-db/

        # handle fasta file if it is a URL or local file
        if [[ "${host_fasta}" == http* ]]; then
            # extract filename from URL, strip .gz suffix
            fasta_filename=\$(basename "${host_fasta}" .gz)
            wget -q ${host_fasta} -O \${fasta_filename}.gz
            gunzip \${fasta_filename}.gz
        else
            # check if local file is gzipped
            if [[ "${host_fasta}" == *.gz ]]; then
                fasta_filename=\$(basename "${host_fasta}" .gz)
                gunzip -c ${host_fasta} > \${fasta_filename}
            else
                # local uncompressed can be used directly
                fasta_filename=\$(basename "${host_fasta}")
            fi
        fi

        # add sequence to database's genomic library
        k2 add-to-library --db kraken2-${host_name}-db/ --threads ${task.cpus} \
                          --files \${fasta_filename} --no-masking

        # build the kraken2 database
        k2 build --db kraken2-${host_name}-db/ --threads ${task.cpus} \
                 --kmer-len 35 --minimizer-len 31

        # remove intermediate files
        k2 clean --db kraken2-${host_name}-db/ --log clean.log

        echo "Kraken2 \$(kraken2 -version | head -n 1 | awk '{print \$3}')" >> versions.txt            
        """
    else if (host_name)
        """
        echo "Download and build kraken reference for named host: ${host_name}"

        # download genomic sequences
        k2 download-library --db kraken2-${host_name}-db/ --threads ${task.cpus} \
                            --library ${host_name} --no-masking

        # install taxonomy
        k2 download-taxonomy --db kraken2-${host_name}-db/

        # build the kraken2 database
        k2 build --db kraken2-${host_name}-db/ --threads ${task.cpus} \
                 --kmer-len 35 --minimizer-len 31

        # copy log files needed to extract relevant genomic info
        cp kraken2-${host_name}-db/library/human/assembly_summary.txt kraken2-${host_name}-db/
        cp kraken2-${host_name}-db/library/human/manifest.txt kraken2-${host_name}-db/

        # remove intermediate files
        k2 clean --db kraken2-${host_name}-db/ --log clean.log
            
        echo "Kraken2 \$(kraken2 -version | head -n 1 | awk '{print \$3}')" >> versions.txt            
        """
    else
        error "Input error, host_name, host_fasta, and db_url are all missing. Please supply at least one valid parameter for database creation"

}