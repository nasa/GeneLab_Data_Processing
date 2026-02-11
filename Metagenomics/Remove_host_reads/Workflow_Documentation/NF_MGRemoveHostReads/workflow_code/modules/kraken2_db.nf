process KRAKEN2_DB {
    tag "Downloading host reads database to ${params.ref_dbs_Dir}"
    publishDir "${params.ref_dbs_Dir}", mode: 'copy'
        
    input:
    tuple val(host), val(host_id), val(fasta_url)

    output:
    path "kraken2-${host_id}-db/"

    script:
    """
    k2 download-taxonomy --db kraken2-${host_id}-db/

    # Download FASTA file and uncompress it
    wget -q ${fasta_url} -O host_assembly.fasta.gz
    gunzip -c host_assembly.fasta.gz > host_assembly.fasta

    kraken2-build --add-to-library host_assembly.fasta --db kraken2-${host_id}-db/ --threads ${task.cpus} --no-masking

    kraken2-build --build --db kraken2-${host_id}-db/ --threads ${task.cpus}

    kraken2-build --clean --db kraken2-${host_id}-db/
            
    """
}