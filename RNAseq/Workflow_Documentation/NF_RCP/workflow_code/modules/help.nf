// modules/help.nf

def showHelp(workflow) {
    log.info """
    ┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅
            ${workflow.manifest.name}
                Version: ${workflow.manifest.version}
    ┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅┅

    Usage Examples:
    ----------------

    **Example 1:** Processing GLDS Datasets Using Genome FASTA and GTF from Ensembl
    ----------------------------------------------------------------------------
    ```
    nextflow run ./main.nf --accession GLDS-194
    ```

    **Example 2:** Processing GLDS Datasets Using Local Genome FASTA and GTF
    ---------------------------------------------------------------------
    ```
    nextflow run ./main.nf --accession GLDS-194 \\
        --reference_fasta </path/to/fasta> --reference_gtf </path/to/gtf>
    ```

    **Example 3:** Processing Other Datasets with a User-Created Runsheet
    -----------------------------------------------------------------
    ```
    nextflow run ./main.nf --runsheet_path </path/to/runsheet>
    ```

    Parameters:
    ------------

    --help
        Show this help message and exit.

    --mode
        Workflow type (default: 'default'). Set to 'microbes' for processing microbes using Bowtie2.

    --accession
        GLDS accession number to process (e.g., 'GLDS-194').

    --runsheet_path
        Path to a local runsheet instead of one automatically generated from a GLDS ISA archive.

    --isa_archive_path
        Path to a local ISA archive instead of retrieving from OSDR.

    --technical_replicates
        Path to a 2-column CSV file grouping duplicate samples. Example row: SampleA1,SampleA

    --reference_table
        URL or path to the reference table CSV file.

    --reference_store_path
        Directory where fetched reference files are downloaded. Default: './References'.

    --derived_store_path
        Directory where derived reference files are saved. Default: './DerivedReferences'.

    --reference_source
        Source of reference files (e.g., 'ensembl', 'ensembl_bacteria', 'ensembl_plants', 'ncbi').

    --reference_version
        Reference version (e.g., '96').

    --reference_fasta
        Path to a local reference FASTA file.

    --reference_gtf
        Path to a local reference GTF file.

    --validate_params
        Enable parameter validation. Default: true.

    --skip_vv
        Skip automated V&V. Default: false.

    --max_flag_code
        Maximum flag code. Default: 80.

    --rseqc_sample_count
        Number of reads to sample for RSeQC. Default: 15,000,000.

    --outdir
        Directory to save output files.

    --publish_dir_mode
        Publish directory mode. Default: 'link'.

    --email
        Email address for notifications.

    --version
        Show pipeline version and exit.

    --stage_local
        Whether to download the raw reads to the local filesystem. Default: true.

    Debug Options:
    --------------

    --limit_samples_to
        Limit processing to the specified number of samples.

    --truncate_to
        Subsample raw reads files to the specified number of reads for each raw reads file.

    --force_single_end
        Force analysis to use single-end processing. For paired-end datasets, only R1 is used. Default: false.

    --genome_subsample
        Subsample the reference genome to a specific chromosome.

    --use_dummy_gene_counts
        Use random gene counts during DESeq2 (for testing purposes). Default: false.

    For any issues, please open a ticket on the GitHub repository: https://github.com/nasa/GeneLab_Data_Processing/issues
    """.stripIndent()
    exit 0
}
