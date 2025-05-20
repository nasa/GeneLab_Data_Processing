process BUILD_STAR_INDEX {
  // Builds STAR index, this is ercc-spike-in, organism, read length and ensembl version specific
  tag "Refs: ${ genome_fasta }, ${ genome_gtf }, Source: ${reference_source}${reference_source.toLowerCase().contains('ensembl') ? ', Version: ' + reference_version : ''}, MaxReadLength: ${ max_read_length }${ params.genome_subsample ? ', GenomeSubsample: ' + params.genome_subsample : ''}"
  storeDir "${ derived_store_path }/STAR_Indices/${ reference_source }/${reference_source.toLowerCase().contains('ensembl') ? reference_version + '/' : ''}${ meta.organism_sci }/RL-${ max_read_length.toInteger() }"

  input:
    val(derived_store_path)
    val(organism_sci)
    val(reference_source)
    val(reference_version)
    tuple path(genome_fasta), path(genome_gtf)
    val(meta)
    val(max_read_length) // Based on fastQC report for all samples


  output:
    path("${ genome_fasta.baseName }_RL-${ max_read_length.toInteger() }"), emit: index_dir
    path("${ genome_fasta.baseName }_RL-${ max_read_length.toInteger() }/genomeParameters.txt") // Check for completion, only successful builds should generate this file, this is required as the process error is NOT currently used to raised an exception in the python wrapper.

  script:
    """
    #################################################
    # Determine genomeSAindexNbases as per manual:
    # 
    # Manual Excerpt: 
    # genomeSAindexNbases         14
    # int: length (bases) of the SA pre-indexing string. Typically between 10 and 15. Longer strings will use much more memory, but allow faster searches. 
    # For small genomes, the parameter --genomeSAindexNbases must be scaled down to min(14, log2(GenomeLength)/2 - 1). 
    #################################################
    min() {
        printf "%s\n" "\${@:2}" | sort "\$1" | head -n1
    }

    # Get bases in the fasta file
    # We substract the number of lines to remove newline characters from the count
    NUM_BASES=\$(( \$(grep -v '>' ${ genome_fasta } | wc -c) - \$(grep -v '>' ${ genome_fasta } | wc -l)))
    echo NUM_BASES=\$NUM_BASES

    # Compute parameter using formula: min(14, log2(GenomeLength)/2 - 1)
    COMPUTED_GenomeSAindexNbases=\$(awk -v a=\$NUM_BASES 'BEGIN { print " ", int(((log(a)/log(2))/2)-1) }')
    echo "log2(GenomeLength)/2 - 1 = \${COMPUTED_GenomeSAindexNbases}"
    COMPUTED_GenomeSAindexNbases=\$(min -g 14 \$COMPUTED_GenomeSAindexNbases)


    STAR --runThreadN ${task.cpus} \
    --runMode genomeGenerate \
    --limitGenomeGenerateRAM ${ task.memory.toBytes() - 100000000 } \
    --genomeSAindexNbases \$COMPUTED_GenomeSAindexNbases \
    --genomeDir ${ genome_fasta.baseName }_RL-${ max_read_length.toInteger() } \
    --genomeFastaFiles ${ genome_fasta } \
    --sjdbGTFfile ${ genome_gtf } \
    --sjdbOverhang ${ max_read_length.toInteger() - 1 }
    """
}