process ALIGN_STAR {
  // Aligns reads against STAR index
  tag "Sample: ${ meta.id }"

  input:
    tuple val( meta ), path( reads )
    path(star_index_dir)

  output:
    path("${ meta.id }/${ meta.id }*"), emit: publishables // used to ensure direct files are available for publishing directive
    path("${ meta.id }/${ meta.id}_Log.final.out"), emit: alignment_logs
    tuple val(meta), path("${ meta.id }/${ meta.id }_Aligned.sortedByCoord.out.bam"), emit: bam_by_coord
    tuple val(meta), path("${ meta.id }/${ meta.id }_Aligned.toTranscriptome.out.bam"), emit: bam_to_transcriptome
    path("${ meta.id }/${ meta.id }_ReadsPerGene.out.tab"), emit: reads_per_gene
    path("${ meta.id }/${ meta.id }_Unmapped.out.*"), emit: bam_unaligned 
    path("versions.yml"), emit: versions

  script:
    """
    STAR --twopassMode Basic \
    --limitBAMsortRAM ${ (task.memory.toBytes() * 0.8).round() } \
    --outFilterType BySJout \
    --outSAMunmapped Within \
    --genomeDir ${ star_index_dir } \
    --outSAMattributes NH HI AS NM MD MC \
    --outFilterMismatchNoverReadLmax 0.04 \
    --outFilterMismatchNmax 999 \
    --outFilterMultimapNmax 20 \
    --alignIntronMin 20 \
    --alignSJoverhangMin 8 \
    ${ meta.paired_end ? "--alignMatesGapMax 1000000" : "" } \
    --alignIntronMax 1000000 \
    --alignSJDBoverhangMin 1 \
    --sjdbScore 1 \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMheaderHD @HD VN:1.4 SO:coordinate \
    --runThreadN ${ task.cpus } \
    --readFilesCommand zcat \
    --quantMode TranscriptomeSAM GeneCounts\
    --outFileNamePrefix '${ meta.id }/${ meta.id }_' \
    --outReadsUnmapped Fastx \
    --readFilesIn ${ reads }

    echo '"${task.process}":' > versions.yml
    echo "    star: \$(STAR --version | sed -e "s/STAR_//g")" >> versions.yml
    """
}