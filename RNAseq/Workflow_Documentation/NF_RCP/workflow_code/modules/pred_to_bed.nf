process PRED_TO_BED {
    // Converts reference genePred into Bed format
    storeDir "${ derived_store_path }/Genome_GTF_BED_Files/${reference_source}/${reference_source.toLowerCase().contains('ensembl') ? reference_version + '/' : ''}${organism_sci}" 
            
    input:
        val(derived_store_path)
        val(organism_sci)
        val(reference_source) // Used for defining storage location 
        val(reference_version) // Used for defining storage location 
        path(genome_pred)

    output:
        path("${ genome_pred.baseName }.bed"), emit: genome_bed
        path("versions.yml")            , emit: versions

    script:
    """
    genePredToBed ${ genome_pred } ${ genome_pred.baseName }.bed
    
    echo '"${task.process}":' > versions.yml
    echo "    genePredToBed: 469" >> versions.yml
    """
}