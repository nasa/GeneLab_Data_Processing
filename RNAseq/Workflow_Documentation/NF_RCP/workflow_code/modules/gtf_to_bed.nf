process GTF_TO_BED {
    // Converts reference gtf into bed 
    storeDir "${ derived_store_path }/Genome_GTF_BED_Files/${reference_source}/${reference_source.toLowerCase().contains('ensembl') ? reference_version + '/' : ''}${organism_sci}/microbes"
    
    input:
        val(derived_store_path)
        val(organism_sci)
        val(reference_source) // Used for defining storage location 
        val(reference_version) // Used for defining storage location 
        path(genome_gtf)

    output:
        path("${ genome_gtf.baseName }.bed"), emit: genome_bed

    script:
    """
    gtf_to_bed.py ${ genome_gtf } ${ genome_gtf.baseName }.bed
    """
}
