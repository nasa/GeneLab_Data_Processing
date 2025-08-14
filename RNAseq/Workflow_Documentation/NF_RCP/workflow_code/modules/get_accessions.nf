// Take the accession in the format 'OSD-###' or 'GLDS-###'
// Query the API to Find all GLDS accessions associated with the OSD accessions
// Create a dictionary with keys 'osd_accession' values 'glds_accession'. 
// Then return the osd_accession and glds_accession.

process GET_ACCESSIONS {
    tag "${accession}"
    
    input:
        val accession
        val api_url

    output:
        path "accessions.txt", emit: accessions_txt

    script:
    """
    get_accessions.py --accession "${accession}" --api_url "${api_url}" > accessions.txt
    """
}