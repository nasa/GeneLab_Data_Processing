#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

//params.accession = "GLDS-487"
//params.target_region = "16S"

process GET_RUNSHEET {

    beforeScript "chmod +x ${projectDir}/bin/create_runsheet.py"
    tag "Retrieving raw sequences and metadata for ${accession}..."
    input:
        tuple val(accession), val(target_region)
        
    output:
        path("*_runsheet.csv"), emit: runsheet
        path("*.zip"), emit: zip
        path("GLparams_file.csv"), emit: params_file
        path("GLfile.csv"), emit: input_file
        path("versions.txt"), emit: version

    script:
        """
        create_runsheet.py --OSD ${accession} --target ${target_region}
        GL-version | grep "GeneLab utils"| sed -E 's/^\\s+//' > versions.txt
        echo "dptools v1.3.4" >> versions.txt
        python --version >> versions.txt 
        """
}


workflow {

    values = Channel.of([params.accession, params.target_region])
    GET_RUNSHEET(values)
    file_ch = GET_RUNSHEET.out.input_file
                     .splitCsv(header:true)

     params_ch = GET_RUNSHEET.out.params_file
                     .splitCsv(header:true)

}
