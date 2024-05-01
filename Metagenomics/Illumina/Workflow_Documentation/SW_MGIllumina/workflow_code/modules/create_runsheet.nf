#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.GLDS_accession = "OSS-466"

process GET_RUNSHEET {

    beforeScript "chmod +x ${baseDir}/bin/create_runsheet.py"

    output:
        path("*_runsheet.csv"), emit: runsheet
        path("*.zip"), emit: zip
        path("GLparams_file.csv"), emit: params_file
        path("GLfile.csv"), emit: input_file

    script:
        """
        create_runsheet.py --OSD ${params.GLDS_accession}
        """
}


workflow {

    GET_RUNSHEET()
    file_ch = GET_RUNSHEET.out.input_file
                     .splitCsv()

     params_ch = GET_RUNSHEET.out.params_file
                     .splitCsv(header:true)

}
