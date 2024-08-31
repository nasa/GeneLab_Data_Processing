#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

//params.GLDS_accession = "OSD-72"
//params.RawFilePattern = null // Pattern of files on OSDR for the OSD accession you want to process 

process GET_RUNSHEET {

    beforeScript "chmod +x ${baseDir}/bin/*" 

    input:
        val(GLDS_accession)
        val(target_region)
    output:
        tuple path("a_*amplicon*.txt"), path("target_assay_table.txt"), path("runsheet.csv"), emit: tables
        path("*.zip"), emit: zip
        path("GLfile.csv"), emit: input_file
        path("GLparams_file.csv"), emit: params_file
        path("versions.txt"), emit: version
    script:
        """
        # Download ISA zip file for the GLDS_accession then unzip it
        GL-download-GLDS-data -g ${GLDS_accession} -p ISA -f && unzip *-ISA.zip

        (head -n1  a_*amplicon*.txt ; \\
         grep "${target_region}" a_*amplicon*.txt) > target_assay_table.txt

        if [ ${params.RawFilePattern} == null ];then
        
            # Attempt to download the sequences using the assay table, if that fails then
            # attempt retrieving all fastq.gz files
            GL-download-GLDS-data -f -g ${GLDS_accession} -a target_assay_table.txt -o Raw_Sequence_Data || \\
            GL-download-GLDS-data -f -g ${GLDS_accession} -p ".fastq.gz" -o Raw_Sequence_Data
        
        else

        
            GL-download-GLDS-data -f -g ${GLDS_accession} -p  ${params.RawFilePattern} -o Raw_Sequence_Data

        fi

       # Handle case where URLs contain the "+" sign and replaces it with %2B       
       if grep -q '+' *wanted-file-download-commands.sh;then
           grep '+' *wanted-file-download-commands.sh | \\
           sort -u | \\
           awk '{gsub(/\\+/,"%2B", \$NF);print}' \\
           > plus_containing_${GLDS_accession}-wanted-file-download-commands.sh
           cat plus_containing_${GLDS_accession}-wanted-file-download-commands.sh | parallel -j $task.cpus
       fi
        
        # Create runsheet, input and parameter files from the target assay table
        create_runsheet.sh target_assay_table.txt > runsheet.csv
        cut -d "," -f1-2 runsheet.csv > GLfile.csv
        cut -d "," -f3- runsheet.csv | uniq > GLparams_file.csv
        
        GL-version 2>&1 | grep "GeneLab utils"| sed -E 's/^\\s+//' > versions.txt
        """
}


workflow {

    GET_RUNSHEET(params.GLDS_accession, params.target_region)
    file_ch = GET_RUNSHEET.out.input_file
                     .splitCsv(header:true)

    params_ch = GET_RUNSHEET.out.params_file
                     .splitCsv(header:true)

}
