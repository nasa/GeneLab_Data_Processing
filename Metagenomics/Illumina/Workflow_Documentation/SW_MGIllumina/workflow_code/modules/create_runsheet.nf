#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.GLDS_accession = "OSD-574"
params.RawFilePattern = null // Pattern of files on OSDR for the OSD accession you want to process 

process GET_RUNSHEET {

    beforeScript "chmod +x ${baseDir}/bin/create_runsheet.sh" 

    input:
        val(GLDS_accession)
    output:
        path("a_*metagenomic*.txt"), emit: assay_TABLE
        path("*.zip"), emit: zip
        path("GLfile.csv"), emit: input_file
        path("versions.txt"), emit: version
    script:
        """
        # Download ISA zip file for the GLDS_accession then unzip it
        GL-download-GLDS-data -g ${GLDS_accession} -p ISA -f && unzip *-ISA.zip

        if [ ${params.RawFilePattern} == null ];then
        
            # Attempt to download the sequences using the assay table, if that fails then
            # attempt retrieving all fastq.gz files
            GL-download-GLDS-data -f -g ${GLDS_accession} -a a_*metagenomic*.txt -o Raw_Sequence_Data || \\
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
        
        # Create runsheet from the assay table
        create_runsheet.sh a_*metagenomic*.txt > GLfile.csv
        GL-version | grep "GeneLab utils"| sed -E 's/^\\s+//' > versions.txt
        """
}


workflow {

    GET_RUNSHEET(params.GLDS_accession)
    file_ch = GET_RUNSHEET.out.input_file
                     .splitCsv(header:true)

}
