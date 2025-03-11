#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/**************************************************************************************** 
*********************  Adapter and primer trimming with cutadapt *************************
****************************************************************************************/


// cutadapt process for paired-end and single-end data
// this process runs cutadapt. It is only execute only if anchored_primers == true
process CUTADAPT_ANCHORED {

    tag "Trimming off primers for ${sample_id} using cutadapt..."
    beforeScript "chmod +x ${projectDir}/bin/*"

    input:
        tuple val(sample_id), path(reads), val(isPaired)
        tuple val(F_primer), val(R_primer)
    output:
        tuple val(sample_id), path("*${params.primer_trimmed_R1_suffix[-5..-1]}"), val(isPaired), emit: reads
        tuple val(sample_id),  path("${sample_id}-cutadapt.log"), emit: logs
        tuple val(sample_id),  path("${sample_id}-trimmed-counts.tsv"), emit: trim_counts
        path("versions.txt"), emit: version
    script:
    """
    R_primer_comp=`echo ${R_primer} |tr ATGCRYSWKMBVDHN  TACGYRSWMKVBHDN |rev`
    F_primer_comp=`echo ${F_primer} |tr ATGCRYSWKMBVDHN  TACGYRSWMKVBHDN |rev`

    F_linked_primer="^${F_primer}...\${R_primer_comp}"
    R_linked_primer="^${R_primer}...\${F_primer_comp}"
    if [ ${isPaired} == 'true' ]; then
        
        # Command depends on if primers are linked or not
        if [ ${params.primers_linked} == "TRUE" ]; then

            if [ ${params.discard_untrimmed} == "TRUE" ]; then
                cutadapt -a \${F_linked_primer} \\
                         -A \${R_linked_primer} \\
                         -o ${sample_id}${params.primer_trimmed_R1_suffix} \\
                         -p ${sample_id}${params.primer_trimmed_R2_suffix}  \\
                         --discard-untrimmed \\
                         -m ${params.min_cutadapt_len} \\
                         ${reads[0]} ${reads[1]}  > ${sample_id}-cutadapt.log 2>&1
            else
                cutadapt -a \${F_linked_primer} \\
                         -A \${R_linked_primer} \\
                         -o ${sample_id}${params.primer_trimmed_R1_suffix} \\
                         -p ${sample_id}${params.primer_trimmed_R2_suffix} \\
                         -m ${params.min_cutadapt_len} \\
                         ${reads[0]} ${reads[1]} > ${sample_id}-cutadapt.log 2>&1
            fi

        else

            if [ ${params.discard_untrimmed} == "TRUE" ]; then

                cutadapt -g ^${F_primer} \\
                         -G ^${R_primer} \\
                         -o ${sample_id}${params.primer_trimmed_R1_suffix} \\
                         -p ${sample_id}${params.primer_trimmed_R2_suffix} \\
                         --discard-untrimmed \\
                         -m ${params.min_cutadapt_len} \\
                         ${reads[0]} ${reads[1]} > ${sample_id}-cutadapt.log 2>&1
          
            else
 
                cutadapt -g ^${F_primer} \\
                         -G ^${R_primer} \\
                         -o ${sample_id}${params.primer_trimmed_R1_suffix} \\
                         -p ${sample_id}${params.primer_trimmed_R2_suffix} \\
                         -m ${params.min_cutadapt_len} \\
                         ${reads[0]} ${reads[1]} > ${sample_id}-cutadapt.log 2>&1
 
            fi

        fi

        paste <( printf "${sample_id}" ) \\
              <( grep "read pairs processed"  ${sample_id}-cutadapt.log | tr -s " " "\\t" | cut -f 5 | tr -d "," ) \\
              <( grep "Pairs written"  ${sample_id}-cutadapt.log | tr -s " " "\\t" | cut -f 5 | tr -d "," ) \\
              > ${sample_id}-trimmed-counts.tsv
    # Single-end
    else

            # Command depends on if primers are linked or not
            if [ ${params.primers_linked} == "TRUE" ]; then

                if [ ${params.discard_untrimmed} == "TRUE" ]; then

                    cutadapt -a \${F_linked_primer} \\
                             -o ${sample_id}${params.primer_trimmed_R1_suffix} \\
                             --discard-untrimmed \\
                             -m ${params.min_cutadapt_len} \\
                             ${reads[0]} > ${sample_id}-cutadapt.log 2>&1
                else

                    cutadapt -a \${F_linked_primer} \\
                             -o ${sample_id}${params.primer_trimmed_R1_suffix} \\
                             -m ${params.min_cutadapt_len} \\
                             ${reads[0]}> ${sample_id}-cutadapt.log 2>&1
                
                fi

            else

                if [ ${params.discard_untrimmed} == "TRUE" ]; then

                    cutadapt -g ^${F_primer} \\
                             -o ${sample_id}${params.primer_trimmed_R1_suffix} \\
                             --discard-untrimmed \\
                             -m ${params.min_cutadapt_len} \\
                             ${reads[0]} > ${sample_id}-cutadapt.log 2>&1
                
                else

                    cutadapt -g ^${F_primer} \\
                             -o ${sample_id}${params.primer_trimmed_R1_suffix} \\
                             -m ${params.min_cutadapt_len} \\
                             ${reads[0]} > ${sample_id}-cutadapt.log 2>&1
                
                fi

            fi

            paste <( printf "${sample_id}" ) \\
                  <( grep "reads processed" > ${sample_id}-cutadapt.log | tr -s " " "\\t" | cut -f 4 | tr -d "," ) \\
                  <( grep "Reads written" ${sample_id}-cutadapt.log | tr -s " " "\\t" | cut -f 5 | tr -d "," ) \\
                  > ${sample_id}-trimmed-counts.tsv
    fi
    
    VERSION=`cutadapt --version`
    echo "cutadapt \${VERSION}" > versions.txt
    """
}


// cutadapt process for paired-end and single-end data
// this process runs cutadapt. It is only execute only if anchored_primers == false
process CUTADAPT_UNANCHORED {

    tag "Trimming off primers for ${sample_id} using cutadapt..."
    beforeScript "chmod +x ${projectDir}/bin/*"

    input:
        tuple val(sample_id), path(reads), val(isPaired)
        tuple val(F_primer), val(R_primer)
    output:
        tuple val(sample_id), path("*${params.primer_trimmed_R1_suffix[-5..-1]}"), val(isPaired), emit: reads
        tuple val(sample_id),  path("${sample_id}-cutadapt.log"), emit: logs
        tuple val(sample_id),  path("${sample_id}-trimmed-counts.tsv"), emit: trim_counts
        path("versions.txt"), emit: version
    script:
    """
    R_primer_comp=`echo ${R_primer} |tr ATGCRYSWKMBVDHN  TACGYRSWMKVBHDN |rev`
    F_primer_comp=`echo ${F_primer} |tr ATGCRYSWKMBVDHN  TACGYRSWMKVBHDN |rev`

    F_linked_primer="${F_primer}...\${R_primer_comp}"
    R_linked_primer="${R_primer}...\${F_primer_comp}"
    if [ ${isPaired} == 'true' ]; then
        
        # Command depends on if primers are linked or not
        if [ ${params.primers_linked} == "TRUE" ]; then

            if [ ${params.discard_untrimmed} == "TRUE" ]; then
                cutadapt -a \${F_linked_primer} \\
                         -A \${R_linked_primer} \\
                         -o ${sample_id}${params.primer_trimmed_R1_suffix} \\
                         -p ${sample_id}${params.primer_trimmed_R2_suffix}  \\
                         --discard-untrimmed \\
                         -m ${params.min_cutadapt_len} \\
                         ${reads[0]} ${reads[1]}  > ${sample_id}-cutadapt.log 2>&1
            else
                cutadapt -a \${F_linked_primer} \\
                         -A \${R_linked_primer} \\
                         -o ${sample_id}${params.primer_trimmed_R1_suffix} \\
                         -p ${sample_id}${params.primer_trimmed_R2_suffix} \\
                         -m ${params.min_cutadapt_len} \\
                         ${reads[0]} ${reads[1]} > ${sample_id}-cutadapt.log 2>&1
            fi

        else

            if [ ${params.discard_untrimmed} == "TRUE" ]; then

                cutadapt -g ${F_primer} \\
                         -G ${R_primer} \\
                         -o ${sample_id}${params.primer_trimmed_R1_suffix} \\
                         -p ${sample_id}${params.primer_trimmed_R2_suffix} \\
                         --discard-untrimmed \\
                         -m ${params.min_cutadapt_len} \\
                         ${reads[0]} ${reads[1]} > ${sample_id}-cutadapt.log 2>&1
          
            else
 
                cutadapt -g ${F_primer} \\
                         -G ${R_primer} \\
                         -o ${sample_id}${params.primer_trimmed_R1_suffix} \\
                         -p ${sample_id}${params.primer_trimmed_R2_suffix} \\
                         -m ${params.min_cutadapt_len} \\
                         ${reads[0]} ${reads[1]} > ${sample_id}-cutadapt.log 2>&1
 
            fi

        fi

        paste <( printf "${sample_id}" ) \\
              <( grep "read pairs processed"  ${sample_id}-cutadapt.log | tr -s " " "\\t" | cut -f 5 | tr -d "," ) \\
              <( grep "Pairs written"  ${sample_id}-cutadapt.log | tr -s " " "\\t" | cut -f 5 | tr -d "," ) \\
              > ${sample_id}-trimmed-counts.tsv
    # Single-end
    else

            # Command depends on if primers are linked or not
            if [ ${params.primers_linked} == "TRUE" ]; then

                if [ ${params.discard_untrimmed} == "TRUE" ]; then

                    cutadapt -a \${F_linked_primer} \\
                             -o ${sample_id}${params.primer_trimmed_R1_suffix} \\
                             --discard-untrimmed \\
                             -m ${params.min_cutadapt_len} \\
                             ${reads[0]} > ${sample_id}-cutadapt.log 2>&1
                else

                    cutadapt -a \${F_linked_primer} \\
                             -o ${sample_id}${params.primer_trimmed_R1_suffix} \\
                             -m ${params.min_cutadapt_len} \\
                             ${reads[0]}> ${sample_id}-cutadapt.log 2>&1
                
                fi

            else

                if [ ${params.discard_untrimmed} == "TRUE" ]; then

                    cutadapt -g ${F_primer} \\
                             -o ${sample_id}${params.primer_trimmed_R1_suffix} \\
                             --discard-untrimmed \\
                             -m ${params.min_cutadapt_len} \\
                             ${reads[0]} > ${sample_id}-cutadapt.log 2>&1
                
                else

                    cutadapt -g ${F_primer} \\
                             -o ${sample_id}${params.primer_trimmed_R1_suffix} \\
                             -m ${params.min_cutadapt_len} \\
                             ${reads[0]} > ${sample_id}-cutadapt.log 2>&1
                
                fi

            fi

            paste <( printf "${sample_id}" ) \\
                  <( grep "reads processed" > ${sample_id}-cutadapt.log | tr -s " " "\\t" | cut -f 4 | tr -d "," ) \\
                  <( grep "Reads written" ${sample_id}-cutadapt.log | tr -s " " "\\t" | cut -f 5 | tr -d "," ) \\
                  > ${sample_id}-trimmed-counts.tsv
    fi
    
    VERSION=`cutadapt --version`
    echo "cutadapt \${VERSION}" > versions.txt
    """
}



// This process combines the cutadapt logs and summarizes them. 
process COMBINE_CUTADAPT_LOGS_AND_SUMMARIZE {

    tag "Combining the logs generated by cutadapt..."

    input:
        path(counts)
        path(logs)
    output:
        path("${params.output_prefix}cutadapt${params.assay_suffix}.log"), emit: logs
        path("${params.output_prefix}trimmed-read-counts${params.assay_suffix}.tsv"), emit: counts
    script:
        """
        cat ${logs} > ${params.output_prefix}cutadapt${params.assay_suffix}.log
        
        cat <( printf "sample\\traw_reads\\tcutadapt_trimmed\\n" ) \\
            <( cat ${counts} ) > ${params.output_prefix}trimmed-read-counts${params.assay_suffix}.tsv
        """
}

