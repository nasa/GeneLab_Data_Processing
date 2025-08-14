#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/**************************************************************************************** 
*********************  Assembly binning *************************************************
****************************************************************************************/



// This process runs metabat2 for binning contigs.
process METABAT_BINNING {

    tag "Binning ${sample_id}-s contigs with metabat2..."

    input:
        tuple val(sample_id), path(assembly), path(bam) 
    output:
        tuple val(sample_id), path("${sample_id}-metabat-assembly-depth.tsv"), emit: depth
        tuple val(sample_id), path("${sample_id}-bin*"), emit: bins, optional: true
        path("versions.txt"), emit: version
    script:
        """
        # Only running if the assembly produced anything
        if [ -s ${assembly} ]; then

            jgi_summarize_bam_contig_depths \\
                    --outputDepth ${sample_id}-metabat-assembly-depth.tsv \\
                    --percentIdentity 97 \\
                    --minContigLength 1000 \\
                    --minContigDepth 1.0  \\
                    --referenceFasta ${assembly} ${bam} 

            # only running if there are contigs with coverage 
            # information in the coverage file we just generated
            if [ `wc -l ${sample_id}-metabat-assembly-depth.tsv | sed 's/^ *//' | cut -f 1 -d " "` -gt 1 ]; then 

                metabat2  \\
                    --inFile ${assembly} \\
                    --outFile ${sample_id}-bin \\
                    --abdFile ${sample_id}-metabat-assembly-depth.tsv \\
                    -t ${task.cpus}

            else

                printf "\\n\\nThere was no coverage info generated in ${sample_id}-metabat-assembly-depth.tsv, so no binning with metabat was performed.\\n\\n" 

            fi

            # changing extensions from .fa to .fasta to match nt fasta extension elsewhere in GeneLab
            find . -name '${sample_id}*.fa' > ${sample_id}-bin-files.tmp
            
            if [ -s ${sample_id}-bin-files.tmp ]; then
                paste -d " " <( sed 's/^/mv /' ${sample_id}-bin-files.tmp ) \\
                             <( sed 's/.fa/.fasta/' ${sample_id}-bin-files.tmp ) \\
                             > ${sample_id}-rename.tmp
                bash ${sample_id}-rename.tmp
            fi

            rm -rf ${sample_id}-bin-files.tmp ${sample_id}-rename.tmp

        else

            touch ${sample_id}-metabat-assembly-depth.tsv
            printf "Binning not performed because the assembly didn't produce anything.\\n" 
        fi
        echo metabat2 \$(metabat2 --help 2>&1 | head -n 2 | tail -n 1| sed 's/.*\\:\\([0-9]*\\.[0-9]*\\).*/\\1/') > versions.txt
        """
}



workflow binning {

    take:
        assembly_ch
        read_mapping_ch


    main:
        binning_ch = METABAT_BINNING(assembly_ch.join(read_mapping_ch))


    emit:
    binning_results = binning_ch.out.bins

}

