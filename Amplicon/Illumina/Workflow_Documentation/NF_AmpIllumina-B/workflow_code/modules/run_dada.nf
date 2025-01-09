#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


/**************************************************************************************** 
*********************  ASV generation using Dada2 ***************************************
****************************************************************************************/

process RUN_R_TRIM {

    tag "Running dada2 on the trimmed reads..."

    input:
        tuple path(sample_IDs_file), val(isPaired)
        path(trimmed_reads)
        path(trimmed_read_counts)
    output:
        path("Filtered_Sequence_Data/*${params.filtered_R1_suffix[-5..-1]}"), emit: reads
        path("Filtered_Sequence_Data/filtered-read-counts${params.assay_suffix}.tsv"), emit: filtered_count
        path("final_outputs/taxonomy${params.assay_suffix}.tsv"), emit: taxonomy
        path("final_outputs/taxonomy-and-counts${params.assay_suffix}.biom"), emit: biom
        path("final_outputs/ASVs${params.assay_suffix}.fasta"), emit: fasta
        path("final_outputs/read-count-tracking${params.assay_suffix}.tsv"), emit: read_count
        path("final_outputs/counts${params.assay_suffix}.tsv"), emit: counts
        path("final_outputs/taxonomy-and-counts${params.assay_suffix}.tsv"), emit: taxonomy_count
        path("versions.txt"), emit: version
    script:
        
        """
        if [ ${isPaired} == true ];then

            mkdir -p  Trimmed_Sequence_Data/ && \\
            mv ${trimmed_read_counts} \\
                *${params.primer_trimmed_R1_suffix} \\
                *${params.primer_trimmed_R2_suffix} \\
                Trimmed_Sequence_Data/

            mkdir -p Filtered_Sequence_Data/ final_outputs/
            
            TRIM_PRIMERS='TRUE'
            Illumina-PE-R-processing.R \
                           "${params.left_trunc}" \
                           "${params.right_trunc}" \
                           "${params.left_maxEE}" \
                           "${params.right_maxEE}" \
                           "\${TRIM_PRIMERS}" \
                           "${sample_IDs_file}" \
                           "Trimmed_Sequence_Data/" \
                           "Filtered_Sequence_Data/" \
                           "${params.primer_trimmed_R1_suffix}" \
                           "${params.primer_trimmed_R2_suffix}" \
                           "${params.filtered_R1_suffix}" \
                           "${params.filtered_R2_suffix}" \
                           "final_outputs/" \
                           "${params.output_prefix}" \
                           "${params.target_region}" \
                           "${params.concatenate_reads_only}" \
                           "${params.assay_suffix}"
        
        else
        
            mkdir -p  Trimmed_Sequence_Data/ && \\
            mv ${trimmed_read_counts} *${params.primer_trimmed_R1_suffix}  Trimmed_Sequence_Data/

            mkdir -p Filtered_Sequence_Data/ final_outputs/

            TRIM_PRIMERS='TRUE'
            Illumina-SE-R-processing.R \
                          "${params.left_trunc}" \
                          "${params.left_maxEE}" \
                          "\${TRIM_PRIMERS}" \
                          "${sample_IDs_file}" \
                          "Trimmed_Sequence_Data/" \
                          "Filtered_Sequence_Data/" \
                          "${params.primer_trimmed_R1_suffix}" \
                          "${params.filtered_R1_suffix}" \
                          "final_outputs/" \
                          "${params.output_prefix}" \
                          "${params.target_region}" \
                          "${params.assay_suffix}"
                          
        fi
        # Sort the taxonomy count by ASV id               
        (head -n 1 final_outputs/taxonomy-and-counts${params.assay_suffix}.tsv; \\
            awk 'NR>1{print}' final_outputs/taxonomy-and-counts${params.assay_suffix}.tsv | sort -V -k1) \\
            > temp_tax_cont.tsv && mv temp_tax_cont.tsv final_outputs/taxonomy-and-counts${params.assay_suffix}.tsv

        R --vanilla --version  |grep "R version" > versions.txt
        get_R_package_version.R
        """
   
}



process RUN_R_NOTRIM {

    tag "Running dada2 on the raw reads..."
    beforeScript "chmod +x ${projectDir}/bin/*"

    input:
        tuple path(sample_IDs_file), val(isPaired)
        path(raw_reads)
        val(raw_read_suffix) //[R1,R2] or [R1]
    output:
        path("Filtered_Sequence_Data/*${params.filtered_R1_suffix[-5..-1]}"), emit: reads
        path("Filtered_Sequence_Data/filtered-read-counts${params.assay_suffix}.tsv"), emit: filtered_count
        path("final_outputs/taxonomy${params.assay_suffix}.tsv"), emit: taxonomy
        path("final_outputs/taxonomy-and-counts${params.assay_suffix}.biom"), emit: biom
        path("final_outputs/ASVs${params.assay_suffix}.fasta"), emit: fasta
        path("final_outputs/read-count-tracking${params.assay_suffix}.tsv"), emit: read_count
        path("final_outputs/counts${params.assay_suffix}.tsv"), emit: counts
        path("final_outputs/taxonomy-and-counts${params.assay_suffix}.tsv"), emit: taxonomy_count
        path("versions.txt"), emit: version  
    script:
        """
        if [ ${isPaired} == true ]; then

            mkdir -p raw_reads/ && \\
            mv *${raw_read_suffix[0]} *${raw_read_suffix[1]} raw_reads/
    
            mkdir -p Filtered_Sequence_Data/ final_outputs/
            
            TRIM_PRIMERS='FALSE'
            Illumina-PE-R-processing.R \
                            "${params.left_trunc}" \
                            "${params.right_trunc}" \
                            "${params.left_maxEE}" \
                            "${params.right_maxEE}" \
                            "\${TRIM_PRIMERS}" \
                            "${sample_IDs_file}" \
                            "raw_reads/" \
                            "Filtered_Sequence_Data/" \
                            "${raw_read_suffix[0]}" \
                            "${raw_read_suffix[1]}" \
                            "${params.filtered_R1_suffix}" \
                            "${params.filtered_R2_suffix}" \
                            "final_outputs/" \
                            "${params.output_prefix}" \
                            "${params.target_region}" \
                            "${params.concatenate_reads_only}" \
                            "${params.assay_suffix}"
            
        else
        
            mkdir -p raw_reads/ && mv *${raw_read_suffix[0]} raw_reads/
            mkdir -p Filtered_Sequence_Data/ final_outputs/
            
            TRIM_PRIMERS='FALSE'
            Illumina-SE-R-processing.R \
                           "${params.left_trunc}" \
                           "${params.left_maxEE}" \
                           "\${TRIM_PRIMERS}" \
                           "${sample_IDs_file}" \
                           "raw_reads/" \
                           "Filtered_Sequence_Data/" \
                           "${raw_read_suffix[0]}" \
                           "${params.filtered_R1_suffix}" \
                           "final_outputs/" \
                           "${params.output_prefix}" \
                           "${params.target_region}" \
                           "${params.assay_suffix}"

        fi
        # Sort the taxonomy count by ASV id               
        (head -n 1 final_outputs/taxonomy-and-counts${params.assay_suffix}.tsv; \\
            awk 'NR>1{print}' final_outputs/taxonomy-and-counts${params.assay_suffix}.tsv | sort -V -k1) \\
            > temp_tax_cont.tsv && mv temp_tax_cont.tsv final_outputs/taxonomy-and-counts${params.assay_suffix}.tsv
        
         R --vanilla --version  |grep "R version" > versions.txt
         get_R_package_version.R
        """
}
