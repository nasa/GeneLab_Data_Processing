#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process RUN_R {

    tag "Assigning taxonomy to OTUs using decipher..."
   
    input:
        path(otus) // config["final_outputs_dir"] + config["output_prefix"] + "OTUs.fasta",
        path(counts) // config["final_outputs_dir"] + config["output_prefix"] + "counts.tsv"
        path(trimmed_read_counts) //config["trimmed_reads_dir"] + config["output_prefix"] + "trimmed-read-counts.tsv"
        path(filtered_read_counts)  // config["filtered_reads_dir"] + config["output_prefix"] +  filtered-read-counts.tsv
    output:
        path("Final_Outputs/${params.output_prefix}taxonomy${params.assay_suffix}.tsv"), emit: taxonomy
        path("Final_Outputs/${params.output_prefix}taxonomy-and-counts${params.assay_suffix}.biom"), emit: biom
        path("Final_Outputs/${params.output_prefix}taxonomy-and-counts${params.assay_suffix}.tsv"), emit: tsv
        path("Final_Outputs/${params.output_prefix}read-count-tracking${params.assay_suffix}.tsv"), emit: read_count
    script:
        """
        mkdir Trimmed_Sequence_Data/ && mv ${trimmed_read_counts}  Trimmed_Sequence_Data/
        mkdir Filtered_Sequence_Data/ && mv ${filtered_read_counts} Filtered_Sequence_Data/
        mkdir Final_Outputs/ && \\
          cp ${otus} Final_Outputs/ && \\
          mv ${counts} Final_Outputs/

        454-IonTorrent-R-processing.R \
                         "${otus}" \
                         "Trimmed_Sequence_Data/" \
                         "Filtered_Sequence_Data/" \
                         "Final_Outputs/" \
                         "${params.output_prefix}" \
                         "${params.target_region}" \
                         "${params.assay_suffix}"
                         
        # Sort the taxonomy count by ASV id
        (head -n 1 "Final_Outputs/${params.output_prefix}taxonomy-and-counts${params.assay_suffix}.tsv"; \\
            awk 'NR>1{print}' "Final_Outputs/${params.output_prefix}taxonomy-and-counts${params.assay_suffix}.tsv" | sort -V -k1) \\
            > temp_tax_cont.tsv && mv temp_tax_cont.tsv "Final_Outputs/${params.output_prefix}taxonomy-and-counts${params.assay_suffix}.tsv"

        """
}
