#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process VSEARCH_DEREP_SAMPLE {
    
    tag "dereplicating ${sample_id}s sequences..."
    label "vsearch"

    input:
        tuple val(sample_id), path(reads)
    output:
        path("${sample_id}-derep.fa.tmp"), emit: derep_temp
        path("versions.txt"), emit: version
    script:
        """
        vsearch --derep_fulllength ${reads} \\
                --strand both \\
                --output "${sample_id}-derep.fa.tmp" \\
                --sizeout \\
                --relabel "sample=${sample_id};seq_" > /dev/null

        vsearch --version 2>&1 |grep "vsearch" |head -n1 | sed -E 's/(vsearch v.+?)_linux.+/\\1/' > versions.txt
        """
}


process VSEARCH_COMBINE_DEREPD_SAMPLES {

    tag "Combining all dereplicated samples..."
    label "vsearch"

    input:
        path(derepd_reads)
    output:
        path("all-samples.fa.tmp")
    script:
        """
        cat ${derepd_reads} > all-samples.fa.tmp
        """
}

 process VSEARCH_PROCESS_ALL {

    tag "Clustering your sequences to OTUs using vsearch..."
    label "vsearch"

    input:
        path(all_samples_fasta)
    output:
        path("${params.output_prefix}OTUs.fasta"), emit: fasta
        path("${params.output_prefix}counts${params.assay_suffix}.tsv"), emit: counts
        path("versions.txt"), emit: version
    script:
        """
        # Dereplicate all
        vsearch --derep_fulllength ${all_samples_fasta} \\
                --strand both \\
                --output all-samples_derep.fa.tmp \\
                --sizein --sizeout

        # Clustering to get rep seqs
        vsearch --cluster_size all-samples_derep.fa.tmp \\
                --id 0.97 \\
                --strand both \\
                --sizein \\
                --sizeout \\
                --relabel "OTU_" \\
                --centroids rep-seqs.fa.tmp

        # Removing singletons
        vsearch --sortbysize rep-seqs.fa.tmp \\
                --minsize 2 \\
                --output rep-seqs-no-singletons.fa.tmp 

        # Chimera check and removal
        vsearch --uchime_denovo rep-seqs-no-singletons.fa.tmp \\
                --sizein \\
                --nonchimeras ${params.output_prefix}OTUs.fasta \\
                --relabel "OTU_"

        # Mapping seqs to OTUs to get OTU abundances per sample
        vsearch --usearch_global ${all_samples_fasta} \\
                -db ${params.output_prefix}OTUs.fasta \\
                --sizein --id 0.97 \\
                --otutabout counts.tmp
                
        sed 's/^#OTU ID/OTU_ID/' counts.tmp \\
            > ${params.output_prefix}counts${params.assay_suffix}.tsv
       
        vsearch --version 2>&1 |grep "vsearch" |head -n1 | sed -E 's/(vsearch v.+?)_linux.+/\\1/' > versions.txt
        """
 }


process REMOVE_LINE_WRAPS {

    tag "Removing line wraps from OTU fasta file..."

    input:
        path(temp_fasta)
    output:
        path("${params.output_prefix}OTUs${params.assay_suffix}.fasta"), emit: fasta
        path("versions.txt"), emit: version
    script:
        """
        # Removing line wraps from fasta file
        bit-remove-wraps ${temp_fasta} \\
            > ${params.output_prefix}OTUs${params.assay_suffix}.fasta.tmp && \\
        mv ${params.output_prefix}OTUs${params.assay_suffix}.fasta.tmp \\
           ${params.output_prefix}OTUs${params.assay_suffix}.fasta

        bit-version 2>&1 |grep "Bioinformatics Tools"|sed -E 's/^\\s+//' > versions.txt      
        """
}


workflow pick_otus {

    take:
        reads_ch

    main:
        VSEARCH_DEREP_SAMPLE(reads_ch)
        derep_temp_ch = VSEARCH_DEREP_SAMPLE.out.derep_temp.collect()
        VSEARCH_COMBINE_DEREPD_SAMPLES(derep_temp_ch) | VSEARCH_PROCESS_ALL
        REMOVE_LINE_WRAPS(VSEARCH_PROCESS_ALL.out.fasta)

        // capture software versions
        software_versions_ch = Channel.empty()
        VSEARCH_DEREP_SAMPLE.out.version | mix(software_versions_ch) | set{software_versions_ch}
        VSEARCH_PROCESS_ALL.out.version | mix(software_versions_ch) | set{software_versions_ch}
        REMOVE_LINE_WRAPS.out.version | mix(software_versions_ch) | set{software_versions_ch}

    emit:
        otus = REMOVE_LINE_WRAPS.out.fasta
        counts = VSEARCH_PROCESS_ALL.out.counts
        versions = software_versions_ch

}
