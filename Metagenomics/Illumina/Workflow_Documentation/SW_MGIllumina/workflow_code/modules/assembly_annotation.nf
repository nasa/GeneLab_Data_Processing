#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
params.cat_db = "/mnt/c/Users/olabi/Documents/bioinformatics/test/processing_info/ref-dbs/CAT_prepare_20210107/2021-01-07_CAT_database/"
params.ko_db_dir = "/mnt/c/Users/olabi/Documents/bioinformatics/test/processing_info/ref-dbs/kofamscan_db/"
params.block_size = 4

/**************************************************************************************** 
**************************  Sequence Assembly Annotation *******************************
****************************************************************************************/
// This process calls genes on each assembly file.
process CALL_GENES {

    tag "Predicting genes for ${sample_id}-s assembly"
    label "call_genes"
    
    input:
        tuple val(sample_id), path(assembly) 
    output:
        // Amino acids, nucleotides and gff
        tuple val(sample_id), path("${sample_id}-genes.faa"), path("${sample_id}-genes.fasta"), path("${sample_id}-genes.gff"), emit: genes
        path("versions.txt"), emit: version
    script:
        """
        # Only running if assembly produced any contigs
        if [ -s ${assembly} ]; then

            prodigal -q -c -p meta -a ${sample_id}-genes.faa \\
                     -d ${sample_id}-genes.fasta \\
                     -f gff -o ${sample_id}-genes.gff \\
                     -i ${assembly} 
        else

            touch ${sample_id}-genes.faa ${sample_id}-genes.fasta ${sample_id}-genes.gff
            printf "Gene-calling not performed because the assembly didn't produce anything.\\n"

        fi
        prodigal -v 2>&1 | grep Prodigal > versions.txt
        """
}        

// Removing line-wraps using bit
process REMOVE_LINEWRAPS {

    tag "Remove line wraps in ${sample_id}-s nucleotide and amino acid file..."
    label "call_genes"
    label "bit"

    input:
        tuple val(sample_id), path(aa), path(nt), path(gff)
   
    output:
        tuple val(sample_id), path("${sample_id}-genes.faa"), path("${sample_id}-genes.fasta"), emit: genes
        path("versions.txt"), emit: version
    script:
        """
         if [ -s ${aa} ] && [ -s ${nt} ]; then
            # Removing line-wraps
            bit-remove-wraps ${aa} > ${sample_id}-genes.faa.tmp 2> /dev/null && \\
            mv ${sample_id}-genes.faa.tmp ${sample_id}-genes.faa
            
            bit-remove-wraps ${nt} > ${sample_id}-genes.fasta.tmp 2> /dev/null && \\
            mv ${sample_id}-genes.fasta.tmp ${sample_id}-genes.fasta
        else

            touch ${sample_id}-genes.faa ${sample_id}-genes.fasta
            printf "Line wrapping not performed because gene-calling wasn't performed on ${sample_id}.\\n"
        fi 
        bit-version |grep "Bioinformatics Tools"|sed -E 's/^\\s+//' > versions.txt
        """
}


// This process runs the gene-level (KO) functional annotation for each sample.
// KO annotatiuon of the predicted amino acids
process KO_ANNOTATION {

    tag "Running KO annotation of ${sample_id}-s predicted amino acids.."
    //label "contig_annotation"
    
    input:
       tuple val(sample_id), path(assembly), path(aa), path(nt)
       path(ko_db_dir)
    output:
        tuple val(sample_id), path("${sample_id}-KO-tab.tmp"), emit: temp_table
        path("versions.txt"), emit: version
    script:
        """
        # only running if assembly produced any contigs and genes were identified (they are required for this)
        if [ -s ${assembly} ] && [ -s ${aa} ]; then

            exec_annotation -p ${ko_db_dir}/profiles/ \\
                            -k ${ko_db_dir}/ko_list \\
                            --cpu ${task.cpus} -f detail-tsv \\
                            -o ${sample_id}-KO-tab.tmp --tmp-dir ${sample_id}-tmp-KO-dir \\
                            --report-unannotated ${aa}

        else

            touch ${sample_id}-KO-tab.tmp
            printf "Functional annotations not performed because the assembly didn't produce anything and/or no genes were identified.\\n"

        fi
        exec_annotation -v > versions.txt
        """
}


process FILTER_KFAMSCAN {

    tag "Filtering ${sample_id}-s KO annotation results..."
    label "bit"
    label "contig_annotation"

    input:
       tuple val(sample_id), path(KO_tab_tmp)   
    output:
        tuple val(sample_id), path("${sample_id}-annotations.tsv"), emit: ko_annotation
        path("versions.txt"), emit: version
    script:
        """
        if [ -s ${KO_tab_tmp} ]; then

            bit-filter-KOFamScan-results -i ${KO_tab_tmp} -o ${sample_id}-annotations.tsv

        else 

            touch ${sample_id}-annotations.tsv
            printf "Nothing to filter since functional annotation was not performed.\\n"

        fi
        bit-version |grep "Bioinformatics Tools"|sed -E 's/^\\s+//' > versions.txt
        """

}

// This process runs the gene- and contig-level taxonomic classifications for each assembly.
process TAX_CLASSIFICATION {

    tag "Taxonomy classification of ${sample_id}-s "
    label "contig_annotation"

    input:
       tuple val(sample_id), path(assembly), path(aa), path(nt)
       path(cat_db)
    output:
        // Gene and contig taxonomy
        tuple val(sample_id), path("${sample_id}-gene-tax.tsv"), path("${sample_id}-contig-tax.tsv"), emit: taxonomy
        path("versions.txt"), emit: version
    script:
        """
        # Only running if assembly produced any contigs and 
        # genes were identified (they are required for this)
        if [ -s ${assembly} ] && [ -s ${aa} ]; then

            CAT contigs -d ${cat_db}/${params.cat_db_sub_dir} -t ${cat_db}/${params.cat_taxonomy_dir} \\
                        -n ${task.cpus} -r 3 --top 4 \\
                        --I_know_what_Im_doing -c ${assembly} \\
                        -p ${aa} -o ${sample_id}-tax-out.tmp \\
                        --no_stars --block_size ${params.block_size} \\
                        --index_chunks 2 --force 

            # Adding names to gene classifications
            CAT add_names -i ${sample_id}-tax-out.tmp.ORF2LCA.txt \\
                          -o ${sample_id}-gene-tax.tmp -t ${cat_db}/${params.cat_taxonomy_dir} \\
                          --only_official --exclude_scores

            # Formatting gene classifications
            bash format-gene-tax-classifications.sh \\
                       ${sample_id}-gene-tax.tmp ${sample_id}-gene-tax.tsv

            # Adding names to contig classifications
            CAT add_names -i ${sample_id}-tax-out.tmp.contig2classification.txt \\
                          -o ${sample_id}-contig-tax.tmp -t ${cat_db}/${params.cat_taxonomy_dir} \\
                          --only_official --exclude_scores

            # Formatting contig classifications
            bash format-contig-tax-classifications.sh \\
                  ${sample_id}-contig-tax.tmp ${sample_id}-contig-tax.tsv

        else

            touch ${sample_id}-gene-tax.tsv ${sample_id}-contig-tax.tsv
            printf "Assembly-based taxonomic classification not performed because the assembly didn't produce anything and/or no genes were identified.\\n" 

        fi
        CAT --version | sed -E 's/(CAT v.+)\\s\\(.+/\\1/'  > versions.txt
        """
}

workflow annotate_assembly {
    take:
        assembly_ch
        ko_db_dir
        cat_db
 
    main:
        CALL_GENES(assembly_ch)
        genes_ch = CALL_GENES.out.genes | REMOVE_LINEWRAPS.out.genes
        KO_ANNOTATION(assembly_ch.join(genes_ch) ko_db_dir)
        KO_ANNOTATION.out.temp_table | FILTER_KFAMSCAN.out.ko_annotation
        TAX_CLASSIFICATION(assembly_ch, genes_ch, cat_db)

}
