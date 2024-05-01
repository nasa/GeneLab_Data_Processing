#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.reduced_tree = "True"

/**************************************************************************************** 
*********************  Bin check and summary ********************************************
****************************************************************************************/

// Summarize bin assemblies
process SUMMARIZE_BIN_ASSEMBLIES {
 
    tag "Getting a summary of the recovered bins..."
    label "bins"
    label "bit"
      
    input:
        path(bins)
    output:
        path("${params.additional_filename_prefix}bin-assembly-summaries.tsv")

    script:
        """
        # Only running if any bins were recovered
        if [ `find -L . -name '*.fasta' | wc -l | sed 's/^ *//'` -gt 0 ]; then

            bit-summarize-assembly *.fasta -o bin-summaries.tmp -t

            # Slimming down the output
            cut -f 1,2,3,5,6,8,11,18,19,20 bin-summaries.tmp \\
                > ${params.additional_filename_prefix}bin-assembly-summaries.tsv

        else

            printf "There were no bins recovered.\\n" \\
              > ${params.additional_filename_prefix}bin-assembly-summaries.tsv

        fi
        """
}


// Runs checkm on recovered bins
process CHECKM_ON_BINS {

    tag "Running checkm on the recovered bins..."
    label "bins"
    
    input:
        path(bins) 
    output:
        path("${params.additional_filename_prefix}bins-checkm-out.tsv")

    script:
        """
        # only running if there were bins recovered
        if [ `find -L . -name '*fasta' | wc -l | sed 's/^ *//'` -gt 0 ]; then

            mkdir -p checkm-working-tmp/

            if [ ${params.reduced_tree} == "True" ]; then

                checkm lineage_wf \\
                       -f ${params.additional_filename_prefix}bins-checkm-out.tsv \\
                       --tab_table \\
                       -t ${task.cpus} \\
                       --reduced_tree \\
                       --pplacer_threads 1 \\
                       -x fasta . checkm-out-tmp/ \\
                       --tmpdir checkm-working-tmp/

            else

                checkm lineage_wf \\
                       -f ${params.additional_filename_prefix}bins-checkm-out.tsv \\
                       --tab_table \\
                       -t ${task.cpus} \\
                       --pplacer_threads 1 \\
                       -x fasta . checkm-out-tmp/ \\
                       --tmpdir checkm-working-tmp/

            fi

        else

            printf "There were no bins recovered, so checkm was not run.\\n" \\
              > ${params.additional_filename_prefix}bins-checkm-out.tsv

        fi
        """
}

process GENERATE_BINS_OVERVIEW_TABLE {

    tag "Generating an overall overview of the recovered bins..."
    label "bins"
    label "bit"

    input:
        path(bin_assembly_summaries) 
        path(bins_checkm_results)
        path(bins) 
    output:
        path("${params.additional_filename_prefix}bins-overview${params.assay_suffix}.tsv")

    script:
        """
        # Only running if there were bins recovered
        if [ `find -L . -name '*.fasta' | wc -l | sed 's/^ *//'` -gt 0 ]; then

            # Making sure none of the intermediate files exist already
            [ -f checkm-estimates.tmp ] && rm -rf checkm-estimates.tmp 
            [ -f checkm-estimates-with-headers.tmp ] && rm -rf checkm-estimates-with-headers.tmp

            for bin in `cut -f 1 ${bin_assembly_summaries} | tail -n +2`
            do

                grep -w -m 1 "^\$bin" ${bins_checkm_results} | \\
                cut -f 12,13,14 >> checkm-estimates.tmp

            done

            # Adding header
            cat <(printf "est. completeness\\test. redundancy\\test. strain heterogeneity\\n") \\
                 checkm-estimates.tmp > checkm-estimates-with-headers.tmp

            # Combining
            paste ${bin_assembly_summaries} checkm-estimates-with-headers.tmp \\
               > ${params.additional_filename_prefix}bins-overview${params.assay_suffix}.tsv

        else
        
            printf "There were no bins recovered.\\n" \\
               > ${params.additional_filename_prefix}bins-overview${params.assay_suffix}.tsv

        fi
        """
}


workflow summarize_bins {

    take:
        binning_ch

    main:
        bins = binning_ch.map{ sample_id, depth, bins -> bins instanceof List ? bins.each{it}: bins }.flatten().collect()
        bin_assembly_summaries_ch = SUMMARIZE_BIN_ASSEMBLIES(bins)
        bins_checkm_results_ch  = CHECKM_ON_BINS(bins)
        table = GENERATE_BINS_OVERVIEW_TABLE(bin_assembly_summaries_ch, bins_checkm_results_ch, bins)
    
    emit:
        bins_checkm_results = bins_checkm_results_ch
        overview_table = table
}
