#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/**************************************************************************************** 
*********************  Summarize Meta assembled genomes (MAGs) **************************
****************************************************************************************/


params.min_est_comp = 90
params.max_est_redund = 10
params.max_est_strain_het = 50

/*
Scratch directory for gtdb-tk, if wanting to use disk space instead of RAM, can be memory intensive;
see https://ecogenomics.github.io/GTDBTk/faq.html#gtdb-tk-reaches-the-memory-limit-pplacer-crashes
leave empty if wanting to use memory, the default, put in quotes the path to a directory that 
already exists if wanting to use disk space
*/

params.gtdb_tk_scratch_location = ""

/* 
    Filters checkm results based on estimate completion, redundancy, and 
    strain heterogeneity. Defaults are conservatively 90, 10, and 50  
*/

process FILTER_CHECKM_RESULTS_AND_COPY_MAGS {

    tag "Filtering checkm-s results..."
    label "mags"
    label "bit"

    input:
        path(bins_checkm_results)
        path(bins) 
    output:
        path("${params.additional_filename_prefix}MAGs-checkm-out.tsv"), emit: MAGs_checkm_out
        path("MAGs_dir/"), emit: MAGs_dir 
    script:
        """        
        # Only running if there were bins recovered
        if [ `find -L . -name '*.fasta' | wc -l | sed 's/^ *//'` -gt 0 ]; then

            cat <( printf "Bin Id\\tMarker lineage\\t# genomes\\t# markers\\t# marker sets\\t0\\t1\\t2\\t3\\t4\\t5+\\tCompleteness\\tContamination\\tStrain heterogeneity\\n" ) \\
                <( awk -F '\\t' ' \$12 >= ${params.min_est_comp} && \$13 <= ${params.max_est_redund} && \$14 <= ${params.max_est_strain_het} ' ${bins_checkm_results} ) \\
                > MAGs-checkm-out.tmp

            sed 's/-bin\\./-MAG-/' MAGs-checkm-out.tmp > ${params.additional_filename_prefix}MAGs-checkm-out.tsv
            
            [ -d MAGs_dir/ ] || mkdir MAGs_dir/
            for MAG in `cut -f 1 MAGs-checkm-out.tmp | tail -n +2`
            do
                new_ID=`echo \$MAG | sed 's/-bin\\./-MAG-/'`
                cp \$MAG.fasta MAGs_dir/\${new_ID}.fasta
            done

        else

            printf "There were no MAGs recovered.\\n" > ${params.additional_filename_prefix}MAGs-checkm-out.tsv

        fi
        """
}



// Assign taxonomy to MAGs with gtdb-tk

process  GTDBTK_ON_MAGS {
   
    tag "Assigning taxonomy to your MAGs with gtdb-tk..."  
    label "mags"

    input:
        path(MAGs_checkm_out)
        path(MAGs_dir)
        path(gtdbtk_db_dir)
        val(use_gtdbtk_scratch_location)
        env(GTDBTK_DATA_PATH)
           
    output:
        path("gtdbtk-out/")
    script:
        """
        # Only running if any MAGs were recovered
        if [ `find -L ${MAGs_dir} -name '*.fasta' | wc -l | sed 's/^ *//'` -gt 0 ]; then

            if [ ${use_gtdbtk_scratch_location} == 'true' ]; then

              [ -d gtdbtk_scratch_location/ ] || mkdir gtdbtk_scratch_location/

                gtdbtk classify_wf \\
                        --scratch_dir gtdbtk_scratch_location/ \\
                        --genome_dir ${MAGs_dir} \\
                        -x fasta \\
                        --out_dir gtdbtk-out/ \\
                        --cpus ${task.cpus} \\
                        --pplacer_cpus 1

            else

                gtdbtk classify_wf \\
                       --genome_dir ${MAGs_dir} \\
                       -x fasta \\
                       --out_dir gtdbtk-out/ \\
                       --cpus ${task.cpus} \\
                       --pplacer_cpus 1

            fi

        else

            mkdir -p gtdbtk-out/
            printf "There were no MAGs recovered.\\n" \\
                   > gtdbtk-out/No-MAGs-recovered.txt
                   
            printf "\\n\\nThere were no MAGs recovered, so GTDB-tk was not run.\\n\\n"

        fi
        """
}



// Summarize MAG assemblies
process  SUMMARIZE_MAG_ASSEMBLIES {

    tag "Summarizing MAG assemblies..."
    label "mags"
    label "bit"

    input:
        path(MAGs_dir)
    output:
        path("${params.additional_filename_prefix}MAG-assembly-summaries.tsv")
    script:
        """
        # Only running if any MAGs were recovered
        if [ `find -L ${MAGs_dir} -name '*.fasta' | wc -l | sed 's/^ *//'` -gt 0 ]; then

            bit-summarize-assembly ${MAGs_dir}/*.fasta -o MAG-summaries.tmp -t

            # Slimming down the output
            cut -f 1,2,3,5,6,8,11,18,19,20 MAG-summaries.tmp \\
               > ${params.additional_filename_prefix}MAG-assembly-summaries.tsv

        else

            printf "There were no MAGs recovered.\\n" \\
               > ${params.additional_filename_prefix}MAG-assembly-summaries.tsv

        fi
        """
}

process  GENERATE_MAGS_OVERVIEW_TABLE {

    tag "Generating an overview table of all MAGs..."
    label "mags"
    label "bit"

    input:
        path(MAG_assembly_summaries)
        path(MAGs_checkm_out)
        path(gtdbtk_out)
        path(MAGs_dir)
    output:
        path("${params.additional_filename_prefix}MAGs-overview${params.assay_suffix}.tsv")

    script:
        """
        # Only running if any MAGs were recovered
        if [ `find -L ${MAGs_dir}  -name '*.fasta' | wc -l | sed 's/^ *//'` -gt 0 ]; then

        #--------------------------- get_MAGs_estimates_and_taxonomy.sh ------------------------------------#
        get_MAGs_estimates_and_taxonomy.sh ${MAGs_dir} ${MAG_assembly_summaries} ${MAGs_checkm_out} ${gtdbtk_out}
        #----------------------------------------------------------------------------------------------------#

            # Adding headers
            cat <(printf "est. completeness\\test. redundancy\\test. strain heterogeneity\\n") \\
                checkm-estimates.tmp > checkm-estimates-with-headers.tmp

            cat <(printf "domain\\tphylum\\tclass\\torder\\tfamily\\tgenus\\tspecies\\n") \\
                gtdb-taxonomies.tmp > gtdb-taxonomies-with-headers.tmp

            paste ${MAG_assembly_summaries} \\
                  checkm-estimates-with-headers.tmp \\
                  gtdb-taxonomies-with-headers.tmp \\
                  > MAGs-overview.tmp

            # Ordering by taxonomy
            head -n 1 MAGs-overview.tmp > MAGs-overview-header.tmp

            tail -n +2 MAGs-overview.tmp | \\
            sort -t \$'\\t' -k 14,20 > MAGs-overview-sorted.tmp

            cat MAGs-overview-header.tmp MAGs-overview-sorted.tmp \\
               > ${params.additional_filename_prefix}MAGs-overview${params.assay_suffix}.tsv

        else

            printf "There were no MAGs recovered.\\n" \\
              > ${params.additional_filename_prefix}MAGs-overview${params.assay_suffix}.tsv

        fi
        """
}


process SUMMARIZE_MAG_LEVEL_KO_ANNOTATIONS {

    tag "Parsing MAG KO annotations..."
    label "mags"
    label "bit"

    input:
        path(MAGs_overview)
        path(gene_coverage_annotation_and_tax_files)
        path(MAGs_dir)
    output:
        path("${params.additional_filename_prefix}MAG-level-KO-annotations${params.assay_suffix}.tsv")

    script:
        """
        # Only running if any MAGs were recovered
        if [ `find -L ${MAGs_dir} -name '*.fasta' | wc -l | sed 's/^ *//'` -gt 0 ]; then

            for MAG in `cut -f 1 ${MAGs_overview} | tail -n +2`
            do

                sample_ID=`echo \$MAG | sed 's/-MAG-[0-9]*\$//'`
                grep "^>" ${MAGs_dir}/\$MAG.fasta | tr -d ">" > curr-contig-ids.tmp

                parse-MAG-annots.py \\
                        -i \${sample_ID}-gene-coverage-annotation-and-tax.tsv \\
                        -w curr-contig-ids.tmp \\
                        -M \$MAG \\
                        -o ${params.additional_filename_prefix}MAG-level-KO-annotations${params.assay_suffix}.tsv

            done

        else

            printf "There were no MAGs recovered.\\n" \\
               > ${params.additional_filename_prefix}MAG-level-KO-annotations${params.assay_suffix}.tsv

        fi
        """
}


process SUMMARIZE_MAG_KO_ANNOTS_WITH_KEGG_DECODER {

    tag "Summarizing MAG KO annotations using kegg decoder..."
    label "mags"


    input:
        path(MAG_level_KO_annotations)
        path(MAGs_dir)
    output:
        path("${params.additional_filename_prefix}MAG-KEGG-Decoder-out${params.assay_suffix}.tsv")

    script:
        """
        # Getting number of MAGs recovered
        num_mags_recovered=`find -L ${MAGs_dir}/ -name '*.fasta' | wc -l | sed 's/^ *//'`
        # Only running if any MAGs were recovered
        if [ \$num_mags_recovered -gt 0 ]; then

            # KEGGDecoder splits on the first underscore to identify unique genome/MAG IDs
            # this can be problematic with how things are named, so we are swapping them all to not have
            # any "_" first, then afterwards we are changing the output table back to the original names so 
            # they match elsewhere (they will still be slightly different in the html output, but that is
            # only manually explored anyway)

            # Making version of input for KEGGDecoder with no underscores
            tr "_" "-" < ${MAG_level_KO_annotations} > mod-MAG-level-KO-annotations.tmp

            # Making mapping file
            paste <( cut -f 1 ${MAG_level_KO_annotations} ) \\
                  <( cut -f 1 mod-MAG-level-KO-annotations.tmp ) \\
                  > MAG-ID-map.tmp

            # Running KEGGDecoder
            # can only create html output if there are more than 1
            if [ \$num_mags_recovered -gt 1 ]; then
                KEGG-decoder -v interactive -i mod-MAG-level-KO-annotations.tmp -o MAG-KEGG-Decoder-out.tmp
            else
                KEGG-decoder -i mod-MAG-level-KO-annotations.tmp -o MAG-KEGG-Decoder-out.tmp
            fi

            # Swapping MAG IDs back in output tsv from KEGGDecoder
            swap-MAG-IDs.py -i MAG-KEGG-Decoder-out.tmp -m MAG-ID-map.tmp -o MAG-KEGG-Decoder-out.tsv && \\
            mv MAG-KEGG-Decoder-out.tsv \\
              ${params.additional_filename_prefix}MAG-KEGG-Decoder-out${params.assay_suffix}.tsv


        else

            printf "There were no MAGs recovered.\\n" \\
               >  ${params.additional_filename_prefix}MAG-KEGG-Decoder-out${params.assay_suffix}.tsv

        fi
        """
}



workflow summarize_mags {
    take:
        bins_checkm_results_ch
        bins_ch
        gtdbtk_db_dir
        use_gtdbtk_scratch_location
        gene_coverage_annotation_and_tax_files_ch 


    main:
        FILTER_CHECKM_RESULTS_AND_COPY_MAGS(bins_checkm_results_ch, bins_ch) 
        MAGs_checkm_out_ch = FILTER_CHECKM_RESULTS_AND_COPY_MAGS.out.MAGs_checkm_out
        MAGs_dir_ch = FILTER_CHECKM_RESULTS_AND_COPY_MAGS.out.MAGs_dir

        gtdbtk_out_ch = GTDBTK_ON_MAGS(MAGs_checkm_out_ch, MAGs_dir_ch, gtdbtk_db_dir, use_gtdbtk_scratch_location, gtdbtk_db_dir)

        MAG_assembly_summaries_ch = SUMMARIZE_MAG_ASSEMBLIES(MAGs_dir_ch)

        MAGs_overview_ch = GENERATE_MAGS_OVERVIEW_TABLE(MAG_assembly_summaries_ch,
                                                        MAGs_checkm_out_ch,
                                                        gtdbtk_out_ch,
                                                        MAGs_dir_ch)

        MAG_level_KO_annotations_ch = SUMMARIZE_MAG_LEVEL_KO_ANNOTATIONS(MAGs_overview_ch, 
                                           gene_coverage_annotation_and_tax_files_ch, 
                                           MAGs_dir_ch)

        SUMMARIZE_MAG_KO_ANNOTS_WITH_KEGG_DECODER(MAG_level_KO_annotations_ch, MAGs_dir_ch)

    emit:
        MAGs_overview = MAGs_overview_ch
        MAGs_dir = MAGs_dir_ch
    
}
