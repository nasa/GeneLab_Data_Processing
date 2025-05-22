process GENERATE_PROTOCOL {
    publishDir "${ch_outdir}/GeneLab",
        mode: params.publish_dir_mode,
        pattern: "*.txt"

    publishDir "${ch_outdir}/Metadata",
        mode: params.publish_dir_mode,
        pattern: "runsheet/*",
        saveAs: { filename ->
            if (filename.startsWith("runsheet/")) return filename.replace("runsheet/", "")
            else return filename
        }

    input:
        val(ch_outdir)
        val(ch_meta)
        val(strandedness)
        path(software_versions_yaml)
        val(reference_source)
        val(reference_version)
        tuple(path(reference_fasta), path(reference_gtf))
        path(runsheet)

    output:
        path("protocol${params.assay_suffix}.txt")
        path("runsheet/${runsheet}")

    script:
        def mode = params.mode == 'microbes' ? '--mode microbes' : ''
        def ref_source = reference_source ? "--reference_source ${reference_source}" : ''
        def ref_version = reference_version ? "--reference_version ${reference_version}" : ''

        """
        generate_protocol.py \
        ${mode} \
        --outdir . \
        --software_table ${software_versions_yaml} \
        --assay_suffix ${params.assay_suffix} \
        --paired_end ${ch_meta.paired_end} \
        --has_ercc ${ch_meta.has_ercc} \
        --workflow_version ${workflow.manifest.version} \
        --strandedness ${strandedness} \
        --organism "${ch_meta.organism_sci}" \
        ${ref_source} \
        ${ref_version} \
        --reference_fasta ${reference_fasta} \
        --reference_gtf ${reference_gtf} \
        --runsheet ${runsheet}

        mkdir -p runsheet
        mv ${runsheet} runsheet/
        """
}