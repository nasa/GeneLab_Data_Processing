process GENERATE_MD5SUMS {
  // Generates tabular data for GeneLab raw and processed data
    publishDir "${ch_outdir}/Metadata",
        mode: params.publish_dir_mode,
        pattern: "*.csv"

    input:
        path(ch_outdir)
    output:
        path("*md5sum*")

    script:
        """
        generate_md5sum_files.py  --outdir ${ch_outdir} --assay-suffix ${params.assay_suffix} \\
        """
}