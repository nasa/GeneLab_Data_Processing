// Include the COPY_READS module
include { COPY_READS } from '../modules/copy_reads.nf'

workflow STAGE_RAW_READS {
    take:
        ch_samples

    main:
        truncate_to = params.truncate_to
        if ( truncate_to ) {
        // Download truncated raw reads
            ch_samples | map { it -> it[0].paired_end ? [it[0], it[1][0], it[1][1]] : [it[0], it[1][0]]}
                 | branch {
                   paired: it.size() == 3
                   single: it.size() == 2
                 }
                 | set{ ch_raw_read_pointers }
            
            
            // TO DO: Move the two splitFastq truncation steps into processes that run in parallel. Low priority since this is just for debugging.

            // PAIRED END
            // Only difference is the splitFastq arg 'pe'
            ch_raw_read_pointers.paired | splitFastq(pe: true, decompress: true, compress: true, limit: truncate_to, by: truncate_to, file: true)
                                        | map { it -> [ it[0], [ it[1], it[2] ] ]}
                                        // | view { it -> "TRUNCATED PAIRED READS ($truncate_to): $it[0]"}
                                        | set { ch_raw_reads }
            // SINGLE END
            // Only difference is the splitFastq arg 'pe'
            ch_raw_read_pointers.single | splitFastq(decompress: true, compress: true, limit: truncate_to, by: truncate_to, file: true)
                                        | map { it -> [ it[0], [ it[1] ] ]}
                                        // | view { it -> "TRUNCATED SINGLE READS ($truncate_to): $it[0]"}
                                        | mix( ch_raw_reads )
                                        | set { ch_raw_reads }

            // Moves the truncated files to expected raw read locations as per samplesheet
            COPY_READS(ch_raw_reads)
            // Collect sample IDs into a file
            COPY_READS.out.raw_reads | map{ it -> it[1] } | collect | set { ch_all_raw_reads }
            COPY_READS.out.raw_reads | map { it[0].id }
                            | collectFile(name: "samples.txt", sort: true, newLine: true)
                            | set { samples_txt }
        } else {
        // download full raw reads
            ch_samples | map { it -> it[0].paired_end ? [it[0], [ it[1][0], it[1][1] ]] : [it[0], [it[1][0]]]}
                         | set { ch_raw_reads }

            // Download the raw reads and publish them to expected raw read locations as per samplesheet
            COPY_READS(ch_raw_reads)
            // Collect sample IDs into a file
            COPY_READS.out.raw_reads | map{ it -> it[1] } | collect | set { ch_all_raw_reads }
            COPY_READS.out.raw_reads | map { it[0].id }
                            | collectFile(name: "samples.txt", sort: true, newLine: true)
                            | set { samples_txt }
        }


    emit:
        raw_reads = COPY_READS.out.raw_reads 
        ch_all_raw_reads = ch_all_raw_reads
        samples_txt =  samples_txt
}

