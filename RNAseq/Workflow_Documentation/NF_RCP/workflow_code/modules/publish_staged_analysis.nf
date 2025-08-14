/*
    This process publishes the staged runsheet and raw reads to the output directory.
    This is only used for the stage-only workflow.
*/

process PUBLISH_STAGED_ANALYSIS {
    publishDir "${ch_outdir}",
    pattern: '{00-RawData/**,Metadata/**}',
    mode: params.publish_dir_mode

    input:
        val(ch_outdir)
        path(runsheet)
        path(raw_reads)
    
    output:
        path("Metadata/*"), emit: metadata
        path("00-RawData/Fastq/*"), emit: raw_reads
    
    script:
        """
        # Create directory structure
        mkdir -p Metadata
        mkdir -p 00-RawData/Fastq
        
        # Move runsheet to Metadata directory
        mv ${runsheet} Metadata/
        
        # Move raw reads to 00-RawData/Fastq directory
        mv ${raw_reads} 00-RawData/Fastq/

        """
} 