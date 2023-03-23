/* Processes dealing with fastQC files
*/

process GET_MAX_READ_LENGTH {
  // Returns max read length based on all supplied fastqc data files
  input:
    path(fastqc_datazip)

  output:
    env(RESULT), emit: length

  script:
    """
    # unzip fastqc datazip
    unzip ${ fastqc_datazip }

    RESULT=\$(get_max_readlength_from_fastqc.py "${ fastqc_datazip.baseName }/fastqc_data.txt")
    """
}
