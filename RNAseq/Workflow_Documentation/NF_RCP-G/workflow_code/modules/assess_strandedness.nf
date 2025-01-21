process ASSESS_STRANDEDNESS {
  tag "Dataset-wide"
  input:
    path("infer_out/*") // a collection of infer_experiment stdout files

  output:
    path("result.txt")

  stub:
    """
    assess_strandedness.py infer_out
    echo "unstranded:0.48595" > result.txt # override original results, this is because heavy truncation and genome subsampling can result in an ambigious strand assignment, which normally is an issue, but should be ignore for stubruns 
    """

  script:
    """
    assess_strandedness.py infer_out
    """
}