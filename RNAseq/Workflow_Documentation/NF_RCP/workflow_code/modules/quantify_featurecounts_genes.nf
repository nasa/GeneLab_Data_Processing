process QUANTIFY_FEATURECOUNTS_GENES {
  // tag "Dataset-wide"
  // An R script that gets the number of non-zero genes per sample in a featureCounts output table
  input:
    path("samples.txt")
    path("03-FeatureCounts/*")

  output:
    path("NumNonZeroGenes${params.assay_suffix}.csv"), emit: num_non_zero_genes

  script:
    """
    Quantitate_non-zero_genes_per_sample_Featurecounts.R
    """

}