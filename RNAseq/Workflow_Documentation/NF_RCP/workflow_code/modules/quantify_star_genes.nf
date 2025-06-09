process QUANTIFY_STAR_GENES {
  // tag "Dataset-wide"
  
  input:
    path("samples.txt")
    path("02-STAR_Alignment/*")
    val(strandedness)

  output:
    tuple path("STAR_Unnormalized_Counts${params.assay_suffix}.csv"), path("STAR_NumNonZeroGenes${params.assay_suffix}.csv"), emit: publishables

  script:
    """
    Quantitate_non-zero_genes_per_sample_STAR.R ${strandedness}
    """

}