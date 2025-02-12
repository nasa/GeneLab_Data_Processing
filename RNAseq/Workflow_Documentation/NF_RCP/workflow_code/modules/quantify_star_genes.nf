process QUANTIFY_STAR_GENES {
  // tag "Dataset-wide"
  
  input:
    path("samples.txt")
    path("02-STAR_Alignment/*")
    val(strandedness)

  output:
    tuple path("STAR_Unnormalized_Counts_GLbulkRNAseq.csv"), path("STAR_NumNonZeroGenes_GLbulkRNAseq.csv"), emit: publishables

  script:
    """
    Quantitate_non-zero_genes_per_sample_STAR.R ${strandedness}
    """

}