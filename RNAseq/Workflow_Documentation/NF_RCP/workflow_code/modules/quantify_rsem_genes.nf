process QUANTIFY_RSEM_GENES {
  // tag "Dataset-wide"
  // An R script that extracts gene counts by sample to a table

  input:
    path("samples.txt")
    path("03-RSEM_Counts/*")

  output:
    tuple path("RSEM_Unnormalized_Counts_GLbulkRNAseq.csv"), path("RSEM_NumNonZeroGenes_GLbulkRNAseq.csv"), emit: publishables

  script:
    """
    Quantitate_non-zero_genes_per_sample_RSEM.R
    """

}