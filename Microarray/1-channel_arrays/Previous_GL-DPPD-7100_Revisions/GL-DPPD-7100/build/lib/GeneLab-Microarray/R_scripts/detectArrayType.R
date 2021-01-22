#!/usr/bin/env Rscript

# source("http://bioconductor.org/biocLite.R")
# biocLite("Risa")

require(Risa)

rsaIN = readISAtab("~/Documents/genelab/rot1/GLDS-4/metadata/GLDS-4_metadata_GSE18388-ISA/")

sampFactors = rsaIN@study.files[[1]]
assayFactors = rsaIN@assay.files[[1]]

# May not be feasible tro detect array type but above command will extract factors