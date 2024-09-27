#!/usr/bin/env Rscript

# Get versions
VERSIONS <-  sprintf("vegan %s\ntidyverse %s\ndendextend %s\nphyloseq %s\nDESeq2 %s\nggrepel %s\ndplyr %s\nRColorBrewer %s\ngrid %s\n", 
		     packageVersion("vegan"), 
		     packageVersion("tidyverse"), 
		     packageVersion("dendextend"),
                     packageVersion("phyloseq"),
                     packageVersion("DESeq2"),
                     packageVersion("ggrepel"),
                     packageVersion("dplyr"),
                     packageVersion("RColorBrewer"),
                     packageVersion("grid"))

# Write versions to file

write(x= VERSIONS, file="versions.txt", append=TRUE)
