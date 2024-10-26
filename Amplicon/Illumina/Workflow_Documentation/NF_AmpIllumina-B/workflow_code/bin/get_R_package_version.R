#!/usr/bin/env Rscript

# Get versions
VERSIONS <-  sprintf("dada2 %s\nDECIPHER %s\nbiomformat %s\n", 
		     packageVersion("dada2"), 
		     packageVersion("DECIPHER"), 
		     packageVersion("biomformat"))

# Write versions to file

write(x= VERSIONS, file="versions.txt", append=TRUE)
