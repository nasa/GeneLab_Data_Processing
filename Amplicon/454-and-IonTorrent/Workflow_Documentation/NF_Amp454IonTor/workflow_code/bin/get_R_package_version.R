#!/usr/bin/env Rscript

# Get versions
VERSIONS <-  sprintf("DECIPHER %s\nbiomformat %s\n", 
		     packageVersion("DECIPHER"), 
		     packageVersion("biomformat"))

# Write versions to file

write(x= VERSIONS, file="versions.txt", append=TRUE)
