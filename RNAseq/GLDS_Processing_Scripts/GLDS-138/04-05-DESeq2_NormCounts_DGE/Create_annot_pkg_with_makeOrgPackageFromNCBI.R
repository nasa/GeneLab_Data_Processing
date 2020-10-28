## Modify makeOrgPackageFromNCBI to work correctly and create an annotation package ##

# Load the AnnotationForge package
library("AnnotationForge")

# Install and load a "GO.db" package
BiocManager::install("GO.db")
library("GO.db")

# Open up the loaded function to make necessary changes. 
# In this script we remove the gene2unigene parts of the script because it points to a gz file on a website that no longer exists.
trace(".primaryFiles", edit=TRUE, where = makeOrgPackageFromNCBI)
# click "Save" on that editor when complete and it will close

# We also have to replace "entrezgene" to "entrezgene_id" in the getEnsemblData portion of the script
trace(".getEnsemblData", edit=TRUE, where = makeOrgPackageFromNCBI)
# Then click "Save"

# One more place where the script calls for unigene that needs to be removed
trace("prepareDataFromNCBI", edit=TRUE, where = makeOrgPackageFromNCBI)
# Then click "Save"

# And now you can run the fulction. Fill in the authro and maintainer info before running. Note that the version field has to be #.#
makeOrgPackageFromNCBI(version="1.0",
                       author = "FirstName LastName <username@domain.com>",
                       maintainer="FirstName LastName <username@domain.com>",
                       outputDir = ".",
                       tax_id = "224308",
                       genus = "Bacillus",
                       species = "subtilis")


## Print session info
sessionInfo()
