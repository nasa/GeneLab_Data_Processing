# install-org-db.R

# Function: Get annotations db from ref table. If no annotations db is defined, create the package name from genus, species, (and strain for microbes), 
# Try to Bioconductor install annotations db. If fail then build the package using AnnotationForge, install it into the current directory.
# Requires ~80GB for NCBIFilesDir file caching
install_annotations <- function(target_organism, refTablePath) {
    if (!file.exists(refTablePath)) {
        stop("Reference table file does not exist at the specified path: ", refTablePath)
    }

    ref_table <- read.csv(refTablePath)
    target_taxid <- ref_table %>%
        filter(species == target_organism) %>%
        pull(taxon)
    
    # Parse organism's name in the reference table to create the org.db name (target_org_db)
    target_species_designation <- ref_table %>%
        filter(species == target_organism) %>%
        pull(species) %>%
        gsub("\\s+", " ", .) %>%
        gsub("[^A-Za-z0-9 ]", "", .)
    
    genus_species <- strsplit(target_species_designation, " ")[[1]]
    if (length(genus_species) < 1) {
        stop("Species designation is not correctly formatted: ", target_species_designation)
    }
    
    genus <- genus_species[1]
    species <- ifelse(length(genus_species) > 1, genus_species[2], "")
    strain <- ref_table %>%
        filter(species == target_organism) %>%
        pull(strain) %>%
        gsub("[^A-Za-z0-9]", "", .)
    
    if (!is.na(strain) && strain != "") {
        species <- paste0(species, strain)
    }
    
    # Get package name or build it if not provided
    target_org_db <- ref_table %>%
        filter(species == target_organism) %>%
        pull(annotations)
    
    if (is.na(target_org_db) || target_org_db == "") {
        cat("\nNo annotation database specified. Constructing package name...\n")
        target_org_db <- paste0("org.", substr(genus, 1, 1), species, ".eg.db")
    }
    
    cat(paste0("\nChecking Bioconductor for '", target_org_db, "'...\n"))
    if (requireNamespace(target_org_db, quietly = TRUE)) {
        cat(paste0("'", target_org_db, "' is already installed.\n"))
    } else {
        cat(paste0("\nAttempting to install '", target_org_db, "' from Bioconductor...\n"))
        BiocManager::install(target_org_db, ask = FALSE)
        if (requireNamespace(target_org_db, quietly = TRUE)) {
            cat(paste0("'", target_org_db, "' has been successfully installed from Bioconductor.\n"))
        } else {
            cat(paste0("\nInstallation from Bioconductor failed, attempting to build '", target_org_db, "'...\n"))
            if (!dir.exists(target_org_db)) {
                tryCatch({
                    BiocManager::install(c("AnnotationForge", "biomaRt", "GO.db"), ask = FALSE)
                    library(AnnotationForge)
                    makeOrgPackageFromNCBI(
                        version = "0.1",
                        author = "Your Name <your.email@example.com>",
                        maintainer = "Your Name <your.email@example.com>",
                        outputDir = "./",
                        tax_id = target_taxid,
                        genus = genus,
                        species = species
                    )
                    install.packages(file.path("./", target_org_db), repos = NULL, type = "source", quiet = TRUE)
                    cat(paste0("'", target_org_db, "' has been successfully built and installed.\n"))
                }, error = function(e) {
                    stop("Failed to build and load the package: ", target_org_db, "\nError: ", e$message)
                })
            } else {
                cat(paste0("Local annotation package ", target_org_db, " already exists. This local package will be installed.\n"))
                install.packages(file.path("./", target_org_db), repos = NULL, type = "source", quiet = TRUE)
            }
        }
    }

    library(target_org_db, character.only = TRUE)
    cat(paste0("Using Annotation Database '", target_org_db, "'.\n"))
    return(target_org_db)
}
