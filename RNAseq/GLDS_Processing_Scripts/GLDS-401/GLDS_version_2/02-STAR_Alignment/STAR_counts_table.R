print("Make STAR counts table")
print("")

work_dir="/GLDS-401/processing_scripts/02-STAR_Alignment/"
align_dir="/GLDS-401/02-STAR_Alignment"

setwd(file.path(work_dir))

### Pull in sample names ###
study <- read.csv(Sys.glob(file.path(work_dir,"samples.txt")), header = FALSE, row.names = 1, stringsAsFactors = TRUE)

##### Import Data
ff <- list.files(file.path(align_dir), pattern = "ReadsPerGene.out.tab", recursive=TRUE, full.names = TRUE)

## Reorder the *genes.results files to match the ordering of the ISA samples
ff <- ff[sapply(rownames(study), function(x)grep(paste0(x,'_ReadsPerGene.out.tab$'), ff, value=FALSE))]

# Remove the first 4 lines
counts.files <- lapply( ff, read.table, skip = 4 )

# Get counts aligned to either strand for undtranded data by selecting col 2, to the first (forward) strand by selecting col 3 or to the second (reverse) strand by selecting col 4
counts <- as.data.frame( sapply( counts.files, function(x) x[ , 4 ] ) )

# Add column and row names
colnames(counts) <- rownames(study)
row.names(counts) <- counts.files[[1]]$V1


##### Export unnormalized counts table
setwd(file.path(align_dir))
write.csv(counts,file='STAR_Unnormalized_Counts.csv')


## print session info ##
print("Session Info below: ")
print("")
sessionInfo()
