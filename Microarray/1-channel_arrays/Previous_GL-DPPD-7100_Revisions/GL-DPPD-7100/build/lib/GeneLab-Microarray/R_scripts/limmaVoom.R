#!/usr/bin/env Rscript

# install.packages("optparse")
# source("http://bioconductor.org/biocLite.R")
# biocLite("limma")

suppressPackageStartupMessages(library("optparse"))

# Read options
option_list = list(
  make_option(
    c("-m", "--countsMatrix"), 
    type = "character", 
    help = "Path to a tab-delimited counts matrix .txt file with the sample names in the column headers"
  ),
  make_option(
    c("-f", "--countsFiles"), 
    type = "character", 
    help = "A comma-separated list of tab-delimited counts .txt files to read in and concatenate into a dataframe within R"
  ),
  make_option(
    c("-d", "--inputDirectory"), 
    type = "character", 
    help = "A path to a directory containing a set of tab-delimited counts .txt files to read in and concatenate into a dataframe within R"
  ),
  make_option(
    c("-i", "--ISApath"), 
    type = "character", 
    help = "Path to the file containing the sample-level metadata"
  ),
  make_option(
    "--group1", 
    type = "character", 
    help = "'_'delimited list of factors to select samples for group 1 [ex: flight_geneKO]"
  ),
  make_option(
    "--group2", 
    type = "character", 
    help = "'_'delimited list of factors to select samples for group 2 [ex: ground_geneKO]"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = "DGE.txt",
    help = "Name of (or path to) file to write results to (default: DGE.txt)"
  ),
  make_option(
    c("-r", "--rmOutliers"), 
    type = "character", 
    help = "Underscore-delimited list of samples to exclude as outliers from differential expression analysis, matching the sample names in the metadata [ex: GSM1234_GSM1235]"
  ),
  make_option(
    "--robust",
    type = "logical",
    default = "TRUE",
    help = "Use of robust setting to protect against outlying genes (Recommended, default: TRUE)"
  ),
  make_option(
    c("-n", "--normalization"), 
    type = "character",
    default = "TMM",
    help = "Normalization technique [TMM (default), none (no normalization), upperquartile, or RLE (relative log expression)] (Recommended, default: TRUE)"
  )
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

addSlash = function(string) {
  # Adds a trailing forward slash to the end of a string (ex path to a driectory) if it is not present
  if (substr(x = string, start = nchar(string), stop = nchar(string)) != "/") {
    string = paste(string, "/", sep = "")
  }
  return(string)
}

if (is.null(opt$countsMatrix) & is.null(opt$countsFiles) & is.null(opt$inputDirectory)) {
  print_help(opt_parser)
  stop("No counts data provided", call. = FALSE)
} else {
  if (!is.null(opt$countsMatrix)) {
    inFH = opt$countsMatrix
  } else if (!is.null(opt$countsFiles)) {
    inFH = opt$countsFiles
  } else {
    inFH = opt$inputDirectory
  }
}

if(sum(c(!is.null(opt$countsMatrix), !is.null(opt$countsFiles), !is.null(opt$inputDirectory))) > 1){
  stop("More than one input option used, please review available input options and use only one", call. = FALSE)
}

if (is.null(opt$ISApath)) {
  print_help(opt_parser)
  stop("No ISA file provided", call. = FALSE)
} else {
  isaFH = opt$ISApath
}

# Read in underscore-delimited factor levels and split into lists
if (!is.null(opt$group1) & !is.null(opt$group2)){
  fact1 = opt$group1
  fact1 = strsplit(fact1, split = "_")[[1]]
  fact2 = opt$group2
  fact2 = strsplit(fact2, split = "_")[[1]]
} else{
  print_help(opt_parser)
  stop("Factor levels not provided or improperly formated", call.=FALSE)
}

outFH = opt$output

suppressPackageStartupMessages(library("limma"))
suppressPackageStartupMessages(library("edgeR"))

# Read in ISA tab file and extract sample file
# isaFH = "../../GLDS-101_metadata_ISA/s_BIOBANK.txt"
tryCatch({
  studyFactors = read.delim(isaFH, header=T, sep="",stringsAsFactors = F)
}, error=function(e){ 
  stop("ISA sample file could not be read by parsing tab-delimited files", call. = F)
})

# From sample file, extract column containing 'Factor Value'
tryCatch({
  factorValues = as.data.frame(studyFactors[, grepl("Factor.Value", colnames(studyFactors))])
  row.names(factorValues) = studyFactors[, grepl("^Sample.Name$", colnames(studyFactors))]
  {
    # Match sample names in file name
    ## [replace('_','-').replace('(','-').replace(')','-').replace(' ','-').strip('-')]
    replaceWithHyphen = c("_", "\\(", "\\)", " ", "\\.")
    removeList = c("^-", "-$")
    for (i in 1:length(replaceWithHyphen)) {
      # Replace other characters with hyphens
      row.names(factorValues) =  gsub(replaceWithHyphen[i], '-', row.names(factorValues))
    }
    for (i in 1:length(removeList)) {
      # Remove leading/trailing hypens
      row.names(factorValues) =  gsub(removeList[i], '', row.names(factorValues))
    }
    
  }
}, error = function(e) {
  stop("Unable to pull factors or sample names from the study level metadata",
       call. = F)
})

# Read counts data
# wd = "/Users/dmattox/Documents/genelab/RNAseq/GLDS-101/RNAseq/Feature counts/"
# inFH = dir(wd)
if (!is.null(opt$inputDirectory)) {
  inFH = dir(inFH)
  inFH = inFH[grepl("\\.tabular", inFH)]
  
  sampNames = gsub("_transcriptomics_.*", "", inFH)
  sampNames = gsub("\\.txt", "", sampNames)
  sampNames = gsub("\\.tab", "", sampNames)
  sampNames = gsub("\\.tabular", "", sampNames)
  sampNames = gsub(".*/", "", sampNames)
  sampNames = gsub("GLDS-\\d*_", "", sampNames) # Extract sample names from the list of input files
  
  cat("Reading in:", inFH[1],"\n")
  tmp = read.delim(inFH[1], header=T, sep = "\t", stringsAsFactors = F)
  noID = grepl("^$",x = tmp[,1]) # Tag missing gene IDs ("")
  if (any(noID)) {
    tmp = tmp[!noID,]
  }
  cnts = as.data.frame(matrix(nrow = nrow(tmp), ncol = length(inFH)),row.names = tmp[,1])
  cnts[,1] = tmp[,2]
  colnames(cnts) = sampNames
  for (i in 2:length(inFH)) {
    cat("Reading in:", inFH[i],"\n")
    tmp = read.delim(inFH[i], header=T, sep = "\t", stringsAsFactors = F)
    if (any(noID)) {
      tmp = tmp[!noID,] # Remove the same row from every separate counts file
    }
    cnts[,i] = tmp[,2]
  }
  cat("\nWarning:",sum(noID),"gene ID(s) were found missing and removed\n")
} else if (!is.null(opt$countsFiles)) {
  # Read in and concatenate counts data from a list of individual files
  # Assumes each file for a study has the same number of rows and they are in the same order
  inFH = strsplit(inFH, split = ",")[[1]]
  
  sampNames = gsub("_transcriptomics_.*", "", inFH)
  sampNames = gsub("\\.txt", "", sampNames)
  sampNames = gsub("\\.tab", "", sampNames)
  sampNames = gsub("\\.tabular", "", sampNames)
  sampNames = gsub(".*/", "", sampNames)
  sampNames = gsub("GLDS-\\d*_", "", sampNames) # Extract sample names from the list of input files
  
  cat("Reading in:", inFH[1],"\n")
  tmp = read.delim(inFH[1], header=T, sep = "\t", stringsAsFactors = F)
  noID = grepl("^$",x = tmp[,1]) # Tag missing gene IDs ("")
  if (any(noID)) {
    tmp = tmp[!noID,]
  }
  cnts = as.data.frame(matrix(nrow = nrow(tmp), ncol = length(inFH)),row.names = tmp[,1])
  cnts[,1] = tmp[,2]
  colnames(cnts) = sampNames
  for (i in 2:length(inFH)) {
    cat("Reading in:", inFH[i],"\n")
    tmp = read.delim(inFH[i], header=T, sep = "\t", stringsAsFactors = F)
    if (any(noID)) {
      tmp = tmp[!noID,] # Remove the same row from every separate counts file
    }
    cnts[,i] = tmp[,2]
  }
  cat("\nWarning:",sum(noID),"gene ID(s) were found missing and removed\n")
} else {
  # If data read in as a counts matrix
  cat("Reading in:", inFH,"\n")
  cnts = read.delim(inFH, header = T, sep = "\t", stringsAsFactors = F)
  row.names(cnts) = cnts[,1]
  cnts = cnts[,-1]
  sampNames = colnames(cnts)
}

# Normalize data and create DGEList object
dge = DGEList(cnts)
if (opt$normalization %in% c("TMM", "RLE", "upperquartile", "none")){
  dge = calcNormFactors(dge, method = opt$normalization)
  cat("Normalized with method:",opt$normalization,"\n\n")
} else{
  warning("Normalization method not recognized, data was not normalized", call. = F)
}



# Link metadata to samples
if(nrow(factorValues) != ncol(cnts)){
  cat("\nWarning: Number of samples in the expression set not equal to the number of samples in the metadata\n")
}

if ( all(sampNames %in% row.names(factorValues)) ){
  # Sample names match exactly between file names and metadata and can be used to order the factors
  factorValues = as.data.frame(factorValues[sampNames,])
  row.names(factorValues) = sampNames # Reset the row names of the factorValues object in the case where there is only one factor and the row names are lost
} else {
  # Match by non-case-sensitive pattern matching
  newOrder = rep(0,ncol(cnts))
  for(i in 1:ncol(cnts)){ # Reorder the factorValues dataframe to match the order of sample names in the expression set
    newOrder[i] = grep(pattern = sampNames[i], x = row.names(factorValues),ignore.case = T)
  }
  factorValues = as.data.frame(factorValues[newOrder,])
  row.names(factorValues) = sampNames # Reset the row names of the factorValues object in the case where there is only one factor and the row names are lost
}

if (!is.null(opt$rmOutliers)){
  outliers = opt$rmOutliers
  outliers = strsplit(outliers, split = "_")[[1]]
}else outliers = list()

group <- rep(0,ncol(cnts)) # Create list to hold group assignments for all
for(i in 1:nrow(factorValues)){ # Assign each sample to a group [3 = both groups 1 & 2, 0 = neither group, 4 = designated as an outlier]
  if(row.names(factorValues)[i] %in% outliers){
    group[i] = 4
  }else if(all(fact1 %in% factorValues[i,]) & all(fact2 %in% factorValues[i,])){
    group[i] = 3
  }else if(all(fact1 %in% factorValues[i,])){
    group[i] = 1
  }else if(all(fact2 %in% factorValues[i,])){
    group[i] = 2
  }
}
# Error handling
if(sum(group == 4) > 0){
  cat("The following samples were indicated to be outliers and were removed from further analysis:\n",rownames(factorValues)[group == 4],"\n")
}
if(sum(group == 3) > 0){
  cat("The following samples belonged to both groups and were removed from further analysis:\n",rownames(factorValues)[group == 3],"\n")
}
if(sum(group == 3) == nrow(factorValues)){
  stop("All of the samples belonged to both groups! Exiting.", call.=F)
}
cat(sum(group == 1),"sample(s) found in group 1:\n",rownames(factorValues)[group == 1],"\n")
cat(sum(group == 2),"sample(s) found in group 2:\n",rownames(factorValues)[group == 2],"\n")
if(sum(group == 0) > 0){
  cat("Warning:",sum(group == 0),"sample(s) not found in either group:\n",rownames(factorValues)[group == 0],"\nIf this is not expected, please ensure the provided factor levels match the factor levels in the study-level metadata exactly\n\n")
}

if(sum(group == 1) == 0 | sum(group == 2) == 0){
  stop("One or both comparison groups were found to be empty. Exiting...", call.=F)
}

dge = dge[,!(group == 4)] # Remove outliers before fitting the linear model
group = group[!(group == 4)]

#Create a design matrix
group = as.factor(group)
design = model.matrix( ~ 0 + group)

# Create voom object
cat("\nApplying voom transformation and performing differential expression analysis...\n")
v = voom(dge, design, plot = F)

# Standard limma differential gene expression analysis
fit = lmFit(v, design)
contrast.matrix <- makeContrasts(group1-group2,levels=design)
fit2 = contrasts.fit(fit, contrast.matrix)
if (opt$robust) {
  fit2 = eBayes(fit2, robust = TRUE)
} else {
  fit2 = eBayes(fit2, robust = FALSE)
}

table = data.frame(topTable(fit2, coef=1, n=Inf, adjust="BH"))
write.table(table,file= outFH,sep="\t", quote = F)
cat("\nAll done! Differential expression information saved to:",outFH,"\n\n")











