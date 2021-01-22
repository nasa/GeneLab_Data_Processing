#!/usr/bin/env Rscript

# install.packages("optparse")
# source("http://bioconductor.org/biocLite.R")
# biocLite("limma")

suppressPackageStartupMessages(library("optparse"))

# Read options
option_list = list(
  make_option(
    c("-d", "--exprData"), 
    type = "character", 
    help = "Name of (or path to) the input file (tab delimited .txt file or binary RData object)"
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
    c("-t", "--trend"),
    type = "logical",
    default = "TRUE",
    help = "Logical option to use limma-trend, setting the trend argument for the eBayes function  (default: TRUE)"
  ),
  make_option(
    c("-r", "--rmOutliers"), 
    type = "character", 
    help = "Underscore-delimited list of samples to exclude as outliers from differential expression analysis, matching the sample names in the metadata [ex: GSM1234_GSM1235]"
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

if (is.null(opt$exprData)) {
  print_help(opt_parser)
  stop("No expression data provided", call. = FALSE)
} else
  inFH = opt$exprData

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

suppressPackageStartupMessages(library("limma"))

# Read in ISA tab file and extract sample file
# isaFH = "../metadata/GLDS-4_metadata_GSE18388-ISA/s_GSE18388.txt"
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
  stop("Unable to pull sample names from the study level metadata",
       call. = F)
})


# Read in an expression value txt file
# inFH = "annotExpValues.txt"
if(grepl(".txt$",x = inFH) == TRUE){
  eset = read.delim(inFH,header=T,sep = "\t",stringsAsFactors = F)
  colnames(eset) = gsub("\\.","-",colnames(eset)) # Keep the sample names standardized (if data read in as a text file, hyphens are swapped for periods)
  rownames(eset) = eset[,1]
  eset[,1] = NULL
}else{
  load(inFH)
}

if(nrow(factorValues) != ncol(eset)){
  cat("\nWarning: Number of samples in the expression set not equal to the number of samples in the metadata\n")
}

#From the eset matrix, determine which columns correspond to which factor values
esetSampNames <- colnames(eset)
if ( all(esetSampNames %in% row.names(factorValues)) ){
  # Sample names match exactly between file names and metadata and can be used to order the factors
  factorValues = as.data.frame(factorValues[esetSampNames,])
  row.names(factorValues) = esetSampNames # Reset the row names of the factorValues object in the case where there is only one factor and the row names are lost
} else {
  # Match by non-case-sensitive pattern matching
  newOrder = rep(0,ncol(eset))
  for(i in 1:ncol(eset)){ # Reorder the factorValues dataframe to match the order of sample names in the expression set
    newOrder[i] = grep(pattern = esetSampNames[i], x = row.names(factorValues),ignore.case = T)
  }
  factorValues = as.data.frame(factorValues[newOrder,])
  row.names(factorValues) = esetSampNames # Reset the row names of the factorValues object in the case where there is only one factor and the row names are lost
}

if (!is.null(opt$rmOutliers)){
  outliers = opt$rmOutliers
  outliers = strsplit(outliers, split = "_")[[1]]
}else outliers = list()

group <- rep(0,ncol(eset)) # Create list to hold group assignments for all
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
  cat("The following samples belonged to both groups and not considered in the diffrential expression analysis:\n",rownames(factorValues)[group == 3],"\n")
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

eset = eset[,!(group == 4)] # Remove outliers before fitting the linear model
group = group[!(group == 4)]

# # Troubleshooting print statements
# cat("\nGroup1:\n",colnames(eset)[group == 1],"\n")
# cat("\nGroup2:\n",colnames(eset)[group == 2],"\n")

#Create a design matrix based on the ordering of the columns within eset
group = as.factor(group)
design = model.matrix( ~ 0 + group)

#This part of the script is straight from Limma documentation
fit <- lmFit(eset, design)
contrast.matrix <- makeContrasts(group1-group2,levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
if (opt$trend) {
  fit2 <- eBayes(fit2, trend = TRUE)
} else {
  fit2 <- eBayes(fit2, trend = FALSE)
}


#Here we write the results to a tab delimited text file that is ordered by adjusted p-value
#coef refers to which column is of interest (1 is log2FC), adjust refers to multiple hypothesis testing method ("BH" = Benjamini & Hochberg)
outFH = opt$output
table <- data.frame(topTable(fit2, coef=1, n=Inf, adjust="BH"))
write.table(table,file= outFH,sep="\t", quote = F)
cat("All done! Differential expression information saved to:",outFH,"\n")

