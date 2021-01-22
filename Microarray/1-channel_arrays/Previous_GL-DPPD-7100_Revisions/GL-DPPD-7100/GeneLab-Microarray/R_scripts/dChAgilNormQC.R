#!/usr/bin/env Rscript

### Dual channel Agilent microarray normalization and quality control

# install.packages("optparse")
# source("http://bioconductor.org/biocLite.R")
# biocLite("limma")
# biocLite("arrayQualityMetrics")
# biocLite("oligo")

suppressPackageStartupMessages(library("optparse"))

# Read options
option_list = list(
  make_option(
    c("-i", "--input"), 
    type = "character", 
    help = "Path to directory containing input raw array files"
  ),
  # make_option(
  #   c("-n", "--normalization"),
  #   type = "character",
  #   default = "normexp",
  #   help = "Normalization method [rma (default, full rma), quantile (no background correction), background (no quant. normalization), log2 (no quant. norm. or background correction)"
  # ),
  make_option(
    c("-o", "--outFile"),
    type = "character",
    default = "normMA",
    help = "Name of the output file [without extension!] (default: normMA)"
  ),
  make_option(
    "--QCoutput",
    type = "logical",
    default = TRUE,
    help = "Output QC_reporting directory of QC plots (default = TRUE)"
  ),
  make_option(
    "--QCpackage",
    type = "character",
    default = "aqm",
    help = "Package used to generate QC plots: aqm (generate html report with arrayQualityMetrics package,default), R (use standard R plotting options to generate QC plots as .png files. Figures may not maintain proper formatting for datasets with many samples)"
  ),
  make_option(
    "--QCDir",
    type = "character",
    default = "./QC_reporting/",
    help = "Path to directory for storing QC output, including a terminal forward slash. Will be created if it does not exist yet (default = './QC_reporting/')"
  ),
  make_option(
    "--pullIDs", 
    type = "logical", 
    default = FALSE, 
    help = "Logical option to try and extract RefSeq gene IDs from a raw text file. It will write a tab-delimited GPL file into the directory containing the raw files (default: false)"
  ),
  make_option(
    "--GLDS", 
    type = "character",
    help = "Full accession number for QC outputs (ie 'GLDS-21' for GLDS-21)"
  )
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

outFH = opt$outFile
QCout = opt$QCoutput
QCpack = opt$QCpackage


addSlash = function(string) {
  # Adds a trailing forward slash to the end of a string (ex path to a driectory) if it is not present
  if (substr(x = string,
             start = nchar(string),
             stop = nchar(string)) != "/") {
    string = paste(string, "/", sep = "")
  }
  return(string)
}

detach_package = function(pkg, character.only = FALSE) {
  if (!character.only)
  {
    pkg <- deparse(substitute(pkg))
  }
  search_item <- paste("package", pkg, sep = ":")
  while (search_item %in% search())
  {
    detach(search_item,
           unload = TRUE,
           character.only = TRUE)
  }
}

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("No path to input directory provided. Please look over the available options", call. = F)
}else{
  inPath = addSlash(opt$input)
  # inPath = "~/Documents/genelab/rot1/GLDS-28/microarray/"
  # setwd(inPath)
}

if (is.null(opt$GLDS)) {
  # Include GLDS accession number in outputs if provided
  cat("Warning: No GLDS accession number provided\n")
  if (grepl("GLDS-[0-9]+", inPath)) {
    glAn = regmatches(inPath, regexpr("GLDS-[0-9]+", inPath)) # Attempt to extract the GLDS accession number from the input path
    cat("Try to guess GLDS accession number... ", glAn,"\n")
  } else{
    glAn = FALSE
  }
} else{
  glAn = opt$GLDS
}

# Load packages
suppressPackageStartupMessages(library("limma"))

inFiles = dir(inPath)
inFiles = inFiles[grepl("_raw.txt$",inFiles)]

if (length(inFiles) > 0){
  cat("Detected raw files:\n")
  for (i in 1:length(inFiles)) {
    cat("\t",inFiles[i],"\n")
  }
  cat("\n")
} else {
  stop("No raw files detected in the current directory",
       call. = F)
}

sampNames = gsub("microarray","",inFiles)
sampNames = gsub("\\..*","",sampNames)
sampNames = gsub(".*/", "", sampNames)
sampNames = gsub("GLDS-\\d*_", "", sampNames)
sampNames = gsub("raw", "", sampNames)
sampNames = gsub("_", "", sampNames) # Extract sample names from the list of filenames

RG = read.maimages(
  files = inFiles,
  source = "agilent.median",
  path = inPath,
  columns = list(
    G = "gMedianSignal",
    Gb = "gBGMedianSignal",
    R = "rMedianSignal",
    Rb = "rBGMedianSignal"
  ),
  annotation = "FeatureNum",
  names = sampNames
)

# Pre-normalization QC step
qcDir = addSlash(opt$QCDir)
if (!file.exists(qcDir)){ # Create QC directory if it does not exist yet
  dir.create(qcDir)
}

if(QCout == T) {
  cat("Performing intial QC\n")
  if (QCpack == "aqm"){
    
    # Load in QC package
    suppressPackageStartupMessages(require(arrayQualityMetrics))
    suppressPackageStartupMessages(require(oligo))
    
    suppressWarnings(
      arrayQualityMetrics(
        expressionset = RG,
        outdir = paste(qcDir, "raw_report", sep = ""),
        force = T
      )
    )
    detach_package(oligo)
  } else {
    warning("\n--QCpackage option was not recognized and quality control is not being performed.\n", call. = F)
  }
}

# Normalization
cat("Normalizing two channel Agilent microarray data...\n")
cat("\tBackground correcting\n")

RGb = backgroundCorrect(RG, method = "normexp") # normexp method is based on the same normal plus exponential convolution model which has previously been used to background correct Affymetrix data as part of the popular RMA algorithm [Ritchie et al. Bioinformatics 2007]

cat("\tNormalizing within arrays\n")

MA = normalizeWithinArrays(RG, method = "loess", weights=NULL) # Agilent specific global loess normalization method
## Loess normalization assumes that the bulk of the probes on the array are not differentially expressed

cat("\tNormalizing between arrays\n")

MA = normalizeBetweenArrays(MA, method = "Aquantile") # Normalize the average intensities ("A") between arrays

# Saving the normalized data
# Create output directory if it does not exist yet
outDir = gsub("(/)[^/]*$", "\\1", outFH) # Strip the filename away from the directory path to the input file
if (grepl("/", outDir)) {
  if (!file.exists(outDir)) {
    # Create the output directory if it does not exist yet
    dir.create(outDir)
  }
}
# outFH = "normMA"
save(MA, file = paste(outFH, ".rda", sep = ""))

cat("Success! Normalized data saved to", outFH, "as a .RData file\n\n")

# Post-normalization QC step
if(QCout == T) {
  cat("Performing post-normalization QC\n")
  if (QCpack == "aqm"){
    suppressPackageStartupMessages(require(oligo))
    
    suppressWarnings(
      arrayQualityMetrics(
        expressionset = MA,
        outdir = "normalized_report",
        force = T
      )
    )
  }
}

if (opt$pullIDs == TRUE) {
  cat("Attempting to extract RefSeq gene IDs from raw microarray file", inFiles[1],"\n")
  tryCatch({
    txt = read.delim(
      paste(inPath, inFiles[1], sep = ""),
      header = T,
      stringsAsFactors = F,
      blank.lines.skip = F
    )
    skipCnt = max(grep("\\*", txt[, 1]))
    txt = read.delim(
      paste(inPath, inFiles[1], sep = ""),
      header = T,
      stringsAsFactors = F,
      skip = (skipCnt + 3),
      blank.lines.skip = T
    )
  }, error = function(e) {
    stop(
      "Problem reading and parsing raw text file. It may not fit standard Agilent formatting conventions\n"
    )
  })
  tryCatch({
    ID = txt[,2]
    refInd = grep("^SystematicName$", colnames(txt))
    GB_ACC = txt[,refInd]
    GPL = cbind(ID,GB_ACC)
  }, error = function(e) {
    stop("Problem recognizing column names in the raw text file. Unable to extract annotation information\n")
  })
  if (glAn != F) {
    gplFH = paste(inPath,glAn,"_GPL.txt",sep="")
  } else {
    gplFH = paste(inPath,"GPL.txt",sep="")
  }
  write.table(
    GPL,
    row.names = F,
    file = gplFH,
    quote = F,
    sep = "\t"
  )
  cat("Success! Extracted annotation information saved to",gplFH,"\n")
}




















