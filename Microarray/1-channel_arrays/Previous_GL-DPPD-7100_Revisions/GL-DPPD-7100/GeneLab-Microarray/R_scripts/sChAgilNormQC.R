#!/usr/bin/env Rscript

### Single channel Agilent microarray normalization and quality control

# install.packages("optparse")
# source("http://bioconductor.org/biocLite.R")
# biocLite("limma")
## install.packages("statmod")
# biocLite("arrayQualityMetrics")
# biocLite("oligo")

suppressPackageStartupMessages(library("optparse")) # Load optparse package to read in arguments

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
    default = "expValues",
    help = "Name of the output file [without extension!] (default: expValues)"
  ),
  make_option(
    c("-t", "--outType"),
    type = "character",
    default = "both",
    help = "Format of output data: R (Rdata object), txt (tab delimited file with identifiers and sample names), both (default)"
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
opt = parse_args(opt_parser) # Parse the arguments into a list

#norm = opt$normalization
outFH = opt$outFile
QCout = opt$QCoutput
QCpack = opt$QCpackage # QCpack = "R"


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
  # More robust function to detach loaded packages
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
}

if (is.null(opt$GLDS)) {
  # Include GLDS accession number in outputs if provided
  cat("Warning: No GLDS accession number provided\n")
  if (grepl("GLDS-[0-9]+", inPath)) { # Check if the GLDS number is in the input directory
    glAn = regmatches(inPath, regexpr("GLDS-[0-9]+", inPath)) # Attempt to extract the GLDS accession number from the input path
    cat("Try to guess GLDS accession number... ", glAn,"\n")
  } else{
    glAn = FALSE # Set to false if not provided and not in the input directory
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


raw = read.maimages(
  files = inFiles,
  source = "agilent.median",
  green.only = T,
  path = inPath,
  columns = list(
    G = "gMedianSignal",
    Gb = "gBGMedianSignal"
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
    
    rawData = new("ExpressionSet", exprs = as.matrix(raw)) # Generate a temprory expression set object for performing quality control
    suppressWarnings(
      arrayQualityMetrics(
        expressionset = rawData,
        outdir = paste(qcDir, "raw_report", sep = ""),
        force = T,
        do.logtransform = T
      )
    )
    rm(rawData)
    detach_package(oligo)
  } else if (QCpack == "R") {
    # Create a directory to hold figures from the raw QC step
    rawDir = paste(qcDir, "raw_report/", sep = "")
    if (!file.exists(rawDir)){ # Create a raw report directory within qcDir if it does not exist yet
      dir.create(rawDir)
    }
    
    # Prepare plotting options
    toMatch = c(8,183,31,45,51,100,101,118,128,139,147,183,254,421,467,477,
                483,493,498,503,508,535,552,575,635,655)
    color = grDevices::colors()[rep(toMatch,10)] # Create a library of colors for plotting
    
    # Intensity distributions of the pm probes from each microarray on the same graph
    cat("\tGenerating initial distribution plots")
    mypms = as.matrix(raw)
    png(
      paste(rawDir, glAn, '_rawDensityDistributions.png', sep = ''),
      width = 800,
      height = 800
    )
    ylims = c(0, .8)
    xlims = c(0, 20)
    for (i in 1:ncol(mypms)) {
      cat(".")
      if (i == 1) {
        plot(
          density(log2(mypms[, i])),
          ylim = ylims,
          xlim = xlims,
          xlab = 'log2(Raw Intensities)',
          main = paste(glAn, ' Raw intensity distributions', sep = ''),
          col = color[i]
        )
        par(new = T)
      } else{
        plot(
          density(log2(mypms[, i])),
          ylim = ylims,
          xlim = xlims,
          axes = F,
          xlab = '',
          ylab = '',
          main = '',
          col = color[i]
        )
        par(new = T)
      }
    }
    legend(
      16,
      0.8,
      col = color[1:length(inFiles)],
      legend = sampNames
      ,
      pch = 15,
      bty = "n",
      cex = 0.9,
      pt.cex = 0.8,
      y.intersp = 0.8
    )
    garbage <- dev.off()
    cat("\n")
    
    # Boxplots
    png(paste(rawDir, glAn, '_rawBoxplot.png', sep = ''),
        width = 800,
        height = 400)
    par(mar = c(7, 5, 1, 1))
    boxplot(
      log2(mypms),
      las = 2,
      outline = FALSE,
      col = color[1:length(inFiles)],
      main = paste(glAn, " Raw intensities", sep = ""),
      names = sampNames
    )
    mtext(
      text = "log2 Intensity",
      side = 2,
      line = 2.5,
      las = 0
    )
    garbage <- dev.off()
    
    # PCA
    cat("\tPerforming PCA of raw data...\n")
    rawPCA = prcomp(mypms)
    png(paste(rawDir, glAn, '_rawPCA.png', sep = ''),
        width = 800,
        height = 800)
    plot(
      rawPCA$rotation[, 1],
      rawPCA$rotation[, 2],
      col = color[1:length(inFiles)],
      pch = 16,
      xlab = paste(
        "PC1, ",
        round(summary(rawPCA)$importance["Proportion of Variance", 1] * 100, digits = 1),
        "% of variance",
        sep = ""
      ),
      ylab = paste(
        "PC2, ",
        round(summary(rawPCA)$importance["Proportion of Variance", 2] * 100, digits = 1),
        "% of variance",
        sep = ""
      ),
      main = paste(glAn, " PCA of raw data", sep = "")
    )
    text(
      rawPCA$rotation[, 1],
      rawPCA$rotation[, 2],
      labels = sampNames,
      cex = 1,
      pos = 3
    )
    garbage <- dev.off()
  } else {
    warning("\n--QCpackage option was not recognized and quality control is not being performed.\n", call. = F)
  }
  
}

# Normalize data
cat("Normalizing single channel Agilent microarray data...\n")
cat("\tBackground correcting\n")

normVals = backgroundCorrect(raw, method = "normexp", verbose = F)

cat("\tNormalizing between arrays\n")

normVals = normalizeBetweenArrays(normVals, method = "quantile")

# Saving the normalized data
eset = normVals$E
row.names(eset) = normVals$genes[[1]]

# Create output directory if it does not exist yet
outDir = gsub("(/)[^/]*$", "\\1", outFH) # Strip the filename away from the directory path to the input file
if (grepl("/", outDir)) {
  if (!file.exists(outDir)) {
    # Create the output directory if it does not exist yet
    dir.create(outDir)
  }
}

outType = opt$outType
if (outType == "both") {
  save(eset, file = paste(outFH, ".rda", sep = ""))
  write.table(
    data.frame("ID" = row.names(eset),eset), # provides the rownames as a labeled column in the saved output
    row.names = F,
    file = paste(outFH, ".txt", sep = ""),
    sep = "\t",
    quote = F
  )
  cat("Success! Normalized data saved to", paste(outFH, sep=""), "as both a .txt and a .RData file\n\n")
} else if (outType == "R") {
  save(eset, file = paste(outFH, ".rda", sep = ""))
  cat("Success! Normalized data saved to", paste(outDir, outFH, sep=""), "as a .RData file\n\n")
} else if (outType == "txt") {
  write.table(
    data.frame("ID" = row.names(eset),eset), # provides the rownames as a labeled column in the saved output
    row.names = F,
    file = paste(outFH, ".txt", sep = ""),
    sep = "\t",
    quote = F
  )
  cat("Success! Normalized data saved to", paste(outDir, outFH, sep=""), "as a .txt file\n\n")
} else{
  print_help(opt_parser)
  stop("Help, I don't know how to save this data!", call. = F)
}

# Post-normalization QC step
if(QCout == T) {
  cat("Performing post-normalization QC\n")
  if (QCpack == "aqm"){
    suppressPackageStartupMessages(require(oligo))
    
    normalizedData = new("ExpressionSet", exprs = as.matrix(eset)) # Generate a temprory expression set object for performing quality control
    suppressWarnings(
      arrayQualityMetrics(
        expressionset = normalizedData,
        outdir = paste(qcDir, "normalized_report", sep = ""),
        force = T,
        do.logtransform = T
      )
    )
    rm(normalizedData)
  } else if (QCpack == "R") {
    # Create a directory to hold figures from the normalized QC step
    normDir = paste(qcDir, "normalized_report/", sep = "")
    if (!file.exists(normDir)){ # Create a raw report directory within qcDir if it does not exist yet
      dir.create(normDir)
    }
    
    # Density distributions
    png(
      paste(normDir, glAn, '_normDensityDistributions.png', sep = ''),
      width = 800,
      height = 800
    )
    ylims = c(0, .8)
    xlims = c(0, 20)
    for (i in 1:ncol(eset)) {
      if (i == 1) {
        plot(
          density(eset[, i]),
          ylim = ylims,
          xlim = xlims,
          xlab = 'Normalized expression values[log2]',
          main = paste(glAn, ' Normalized expression distributions', sep = ''),
          col = color[i]
        )
        par(new = T)
      } else{
        plot(
          density(eset[, i]),
          ylim = ylims,
          xlim = xlims,
          axes = F,
          xlab = '',
          ylab = '',
          main = '',
          col = color[i]
        )
        par(new = T)
      }
    }
    legend(
      13,
      0.8,
      col = color[1:length(inFiles)],
      legend = sampNames
      ,
      pch = 15,
      bty = "n",
      cex = 0.9,
      pt.cex = 0.8,
      y.intersp = 0.8
    )
    garbage <- dev.off()
    
    # Boxplots
    png(
      paste(normDir, glAn, '_normBoxplot.png', sep = ''),
      width = 800,
      height = 400
    )
    par(mar = c(7, 5, 1, 1))
    boxplot(
      eset,
      las = 2,
      outline = FALSE,
      col = color[1:length(inFiles)],
      main = paste(glAn, " Normalized intensities", sep = ""),
      names = sampNames
    )
    mtext(
      text = "log2 Intensity",
      side = 2,
      line = 2.5,
      las = 0
    )
    garbage <- dev.off()
    
    # PCA
    cat("\tPerforming PCA of normalized data...\n")
    normPCA = prcomp(eset)
    png(paste(normDir, glAn, '_normPCA.png', sep = ''),
        width = 800,
        height = 800)
    plot(
      normPCA$rotation[, 1],
      normPCA$rotation[, 2],
      col = color[1:length(inFiles)],
      pch = 16,
      xlab = paste(
        "PC1, ",
        round(summary(normPCA)$importance["Proportion of Variance", 1] * 100, digits = 1),
        "% of variance",
        sep = ""
      ),
      ylab = paste(
        "PC2, ",
        round(summary(normPCA)$importance["Proportion of Variance", 2] * 100, digits = 1),
        "% of variance",
        sep = ""
      ),
      main = paste(glAn, " PCA of normalized data", sep = "")
    )
    text(
      normPCA$rotation[, 1],
      normPCA$rotation[, 2],
      labels = sampNames,
      cex = 1,
      pos = 3
    )
    garbage <- dev.off()
  } else {
    
  }
}

if (opt$pullIDs == TRUE) {
  cat("Attempting to extract RefSeq gene IDs from raw microarray file",
      inFiles[1],
      "\n")
  tryCatch({
    txt = read.delim(
      paste(inPath, inFiles[1], sep = ""),
      header = F,
      stringsAsFactors = F,
      blank.lines.skip = F
    )
    skipCnt = max(grep("\\*", txt[, 1]))
    txt = read.delim(
      paste(inPath, inFiles[1], sep = ""),
      header = T,
      stringsAsFactors = F,
      skip = (skipCnt + 1)
    )
  }, error = function(e) {
    stop(
      "Problem reading and parsing raw text file. It may not fit standard Agilent formatting conventions\n"
    )
  })
  tryCatch({
    ID = txt[, 2]
    refInd = grep("^SystematicName$", colnames(txt))
    GB_ACC = txt[, refInd]
    GPL = cbind(ID, GB_ACC)
  }, error = function(e) {
    stop(
      "Problem recognizing column names in the raw text file. Unable to extract annotation information\n"
    )
  })
  if (glAn != F) {
    gplFH = paste(inPath, glAn, "_GPL.txt", sep = "")
  } else {
    gplFH = paste(inPath, "GPL.txt", sep = "")
  }
  write.table(
    GPL,
    row.names = F,
    file = gplFH,
    quote = F,
    sep = "\t"
  )
  cat("Success! Extracted annotation information saved to",
      gplFH,
      "\n")
}
