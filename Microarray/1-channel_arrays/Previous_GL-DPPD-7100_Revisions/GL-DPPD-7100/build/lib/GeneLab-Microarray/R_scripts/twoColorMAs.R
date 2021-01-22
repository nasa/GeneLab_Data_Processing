#!/usr/bin/env Rscript

# install.packages("optparse")
# source("http://bioconductor.org/biocLite.R")
# biocLite("limma")
## install.packages("statmod")
# biocLite("arrayQualityMetrics")

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
    default = "expValues",
    help = "Name of the output file [without extension!] (default: not_sure_yet)"
  ),
  make_option(
    c("-i", "--ISApath"), 
    type = "character", 
    help = "Path to the directory containing the dataset metadata"
  ),
  make_option(
    "--QCDir",
    type = "character",
    default = "./QC_reporting/",
    help = "Path to directory for storing QC output, including a terminal forward slash. Will be created if it does not exist yet (default = './QC_reporting/')"
  ),
  make_option(
    c("-g", "--gpl"), 
    type = "character", 
    help = "Path to the file containing custom array annotation information"
  ),
  make_option(
    "--GLDS", 
    type = "character", 
    help = "Full accession number for QC outputs (ie 'GLDS-21' for GLDS-21)"
  )
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

#norm = opt$normalization
outFH = opt$outFile
if (!is.null(opt$gpl)){
  annotFH = opt$gpl # annotFH = "GPL15420_features_probes.txt"
}

addSlash = function(string) {
  # Adds a trailing forward slash to the end of a string (ex path to a driectory) if it is not present
  if (substr(x = string,
             start = nchar(string),
             stop = nchar(string)) != "/") {
    string = paste(string, "/", sep = "")
  }
  return(string)
}

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("No path to input directory provided. Please look over the available options", call. = F)
}else{
  inPath = opt$input
  #setwd(inPath)
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

#relDir = getwd() # Save the path to the original directory (for relative path handling)

if (is.null(opt$ISApath)) {
  # Check if metadata is provided
  print_help(opt_parser)
  stop("No path to metadata directory provided. Please look over the available options",
       call. = F)
} else{
  isaPath = opt$ISApath
  isaPath = addSlash(isaPath)
  tryCatch({
    isaFiles = dir(isaPath) # Identify files in metadata directory
    sFile = isaFiles[grep("^s_*", isaFiles)] # 
    aFile = isaFiles[grep("^a_*", isaFiles)]
    sampFactors = read.delim(
      paste(isaPath, sFile, sep = ""),
      header = T,
      sep = "",
      stringsAsFactors = F
    )
    assFactors = read.delim(
      paste(isaPath, aFile, sep = ""),
      header = T,
      sep = "",
      stringsAsFactors = F
    )
  }, error = function(e) {
    stop("ISA files could not be read by parsing the tab-delimited files",
         call. = F)
  })
}

# Link samples to factors and labels
factors = merge(x = sampFactors, y = assFactors, by = "Sample.Name") # Join the information from the ISA tab delimited assay and sample files

suppressPackageStartupMessages(library("limma"))
suppressPackageStartupMessages(library("arrayQualityMetrics"))

# inPath = "~/Documents/genelab/rot1/GLDS-28/microarray/"
# setwd(inPath)

inFiles = dir(inPath)
inFiles = inFiles[grepl("_raw.txt$",inFiles)]

targets = data.frame(FileName = inFiles)
row.names(targets) = gsub(".*_microarray_","",targets$FileName)
row.names(targets) = gsub("_raw.*","",row.names(targets))

# Pull out sample descriptions, names, and dyes from the metadata
factNames = factors[ , grep("^Sample.Name$",colnames(factors),ignore.case = T)]
labs = factors[ , grep("^Label$",colnames(factors),ignore.case = T)]
facts = factors[ , grep("Comment..Sample_source_name.",colnames(factors),ignore.case = T)]

Cy3 = c("space","mars") # need to link file names to sample names in metadata
Cy5 = c("earth","earth")
targets = cbind(targets, Cy3, Cy5)

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
  names = row.names(targets)
)

# QC step
qcDir = opt$QCDir
suppressWarnings(
  arrayQualityMetrics(expressionset = RG,
                      outdir = paste(qcDir, "raw_report", sep = ""),
                      force = T)
  
)

# Normalization
RGb = backgroundCorrect(RG, method = "normexp") # normexp method is based on the same normal plus exponential convolution model which has previously been used to background correct Affymetrix data as part of the popular RMA algorithm [Ritchie et al. Bioinformatics 2007]
MA = normalizeWithinArrays(RG, method = "loess", weights=NULL) # Agilent specific global loess normalization method
## Loess normalization assumes that the bulk of the probes on the array are not differentially expressed
MAq = normalizeBetweenArrays(MA, method = "Aquantile") # Normalize the average intensities ("A") between arrays

# QC step
suppressWarnings(
  arrayQualityMetrics(expressionset = RG,
                      outdir = paste(opt$QCDir, "raw_report", sep = ""),
                      force = T)
  
)

# Feature number to BSU locus
annotFun = function(oldLab, newLink, newLab, ...){
  # Function to switch from one probe label (oldLab) to another (newLab), provided
  ## a linking variable (newLink) with common labels to the labels in the oldLab list organized in the same order as the newLab list
  # Returns a translation from the the oldLab to the newLab with "" for unmapped probes, ordered by the oldLab
  newOrder = match(oldLab,newLink)
  out = newLab[newOrder]
  return(out)
}

annot = read.delim(annotFH,header = T, stringsAsFactors = F)
BSUs = annotFun(
  oldLab = MAq$genes[,1], 
  newLink = annot$FeatureNum,
  newLab = annot$BSU
  )
noIDCnt = sum(BSUs == "") # Number of unmapped probes
MAq$genes[,1] = BSUs # Translate probe IDs
MAq = MAq[!MAq$genes[,1] == "", ] # Remove umapped probes
BSUs = BSUs[!BSUs == ""]

# Filter out probes with missing values
naTag = (apply(is.na(MAq$M), 1, any) | apply(is.na(MAq$A), 1, any))
MAq = MAq[!naTag,]
BSUs = BSUs[!naTag]

if (any(opt$dupProbes %in% c("average", "max"))) {
  cat("Filtering out multiple probes per gene ID...\n")
  
  if (opt$dupProbes == "average") {
    # Collapse multiple probes per gene ID by averaging expression values across all samples
    rmRowTag = rep(TRUE, length(MAq$genes[,1])) # Tag rows to drop (set single or averaged probes to FALSE below)
    for (i in 1:length(rmRowTag)) {
      if (sum(BSUs == BSUs[i]) > 1) {
        inds = grep(BSUs[i], BSUs) # List of indices at which a probe for a given gene ID occur
        MAq$M[inds[1], ] = apply(X = MAq$M[inds, ],
                                FUN = mean,
                                MARGIN = 2) # Changes the values of the first occurence of a probe to the sample-specific average of the values from all the probes for that gene ID
        MAq$A[inds[1], ] = apply(X = MAq$A[inds, ],
                                 FUN = mean,
                                 MARGIN = 2) # Changes the values of the first occurence of a probe to the sample-specific average of the values from all the probes for that gene ID
        rmRowTag[inds[1]] = FALSE
      } else
        rmRowTag[i] = FALSE
    }
    nDups = sum(rmRowTag)
    MAq = MAq[!rmRowTag, ]

    cat("\tUnampped probes removed:", noIDCnt, "\n")
    cat("\tProbes with missing values removed:", sum(naTag), "\n")
    cat("\tDuplicated probes removed:", nDups, "\n\n")
    cat("Annotated probes remaining:", length(MAq$genes[,1]), "\n\n")
    if (length(MAq$genes[,1]) > length(unique(MAq$genes[,1]))) {
      cat("\n\tWarning: non-unique probe to ID mappings remain \n")
    }
    
  } else if (opt$dupProbes == "max") {
    # Collapse multiple probes per gene ID by selecting a representative with the highest mean intensity across all samples
    rmRowTag = rep(TRUE, length(MAq$genes[,1])) # Tag rows to drop (set single or highest expressing probes to FALSE below)
    for (i in 1:length(MAq$genes[,1])) {
      if (sum(BSUs == BSUs[i]) > 1) {
        inds = grep(BSUs[i], BSUs)
        top = 0
        keep = 0
        for (j in 1:length(inds)) {
          curr = mean(as.numeric(MAq$A[inds[j], ]))
          if (is.na(curr)){ # Check if at least one of the samples has a missing value for the current probe
            curr = 0
          }
          if (curr > top) {
            top = curr
            keep = inds[j]
          }
        }
        rmRowTag[keep] = FALSE
      } else
        rmRowTag[i] = FALSE
    }
    nDups = sum(rmRowTag)
    MAq = MAq[!rmRowTag, ]

    cat("\tUnampped probes removed:", noIDCnt, "\n")
    cat("\tProbes with missing values removed:", sum(naTag), "\n")
    cat("\tDuplicated probes removed:", nDups, "\n\n")
    cat("Annotated probes remaining:", length(MAq$genes[,1]), "\n\n")
    if (length(MAq$genes[,1]) > length(unique(MAq$genes[,1]))) {
      cat("\n\tWarning: non-unique probe to ID mappings remain \n")
    }
  }
} else{
  stop("Method for dealing with probes mapped to the same gene IDs not recognized\n",
       call. = F)
}



# Differential expression analysis
targets2 = targetsA2C(targets)
u = unique(targets2$Target)
f = factor(targets2$Target, levels = u)
design = model.matrix(~0 + f)
colnames(design) = u

corfit = intraspotCorrelation(MAq, design)
fit = lmscFit(MAq, design, correlation = corfit$consensus.correlation)

cont.matrixMars = makeContrasts(mars-earth, levels = design)
cont.matrixSpace = makeContrasts(space-earth, levels = design)

fit2Mars = contrasts.fit(fit, cont.matrixMars)
fit2Mars = eBayes(fit2Mars)
tableMars = data.frame(topTable(fit2Mars, coef=1, n=Inf, adjust="BH"))
fit2Space = contrasts.fit(fit, cont.matrixSpace)
fit2Space = eBayes(fit2Space)
tableSpace = data.frame(topTable(fit2Space, coef=1, n=Inf, adjust="BH"))


write.table(table,file=opt$output,sep="\t", quote = F)
cat("All done! Differential expression information saved to:",opt$output,"\n")


# design = modelMatrix(targets, ref = "earth")
# fit = lmFit(MA, design)
# fit = eBayes(fit)



# write.table(
#   MAq,
#   file = paste(outFH,".txt", sep = ""),
#   sep = "\t",
#   quote = F,
#   row.names = FALSE
# )

# Potential QC
suppressWarnings(
  arrayQualityMetrics(
    expressionset = MAq,
    outdir = "MA_test_report",
    force = T
  )
)
