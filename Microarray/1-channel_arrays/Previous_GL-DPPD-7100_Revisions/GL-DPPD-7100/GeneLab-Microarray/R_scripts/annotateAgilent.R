#!/usr/bin/env Rscript

# install.packages("optparse")
# source("http://bioconductor.org/biocLite.R")
# biocLite("limma")

suppressPackageStartupMessages(library("optparse"))

# Read options
option_list = list(
  make_option(
    c("-i", "--input"), 
    type = "character", 
    help = "Name of (or path to) the input normalized data (tab delimited .txt file or binary Rdata object)"
  ),
  make_option(
    c("-a", "--annotation"), 
    type = "character", 
    default = "search",
    help = "Set to 'search' to automatically look for a GPL file in the directory specified with the '--gplDir' option, otherwise provide a path to the file containing custom array annotation information, containing column names but no additional header"
  ),
  make_option(
    "--gplDir",
    type = "character",
    help = "Specify the directory to look for a GPL file if the '-a/--annotation' option is set to 'search'. (Default: the same directory as the input file)"
  ),
  make_option(
    c("-p", "--probeIDs"), 
    type = "character", 
    default = "ID",
    help = "Column name in the provided annotation file for the column containing the probe IDs as they are listed in the normalized data (default: 'ID'). If set to 'detect', the script will try to identify the column in the annotation file with the most similarity to the probe IDs in the input data"
  ),
  make_option(
    c("-g", "--geneIDs"), 
    type = "character", 
    default = "GB_ACC",
    help = "Column name in the provided annotation file for the column containing the desired gene IDs to use to link to the probe IDs (default: 'GB_ACC'). Note, if there spaces in the column name, replace them with periods (ie RefSeq Transcript ID --> RefSeq.Transcript.ID. If set to 'detect', the script will try to identify the columns with RefSeq gene IDs or potential alternative gene IDs"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = "annotExpValues",
    help = "Name of (or path to) file to write results to, WITHOUT file extension (default: annotExpValues)"
  ),
  make_option(
    c("-t", "--outType"),
    type = "character",
    default = "both",
    help = "Format of output data: R (Rdata object), txt (tab delimited file with identifiers and sample names [single ch only]), both [single ch only]"
  ),
  make_option(
    "--QCDir",
    type = "character",
    default = "./QC_reporting/",
    help = "Path to directory to store the annotation output. Will be created if it does not exist yet, but it is recommended to use the same directory as was used for QC with the normalization step (default = './QC_reporting/')"
  ),
  make_option(
    c("-d", "--dupProbes"),
    type = "character",
    default = "max",
    help = "Method for handling multiple probes [max (default, probe with the highest mean expression), average (mean of all probes for a gene), topvar (highest variance with nsFilter function)"
  ),
  make_option(
    "--GLDS", 
    type = "character", 
    help = "Full accession number for QC outputs (ie 'GLDS-21' for GLDS-21)"
  )
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

addSlash = function(string) {
  # Adds a trailing forward slash to the end of a string (ex path to a driectory) if it is not present
  if (substr(x = string,
             start = nchar(string),
             stop = nchar(string)) != "/") {
    string = paste(string, "/", sep = "")
  }
  return(string)
}

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call. = FALSE)
}else { inFH = opt$input } 

# Load packages
suppressPackageStartupMessages(require("methods"))

outFH = opt$outFile

qcDir = opt$QCDir

if (!is.null(opt$annotation)){
  if (opt$annotation == "search") {
    if (is.null(opt$gplDir)) {
      cat("Looking in the 'input' directory for a GPL file...\n")
      inDir = gsub("(/)[^/]*$", "\\1", inFH) # Strip the filename away from the directory path to the input file
      if (!grepl("/", inDir)) {
        # If the script was run from the directory containing the raw data and no path was provided, the script looks in that directory
        inDir = "./"
      }
    } else {
      cat("Looking in",opt$gplDir,"for a GPL file...\n")
      inDir = opt$gplDir
      inDir = addSlash(inDir)
    }
    files = dir(inDir)
    files = files[grepl("GPL[[:digit:]]*", files)]
    if (length(files) == 1) {
      annotFH = paste(inDir,files[1],sep="")
      cat("\t",annotFH,"identified as an annotation file\n")
    } else if (length(files) == 0) {
      stop(paste("No GPL file found in ",inDir,"\n",sep=""),
           call. = F)
    } else if (length(files) > 1) {
      stop(paste("Multiple GPL files found in ",inDir,"\nConsider providing the specific file with the '--annotation' option\n",sep=""),
           call. = F)
    }
  } else {
    annotFH = opt$annotation # annotFH = "GPL10094-20413.txt"
  }
} else {
  stop("Annotation file not specified and unable to identify a single GPL annotation file in the input directory\n", call. = F)
}

# Read in annotation file
tryCatch({
  annot = read.delim(annotFH, header = T, stringsAsFactors = F)
  if (grepl("^#",annot[1,1]) | grepl("^!",annot[1,1])) {
    cat("Extra header lines detected. Attempting to read this file in again without these lines but if this fails, consider manually trimming lines above the column names in the GPL file\n")
    skipCnt = 0
    for (i in 1:nrow(annot)) {
      if (grepl("^#",annot[i,1]) | grepl("^!",annot[i,1])) {
        # Looks for the last line that begins with an # or ! and assumes the next line below it contains the column names
        skipCnt = i
      } else {
        break()
      }
    }
    if ((skipCnt > 0) & (skipCnt < nrow(annot))){
      annot = read.delim(annotFH, header = T, stringsAsFactors = F, skip = (skipCnt + 1)) # skipCnt + 1 because of the first line being read as a header the first time
    }
  }
}, error = function(e) {
  stop(paste("Unable to read in ",annotFH, sep = ""),
       call. = F)
})

# Read in normalized data
# inFH = "/Users/dmattox/Documents/genelab/rot1/GLDS-41/microarray/expValues.txt"
tryCatch({
  if (grepl(".txt$", x = inFH) == TRUE) {
    eset = read.delim(
      inFH,
      header = T,
      sep = "\t",
      stringsAsFactors = F
    )
    rownames(eset) = eset[,1]
    eset[,1] = NULL
    cat(inFH,"read in successfully\n")
  } else if (grepl(".rda$", x = inFH) == TRUE) {
    load(inFH)
    cat(inFH,"read in successfully\n")
  } else {
    stop("File extension not recognized\n", call. = F)
  }
}, error = function(e) {
  stop("Input file was not recognized", call. = F)
})

if (exists("MA")) {
  # Determine if working with single channel or two channel microarray data
  dCh = TRUE
} else {
  dCh = FALSE
}

if (is.null(opt$GLDS)) {
  # Include GLDS accession number in outputs if provided
  glAn = ''
  cat("Warning: No GLDS accession number provided\n")
  if (grepl("GLDS-[0-9]+", inFH)) {
    glAn = regmatches(inFH, regexpr("GLDS-[0-9]+", inFH)) # Attempt to extract the GLDS accession number from the input path
    cat("Try to guess GLDS accession number... ", glAn,"\n")
  } else{
    glAn = FALSE
  }
} else{
  glAn = opt$GLDS
}

annotFun = function(oldProbes, newProbes, newIDs){
  # Function to switch from one probe label (oldProbes) to another (newIDs), provided
  ## a linking variable (newProbes) with common labels to the labels in the oldProbes list organized in the same order as the newIDs list
  # Returns a translation from the the oldProbes to the newIDs with "" for unmapped probes, ordered by the oldProbes
  newOrder = match(oldProbes,newProbes)
  out = newIDs[newOrder]
  return(out)
}

# Identify the column in the annotation file containing the most probe IDs
if (opt$probeIDs == 'detect') {
  cat("Looking for a column in the annotation file that matches the probe IDs found in the raw data\n")
  if (dCh == T) {
    probes = MA$genes[,1]
  } else {
    probes = row.names(eset)
  }
  sim = matrix(rep(0,ncol(annot)))
  for (i in 1:ncol(annot)) {
    sim[i] = sum(probes %in% annot[,i])
  }
  topInd = which.max(sim)
  newProbeName = colnames(annot)[topInd]
  cat("Probe ID column in annotation file identified as:",newProbeName,"\n\n")
  if ((sim[topInd] / length(probes)) < 0.5 ) {
    warning(paste("Only ", round(sim[topInd] / length(probes)*100,digits=1), "% of the probe IDs from the raw data are contained in the selected column from the annotation file\n",sep=""))
  }
} else {
  newProbeName = opt$probeIDs # newProbeName = "ID"
}

# Identify the column in the annotation file with columns with RefSeq gene IDs or potential alternative gene IDs
if (opt$geneIDs == 'detect') {
  RefSeq = matrix(rep(0, ncol(annot))) # Looking for a column of RefSeq IDs
  for (i in 1:ncol(annot)) {
    RefSeq[i] = sum(grepl("^[XN][MRPG]_", annot[, i])) / nrow(annot)
  }
  if (any(RefSeq > 0.4)) {
    topInd = which(RefSeq == max(RefSeq)) # Return all columns of annot containing the most RefSeq IDs
    if (length(topInd) > 1) {
      # Tie breaker
      for (i in 1:length(topInd)) {
        j = topInd[i]
        RefSeq[j] = sum(grepl("^NM_", annot[, j])) / nrow(annot) # Recalculate percentages based on a more strict RefSeq regex match
      }
      topInd = which.max(RefSeq) # In the case of a second tie, which.max just takes the lower index
    }
    newIDName = colnames(annot)[topInd]
    cat(
      "Identified",
      newIDName,
      "as a RefSeq containing column in the annotation file with",
      round(RefSeq[topInd] * 100, 2),
      "percent of the rows matching RefSeq formatting\n"
    )
  } else {
    BSU = matrix(rep(0, ncol(annot))) # Looking for a column of Bacillus subtilis gene IDs
    for (i in 1:ncol(annot)) {
      BSU[i] = sum(grepl("^BSU[[:digit:]]{5}", annot[, i])) / nrow(annot)
    }
    if (any(BSU > 0.4)) {
      topInd = which.max(BSU)
      newIDName = colnames(annot)[topInd]
      cat(
        "Identified",
        newIDName,
        "as a BSU ID containing column in the annotation file with",
        round(BSU[topInd] * 100, 2),
        "percent of the rows matching RefSeq formatting\n"
      )
    } else {
      ATMG = matrix(rep(0, ncol(annot))) # Looking for a column of Arabidopsis thaliana gene IDs
      for (i in 1:ncol(annot)) {
        ATMG[i] = sum(grepl("^ATMG[[:digit:]]{5}", annot[, i])) / nrow(annot)
      }
      if (any(ATMG > 0.4)) {
        topInd = which.max(ATMG)
        newIDName = colnames(annot)[topInd]
        cat(
          "Identified",
          newIDName,
          "as an ATMG ID containing column in the annotation file with",
          round(ATMG[topInd] * 100, 2),
          "percent of the rows matching RefSeq formatting\n"
        )
      } else {
        yeast = matrix(rep(0, ncol(annot))) # Looking for a column of yeast gene IDs
        for (i in 1:ncol(annot)) {
          yeast[i] = sum(grepl("^Y[[:upper:]]{2}(\\d)+", annot[, i])) / nrow(annot)
        }
        if (any(yeast > 0.4)) {
          topInd = which.max(yeast)
          newIDName = colnames(annot)[topInd]
          cat(
            "Identified",
            newIDName,
            "as an yeast ID containing column in the annotation file with",
            round(yeast[topInd] * 100, 2),
            "percent of the rows matching RefSeq formatting\n"
          )
        } else{
          stop("No columns with sufficent numbers of recognized gene IDs detected")
        }
      }
    }
  }
} else {
  newIDName = opt$geneIDs # newIDName = "GB_ACC"
}

## Two channel annotation
if (dCh == TRUE) {
  cat("Annotating two channel microarray data from",inFH,"\n\n")
  cat("Matching gene IDs to probe names...\n")
  geneIDs = annotFun(
    oldProbes = MA$genes[,1], 
    newProbes = annot[,newProbeName],
    newIDs = annot[,newIDName]
  )
  
  if (any(grepl("///", geneIDs))) {
    tryCatch({
      for (i in 1:length(geneIDs)) {
        geneIDs[i] = strsplit(geneIDs[i], split = "///")[[1]][1]
      }
      geneIDs = gsub(" ", "", geneIDs)
    }, error = function(e) {
      warning("Error in stripping multiple annotations for the same probe. Gene IDs may need further reformatting\n")
    })
  }
  
  if (any(grepl("(\\|)", geneIDs))) {
    tryCatch({
      for (i in 1:length(geneIDs)) {
        tmp = strsplit(geneIDs[i], split = "\\|")[[1]]
        geneIDs[i] = tmp[grep("^([[:upper:]]){2}_", tmp)][1]
      }
      geneIDs = gsub(" ", "", geneIDs)
    }, error = function(e) {
      warning("Error in stripping multiple annotations for the same probe. Gene IDs may need further reformatting\n")
    })
  }
  
  cat("Removing unlabeled probe names...\n")
  if ( sum(grepl("^([[:upper:]]){2}_",geneIDs))/length(geneIDs) > 0.5 ) {
    # If using primarily RefSeq IDs in the newIDs column, only use rows following RefSeq formatting (ie "^NM_" )
    noIDTag = ( !grepl("^([[:upper:]]){2}_",geneIDs) | geneIDs == "")
  } else if ( sum(grepl("PA([[:digit:]]){4}",geneIDs))/length(geneIDs) > 0.5 ) {
    # If using primarily p. aeruginosa IDs in the newIDs column, only use rows following p. aeruginosa IDs formatting (ie "^PAxxxx" )
    noIDTag = ( !grepl("PA([[:digit:]]){4}",geneIDs) | geneIDs == "")
  } else {
    noIDTag = (geneIDs == "" | grepl("corner", geneIDs, ignore.case = T))
  }
  noIDCnt = sum(noIDTag) # Number of unmapped probes
  MA = MA[!noIDTag, ] # Remove umapped probes
  geneIDs = geneIDs[!noIDTag]
  
  
  
  # Filter out probes with missing values
  naTag = (apply(is.na(MA$M), 1, any) | apply(is.na(MA$A), 1, any))
  MA = MA[!naTag,]
  geneIDs = geneIDs[!naTag]
  naCnt = sum(naTag)
  
  if (any(opt$dupProbes %in% c("average", "max"))) {
    cat("Filtering out multiple probes per gene ID...\n")
    
    if (opt$dupProbes == "average") {
      # Collapse multiple probes per gene ID by averaging expression values across all samples
      rmRowTag = rep(TRUE, length(MA$genes[,1])) # Tag rows to drop (set single or averaged probes to FALSE below)
      for (i in 1:length(rmRowTag)) {
        if (sum(geneIDs == geneIDs[i]) > 1) {
          inds = grep(geneIDs[i], geneIDs) # List of indices at which a probe for a given gene ID occur
          MA$M[inds[1], ] = apply(X = MA$M[inds, ],
                                  FUN = mean,
                                  MARGIN = 2) # Changes the values of the first occurence of a probe to the sample-specific average of the values from all the probes for that gene ID
          MA$A[inds[1], ] = apply(X = MA$A[inds, ],
                                  FUN = mean,
                                  MARGIN = 2) # Changes the values of the first occurence of a probe to the sample-specific average of the values from all the probes for that gene ID
          rmRowTag[inds[1]] = FALSE
        } else
          rmRowTag[i] = FALSE
      }
      nDups = sum(rmRowTag)
      MA = MA[!rmRowTag, ]
      geneIDs = geneIDs[!rmRowTag]
      
      cat("\tUnmapped probes removed:", noIDCnt, "\n")
      cat("\tProbes with missing values removed:", naCnt, "\n")
      cat("\tDuplicated probes removed:", nDups, "\n\n")
      cat("Annotated probes remaining:", length(MA$genes[,1]), "\n\n")
      if (length(geneIDs) > length(unique(geneIDs))) {
        cat("\n\tWarning: non-unique probe to ID mappings remain \n")
      }
      
      MA$genes[,1] = geneIDs # Set probe IDs to gene IDs
      
    } else if (opt$dupProbes == "max") {
      # Collapse multiple probes per gene ID by selecting a representative with the highest mean intensity across all samples
      rmRowTag = rep(TRUE, length(MA$genes[,1])) # Tag rows to drop (set single or highest expressing probes to FALSE below)
      for (i in 1:length(rmRowTag)) {
        if (sum(geneIDs == geneIDs[i]) > 1) {
          inds = grep(geneIDs[i], geneIDs)
          top = 0
          keep = 0
          for (j in 1:length(inds)) {
            curr = mean(as.numeric(MA$A[inds[j], ]))
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
      MA = MA[!rmRowTag, ]
      geneIDs = geneIDs[!rmRowTag]
      
      cat("\tUnmapped probes removed:", noIDCnt, "\n")
      cat("\tProbes with missing values removed:", naCnt, "\n")
      cat("\tDuplicated probes removed:", nDups, "\n\n")
      cat("Annotated probes remaining:", length(MA$genes[,1]), "\n\n")
      if (length(geneIDs) > length(unique(geneIDs))) {
        cat("\n\tWarning: non-unique probe to ID mappings remain \n")
      }
      
      MA$genes[,1] = geneIDs # Set probe IDs to gene IDs
    }
  } else{
    stop("Method for dealing with probes mapped to the same gene IDs not recognized\n",
         call. = F)
  }
  
  # Output annotation report to the specified QC directory
  summDir = paste(qcDir, "summary_report/", sep = "")
  if (!file.exists(summDir)){ # Create a summary report directory within qcDir if it does not exist yet
    dir.create(summDir)
  }
  
  AR = c(
    paste("Unmapped probes removed:", noIDCnt),
    paste("Probes with missing values removed:", naCnt),
    paste("Duplicated probes removed:", nDups),
    paste("Annotated probes remaining:", length(MA$genes[,1]))
  )
  if (glAn != FALSE) {
    write.table(
      AR,
      file = paste(summDir, glAn, "_annotReport.txt", sep = ""),
      quote = F,
      col.names = F,
      row.names = F
    )
    cat("Annotation report generated!",paste(summDir, glAn, "_annotReport.txt", sep = ""),"\n")
  } else {
    write.table(
      AR,
      file = paste(summDir,"annotReport.txt", sep = ""),
      quote = F,
      col.names = F,
      row.names = F
    )
    cat("Annotation report generated!",paste(summDir,"annotReport.txt", sep = ""),"\n")
  }
  
  # Save filtered expression values
  outFH = opt$output
  if (opt$outType == "both") {
    save(MA, file = paste(outFH, ".rda", sep = ""))
    cat("Success! Annotated data saved to", outFH, "as a .RData file \nTwo channel microarray data cannot be saved in a text file yet, sorry!\n")
  } else if (opt$outType == "R") {
    save(MA, file = paste(outFH, ".rda", sep = ""))
    cat("Success! Annotated data saved to", outFH, "as a .RData file \n")
  } else if (opt$outType == "txt") {
    cat("Two channel microarray data cannot be saved in a text file yet, sorry!\nSaving as an Rdata object instead:",paste(outFH, ".rda", sep = ""),"\n")
  } else{
    print_help(opt_parser)
    stop("Help, I don't know how to save this data!\n", call. = F)
  }
  
  
}


## Single channel annotation
if (dCh == FALSE) {
  cat("Annotating single channel microarray data from",inFH,"\n\n")
  cat("Matching gene IDs to probe names...\n")
  geneIDs = annotFun(
    oldProbes = row.names(eset), 
    newProbes = annot[,newProbeName],
    newIDs = annot[,newIDName]
  )
  
  if (any(grepl("///", geneIDs))) {
    tryCatch({
      for (i in 1:length(geneIDs)) {
        geneIDs[i] = strsplit(geneIDs[i], split = "///")[[1]][1]
      }
      geneIDs = gsub(" ", "", geneIDs)
    }, error = function(e) {
      warning("Error in stripping multiple annotations for the same probe. Gene IDs may need further reformatting\n")
    })
  }
  
  if (any(grepl("(\\|)", geneIDs))) {
    tryCatch({
      for (i in 1:length(geneIDs)) {
        tmp = strsplit(geneIDs[i], split = "\\|")[[1]]
        geneIDs[i] = tmp[grep("^([[:upper:]]){2}_", tmp)][1]
      }
      geneIDs = gsub(" ", "", geneIDs)
    }, error = function(e) {
      warning("Error in stripping multiple annotations for the same probe. Gene IDs may need further reformatting\n")
    })
  }
  
  cat("Removing unlabeled probe names...\n")
  if ( sum(grepl("^([[:upper:]]){2}_",geneIDs))/length(geneIDs) > 0.5 ) {
    # If using primarily RefSeq IDs in the newIDs column, only use rows following RefSeq formatting (ie "^NM_" )
    noIDTag = ( !grepl("^([[:upper:]]){2}_",geneIDs) | geneIDs == "")
  } else if ( sum(grepl("PA([[:digit:]]){4}",geneIDs))/length(geneIDs) > 0.5 ) {
    # If using primarily p. aeruginosa IDs in the newIDs column, only use rows following p. aeruginosa IDs formatting (ie "^PAxxxx" )
    noIDTag = ( !grepl("PA([[:digit:]]){4}",geneIDs) | geneIDs == "")
  } else {
    noIDTag = (geneIDs == "" | grepl("corner", geneIDs, ignore.case = T))
  }
  noIDCnt = sum(noIDTag) # Number of unmapped probes
  eset = eset[!noIDTag, ] # Remove umapped probes
  geneIDs = geneIDs[!noIDTag]
  
  
  
  # Filter out probes with missing values
  naTag = (apply(is.na(eset), 1, any))
  eset = eset[!naTag,]
  geneIDs = geneIDs[!naTag]
  naCnt = sum(naTag)
  
  if (any(opt$dupProbes %in% c("average", "max"))) {
    cat("Filtering out multiple probes per gene ID...\n")
    
    if (opt$dupProbes == "average") {
      # Collapse multiple probes per gene ID by averaging expression values across all samples
      rmRowTag = rep(TRUE, nrow(eset)) # Tag rows to drop (set single or averaged probes to FALSE below)
      for (i in 1:length(rmRowTag)) {
        if (sum(geneIDs == geneIDs[i]) > 1) {
          inds = grep(geneIDs[i], geneIDs) # List of indices at which a probe for a given gene ID occur
          eset[inds[1], ] = apply(X = eset[inds, ],
                                  FUN = mean,
                                  MARGIN = 2) # Changes the values of the first occurence of a probe to the sample-specific average of the values from all the probes for that gene ID
          rmRowTag[inds[1]] = FALSE
        } else
          rmRowTag[i] = FALSE
      }
      nDups = sum(rmRowTag)
      eset = eset[!rmRowTag, ]
      geneIDs = geneIDs[!rmRowTag]
      
      cat("\tUnmapped probes removed:", noIDCnt, "\n")
      cat("\tProbes with missing values removed:", naCnt, "\n")
      cat("\tDuplicated probes removed:", nDups, "\n\n")
      cat("Annotated probes remaining:", nrow(eset), "\n\n")
      if (length(geneIDs) > length(unique(geneIDs))) {
        cat("\n\tWarning: non-unique probe to ID mappings remain \n")
      }
      
      row.names(eset) = geneIDs # Set probe IDs to gene IDs
      
    } else if (opt$dupProbes == "max") {
      # Collapse multiple probes per gene ID by selecting a representative with the highest mean intensity across all samples
      rmRowTag = rep(TRUE, nrow(eset)) # Tag rows to drop (set single or highest expressing probes to FALSE below)
      for (i in 1:length(rmRowTag)) {
        if (sum(geneIDs == geneIDs[i]) > 1) {
          inds = grep(geneIDs[i], geneIDs)
          top = 0
          keep = 0
          for (j in 1:length(inds)) {
            curr = mean(as.numeric(eset[inds[j], ]))
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
      eset = eset[!rmRowTag, ]
      geneIDs = geneIDs[!rmRowTag]
      
      cat("\tUnmapped probes removed:", noIDCnt, "\n")
      cat("\tProbes with missing values removed:", naCnt, "\n")
      cat("\tDuplicated probes removed:", nDups, "\n\n")
      cat("Annotated probes remaining:", nrow(eset), "\n\n")
      if (nrow(eset) > length(unique(row.names(eset)))) {
        cat("\n\tWarning: non-unique probe to ID mappings remain \n")
      }
      
      row.names(eset) = geneIDs # Set probe IDs to gene IDs
    }
  } else{
    stop("Method for dealing with probes mapped to the same gene IDs not recognized\n",
         call. = F)
  }
  
  # Output annotation report to the specified QC directory
  summDir = paste(qcDir, "summary_report/", sep = "")
  if (!file.exists(summDir)){ # Create a summary report directory within qcDir if it does not exist yet
    dir.create(summDir)
  }
  
  AR = c(
    paste("Unmapped probes removed:", noIDCnt),
    paste("Probes with missing values removed:", naCnt),
    paste("Duplicated probes removed:", nDups),
    paste("Annotated probes remaining:", nrow(eset))
  )
  if (glAn != FALSE) {
    write.table(
      AR,
      file = paste(summDir, glAn, "_annotReport.txt", sep = ""),
      quote = F,
      col.names = F,
      row.names = F
    )
    cat("Annotation report generated!",paste(summDir, glAn, "_annotReport.txt", sep = ""),"\n")
  } else {
    write.table(
      AR,
      file = paste(summDir,"annotReport.txt", sep = ""),
      quote = F,
      col.names = F,
      row.names = F
    )
    cat("Annotation report generated!",paste(summDir,"annotReport.txt", sep = ""),"\n")
  }
  
  # Save filtered expression values
  outFH = opt$output
  colnames(eset) = gsub("\\.","-",colnames(eset)) # Keep the sample names standardized (if data read in as a text file, hyphens are swapped for periods)
  if (opt$outType == "both") {
    save(eset, file = paste(outFH, ".rda", sep = ""))
    write.table(
      data.frame("ID" = row.names(eset),eset), # provides the rownames as a labelled column in the saved output
      row.names = F,
      file = paste(outFH, ".txt", sep = ""),
      sep = "\t",
      quote = F
    )
    cat("Success! Annotated data saved to", outFH, "as both a .txt and a .RData file \n")
  } else if (opt$outType == "R") {
    save(eset, file = paste(outFH, ".rda", sep = ""))
    cat("Success! Annotated data saved to", outFH, "as a .RData file \n")
  } else if (opt$outType == "txt") {
    write.table(
      data.frame("ID" = row.names(eset),eset), # provides the rownames as a labelled column in the saved output
      row.names = F,
      file = paste(outFH, ".txt", sep = ""),
      sep = "\t",
      quote = F
    )
    cat("Success! Annotated data saved to", outFH, "as a .txt file \n")
  } else{
    print_help(opt_parser)
    stop("Help, I don't know how to save this data!\n", call. = F)
  }
}






