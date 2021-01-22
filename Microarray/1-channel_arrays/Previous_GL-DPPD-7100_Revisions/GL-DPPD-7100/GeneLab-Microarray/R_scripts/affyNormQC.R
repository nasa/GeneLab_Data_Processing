#!/usr/bin/env Rscript

### Single channel Affymetrix microarray normalization and quality control


# install.packages("optparse")
# source("http://bioconductor.org/biocLite.R")
# biocLite("affy")
# biocLite("affyPLM")
# biocLite("oligo")
# biocLite("arrayQualityMetrics") # Known dependencies for arrayQualityMetrics listed below
#   biocLite("hexbin")
#   biocLite("jsonlite")
#   biocLite("openssl")
#   biocLite("stringi")
#   biocLite("reshape2")
#   biocLite("Cairo")

suppressPackageStartupMessages(library("optparse")) # Load optparse package to read in arguments

relDir = getwd() # Store the directory from which this script was called for relative paths

# Read options
option_list = list(
  make_option(
    c("-i", "--input"), 
    type = "character", 
    help = "Path to directory containing input .CEL files"
  ),
  make_option(
    c("-n", "--normalization"),
    type = "character",
    default = "rma",
    help = "Normalization method [rma (default, full rma), quantile (no background correction), background (no quant. normalization), log2 (no quant. norm. or background correction)"
  ),
  make_option(
    c("-o", "--outFile"),
    type = "character",
    default = "expValues",
    help = "Name of the output file [without extension!] (default: expValues)"
  ),
  make_option(
    "--outDir",
    type = "character",
    default = "./",
    help = "Path to an output directory, including a terminal forward slash (default: directory this script is called from)"
  ),
  make_option(
    c("-t", "--outType"),
    type = "character",
    default = "both",
    help = "Format of output data: R (Rdata object), txt (tab delimited file with identifiers and sample names), both (default)"
  ),
  make_option(
    "--outputData",
    type = "logical",
    default = TRUE,
    help = "Output data at all (default TRUE)"
  ),
  make_option(
    c("-a", "--arrayInfoOnly"),
    type = "logical",
    default = FALSE,
    help = "Detect-affy-array-only mode. If true, script will exit after outputting the arrayInfo file. (Default: FALSE)"
  ),
  make_option(
    "--QCpackage",
    type = "character",
    default = "aqm",
    help = "Package used to generate QC plots: aqm (generate html report with arrayQualityMetrics package,default), R (use standard R plotting options to generate QC plots as .png files. Figures may not maintain proper formatting for datasets with many samples)"
  ),
  make_option(
    # Note: this option will not work on the DP server when reading in the arrays with the "oligo" package (ie ST type arrays) due to the BLAS multi-threading issue
    "--NUSEplot",
    type = "logical",
    default = FALSE,
    help = "Include a NUSE and RLE plot in the QC output, adds significantly to runtime (default = FALSE). Note: this argument is only considered with 'R' selected for --QCPackage"
  ),
  make_option(
    "--QCoutput",
    type = "logical",
    default = TRUE,
    help = "Output QC_reporting directory of QC plots (default = TRUE)"
  ),
  make_option(
    "--QCDir",
    type = "character",
    default = "./QC_reporting/",
    help = "Path to directory for storing QC output. Will be created if it does not exist yet (default = './QC_reporting/')"
  ),
  make_option(
    "--GLDS", 
    type = "character", 
    help = "Full accession number for QC outputs (ie 'GLDS-21' for GLDS-21). If it's not provided, the script will try to detect it from the path provided to the '-i/--input' option"
  )
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser) # Parse the arguments into a list

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

# Saving options from the list as global variables
norm = opt$normalization 
QCout = opt$QCoutput
QCpack = opt$QCpackage
NUSEplot = opt$NUSEplot

if (is.null(opt$input)) {
  # Check for provided input directory
  print_help(opt_parser)
  stop("No path to input directory provided. Please look over the available options",
       call. = F)
} else{
  inPath = addSlash(opt$input)
  setwd(inPath) # Change the working directory to the directory containing the raw files
}

if (is.null(opt$GLDS)) {
  # Include GLDS accession number in outputs if provided
  cat("Warning: No GLDS accession number provided\n")
  if (grepl("GLDS-[0-9]+", inPath)) { # Check if a GLDS is in the input directory 
    glAn = regmatches(inPath, regexpr("GLDS-[0-9]+", inPath)) # Attempt to extract the GLDS accession number from the input path
    cat("Trying to guess GLDS accession number... ", glAn,"\n")
  } else{
    glAn = FALSE # Set to false if not provided and not in the input directory
  }
} else{
  glAn = opt$GLDS
}

# Load affy package to read in .CEL files
suppressPackageStartupMessages(require(affy)) # Load affy-specific package for reading in arrays

celFiles = list.celfiles(full.names = TRUE) # Pull all .CEL files from the input directory into a list

if (length(celFiles) > 0){
  # If the list of CEL files is populated, list the contents
  cat("Detected .CEL files:\n")
  for (i in 1:length(celFiles)) {
    cat("\t",celFiles[i],"\n")
  }
  cat("\n")
} else {
  stop("No .CEL files detected in the current directory",
       call. = F)
}

# Strip the sample names from the file names
sampNames = gsub("_microarray", "", celFiles)
sampNames = gsub(".CEL", "", sampNames)
sampNames = gsub(".*/", "", sampNames)
sampNames = gsub("GLDS-\\d*_", "", sampNames)
sampNames = gsub("raw", "", sampNames)
sampNames = gsub("_", "", sampNames)

useOligo = tryCatch({ 
  # Tries to read in the .CEL files using the affy package, but if it fails (in the case of some of the newer human microarrays), it switches over to the oligo package. If the .CEL files load with the affy package but are better handled with the oligo package, it will also detatch the affy package and load the oligo package
  suppressWarnings(expr = {
    raw = ReadAffy(filenames = celFiles, # Read the raw files with the affy package
                   sampleNames = sampNames)
  })
  arrInfo = c("Affymetrix",
              as.character(raw@cdfName)) # Store array infromation
  if (grepl("-st-", raw@cdfName, ignore.case = T)) {
    # If the array is an ST type (best handled with the oligo package), remove the affy package and load the oligo package
    detach_package(affy)
    rm(raw)
    suppressPackageStartupMessages(require(oligo))
    raw = read.celfiles(filenames = celFiles, # Read in the raw files with the oligo package
                        sampleNames = sampNames)
    st = T # Set a marker variable to indicate use of oligo package
    
  } else {
    # If the affy package is sufficent, load the probe-level modeling package and set marker variable for using affy
    suppressPackageStartupMessages(require(affyPLM))
    st = F
  }
  cat("") # 
}, error = function(e) {
  # If reading in the raw files with the affy package fails (occurs for newer human arrays), returning TRUE sets useOligo to TRUE
  cat(
    "Could not read in provided .CEL files with affy package. Attempting to read them with oligo package...\n"
  )
  return(TRUE)
})

if (!is.null(useOligo)){
  if (useOligo == TRUE) {
    # If the useOligo variable exists and is TRUE, detach the affy package and load the oligo package
    tryCatch({
      detach_package(affy)
      suppressPackageStartupMessages(require(oligo))
      raw = read.celfiles(filenames = celFiles, # Read in the raw files with the oligo package
                          sampleNames = sampNames)
      st = T # Set the oligo marker varible to TRUE
      ver = raw@annotation # Pull package information from the oligo object
      if (grepl("^pd.", ver)) {
        # Convert from the pd.* annotation package to the standard array version name
        ver = gsub("^pd.", "", ver)
        ver = gsub("\\.", "-", ver)
        ver = gsub("(\\d)(-)(\\d)", "\\1_\\3", ver) # Replace hyphens with underscores
      }
      arrInfo = c("Affymetrix", as.character(ver)) # Store the array information
    }, error = function(e) {
      stop("Unable to read in .CEL files with oligo package or affy package", call. = F)
    })
  }
}


setwd(relDir) # Return the working directory to directory script was called from to enable use of relative paths
# Create QC output directory
qcDir = addSlash(opt$QCDir)
if (!file.exists(qcDir)){ # Create QC directory if it does not exist yet
  dir.create(qcDir)
}

# Output array information to a separate file
summDir = paste(qcDir, "summary_report/", sep = "")
if (!file.exists(summDir)){ # Create a summary report directory within qcDir if it does not exist yet
  dir.create(summDir)
}

if (glAn != FALSE) {
  # Save the array information formatted with the GLDS
  write.table(
    arrInfo,
    file = paste(summDir, glAn, "_arrayInfo.txt", sep = ""),
    quote = F,
    col.names = F,
    row.names = F
  )
} else {
  # Save the array information formatted without the GLDS
  write.table(
    arrInfo,
    file = paste(summDir, "arrayInfo.txt", sep = ""),
    quote = F,
    col.names = F,
    row.names = F
  )
}

cat("Array type detected.\n")

# Exit script if arrayInfoOnly mode is True
if (opt$arrayInfoOnly == TRUE) {
  cat("Detect-array-type-only mode on, exiting.\n")
  quit(save = "no",
       status = 0,
       runLast = FALSE)
}

## Raw QC
if(QCout == T) {
  cat("Performing intial QC\n")
  
  if (QCpack == "aqm"){
    # Load in the arrayQualityMetrics QC package
    suppressPackageStartupMessages(require(arrayQualityMetrics))
    
    suppressWarnings(
      arrayQualityMetrics( # Generate the raw QC report in the generated "raw_report" directory within the QC directory
        expressionset = raw,
        outdir = paste(qcDir, "raw_report", sep = ""),
        force = T,
        do.logtransform = T # Log-transform because the values haven't been normalized yet
      )
    )
  } else if ( QCpack == "R") { # Generate QC plots using base R plotting options
    # Create a directory to hold figures from the raw QC step
    rawDir = paste(qcDir, "raw_report/", sep = "")
    if (!file.exists(rawDir)){ # Create a raw report directory within qcDir if it does not exist yet
      dir.create(rawDir)
    }
    
    # Prepare plotting options
    toMatch = c(8,183,31,45,51,100,101,118,128,139,147,183,254,421,467,477,
                483,493,498,503,508,535,552,575,635,655) # A list of 26 indices for contrasting colors to use for plots
    color = grDevices::colors()[rep(toMatch,10)] # Create a library of colors for plotting, repeating 10 times because who knows how many sample there will be
    
    # Generate pseudo-images
    cat("\tGenerating raw images")
    if (st == T) {
      for (i in 1:length(celFiles)) {
        png(
          paste(rawDir, glAn, '_', sampNames[i], '_image.png', sep = ''),
          width = 800,
          height = 800
        )
        image(raw, which = i)
        garbage <- dev.off()
        cat(".")
      }
    } else{
      nblines = length(celFiles) %/% 4 + as.numeric((length(celFiles) %% 4) != 0)
      png(
        paste(rawDir, glAn, '_images.png', sep = ''),
        width = 800,
        height = 200 * nblines
      )
      par(mfrow = c(nblines, 4))
      image(raw)
      garbage <- dev.off()
    }
    cat("\n")
    
    # Intensity distributions of the pm probes from each microarray on the same graph
    cat("\tGenerating initial distribution plots")
    mypms = pm(raw)
    png(
      paste(rawDir, glAn, '_rawDensityDistributions.png', sep = ''),
      width = 800,
      height = 800
    )
    ylims = c(0, .8)
    xlims = c(0, 16)
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
      13,
      0.8,
      col = color[1:length(celFiles)],
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
    
    # Generate intesity boxplots
    png(paste(rawDir, glAn, '_rawBoxplot.png', sep = ''),
        width = 800,
        height = 400)
    par(mar = c(7, 5, 1, 1))
    if (st == T) {
      invisible(capture.output(
        # Silence oligo::rma printing statement
        boxplot(
          oligo::rma(
            raw,
            background = FALSE,
            normalize = FALSE,
            subset = NULL,
            target = "core"
          ),
          las = 2,
          names = sampNames,
          main = paste(glAn, " Raw intensities", sep = ""),
          col = color[1:length(celFiles)]
        )
      ))
    } else{
      boxplot(
        raw,
        las = 2,
        outline = FALSE,
        col = color[1:length(celFiles)],
        main = paste(glAn, " Raw intensities", sep = ""),
        names = sampNames
      )
    }
    mtext(
      text = "log2 Intensity",
      side = 2,
      line = 2.5,
      las = 0
    )
    garbage <- dev.off()
    
    # Generate PCA plot
    cat("\tPerforming PCA of raw data...\n")
    rawPCA = prcomp(mypms)
    png(paste(rawDir, glAn, '_rawPCA.png', sep = ''),
        width = 800,
        height = 800)
    plot(
      rawPCA$rotation[, 1],
      rawPCA$rotation[, 2],
      col = color[1:length(celFiles)],
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
    
    # If not on the DP sever (or the BLAS issue is resolved), use this option for more thorough QC analysis
    if (NUSEplot == T) {
      cat("\tFitting probe-level model and generating RLE/NUSE plots...\n")
      if (st == T) {
        # oligo package specific code for RLE and NUSE plots
        Pset = fitProbeLevelModel(raw) # Fit a probe-level model to the raw data
        # RLE plot
        png(paste(rawDir, glAn, '_RLE.png', sep = ''),
            width = 800,
            height = 600)
        par(mar = c(7, 5, 1, 1))
        RLE(
          Pset,
          col = color[1:length(sampNames)],
          names = sampNames,
          las = 2,
          main = paste(glAn, " Relative Log Expression (RLE) plot", sep = "")
        )
        abline(h = 0, lty = 1, col = "red")
        garbage <- dev.off()
        # NUSE plot
        png(paste(rawDir, glAn, '_NUSE.png', sep = ''),
            width = 800,
            height = 600)
        par(mar = c(7, 5, 1, 1))
        NUSE(Pset, col = color[1:length(sampNames)], las = 2)
        title(main = paste(glAn, " NUSE plot of microarray experiments", sep =
                             ""))
        abline(h = 1.1, lty = 1, col = "red")
        garbage <- dev.off()
        
      } else{
        # affy package specific code for RLE and NUSE plots
        Pset = fitPLM(raw) # Fit a probe-level model to the raw data
        # RLE plot
        png(paste(rawDir, glAn, '_RLE.png', sep = ''),
            width = 800,
            height = 600)
        par(mar = c(7, 5, 1, 1))
        RLE(
          Pset,
          col = color[1:length(sampNames)],
          names = sampNames,
          las = 2,
          main = "Relative Log Expression (RLE) plot"
        )
        abline(h = 0, lty = 1, col = "red")
        garbage <- dev.off()
        # NUSE plot
        png(paste(rawDir, glAn, '_NUSE.png', sep = ''),
            width = 800,
            height = 600)
        par(mar = c(7, 5, 1, 1))
        NUSE(Pset, col = color[1:length(sampNames)], las = 2)
        title(main = paste(glAn, " NUSE plot of microarray experiments", sep =""))
        abline(h = 1.1, lty = 1, col = "red")
        garbage <- dev.off()
      }
    }
  }else {
    warning("\n--QCpackage option was not recognized and quality control is not being performed.\n", call. = F)
  }
}



outFH = opt$outFile
if (opt$outputData == TRUE) {
  ## Normalize
  cat("\nNormalizing with selected normalization technique...\n")
  if (norm == 'rma') {
    expset = rma(raw) # Perform standard RMA normalization
  } else if (norm == 'quantile') {
    expset = rma(raw, background = F, normalize = T)
  } else if (norm == 'background') {
    expset = rma(raw, background = T, normalize = F)
  } else if (norm == 'log2') {
    expset = rma(raw, background = F, normalize = F)
  } else{
    stop(
      "Normalization did not occur, please examine script inputs and default values\n",
      call. = F
    )
  }
  
  outDir = addSlash(opt$outDir)
  if (grepl("/", outDir)) {
    if (!file.exists(outDir)) {
      # Create the output directory if it does not exist yet
      dir.create(outDir)
    }
  }
  eset = exprs(expset) # Extract the intensity values from the normalized ExpressionSet object
  if (opt$outType == "both") {
    save(eset, file = paste(outDir, outFH, ".rda", sep = "")) # Save as a binary RData file that loads as a variable named eset
    write.table(
      data.frame("ID" = row.names(eset),eset), # Creates a temporary dataframe, providing the row names as a labeled column in the saved output
      row.names = F,
      file = paste(outDir, outFH, ".txt", sep = ""),
      sep = "\t",
      quote = F
    )
    cat("Success! Normalized data saved to", paste(outDir, outFH, sep=""), "as both a .txt and a .RData file\n\n")
  } else if (opt$outType == "R") {
    save(eset, file = paste(outDir, outFH, ".rda", sep = ""))
    cat("Success! Normalized data saved to", paste(outDir, outFH, sep=""), "as a .RData file\n\n")
  } else if (opt$outType == "txt") {
    write.table(
      data.frame("ID" = row.names(eset),eset), # provides the rownames as a labeled column in the saved output
      row.names = F,
      file = paste(outDir, outFH, ".txt", sep = ""),
      sep = "\t",
      quote = F
    )
    cat("Success! Normalized data saved to", paste(outDir, outFH, sep=""), "as a .txt file\n\n")
  } else{
    print_help(opt_parser)
    stop("Help, I don't know how to save this data!", call. = F)
  }
  # Post-normalization QC
  if(QCout == T) {
    cat("Post normalization QC steps...\n")
    if (QCpack == "aqm") {
      suppressWarnings(
        arrayQualityMetrics( # Generates an HTML report with interactive quality control graphs and figures within the QCdir in the "normalized_report" directory
          expressionset = expset,
          outdir = paste(qcDir, "normalized_report", sep = ""),
          force = T
        )
      )
    } else if (QCpack == "R") {
      # Generate QC figures with base R plotting options following normalization
      normDir = paste(qcDir, "normalized_report/", sep = "")
      if (!file.exists(normDir)){ # Create a directory to hold figures from the normalized QC step
        dir.create(normDir)
      }
      
      # Generating overlaying density distributions
      png(
        paste(normDir, glAn, '_normDensityDistributions.png', sep = ''),
        width = 800,
        height = 800
      )
      ylims = c(0, .8)
      xlims = c(0, 16)
      normVals = exprs(expset)
      for (i in 1:ncol(normVals)) {
        if (i == 1) {
          plot(
            density(normVals[, i]),
            ylim = ylims,
            xlim = xlims,
            xlab = 'Normalized expression values[log2]',
            main = paste(glAn, ' Normalized expression distributions', sep = ''),
            col = color[i]
          )
          par(new = T)
        } else{
          plot(
            density(normVals[, i]),
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
        col = color[1:length(celFiles)],
        legend = sampNames
        ,
        pch = 15,
        bty = "n",
        cex = 0.9,
        pt.cex = 0.8,
        y.intersp = 0.8
      )
      garbage <- dev.off()
      
      # Generating boxplots
      png(
        paste(normDir, glAn, '_normBoxplot.png', sep = ''),
        width = 800,
        height = 400
      )
      par(mar = c(7, 5, 1, 1))
      if (st == T) {
        boxplot(
          normVals,
          las = 2,
          outline = FALSE,
          col = color[1:length(celFiles)],
          main = paste(glAn, " Normalized intensities", sep = ""),
          transfo = 'identity',
          names = sampNames
        )
        mtext(
          text = "log2 Intensity",
          side = 2,
          line = 2.5,
          las = 0
        )
      } else{
        boxplot(
          normVals,
          las = 2,
          outline = FALSE,
          col = color[1:length(celFiles)],
          main = paste(glAn, " Normalized intensities", sep = ""),
          names = sampNames
        )
        mtext(
          text = "log2 Intensity",
          side = 2,
          line = 2.5,
          las = 0
        )
      }
      garbage <- dev.off()
      
      # Generating pseudo-MA plots to a median reference array
      cat("\tGenerating MA plots from the normalized data...\n")
      nblines = length(celFiles) %/% 3 + as.numeric((length(celFiles) %% 3) !=
                                                      0)
      png(
        paste(normDir, glAn, '_normPlotMA.png', sep = ''),
        width = 800,
        height = 300 * nblines
      )
      par(mfrow = c(nblines, 3))
      MAplot(expset)
      garbage <- dev.off()
      
      # Generating a PCA plot
      cat("\tPerforming PCA of normalized data...\n")
      normPCA = prcomp(normVals)
      png(paste(normDir, glAn, '_normPCA.png', sep = ''),
          width = 800,
          height = 800)
      plot(
        normPCA$rotation[, 1],
        normPCA$rotation[, 2],
        col = color[1:length(celFiles)],
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
    }
  }
}
