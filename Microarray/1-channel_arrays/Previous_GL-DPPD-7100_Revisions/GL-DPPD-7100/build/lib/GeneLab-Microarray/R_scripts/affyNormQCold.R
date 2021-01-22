#!/usr/bin/env Rscript

# install.packages("optparse")
# source("http://bioconductor.org/biocLite.R")
# biocLite("affy")
# biocLite("affyPLM")
# biocLite("oligo")

suppressPackageStartupMessages(library("optparse"))

relDir = getwd()

# Read options
option_list=list(
  make_option(c("-i","--input"),type="character",help="Path to directory containing input .CEL files"),
  make_option(c("-n","--normalization"),type="character",default="rma",help="Normalization method [rma (default, full rma), quantile (no background correction), background (no quant. normalization), log2 (no quant. norm. or background correction)"),
  make_option(c("-o","--outFile"),type="character",default="expValues",help="Name of the output file [without extension!] (default: expValues)"),
  make_option("--outDir",type="character",default="./",help="Path to an output directory, including a terminal forward slash (default: directory this script is called from)"),
  make_option(c("-t","--outType"),type="character",default="both",help="Format of output data: R (Rdata object), txt (tab delimited file with identifiers and sample names), both(default)"),
  make_option("--outputData",type="logical",default=TRUE,help="Output data at all (default TRUE)"),
  make_option(c("-a","--arrayInfoOnly"),type="logical",default=FALSE,help="Detect-affy-array-only mode. If true, script will exit after outputting the arrayInfo file. (Default: FALSE)"),
  make_option("--QCoutput",type="logical",default=TRUE,help="Output QC_reporting directory of QC plots (default = TRUE)"),
  make_option("--QCDir",type="character",default="./QC_reporting/",help="Path to directory for storing QC output, including a terminal forward slash. Will be created if it does not exist yet (default = './QC_reporting/')"),
  make_option("--NUSEplot",type="logical",default=FALSE,help="Include a NUSE plot in the QC output, adds significantly to runtime (default = FALSE)"),
  make_option("--GLDS",type="character",help="GLDS accession number for plot outputs (ie '21' for GLDS-21)")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

addSlash = function(string){
  # Adds a trailing forward slash to the end of a string (ex path to a driectory) if it is not present
  if(substr(x = string,start = nchar(string), stop = nchar(string)) != "/"){
    string = paste(string,"/",sep="")
  }
  return(string)
}

stopQuiet = function(...) {
  blankMsg = sprintf("\r%s\r", paste(rep(" ", getOption("width")-1L), collapse=" "))
  stop(simpleError(blankMsg))
}

norm = opt$normalization
QCout = opt$QCoutput
NUSEplot = opt$NUSEplot

if (is.null(opt$input)){ # Check for provided input directory
  print_help(opt_parser)
  stop("No path to input directory provided. Please look over the available options", call. = F)
}else{
  inPath = addSlash(opt$input)
  setwd(inPath) # Change the working directory to the directory containing the raw files
}

if (is.null(opt$GLDS)){ # Include GLDS accession number in outputs if provided
  glAn = ''
  cat("Warning: No GLDS accession number provided\n")
  if(grepl("GLDS-[0-9]+",inPath)){
    glAn = regmatches(inPath, regexpr("GLDS-[0-9]+",inPath)) # Attempt to extract the GLDS accession number from the input path
  }
}else{
  glAn = paste('GLDS-',opt$GLDS,sep='')
}

detach_package = function(pkg, character.only = FALSE){
  if(!character.only)
  {
    pkg <- deparse(substitute(pkg))
  }
  search_item <- paste("package", pkg, sep = ":")
  while(search_item %in% search())
  {
    detach(search_item, unload = TRUE, character.only = TRUE)
  }
}

# Load initial libraries
suppressPackageStartupMessages(require(affy))

# setwd("~/Documents/genelab/rot1/GLDS-4/microarray/")
celFiles <- list.celfiles(full.names=TRUE)
sampNames = gsub("_microarray_.*","",celFiles)
sampNames = gsub(".CEL","",sampNames)
sampNames = gsub(".*/","",sampNames)
sampNames = gsub("GLDS-\\d*_","",sampNames)# Extract sample names form the list of .CEL files

tryCatch({suppressWarnings(expr = {raw = ReadAffy()})}, error=function(e){
  stop("No .CEL files detected in the current directory", call. = F)
  })

arrInfo = c("Affymetrix",as.character(raw@cdfName))

if (grepl("-st-",raw@cdfName,ignore.case = T)){
  detach_package(affy)
  rm(raw)
  suppressPackageStartupMessages(require(oligo))
  raw = read.celfiles(celFiles)
  st = T
}else{
  suppressPackageStartupMessages(require(affyPLM))
  st = F
}

setwd(relDir) # Return the working directory to direcotry script was called from to enable use of relative paths
# Create QC output directory
qcDir = addSlash(opt$QCDir)
if(!file.exists(qcDir)) dir.create(qcDir)

# Output array information to a separate file
write.table(arrInfo,file = paste(qcDir,glAn,"_arrayInfo.txt",sep=""),quote = F,
            col.names = F, row.names = F)
cat("Array type detected.\n")

# Exit script if arrayInfoOnly mode is True
if(opt$arrayInfoOnly == TRUE){
  cat("Detect-array-type-only mode on, exiting.\n")
  quit(save = "no", status = 0, runLast = FALSE)
}

## Raw QC
if(QCout == T){
  cat("Performing intial QC\n")
  # Prepare plotting options
  toMatch = c(8,183,31,45,51,100,101,118,128,139,147,183,254,421,467,477,
              483,493,498,503,508,535,552,575,635,655)
  color = grDevices::colors()[rep(toMatch,10)] # Create a library of colors for plotting
  
  #Images
  cat("\tGenerating raw images")
  if(st == T){
    for(i in 1:length(celFiles)){
      png(paste(qcDir,glAn,'_',sampNames[i],'_image.png',sep=''),width=800, height = 800)
      image(raw, which = i)
      garbage <- dev.off()
      cat(".")
    }
  }else{
    nblines=length(celFiles)%/%4 + as.numeric((length(celFiles)%%4)!=0)
    png(paste(qcDir,glAn,'_images.png',sep=''),width=800,height = 200*nblines)
    par(mfrow=c(nblines,4))
    image(raw)
    garbage <- dev.off()
  }
  cat("\n")

  #MA plot
  cat("\tGenerating raw data MA plots...\n")
  nblines=length(celFiles)%/%3 + as.numeric((length(celFiles)%%3)!=0)
  png(paste(qcDir,glAn,'_rawPlotMA.png',sep=''),width=800, height = 300*nblines )
  par(mfrow=c(nblines,3))
  if(st == T){
    MAplot(raw)
  }else{
    MAplot(raw,type="pm")
  }
  garbage <- dev.off()
  
  # Intensity distributions of the pm probes from each microarray on the same graph
  cat("\tGenerating initial distribution plots")
  mypms = pm(raw)
  png(paste(qcDir,glAn,'_rawDensityDistributions.png',sep=''),width=800,height=800 )
  ylims = c(0,.8)
  xlims = c(0,16)
  for(i in 1:ncol(mypms)){
    cat(".")
    if(i == 1){
      plot(density(log2(mypms[,i])),ylim = ylims,xlim=xlims,xlab='log2(Raw Intensities)',main=paste(glAn,' Raw intensity distributions',sep=''),col=color[i])
      par(new=T)
    }else{
      plot(density(log2(mypms[,i])),ylim = ylims,xlim=xlims,axes=F,xlab='',ylab='',main='',col=color[i])
      par(new=T)
    }
  }
  legend(13,0.8,col=color[1:length(celFiles)],legend=sampNames
         ,pch=15,bty = "n",cex = 0.9,pt.cex = 0.8,y.intersp = 0.8)
  garbage <- dev.off()
  cat("\n")
  
  # Boxplots
  png(paste(qcDir,glAn,'_rawBoxplot.png',sep=''),width=800,height = 400)
  par(mar=c(7,5,1,1))
  if(st == T){
    invisible(capture.output( # Silence oligo::rma printing statement
      boxplot(oligo::rma(raw, background=FALSE, normalize=FALSE, subset=NULL, target="core"), las=2,
              names = sampNames, main=paste(glAn," Raw intensities",sep=""),col=color[1:length(celFiles)])
    ))
  }else{
    boxplot(raw,las=2,outline=FALSE,col=color[1:length(celFiles)],main = paste(glAn," Raw intensities",sep=""),names=sampNames)
  }
  mtext(text="log2 Intensity", side=2, line=2.5, las=0)
  garbage <- dev.off()
  
  # PCA
  cat("\tPerforming PCA of raw data...\n")
  rawPCA = prcomp(mypms)
  png(paste(qcDir,glAn,'_rawPCA.png',sep=''),width=800,height = 800)
  plot(rawPCA$rotation[,1],rawPCA$rotation[,2],col=color[1:length(celFiles)],pch=16,
       xlab = paste("PC1, ",round(summary(rawPCA)$importance["Proportion of Variance",1]*100,digits = 1),"% of variance",sep=""),
       ylab = paste("PC2, ",round(summary(rawPCA)$importance["Proportion of Variance",2]*100,digits = 1),"% of variance",sep=""),
       main=paste(glAn," PCA of raw data",sep="")
  )
  text(rawPCA$rotation[,1],rawPCA$rotation[,2],labels = sampNames, cex = 1,pos = 3)
  garbage <- dev.off()
  
  #NUSE plot
  if(NUSEplot == T){
    cat("\tFitting probe-level model and generating RLE/NUSE plots...\n")
    if(st == T){
      Pset = fitProbeLevelModel(raw)
      # RLE plot
      png(paste(qcDir,glAn,'_RLE.png',sep=''),width=800,height = 600)
      par(mar=c(7,5,1,1))
      RLE(Pset, col = color[1:length(sampNames)],
          names = sampNames, las=2, main=paste(glAn," Relative Log Expression (RLE) plot",sep=""))
      abline(h=0,lty=1,col="red")
      garbage <- dev.off()
      # NUSE plot
      png(paste(qcDir,glAn,'_NUSE.png',sep=''),width=800,height = 600)
      par(mar=c(7,5,1,1))
      NUSE(Pset, col = color[1:length(sampNames)], las=2)
      title(main=paste(glAn," NUSE plot of microarray experiments",sep=""))
      abline(h=1.1,lty=1,col="red")
      garbage <- dev.off()
    }else{
      Pset=fitPLM(raw)
      # RLE plot
      png(paste(qcDir,glAn,'_RLE.png',sep=''),width=800,height = 600)
      par(mar=c(7,5,1,1))
      RLE(Pset, col = color[1:length(sampNames)],
          names = sampNames, las=2, main="Relative Log Expression (RLE) plot")
      abline(h=0,lty=1,col="red")
      garbage <- dev.off()
      # NUSE plot
      png(paste(qcDir,glAn,'_NUSE.png',sep=''),width=800,height = 600)
      par(mar=c(7,5,1,1))
      NUSE(Pset,col = color[1:length(sampNames)], las=2)
      title(main=paste(glAn," NUSE plot of microarray experiments",sep=""))
      abline(h=1.1,lty=1,col="red")
      garbage <- dev.off()
    }
  }
}

outFH = opt$outFile
if(opt$outputData == TRUE){
  
  ## Normalize
  cat("\nNormalizing with selected normalization technique...\n")
  if(norm=='rma'){
    eset = rma(raw)
  }else if(norm=='quantile'){
    eset = rma(raw, background = F, normalize = T)
  }else if(norm=='background'){
    eset = rma(raw, background = T, normalize = F)
  }else if(norm=='log2'){
    eset = rma(raw, background = F, normalize = F)
  }else{
    stop("Normalization did not occur, please examine script inputs and default values",call. = F)
  }
  
  outDir = addSlash(opt$outDir)
  if(opt$outType == "both"){
    save(eset,file=paste(outDir,outFH,".rda",sep=""))
    write.table(exprs(eset),file=paste(outDir,outFH,".txt",sep=""),sep="\t",quote = F)
  }else if(opt$outType == "R"){
    save(eset,file=paste(outDir,outFH,".rda",sep=""))
  }else if(opt$outType == "txt"){
    write.table(exprs(eset),file=paste(outDir,outFH,".txt",sep=""),sep="\t",quote = F)
  }else{
    print_help(opt_parser)
    stop("Help, I don't know how to save this data!",call. = F)
  }
}

if(QCout == T){
  cat("Post normalization QC steps...\n")
  # Post-normalization QC
  png(paste(qcDir,glAn,'_normDensityDistributions.png',sep=''),width=800,height=800 )
  ylims = c(0,.8)
  xlims = c(0,16)
  normVals = exprs(eset)
  for(i in 1:ncol(normVals)){
    if(i == 1){
      plot(density(normVals[,i]),ylim = ylims,xlim=xlims,xlab='Normalized expression values[log2]',main=paste(glAn,' Normalized expression distributions',sep=''),col=color[i])
      par(new=T)
    }else{
      plot(density(normVals[,i]),ylim = ylims,xlim=xlims,axes=F,xlab='',ylab='',main='',col=color[i])
      par(new=T)
    }
  }
  legend(13,0.8,col=color[1:length(celFiles)],legend=sampNames
         ,pch=15,bty = "n",cex = 0.9,pt.cex = 0.8,y.intersp = 0.8)
  garbage <- dev.off()
  
  # Boxplots
  png(paste(qcDir,glAn,'_normBoxplot.png',sep=''),width=800,height = 400)
  par(mar=c(7,5,1,1))
  if(st == T){
    boxplot(normVals,las=2,outline=FALSE,col=color[1:length(celFiles)],main=paste(glAn," Normalized intensities",sep=""),transfo='identity',names=sampNames)
    mtext(text="log2 Intensity", side=2, line=2.5, las=0) 
  }else{
    boxplot(normVals,las=2,outline=FALSE,col=color[1:length(celFiles)],main=paste(glAn," Normalized intensities",sep=""),names=sampNames)
    mtext(text="log2 Intensity", side=2, line=2.5, las=0)
  }
  garbage <- dev.off()
  
  #MA plot
  cat("\tGenerating MA plots from the normalized data...\n")
  nblines=length(celFiles)%/%3 + as.numeric((length(celFiles)%%3)!=0)
  png(paste(qcDir,glAn,'_normPlotMA.png',sep=''),width=800, height = 300*nblines )
  par(mfrow=c(nblines,3))
  MAplot(eset)
  garbage <- dev.off()
  
  # PCA
  cat("\tPerforming PCA of normalized data...\n")
  normPCA = prcomp(normVals)
  png(paste(qcDir,glAn,'_normPCA.png',sep=''),width=800,height = 800)
  plot(normPCA$rotation[,1],normPCA$rotation[,2],col=color[1:length(celFiles)],pch=16,
       xlab = paste("PC1, ",round(summary(normPCA)$importance["Proportion of Variance",1]*100,digits = 1),"% of variance",sep=""),
       ylab = paste("PC2, ",round(summary(normPCA)$importance["Proportion of Variance",2]*100,digits = 1),"% of variance",sep=""),
       main=paste(glAn," PCA of normalized data",sep="")
  )
  text(normPCA$rotation[,1],normPCA$rotation[,2],labels = sampNames, cex = 1,pos = 3)
  garbage <- dev.off()
}