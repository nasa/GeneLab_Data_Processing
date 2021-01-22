#!/usr/bin/env Rscript

## A ~*~potentially~*~ generalizable script to extract the raw data from horrible processed-file formats
# Initially built for GLDS-28

# install.packages("optparse")

suppressPackageStartupMessages(library("optparse"))

# Read options
option_list = list(
  make_option(c("-i", "--input"), type = "character", help = "Path to file with buried raw array information [required]"),
  make_option(c("-o", "--output"), type = "character", help = "Path and filename for saving extracted data [required]")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

norm = opt$normalization

if (is.null(opt$input)) {
  # Check and set input file handle
  print_help(opt_parser)
  stop("No path to input file provided. Please look over the available options\n",
       call. = F)
} else
  inFH = opt$input

if (is.null(opt$output)) {
  # Check and set output file handle
  print_help(opt_parser)
  stop("No path to out file provided. Please look over the available options\n",
       call. = F)
} else
  outFH = opt$output

cat("\nReading in file:", inFH, "\n")
test = read.delim(inFH, stringsAsFactors = F, header = F) # Read in processed file

cat("Extracting raw values...\n")
cat("\tGreen foreground\n")
gMedianSignal = test[, grep("gMedianSignal", test)] # Idenitfy median foreground intensity columns
cat("\tRed foreground\n")
rMedianSignal = test[, grep("rMedianSignal", test)]
cat("\tGreen background\n")
gBGMedianSignal = test[, grep("gBGMedianSignal", test)] # Idenitfy median background intensities columns
cat("\tRed background\n")
rBGMedianSignal = test[, grep("rBGMedianSignal", test)]
cat("\tProbe IDs\n")
startInd = grep("gMedianSignal", gMedianSignal) + 1 # Indentify starting row index of raw values
FeatureNum = test[startInd:nrow(test), grep("FeatureNum", test)] # Add unique feature numbers for mapping to genes
gMedianSignal = gMedianSignal[startInd:length(gMedianSignal)] # Extract raw values
rMedianSignal = rMedianSignal[startInd:length(rMedianSignal)]
gBGMedianSignal = gBGMedianSignal[startInd:length(gBGMedianSignal)]
rBGMedianSignal = rBGMedianSignal[startInd:length(rBGMedianSignal)]
raw = cbind(FeatureNum, rMedianSignal, gMedianSignal, rBGMedianSignal, gBGMedianSignal) # Bind raw value vectors into dataframe

cat("Saving... ")
write.table(raw, file = outFH, sep = "\t", quote = F, row.names = FALSE)

cat("Success! Raw values saved to:", outFH, "\n\n")
