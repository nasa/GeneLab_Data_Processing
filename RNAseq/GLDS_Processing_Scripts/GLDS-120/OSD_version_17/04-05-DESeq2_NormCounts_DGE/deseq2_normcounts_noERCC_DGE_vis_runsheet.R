#### Set up your environment #####

## Import libraries (tximport, DESeq2, tidyverse, stringr)

library(tximport)
library(DESeq2)
library(tidyverse)
library(stringr)


## Define organism
organism <- "ARABIDOPSIS"

### Define the location of the input data and where the ouput data will be printed to ###
runsheet_path="/GLDS-120/Metadata/GLDS-120_bulkRNASeq_v1_runsheet.csv"
work_dir="/GLDS-120/processing_scripts/04-05-DESeq2_NormCounts_DGE"
counts_dir="/GLDS-120/03-RSEM_Counts"
norm_output="/GLDS-120/04-DESeq2_NormCounts"
DGE_output="/GLDS-120/05-DESeq2_DGE"

### Pull in the GeneLab annotation table (GL-DPPD-7110_annotations.csv) file ###

org_table_link <- "https://raw.githubusercontent.com/nasa/GeneLab_Data_Processing/master/GeneLab_Reference_Annotations/Pipeline_GL-DPPD-7110_Versions/GL-DPPD-7110/GL-DPPD-7110_annotations.csv"
org_table <- read.table(org_table_link, sep = ",", header = TRUE)

## Define link to the GeneLab annotation table for organism of interest
annotations_link <- org_table[org_table$name == organism, "genelab_annots_link"]


### Set your working directory to the directory where you will execute your DESeq2 script from ###

setwd(file.path(work_dir))


####### Pull all factors for each sample in the study from the runsheet #####

compare_csv_from_runsheet <- function(runsheet_path) {
    df = read.csv(runsheet_path)
    # get only Factor Value columns
    factors = as.data.frame(df[,grep("Factor.Value", colnames(df), ignore.case=TRUE)])
    colnames(factors) = paste("factor",1:dim(factors)[2], sep= "_")
    result = data.frame(sample_id = df[,c("Sample.Name")], factors)	
    return(result)
}


### Load metadata from runsheet csv file ###

compare_csv <- compare_csv_from_runsheet(runsheet_path)


## Create data frame containing all samples and respective factors

study <- as.data.frame(compare_csv[,2:dim(compare_csv)[2]])
colnames(study) <- colnames(compare_csv)[2:dim(compare_csv)[2]]
rownames(study) <- compare_csv[,1]


#### If running from a metadata csv file and not ISA use the command below
#study <- read.csv(Sys.glob(file.path(metadata_dir,"*metadata_mouse.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE)

## Set your working directory to the directory containing the DGE script

setwd(file.path(work_dir))


##### Format groups and indicate the group that each sample belongs to #####

if (dim(study) >= 2){
  group<-apply(study,1,paste,collapse = " & ") # concatenate multiple factors into one condition per sample
} else{
  group<-study[,1]
}
group_names <- paste0("(",group,")",sep = "") # human readable group names
group <- sub("^BLOCKER_", "",  make.names(paste0("BLOCKER_", group))) # group naming compatible with R models, this maintains the default behaviour of make.names with the exception that 'X' is never prepended to group names
names(group) <- group_names
rm(group_names)


##### Format contrasts table, defining pairwise comparisons for all groups #####

contrast.names <- combn(levels(factor(names(group))),2) # generate matrix of pairwise group combinations for comparison
contrasts <- apply(contrast.names, MARGIN=2, function(col) sub("^BLOCKER_", "",  make.names(paste0("BLOCKER_", stringr::str_sub(col, 2, -2))))) # limited make.names call for each group (also removes leading parentheses)
contrast.names <- c(paste(contrast.names[1,],contrast.names[2,],sep = "v"),paste(contrast.names[2,],contrast.names[1,],sep = "v")) # format combinations for output table files names
contrasts <- cbind(contrasts,contrasts[c(2,1),])
colnames(contrasts) <- contrast.names
rm(contrast.names) 


##### Import RSEM raw (gene) count data #####

files <- list.files(file.path(counts_dir),pattern = ".genes.results", full.names = TRUE)

## Reorder the *genes.results files to match the ordering of the samples in the metadata
files <- files[sapply(rownames(study), function(x)grep(paste0(x,".genes.results$"), files, value=FALSE))]

names(files) <- rownames(study)
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)

## Add 1 to genes with lengths of zero - needed to make DESeqDataSet object 
txi.rsem$length[txi.rsem$length == 0] <- 1


## Create data frame defining which group each sample belongs to
sampleTable <- data.frame(condition=factor(group))
rownames(sampleTable) <- colnames(txi.rsem$counts)


##### Make DESeqDataSet object #####

dds <- DESeqDataSetFromTximport(txi.rsem, sampleTable, ~condition)
summary(dds)

##### Filter out genes with counts of less than 10 in all samples #####

keepGenes <- rowSums(counts(dds)) > 10
dds <- dds[keepGenes,]
summary(dds)
dim(dds)

##### Perform DESeq analysis #####

dds_1 <- DESeq(dds)


##### Generate F statistic p-value (similar to ANOVA p-value) using DESeq2 likelihood ratio test (LRT) design #####

dds_1_lrt <- DESeq(dds_1, test = "LRT", reduced = ~ 1)
res_1_lrt <- results(dds_1_lrt)

### Create a data frame containing normalized counts ###

normCounts <- as.data.frame(counts(dds_1, normalized=TRUE))

## Add 1 to all counts to avoid issues with downstream calculations

normCounts <- normCounts +1

### Start the DGE output table with the normalized counts for all samples ###

## reduced output table 1 will be used to generate human-readable DGE table
reduced_output_table_1 <- normCounts

## output table 1 will be used to generate computer-readable DGE table, which is used to create GeneLab visualization plots
output_table_1 <- normCounts

## Iterate through Wald Tests to generate pairwise comparisons of all groups
for (i in 1:dim(contrasts)[2]){
  res_1 <- results(dds_1, contrast=c("condition",contrasts[1,i],contrasts[2,i]))
  res_1 <- as.data.frame(res_1@listData)[,c(2,4,5,6)]
  colnames(res_1) <-c(paste0("Log2fc_",colnames(contrasts)[i]),paste0("Stat_",colnames(contrasts)[i]),paste0("P.value_",colnames(contrasts)[i]),paste0("Adj.p.value_",colnames(contrasts)[i]))
  output_table_1<-cbind(output_table_1,res_1)
  reduced_output_table_1 <- cbind(reduced_output_table_1,res_1)
  rm(res_1)
}


## Generate and add all sample mean column to the normalized counts table
output_table_1$All.mean <- rowMeans(normCounts, na.rm = TRUE, dims = 1)
reduced_output_table_1$All.mean <- rowMeans(normCounts, na.rm = TRUE, dims = 1)

## Generate and add all sample stdev column to the normalized counts table
output_table_1$All.stdev <- rowSds(as.matrix(normCounts), na.rm = TRUE, dims = 1)
reduced_output_table_1$All.stdev <- rowSds(as.matrix(normCounts), na.rm = TRUE, dims = 1)

## Add F statistic p-value (similar to ANOVA p-value) column to the normalized counts table
output_table_1$LRT.p.value <- res_1_lrt@listData$padj
reduced_output_table_1$LRT.p.value <- res_1_lrt@listData$padj

## Generate and add group mean and stdev columns to the normalized counts table
tcounts <- as.data.frame(t(normCounts))
tcounts$group <- group
group_means <- as.data.frame(t(aggregate(. ~ group,data = tcounts,mean)))
group_means <- group_means[-c(1),]
colnames(group_means) <- paste0("Group.Mean_",levels(factor(names(group))))
group_stdev <- as.data.frame(t(aggregate(. ~ group,data = tcounts,sd)))
group_stdev <- group_stdev[-c(1),]
colnames(group_stdev) <- paste0("Group.Stdev_",levels(factor(names(group))))

output_table_1 <- cbind(output_table_1,group_means)
reduced_output_table_1 <- cbind(reduced_output_table_1,group_means)

output_table_1 <- cbind(output_table_1,group_stdev)
reduced_output_table_1 <- cbind(reduced_output_table_1,group_stdev)

rm(group_stdev,group_means,tcounts)


### Add columns needed to generate GeneLab visulaization plots to the normalized counts table

## Add column to indicate the sign (positive/negative) of log2fc for each pairwise comparison
updown_table <- sign(output_table_1[,grep("Log2fc_",colnames(output_table_1))])
colnames(updown_table) <- gsub("Log2fc","Updown",grep("Log2fc_",colnames(output_table_1),value = TRUE))
output_table_1 <- cbind(output_table_1,updown_table)
rm(updown_table)

## Add column to indicate contrast significance with p <= 0.1
sig.1_table <- output_table_1[,grep("P.value_",colnames(output_table_1))]<=.1
colnames(sig.1_table) <- gsub("P.value","Sig.1",grep("P.value_",colnames(output_table_1),value = TRUE))
output_table_1 <- cbind(output_table_1,sig.1_table)
rm(sig.1_table)

## Add column to indicate contrast significance with p <= 0.05
sig.05_table <- output_table_1[,grep("P.value_",colnames(output_table_1))]<=.05
colnames(sig.05_table) <- gsub("P.value","Sig.05",grep("P.value_",colnames(output_table_1),value = TRUE))
output_table_1 <- cbind(output_table_1,sig.05_table)
rm(sig.05_table)

## Add columns for the volcano plot with p-value and adjusted p-value
log_pval_table <- log2(output_table_1[,grep("P.value_",colnames(output_table_1))])
colnames(log_pval_table) <- paste0("Log2_",colnames(log_pval_table))
output_table_1 <- cbind(output_table_1,log_pval_table)
rm(log_pval_table)
log_adj_pval_table <- log2(output_table_1[,grep("Adj.p.value_",colnames(output_table_1))])
colnames(log_adj_pval_table) <- paste0("Log2_",colnames(log_adj_pval_table))
output_table_1 <- cbind(output_table_1,log_adj_pval_table)
rm(log_adj_pval_table)


## Prepare PCA table for GeneLab visualization plots ##

exp_raw <- log2(normCounts)
PCA_raw <- prcomp(t(exp_raw), scale = FALSE)


### Read in GeneLab annotation table for the organism of interest ###

annot <- read.table(annotations_link, sep = "\t", header = TRUE, quote = "", comment.char = "", row.names = 1)

### Combine annotations table and the DGE table

output_table_1 <- merge(annot, output_table_1, by='row.names', all.y=TRUE)
output_table_1 <- output_table_1 %>% 
  rename(
    TAIR = Row.names
  )


reduced_output_table_1 <- merge(annot, reduced_output_table_1, by='row.names', all.y=TRUE)
reduced_output_table_1 <- reduced_output_table_1 %>% 
  rename(
    TAIR = Row.names
  )


### Export unnormalized and normalized counts tables ###

normCounts_exp <- as.data.frame(counts(dds_1, normalized=TRUE))

write.csv(txi.rsem$counts,file.path(norm_output, "RSEM_Unnormalized_Counts.csv"))
write.csv(normCounts_exp,file.path(norm_output, "Normalized_Counts.csv"))


### Export sample grouping and contrasts tables ###

write.csv(sampleTable,file.path(DGE_output, "SampleTable.csv"))
write.csv(contrasts,file.path(DGE_output, "contrasts.csv"))


### Export human-readable DGE table ###

write.csv(reduced_output_table_1,file.path(DGE_output, "differential_expression.csv"), row.names = FALSE)

### Export computer-readable DGE and PCA tables used for GeneLab visualization ###

write.csv(output_table_1,file.path(DGE_output, "visualization_output_table.csv"), row.names = FALSE)
write.csv(PCA_raw$x,file.path(DGE_output, "visualization_PCA_table.csv"), row.names = TRUE)


## print session info ##
print(" ")
print("Session Info below: ")
sessionInfo()

