library(tximport)
library(DESeq2)
library(tidyverse)

## Install annotation package created with the makeOrgPackageFromNCBI function of the AnnotationForge Bioconductor package
package_dir <- "/path/to/annotation_DBI_packages"
install.packages(file.path(package_dir,"org.Bsubtilis.eg.db.tar.gz"), repos=NULL)

organism <- "BACSU"

work_dir="/path/to/GLDS-138/processing_scripts/04-05-DESeq2_NormCounts_DGE"
counts_dir="/path/to/GLDS-138/03-RSEM_Counts"
norm_output="/path/to/GLDS-138/04-DESeq2_NormCounts"
DGE_output="/path/to/GLDS-138/05-DESeq2_DGE"

setwd(file.path(work_dir))

study <- read.csv(Sys.glob(file.path(work_dir,"*metadata.csv")), header = TRUE, row.names = 1, stringsAsFactors = TRUE)

##### Group Formatting
if (dim(study) >= 2){
	group<-apply(study,1,paste,collapse = " & ") # concatenate multiple factors into one condition per sample
} else{
	group<-study[,1]
}
group_names <- paste0("(",group,")",sep = "") # human readable group names
group <- make.names(group) # group naming compatible with R models
names(group) <- group_names
rm(group_names)

##### Contrast Formatting
contrasts <- combn(levels(factor(group)),2) # generate matrix of pairwise group combinations for comparison
contrast.names <- combn(levels(factor(names(group))),2)
contrast.names <- c(paste(contrast.names[1,],contrast.names[2,],sep = "v"),paste(contrast.names[2,],contrast.names[1,],sep = "v")) # format combinations for output table files names
contrasts <- cbind(contrasts,contrasts[c(2,1),])
colnames(contrasts) <- contrast.names
rm(contrast.names) 

##### Import Data
files <- list.files(file.path(counts_dir),pattern = ".genes.results", full.names = TRUE)
names(files) <- rownames(study)
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)

# add 1 to genes with lengths of zero if necessary
txi.rsem$length[txi.rsem$length == 0] <- 1

# make DESeqDataSet object
sampleTable <- data.frame(condition=factor(group))
rownames(sampleTable) <- colnames(txi.rsem$counts)

dds <- DESeqDataSetFromTximport(txi.rsem, sampleTable, ~condition)

# filter out genes with counts of less than 10 in all conditions
keep <- rowSums(counts(dds)) > 10
dds <- dds[keep,]
summary(dds)


#### Perform DESeq analysis
dds_1 <- DESeq(dds)

# export unnormalized, normalized, and ERCC normalized counts
normCounts = as.data.frame(counts(dds_1, normalized=TRUE))
setwd(file.path(norm_output))
write.csv(txi.rsem$counts,file='Unnormalized_Counts.csv')
write.csv(normCounts,file='Normalized_Counts.csv')
write.csv(sampleTable,file='SampleTable.csv')

setwd(file.path(work_dir))

normCounts <- normCounts +1

dds_1_lrt <- DESeq(dds_1, test = "LRT", reduced = ~ 1)
res_1_lrt <- results(dds_1_lrt)

organism_table <- read.csv(file.path(work_dir,"organisms.csv"))

##### Generate annotated DGE tables

assign("has_internet_via_proxy", TRUE, environment(curl::has_internet))
library(STRINGdb) # for String database annotations
library(PANTHER.db) # for GOSLIM annotations

ann.dbi <- organism_table$annotations[organism_table$name == organism] # Organism specific gene annotation database
ann.dbi=as.character(ann.dbi)
if(!require(ann.dbi, character.only=TRUE)) {
	BiocManager::install(ann.dbi, ask = FALSE)
	library(ann.dbi, character.only=TRUE)
}

## for normalized counts
# start output table with normalized sample expression values
output_table_1 <- normCounts
reduced_output_table_1 <- normCounts

##### Iterate through Wald Tests
for (i in 1:dim(contrasts)[2]){
	res_1 <- results(dds_1, contrast=c("condition",contrasts[1,i],contrasts[2,i]))
	res_1 <- as.data.frame(res_1@listData)[,c(2,5,6)]
	colnames(res_1)<-c(paste0("Log2fc_",colnames(contrasts)[i]),paste0("P.value_",colnames(contrasts)[i]),paste0("Adj.p.value_",colnames(contrasts)[i]))
	output_table_1<-cbind(output_table_1,res_1)
	reduced_output_table_1 <- cbind(reduced_output_table_1,res_1)
	rm(res_1)
}

# Gene Annotation columns
keytype = "ALIAS"
annot <- data.frame(rownames(output_table_1), stringsAsFactors = FALSE)
colnames(annot)[1]<-keytype
if ("SYMBOL" %in% columns(eval(parse(text = ann.dbi),env=.GlobalEnv))){
	annot$SYMBOL<-mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = rownames(output_table_1),keytype = keytype, column = "SYMBOL", multiVals = "first")
}
if ("GENENAME" %in% columns(eval(parse(text = ann.dbi),env=.GlobalEnv))){
        annot$GENENAME<-mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = rownames(output_table_1),keytype = keytype, column = "GENENAME", multiVals = "first")
}
if ("ENSEMBL" %in% columns(eval(parse(text = ann.dbi),env=.GlobalEnv))){
        annot$ENSEMBL<-mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = rownames(output_table_1),keytype = keytype, column = "ENSEMBL", multiVals = "first")
}
if ("REFSEQ" %in% columns(eval(parse(text = ann.dbi),env=.GlobalEnv))){
        annot$REFSEQ<-mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = rownames(output_table_1),keytype = keytype, column = "REFSEQ", multiVals = "first")
}
if ("ENTREZID" %in% columns(eval(parse(text = ann.dbi),env=.GlobalEnv))){
        annot$ENTREZID<-mapIds(eval(parse(text = ann.dbi),env=.GlobalEnv),keys = rownames(output_table_1),keytype = keytype, column = "ENTREZID", multiVals = "first")
}

string_db <- STRINGdb$new( version="10", species=organism_table$taxon[organism_table$name == organism],score_threshold=0)
string_map <- string_db$map(annot,"SYMBOL",removeUnmappedRows = FALSE, takeFirst = TRUE)[,c(1,6)]
string_map <- string_map[!duplicated(string_map$SYMBOL),]

annot <- dplyr::left_join(annot,string_map, by = "SYMBOL")
pthOrganisms(PANTHER.db) <- organism
panther <- mapIds(PANTHER.db,keys = annot$ENTREZID,keytype = "ENTREZ",column = "GOSLIM_ID", multiVals = "list")
panther <- na.omit(panther)
annot$GOSLIM_IDS <- panther
rm(string_db,string_map,panther,keytype)

# add all sample mean column
output_table_1$All.mean <- rowMeans(normCounts, na.rm = TRUE, dims = 1)
reduced_output_table_1$All.mean <- rowMeans(normCounts, na.rm = TRUE, dims = 1)

# add all sample stdev column
output_table_1$stdev <- rowSds(as.matrix(normCounts), na.rm = TRUE, dims = 1)
reduced_output_table_1$All.stdev <- rowSds(as.matrix(normCounts), na.rm = TRUE, dims = 1)

# add F statistic p-value (similar to ANOVA p-value) column
output_table_1$LRT.p.value <- res_1_lrt@listData$padj
reduced_output_table_1$LRT.p.value <- res_1_lrt@listData$padj

# add group mean and stdev columns
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

# add updown columns (sign of logfc columns)

updown_table <- sign(output_table_1[,grep("Log2fc_",colnames(output_table_1))])
colnames(updown_table) <- gsub("Log2fc","Updown",grep("Log2fc_",colnames(output_table_1),value = TRUE))
output_table_1 <- cbind(output_table_1,updown_table)
rm(updown_table)

# add contrast significance columns

sig.1_table <- output_table_1[,grep("P.value_",colnames(output_table_1))]<=.1
colnames(sig.1_table) <- gsub("P.value","Sig.1",grep("P.value_",colnames(output_table_1),value = TRUE))
output_table_1 <- cbind(output_table_1,sig.1_table)
rm(sig.1_table)

sig.05_table <- output_table_1[,grep("P.value_",colnames(output_table_1))]<=.05
colnames(sig.05_table) <- gsub("P.value","Sig.05",grep("P.value_",colnames(output_table_1),value = TRUE))
output_table_1 <- cbind(output_table_1,sig.05_table)
rm(sig.05_table)

# add volcano plot columns

log_pval_table <- log2(output_table_1[,grep("P.value_",colnames(output_table_1))])
colnames(log_pval_table) <- paste0("Log2_",colnames(log_pval_table))
output_table_1 <- cbind(output_table_1,log_pval_table)
rm(log_pval_table)
log_adj_pval_table <- log2(output_table_1[,grep("Adj.p.value_",colnames(output_table_1))])
colnames(log_adj_pval_table) <- paste0("Log2_",colnames(log_adj_pval_table))
output_table_1 <- cbind(output_table_1,log_adj_pval_table)
rm(log_adj_pval_table)

# add annotations to table

output_table_1 <- cbind(annot,output_table_1)
reduced_output_table_1 <- cbind(annot,reduced_output_table_1)
rownames(output_table_1) <- NULL
rownames(reduced_output_table_1) <- NULL

output_table_1$GOSLIM_IDS <- vapply(output_table_1$GOSLIM_IDS, paste, collapse = ", ", character(1L))
reduced_output_table_1$GOSLIM_IDS <- vapply(reduced_output_table_1$GOSLIM_IDS, paste, collapse = ", ", character(1L))

write.csv(output_table_1,file.path(DGE_output, "visualization_output_table.csv"), row.names = FALSE)
write.csv(contrasts,file.path(DGE_output, "contrasts.csv"))
write.csv(reduced_output_table_1,file.path(DGE_output, "differential_expression.csv"), row.names = FALSE)

##### Generate PCA
exp_raw <- log2(normCounts)
PCA_raw <- prcomp(t(exp_raw), scale = FALSE)
write.csv(PCA_raw$x,file.path(DGE_output, "visualization_PCA_table.csv"), row.names = TRUE)
rm(exp_raw,PCA_raw)


## print session info ##
print("Session Info below: ")
sessionInfo()
