


buildTargets <- function(opt) {
  cat("Reading ISAtab study metadata file..., ",opt$isafile$name, "\n")
  na_strings = c('NA','na','null','NULL','Null','')
 
  td <- tempdir()
  unlink(list.files(td, full.names = TRUE))
  unzip(opt$isafile$datapath, exdir = td, junkpaths = FALSE)
  isa <- Risa::readISAtab(path = td)

  assay_idx <- which((isa@assay.technology.types == "DNA microarray") & (isa@assay.measurement.types == "transcription profiling"))
  all_targets<-vector("list", length(assay_idx))
  for (which_assay in 1:length(assay_idx)){
    table <- isa@assay.tabs[[assay_idx[which_assay]]]@assay.file
    files<-data.frame(V1=table$`Array Data File`)
    try({files$V2<-table$`Comment[Array Data File Name]`})
    try({files$V3<-table$`Parameter Value[array data file,http://www.ebi.ac.uk/efo/EFO_0004098,EFO]`})
    try({files$V4<-table$`Characteristics[array data file,http://www.ebi.ac.uk/efo/EFO_0004098,EFO]`})
    opt_files <- limma::removeExt(opt$datafiles$name, ".gz")
    files_match<- as.data.frame(sapply(opt_files,function(x){grepl(x,files)}))
    if (dim(files)[2]==1){
      files_match <- t(which(apply(files_match,2,any)))
    }else{
      files_match <- which(apply(files_match,1,any))
    }
    
    files<-files[,files_match]
    # replace <- sapply(opt_files,function(x){grep(x,files)})
    # files[replace]<-opt_files

    target<-list()
    study_samples <- which(table$`Sample Name` %in% isa@study.files[[1]]$`Sample Name`)
    target$labels = length(unique(table$Label))
    labelnames <- unique(table$Label)
    target$seperate_channel_files <- eval((length(unique(files)) == 2*length(unique(table$`Hybridization Assay Name`))))
    factors <- as.data.frame(as.character(isa@factors[[which_assay]][[1]]))
    #factors <- as.data.frame(table[,which(grepl("Factor",colnames(table)))])
    factors <- as.data.frame(factors[study_samples,])
    colnames(factors)<-paste("factor",1:dim(factors)[2], sep = "_")
    is.na(factors) <- Reduce("|", lapply(na_strings, "==", factors)) # avoids whitespace issues for group strings
    group <- apply(factors,1,function(x) {paste(na.omit(x),collapse=" & ")}) # groups are concatenations of non empty factor values per sample
    table$`Hybridization Assay Name` <- sub("\\ \\d+", "", table$`Hybridization Assay Name`)
    isa@study.files[[1]]$`Source Name`[study_samples] <- sub("\\ \\d+", "", isa@study.files[[1]]$`Source Name`[study_samples])
    library(tidyr)
    if ((target$labels == 1) && (target$seperate_channel_files == FALSE)){
      cat("One color microarray design...","\n")
      target$t1 <- data.frame(SampleName=table$`Sample Name`, Group=group, ArrayName=table$`Hybridization Assay Name`, FileName=files)
      target$t1$SourceName <- isa@study.files[[1]]$`Source Name`[study_samples]
      target$t1 <- cbind(target$t1,factors)
      target$paired_samples <- (length(unique(target$t1$SourceName)) < length(unique(target$t1$SampleName)))
    }else if ((target$labels == 2) && (target$seperate_channel_files == FALSE)){
      cat("Two color microarray design...","\n")
      target$t1 <- data.frame(Label=table$Label, Group=group, ArrayName=table$`Hybridization Assay Name`, FileName=files)
      #target$t1 <- data.frame(Label=table$Label, Group=group, ArrayName=files, FileName=files)
      
      target$t1 <- pivot_wider(target$t1, names_from = Label, values_from = Group)
    }else if ((target$labels == 1) && (target$seperate_channel_files == TRUE)){
      cat("One color microarray design with seperate channel files...","\n")
      target$t1 <- data.frame(SampleName=table$`Sample Name`, Group=group, ArrayName=table$`Hybridization Assay Name`, FileName=files)      
      target$t1$SourceName <- isa@study.files[[1]]$`Source Name`[study_samples]
      target$t1 <- target$t1[order(target$t1$SampleName,target$t1$FileName),]
      target$t1$Label <- rep(c("FileName_1","FileName_2"),dim(table)[1]/2)
      target$t1 <- pivot_wider(target$t1, names_from = Label, values_from = FileName)
      target$t1 <- cbind(target$t1,factors)
      target$paired_samples <- (length(unique(target$t1$SourceName)) < length(unique(target$t1$SampleName)))
    }else if ((target$labels == 2) && (target$seperate_channel_files == TRUE)){
      cat("Two color microarray design with seperate channel files...","\n")
      target$t1 <- data.frame(Label=table$Label, Group=group, ArrayName=table$`Hybridization Assay Name`, FileName=files)
      target$t1 <- target$t1[order(target$t1$FileName,target$t1$Label),]
      target$t1$FileLabel <- paste0("FileName",target$t1$Label)
      t1a <- pivot_wider(target$t1[,c("Label","Group","ArrayName")], names_from = Label, values_from = Group)
      t1b <- pivot_wider(target$t1[,c("FileLabel","FileName","ArrayName")], names_from = FileLabel, values_from = FileName)
      target$t1 <- dplyr::left_join(t1a,t1b,by = "ArrayName")
    }else {
      cat("Error detecting labeling pattern", "\n")
    }
    
    target$technical_replicates <- (length(unique(target$t1$ArrayName)) < dim(target$t1)[1])
    
    if (target$labels == 2){
      group_pairs <- target$t1[,labelnames]
      Ref <- as.data.frame(t(apply(group_pairs,1,function(x){unique(group) %in% x})))
      Ref <- apply(Ref,2,all)
      Ref_group <- unique(group)[Ref] # groups common to all arrays
      sample_pairs <- table[,which(colnames(table) %in% c("Label","Sample Name", "Hybridization Assay Name"))]
      sample_pairs <- pivot_wider(sample_pairs, names_from = Label, values_from = `Sample Name`, id_cols = `Hybridization Assay Name`)
      sample_pairs <- sample_pairs[,-c(1)]
      common_samples <- as.data.frame(t(apply(sample_pairs,1,function(x){unique(table$`Sample Name`) %in% x})))
      common_samples <- unique(table$`Sample Name`)[which(apply(common_samples,2,all))]
      
      contrasts <- unique(rbind(group_pairs,group_pairs[,c(2,1)],deparse.level = 0)) # available contrast group pairs
      contrast_space<-t(combn(unique(group),2)) # set of all contrast group pairs

      if(length(Ref_group) == 2){
        cat("Replicate array design ...", "\n")
        target$design <- "Replicate Array"
      }else if (length(common_samples) == 1){
        cat("Commnon reference design ...", "\n")
        target$design <- "Common Reference"
        replace_group <- target$t1[,labelnames]
        replace_group[sample_pairs == common_samples] <- "Ref"
        target$t1[,labelnames] <- replace_group
      }else if ((length(Ref_group) <= 1) && (all(apply(contrast_space,1,function(x){x %in% contrasts})))){
        cat("Direct two color design ...", "\n")
        target$design <- "Direct Two-Color"
      }else if ((length(Ref_group) <= 1) && (!all(apply(contrast_space,1,function(x){x %in% contrasts})))){
        cat("Seperate channel unconnected design ...", "\n")
        target$design <- "Separate Channels"
      }
    }
    
    label_cols <- which(colnames(target$t1) %in% labelnames)
    target$t2 <- target$t1
    try({target$t2$Group <- sapply(target$t1$Group,function(x){paste0("(",x,")")})},silent = TRUE)
    try({target$t2[,label_cols[1]] <- sapply(target$t1[,label_cols[1]],function(x){paste0("(",x,")")})},silent = TRUE)
    try({target$t2[,label_cols[2]] <- sapply(target$t1[,label_cols[2]],function(x){paste0("(",x,")")})},silent = TRUE)
    
    target$t3 <- target$t1
    try({target$t3$Group <- sapply(target$t1$Group,function(x){make.names(x,unique = FALSE, allow_ = TRUE)})},silent = TRUE) #Make syntactically valid group names
    try({target$t3[,label_cols[1]] <- sapply(target$t1[,label_cols[1]],function(x){make.names(x,unique = FALSE, allow_ = TRUE)})},silent = TRUE)
    try({target$t3[,label_cols[2]] <- sapply(target$t1[,label_cols[2]],function(x){make.names(x,unique = FALSE, allow_ = TRUE)})},silent = TRUE)
    
 
    
    all_targets[[which_assay]]<-target
  }
  return(all_targets)
}

gfun <- function(expr, cov, p=0.05) {
  library(genefilter)
  f2 <- ttest(cov, p=p)
  ffun <- filterfun(f2)
  which <- genefilter(expr, ffun)
}


knnCV <- function(EXPR, selectfun, cov, Agg, pselect = 0.01, Scale=FALSE) {
  libary(genefilter)
  nc <- ncol(EXPR)
  outvals <- rep(NA, nc)
  for(i in 1:nc) {
    v1 <- EXPR[,i]
    expr <- EXPR[,-i]
    glist <- selectfun(expr, cov[-i], p=pselect)
    expr <- expr[glist,]
    if( Scale ) {
      expr <- scale(expr)
      v1 <- as.vector(scale(v1[glist]))
    }
    else
      v1 <- v1[glist]
    out <- paste("iter ",i, " num genes= ", sum(glist), sep="")
    print(out)
    Aggregate(row.names(expr), Agg)
    if( length(v1) == 1)
      outvals[i] <- knn(expr, v1, cov[-i], k=5)
    else
      outvals[i] <- knn(t(expr), v1, cov[-i], k=5)
  }
  return(outvals)
}

