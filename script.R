#Author: Ivy Kosater, June 2024
library(dplyr)
library(tidyr)
library(limma)
library(edgeR)
library(org.Hs.eg.db)
library(annotate)
library(ggpubr)
library(DESeq2)

#Main function for running limma-voom 
limma_run <- function(counts_rds, metadata_rds){
  df <- readRDS(counts_rds)
  meta <- readRDS(metadata_rds)
  
  #Change colnames in df to IDs
  colnames(df) =  meta$`RNAseq ID`
  
  #Create dgelist obj
  d0 <- DGEList(df)
  
  #Calculate normalization factors
  d0 <- calcNormFactors(d0)
  
  #Filter lowly-expressed genes (cutoff can be flexible)
  cutoff <- 20
  drop <- which(apply(cpm(d0),1,max) < cutoff)
  d <- d0[-drop,]
  dim(d)
  
  subtype <- meta$Subtypes
  prop.score <- meta$`EXPANSION PROPENSITY SCORE`
  cap <- meta$`CAP100 Score`
  
  #Voom transformation
  mm <- model.matrix(~prop.score + subtype) #Model used for limma. Testing for differential gene expression based on propensity score with
  #cell subtype as covariate 
  
  #Fit the linear model, save DE table
  fit <- lmFit(y, mm)
  fit <- eBayes(fit)
  top <- topTable(fit, coef="prop.score",
                  number = length(rownames(df)))
  
  x <- data.frame(getSYMBOL(rownames(top), data='org.Hs.eg'))
  
  top$symbol <- x$getSYMBOL.rownames.top...data....org.Hs.eg..
  
  write.csv(top, "results.table.csv")
}

#This function runs deseq with the same model
deseq_run <- function(counts_rds, metadata_rds){
  
  df <- readRDS(counts_rds)
  meta <- readRDS(metadata_rds)
  
  #Change colnames in df to IDs
  colnames(df) =  meta$`RNAseq ID`
  
  names(meta)[names(meta) == 'EXPANSION PROPENSITY SCORE'] <- 'prop.score'
  dds <- DESeqDataSetFromMatrix(as.matrix(df), colData = meta, design = ~prop.score + Subtypes)
  dds <- DESeq(dds)
  res <- results(dds, name="prop.score")
  
  x <- data.frame(getSYMBOL(rownames(res.df), data='org.Hs.eg'))
  res.df$symbol <- x$getSYMBOL.rownames.res.df...data....org.Hs.eg..
}

