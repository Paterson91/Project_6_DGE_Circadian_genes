#install.packages("readxl")
#install.packages("factoextra")
#install.packages("viridis")
#install.packages("rafalib")
#source("https://bioconductor.org/biocLite.R")
#biocLite("limma")
#biocLite("gplots")
library(rafalib)
library(readxl)
library("RColorBrewer")
library(limma)
library("gplots")
library(factoextra)
library(dplyr)
library(tibble)
library(tidyverse)
library(viridis)
library(DESeq2)
library(genefilter)
require("biomaRt")
library(org.Mm.eg.db)

setwd("/Users/ap14958/OneDrive - University of Bristol/Genomics Facility Bioinformatics/Project #6 Hugh Piggins/My Analysis/All genes/")
full_table = data.matrix(read.csv("QIAseqUltraplexRNA_2693_counts.csv", header = TRUE, row.names = 1, stringsAsFactors = FALSE))
rownames(full_table) <- gsub("\\..*","",rownames(full_table))

metadata_all = read.table("SampleMetaData_omit.csv", sep = ",", header = TRUE)
metadata_all

ColourStrategy_all = metadata_all$Timepoint
ColourStrategy_metatime = metadata_all$MetaTime

#Isolate samples of interest -ALL
order1_all = as.character(metadata_all$Sample)
iso_table_all = full_table[, match(order1_all, colnames(full_table))]
iso_table_all[,1:ncol(iso_table_all)] <- sapply(round(iso_table_all[,1:ncol(iso_table_all)], 0), as.integer)
head(iso_table_all)

####################################
#### Heatmap Generation
####################################

cols <- palette(magma(256))[as.fumeric(as.character(ColourStrategy_all))]
heatmap.2(iso_table_all, col=viridis(256), Colv=TRUE, 
          scale="row", key=T, keysize=1, symkey=T, dendrogram = 'column',cexRow = 1.3, lhei = c(1.5,5),
          density.info="none", labRow=FALSE, trace="none",margins=c(10,10),
          cexCol=1.3, labCol = ColourStrategy_all, ColSideColors=cols,
          main=paste0("All Clustered"))

dev.copy(png, width = 1000,height = 2000, paste0("Clustered Heatmap All.png"))
dev.off()
dev.off()

heatmap.2(iso_table_all, col=viridis(256), Colv=FALSE, 
          scale="row", key=T, keysize=1, symkey=T, dendrogram = 'column',cexRow = 1.3, lhei = c(1.5,5),
          density.info="none", labRow=FALSE, trace="none",margins=c(10,10),
          cexCol=1.3, labCol = ColourStrategy_all, ColSideColors=cols,
          main=paste0("All Clustered"))
dev.copy(png, width = 1000,height = 2000, paste0("Heatmap All.png"))
dev.off()
dev.off()
gc()

.rs.restartR()
