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

setwd("/Users/ap14958/OneDrive - University of Bristol/Genomics Facility Bioinformatics/Project #6 Hugh Piggins/My Analysis")
full_table = data.matrix(read.csv("QIAseqUltraplexRNA_2693_counts.csv", header = TRUE, row.names = 1, stringsAsFactors = FALSE))
rownames(full_table) <- gsub("\\..*","",rownames(full_table))

goi_full = read.csv("GeneList.csv", header = TRUE)

mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("mmusculus_gene_ensembl", mart)

annotLookup1 <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id", "external_gene_name"),
  filter="external_gene_name",
  values=goi_full$external_gene_name,
  uniqueRows=TRUE)

goi_lookup = dplyr::left_join(goi_full,annotLookup1, by = "external_gene_name")
goi_lookup

metadata = read.table("SampleMetaData.csv", sep = ",", header = TRUE)
metadata

ColourStrategy = metadata$Timepoint

goi_role = "Clock Genes"
goi_role_lookup=goi_lookup[goi_lookup$Role==goi_role,]

#Isolate samples of interest
order1 = as.character(metadata$Sample)
iso_table = full_table[, match(order1, colnames(full_table))]
iso_table[,1:ncol(iso_table)] <- sapply(round(iso_table[,1:ncol(iso_table)], 0), as.integer)
head(iso_table)

#Isolate Genes of interest
order2 = as.character(goi_role_lookup$ensembl_gene_id)
iso_table1 = iso_table[match(order2, rownames(iso_table)),]
iso_table1[,1:ncol(iso_table1)] <- sapply(round(iso_table1[,1:ncol(iso_table1)], 0), as.integer)
head(iso_table1)


####################################
#### Heatmap Generation
####################################

#vst = varianceStabilizingTransformation(iso_table1, blind = TRUE)
#select <- vst[order(rowMeans(vst),decreasing=T),]

#Exchange Ensembl IDs for Gene Names

select1 = dplyr::left_join(rownames_to_column(as.data.frame(iso_table1)), goi_role_lookup, by =c("rowname"="ensembl_gene_id"))
select1 = na.omit(select1)
rownames(select1) = select1$external_gene_name
select1$rowname=NULL
select1$Role=NULL
select1$external_gene_name=NULL
select1 = as.matrix(select1)

cols <- palette(magma(256))[as.fumeric(as.character(ColourStrategy))]
heatmap.2(select1, col=viridis(256), Colv=TRUE, 
          scale="row", key=T, keysize=1, symkey=T, dendrogram = 'column',cexRow = 1.3, lhei = c(1.5,5),
          density.info="none", trace="none",margins=c(10,10),
          cexCol=1.3, labRow=row.names(select1), labCol = ColourStrategy, ColSideColors=cols,
          main=paste0(goi_role, " Clustered"))
dev.copy(png, width = 1000,height = 500, paste0(goi_role, " - Clustered Heatmap.png"))
dev.off()

heatmap.2(select1, col=viridis(256), Colv=FALSE, 
          scale="row", key=T, keysize=1, symkey=T, dendrogram = 'none',cexRow = 1.3, lhei = c(1.5,5),
          density.info="none", trace="none",margins=c(10,10),
          cexCol=1.3, labRow=row.names(select1), labCol = ColourStrategy, ColSideColors=cols,
          main=paste0(goi_role))
dev.copy(png, width = 1000,height = 500, paste0(goi_role, " - Heatmap.png"))
dev.off()




#dev.off()
#gc()

.rs.restartR()
