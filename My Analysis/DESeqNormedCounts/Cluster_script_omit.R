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

setwd("/Users/ap14958/OneDrive - University of Bristol/Genomics Facility Bioinformatics/Project #6 Hugh Piggins/My Analysis/UpsetR Comparison/")
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

metadata_omit = read.table("SampleMetaData_2.csv", sep = ",", header = TRUE)
metadata_omit

ColourStrategy_omit = metadata_omit$Timepoint


goi_role_lookup=goi_lookup

#Isolate samples of interest -Omit
order1_omit = as.character(metadata_omit$Sample)
iso_table_omit = full_table[, match(order1_omit, colnames(full_table))]
iso_table_omit[,1:ncol(iso_table_omit)] <- sapply(round(iso_table_omit[,1:ncol(iso_table_omit)], 0), as.integer)
head(iso_table_omit)
#Isolate Genes of interest
order2_omit= as.character(goi_role_lookup$ensembl_gene_id)
iso_table1_omit = iso_table_omit[match(order2_omit, rownames(iso_table_omit)),]
iso_table1_omit[,1:ncol(iso_table1_omit)] <- sapply(round(iso_table1_omit[,1:ncol(iso_table1_omit)], 0), as.integer)
head(iso_table1_omit)


####################################
#### Heatmap Generation
####################################

#Exchange Ensembl IDs for Gene Names

select1_omit = dplyr::left_join(rownames_to_column(as.data.frame(iso_table1_omit)), goi_role_lookup, by =c("rowname"="ensembl_gene_id"))
select1_omit = na.omit(select1_omit)
rownames(select1_omit) = select1_omit$external_gene_name
select1_omit$rowname=NULL
select1_omit$Role=NULL
select1_omit$external_gene_name=NULL
select1_omit = as.matrix(select1_omit)

cols <- palette(magma(256))[as.fumeric(as.character(ColourStrategy_omit))]
heatmap.2(select1_omit, col=viridis(256), Colv=TRUE, 
          scale="row", key=T, keysize=1, symkey=T, dendrogram = 'column',cexRow = 1.3, lhei = c(1.5,5),
          density.info="none", trace="none",margins=c(10,10),
          cexCol=1.3, labRow=row.names(select1_omit), labCol = ColourStrategy_omit, ColSideColors=cols,
          main=paste0(goi_role, " Clustered"))
dev.copy(png, width = 1000,height = 500, paste0(goi_role, " - Clustered Heatmap Omit.png"))
dev.off()
dev.off()

heatmap.2(select1_omit, col=viridis(256), Colv=FALSE, 
          scale="row", key=T, keysize=1, symkey=T, dendrogram = 'none',cexRow = 1.3, lhei = c(1.5,5),
          density.info="none", trace="none",margins=c(10,10),
          cexCol=1.3, labRow=row.names(select1_omit), labCol = ColourStrategy_omit, ColSideColors=cols,
          main=paste0(goi_role))
dev.copy(png, width = 1000,height = 500, paste0(goi_role, " - Heatmap Omit.png"))
dev.off()
dev.off()


#dev.off()
gc()

.rs.restartR()

