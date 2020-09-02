#install.packages("readxl")
#install.packages("factoextra")
#install.packages("viridis")
#install.packages("rafalib")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install()
#BiocManager::install(c("limma", "gplots","tidyverse","dplyr","biomaRt","org.Mm.eg.db"))
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
library(genefilter)
require("biomaRt")
library(org.Mm.eg.db)

setwd("/Users/ap14958/OneDrive - University of Bristol/Genomics Facility Bioinformatics/Project #6 Hugh Piggins/My Analysis/UpsetR Subset Heatmaps/")
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

metadata_all = read.table("SampleMetaData_2.csv", sep = ",", header = TRUE)
metadata_all

ColourStrategy_all = metadata_all$Timepoint

goi_role = "5"
goi_role_lookup=goi_lookup[goi_lookup$Subset==goi_role,]

#Isolate samples of interest -ALL
order1_all = as.character(metadata_all$Sample)
iso_table_all = full_table[, match(order1_all, colnames(full_table))]
iso_table_all[,1:ncol(iso_table_all)] <- sapply(round(iso_table_all[,1:ncol(iso_table_all)], 0), as.integer)
head(iso_table_all)
#Isolate Genes of interest
order2_all= as.character(goi_role_lookup$ensembl_gene_id)
iso_table1_all = iso_table_all[match(order2_all, rownames(iso_table_all)),]
iso_table1_all[,1:ncol(iso_table1_all)] <- sapply(round(iso_table1_all[,1:ncol(iso_table1_all)], 0), as.integer)
head(iso_table1_all)



####################################
#### Heatmap Generation
####################################

#Exchange Ensembl IDs for Gene Names

select1_all = dplyr::left_join(rownames_to_column(as.data.frame(iso_table1_all)), goi_role_lookup, by =c("rowname"="ensembl_gene_id"))
select1_all = na.omit(select1_all)
rownames(select1_all) = select1_all$external_gene_name
select1_all$rowname=NULL
select1_all$Subset=NULL
select1_all$external_gene_name=NULL
select1_all = as.matrix(select1_all)

cols <- palette(magma(256))[as.fumeric(as.character(ColourStrategy_all))]
heatmap.2(select1_all, col=viridis(256), Colv=TRUE, 
          scale="row", key=T, keysize=1, symkey=T, dendrogram = 'column',cexRow = 1.3, lhei = c(1.5,5),
          density.info="none", trace="none",margins=c(10,10),
          cexCol=1.3, labRow=row.names(select1_all), labCol = ColourStrategy_all, ColSideColors=cols,
          main=paste0("Subset ", goi_role, " Clustered - Omit"))
dev.copy(png, width = 1000,height = 500, paste0("Subset ", goi_role," - Clustered Heatmap All - Omit.png"))
dev.off()
 dev.off()

heatmap.2(select1_all, col=viridis(256), Colv=FALSE, 
          scale="row", key=T, keysize=1, symkey=T, dendrogram = 'none',cexRow = 1.3, lhei = c(1.5,5),
          density.info="none", trace="none",margins=c(10,10),
          cexCol=1.3, labRow=row.names(select1_all), labCol = ColourStrategy_all, ColSideColors=cols,
          main=paste0("Subset ", goi_role, " All - Omit"))
dev.copy(png, width = 1000,height = 500, paste0("Subset ", goi_role," - Heatmap All - Omit.png"))
dev.off()
dev.off()
gc()

.rs.restartR()
