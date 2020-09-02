#install.packages("readxl")
#install.packages("factoextra")
#install.packages("viridis")
#install.packages("rafalib")
#source("https://bioconductor.org/biocLite.R")
#biocLite("limma")
#biocLite("gplots")
#BiocManager::install("org.Mm.eg.db")
library(rafalib)
library(readxl)
library("RColorBrewer")
library(limma)
library("gplots")
library(factoextra)
library(dplyr)
library(tidyverse)
library(viridis)
library(genefilter)
library(reshape2)

library(clusterProfiler)
require("biomaRt")
library(org.Mm.eg.db)

setwd("/Users/ap14958/OneDrive - University of Bristol/Genomics Facility Bioinformatics/Project #6 Hugh Piggins/My Analysis/All Genes/")
full_table = data.matrix(read.csv("QIAseqUltraplexRNA_2693_counts.csv", header = TRUE, row.names = 1, stringsAsFactors = FALSE))

metadata = read.table("SampleMetaData_omit.csv", sep = "\t", header = TRUE)

metadata

#Samples to compare

set1 = "Exp7"
set2 = "Exp6"
numbertoisolate = 1000
#Rank by count
X=15

metadata_set1=metadata[metadata$Condition==set1,]
metadata_set2=metadata[metadata$Condition==set2,]

order1 = as.character(metadata_set1$Original_Sample_Name)
iso_table1 = full_table[, match(order1, colnames(full_table))]
iso_table1[,1:ncol(iso_table1)] <- sapply(round(iso_table1[,1:ncol(iso_table1)], 0), as.integer)
head(iso_table1)

order2 = as.character(metadata_set2$Original_Sample_Name)
iso_table2 = full_table[, match(order2, colnames(full_table))]
iso_table2[,1:ncol(iso_table1)] <- sapply(round(iso_table2[,1:ncol(iso_table2)], 0), as.integer)
head(iso_table2)

select1 <- iso_table1[order(rowMeans(iso_table1),decreasing=T)[1:numbertoisolate],]
select2 <- iso_table2[order(rowMeans(iso_table2),decreasing=T)[1:numbertoisolate],]

select1= data.frame(GeneID = rownames(select1), select1, row.names = NULL)
select2= data.frame(GeneID = rownames(select2), select2, row.names = NULL)

head(select1)
head(select2)

mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("mmusculus_gene_ensembl", mart)

annotLookup1 <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id", "gene_biotype", "external_gene_name","entrezgene"),
  filter="ensembl_gene_id",
  values=select1$GeneID,
  uniqueRows=TRUE)

annotLookup1 <- data.frame(
  select1$GeneID[match(annotLookup1$ensembl_gene_id, select1$GeneID)],
  annotLookup1)

annotLookup2 <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id", "gene_biotype", "external_gene_name","entrezgene"),
  filter="ensembl_gene_id",
  values=select2$GeneID,
  uniqueRows=TRUE)

annotLookup2 <- data.frame(
  select2$GeneID[match(annotLookup2$ensembl_gene_id, select2$GeneID)],
  annotLookup2)

#Gene Ontology

ggo_MF1 <- groupGO(gene     = as.character(annotLookup1$entrezgene),
                keyType = "ENTREZID",
                OrgDb    = org.Mm.eg.db,
                ont      = "MF",
                level    = 3,
                readable = TRUE)
ggo_MF2 <- groupGO(gene     = as.character(annotLookup2$entrezgene),
                keyType = "ENTREZID",
                OrgDb    = org.Mm.eg.db,
                ont      = "MF",
                level    = 3,
                readable = TRUE)
ggo_BP1 <- groupGO(gene     = as.character(annotLookup1$entrezgene),
                   keyType = "ENTREZID",
                   OrgDb    = org.Mm.eg.db,
                   ont      = "BP",
                   level    = 3,
                   readable = TRUE)
ggo_BP2 <- groupGO(gene     = as.character(annotLookup2$entrezgene),
                   keyType = "ENTREZID",
                   OrgDb    = org.Mm.eg.db,
                   ont      = "BP",
                   level    = 3,
                   readable = TRUE)
ggo_CC1 <- groupGO(gene     = as.character(annotLookup1$entrezgene),
                   keyType = "ENTREZID",
                   OrgDb    = org.Mm.eg.db,
                   ont      = "CC",
                   level    = 3,
                   readable = TRUE)
ggo_CC2 <- groupGO(gene     = as.character(annotLookup2$entrezgene),
                   keyType = "ENTREZID",
                   OrgDb    = org.Mm.eg.db,
                   ont      = "CC",
                   level    = 3,
                   readable = TRUE)

#barplot(ggo_MF2, drop=TRUE, showCategory=12, border = c(10,10))
#dotplot(ego2, showCategory=12)

#Plotting

ggo_MF1.1 = as.data.frame(ggo_MF1)
ggo_MF1.2 <- ggo_MF1.1[order(ggo_MF1.1$Count,decreasing=T),]
ggo_MF1.3 = ggo_MF1.2[1:X,]
head(ggo_MF1.3)

ggo_MF2.1 = as.data.frame(ggo_MF2)
ggo_MF2.2 <- ggo_MF2.1[order(ggo_MF2.1$Count,decreasing=T),]
ggo_MF2.3 = ggo_MF2.2[1:X,]
head(ggo_MF2.3)

ggo_BP1.1 = as.data.frame(ggo_BP1)
ggo_BP1.2 <- ggo_BP1.1[order(ggo_BP1.1$Count,decreasing=T),]
ggo_BP1.3 = ggo_BP1.2[1:X,]
head(ggo_BP1.3)

ggo_BP2.1 = as.data.frame(ggo_BP2)
ggo_BP2.2 <- ggo_BP2.1[order(ggo_BP2.1$Count,decreasing=T),]
ggo_BP2.3 = ggo_BP2.2[1:X,]
head(ggo_BP2.3)

ggo_CC1.1 = as.data.frame(ggo_CC1)
ggo_CC1.2 <- ggo_CC1.1[order(ggo_CC1.1$Count,decreasing=T),]
ggo_CC1.3 = ggo_CC1.2[1:X,]
head(ggo_CC1.3)

ggo_CC2.1 = as.data.frame(ggo_CC2)
ggo_CC2.2 <- ggo_CC2.1[order(ggo_CC2.1$Count,decreasing=T),]
ggo_CC2.3 = ggo_CC2.2[1:X,]
head(ggo_MF2.3)

merge = merge(ggo_MF1.3, ggo_MF2.3, by="ID")
merge
compare = merge[,c(2,3,7)]
colnames(compare)=c("Description",set1,set2)
compare.m=melt(compare, id.vars = "Description")
colnames(compare.m) = c("GO Term", "Experiment", "Count")
MF <- ggplot(data=compare.m, aes(x=`GO Term`, y=Count, fill=Experiment)) +
  geom_bar(stat="identity", position=position_dodge()) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
MF

merge = merge(ggo_BP1.3, ggo_BP2.3, by="ID")
merge
compare = merge[,c(2,3,7)]
colnames(compare)=c("Description",set1,set2)
compare.m=melt(compare, id.vars = "Description")
colnames(compare.m) = c("GO Term", "Experiment", "Count")
BP <- ggplot(data=compare.m, aes(x=`GO Term`, y=Count, fill=Experiment)) +
  geom_bar(stat="identity", position=position_dodge()) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
BP

merge = merge(ggo_CC1.3, ggo_CC2.3, by="ID")
merge
compare = merge[,c(2,3,7)]
colnames(compare)=c("Description",set1,set2)
compare.m=melt(compare, id.vars = "Description")
colnames(compare.m) = c("GO Term", "Experiment", "Count")
CC <- ggplot(data=compare.m, aes(x=`GO Term`, y=Count, fill=Experiment)) +
  geom_bar(stat="identity", position=position_dodge()) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
CC


#Gene Enrichment Analysis

#ego1 <- enrichGO(gene          = as.character(annotLookup1$entrezgene),
#                 OrgDb         = org.Mm.eg.db,
#                 keyType = "ENTREZID",
#                 ont           = "MF",
#                 pAdjustMethod = "BH",
#                 pvalueCutoff  = 0.01,
#                 qvalueCutoff  = 0.05,
#                 readable      = TRUE)

#ego2 <- enrichGO(gene          = as.character(annotLookup2$entrezgene),
#                 OrgDb         = org.Mm.eg.db,
#                 keyType = "ENTREZID",
#                 ont           = "MF",
#                 pAdjustMethod = "BH",
#                 pvalueCutoff  = 0.01,
#                 qvalueCutoff  = 0.05,
#                 readable      = TRUE)