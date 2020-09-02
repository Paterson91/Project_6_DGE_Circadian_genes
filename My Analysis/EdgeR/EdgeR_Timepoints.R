setwd("/Users/ap14958/OneDrive - University of Bristol/Genomics Facility Bioinformatics/Project #6 Hugh Piggins/My\ Analysis/EdgeR/")
#Create DGEList data class
library("edgeR")
library("ggplot2")
require("biomaRt")
library(org.Rn.eg.db)

#control = "Day"
#test = "Night"

files <- read.table("QIAseqUltraplexRNA_2693.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE, row.names = 1)
head(files)
rowClean <- function(x){ rownames(x) <- gsub("\\..*", "", rownames(x)); x } 
counts <- rowClean(files)
head(counts)
counts$gene=NULL
counts$gene.id=NULL
counts$strand=NULL
counts$chrom=NULL
counts$loc.5.loc.3=NULL

metadata = read.csv("SampleMetaData_all.csv", header = TRUE)
head(metadata)

order1_all = as.character(metadata$Sample)

iso_table_all = counts[, match(order1_all, colnames(counts))]


x = DGEList(iso_table_all, group = metadata$Timepoint)
x$samples

keep <- filterByExpr(x)
y <- x[keep, , keep.lib.sizes=FALSE]
design = model.matrix(~0+group, data = y$samples)
colnames(design) <- levels(y$samples$group)

#Filter low counts - here we filter CPM = 0.5 (7 Counts in lowest dataset)
keep <- rowSums(cpm(y)>0.5) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]

#Normalisation

z <- calcNormFactors(y)
z$samples

#GLM Quasi-Like F Stat Differential Expression

final <- estimateDisp(z, design)
plotBCV(final)

fit <- glmQLFit(final,design)
plotQLDisp(fit)

CONTRASTS <- makeContrasts( Group1vs2 = Early_Day-Late_Day,
                            Group1vs3 = Early_Day-Mid_Day,
                            Group1vs4 = Early_Day-Mid_Night,
                            Group2vs3 = Late_Day-Mid_Day,
                            Group2vs4 = Late_Day-Mid_Night,
                            Group3vs4 = Mid_Day-Mid_Night,
                            levels = design )

for (i in 1:ncol(CONTRASTS)){
  print(CONTRASTS[,i])
  assign(paste0("glmQLFTest_",i), glmQLFTest(fit, contrast = CONTRASTS[,i]))
}

full1 = topTags(glmQLFTest_1, n = Inf, adjust.method = "BH", sort.by = "PValue")
full2 = topTags(glmQLFTest_2, n = Inf, adjust.method = "BH", sort.by = "PValue")
full3 = topTags(glmQLFTest_3, n = Inf, adjust.method = "BH", sort.by = "PValue")
full4 = topTags(glmQLFTest_4, n = Inf, adjust.method = "BH", sort.by = "PValue")
full5 = topTags(glmQLFTest_5, n = Inf, adjust.method = "BH", sort.by = "PValue")
full6 = topTags(glmQLFTest_6, n = Inf, adjust.method = "BH", sort.by = "PValue")


mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("mmusculus_gene_ensembl", mart)


annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id", "external_gene_name"),
  filter="ensembl_gene_id",
  values=rownames(full1),
  uniqueRows=TRUE)
annotLookup1 <- data.frame(full1[match(annotLookup$ensembl_gene_id, row.names(full1)),],
                           annotLookup)
annotLookup1$ensembl_gene_id=NULL
annotLookup2=annotLookup1[,c("external_gene_name","logFC","logCPM","F","PValue","FDR")]
colnames(annotLookup2)=c("Gene Name","logFC","logCPM","F","PValue","FDR")

write.csv(annotLookup2, file ="Early_Day-Late_Day-final_EdgeR.csv",row.names=TRUE)



annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id", "external_gene_name"),
  filter="ensembl_gene_id",
  values=rownames(full2),
  uniqueRows=TRUE)
annotLookup1 <- data.frame(full2[match(annotLookup$ensembl_gene_id, row.names(full2)),],
                           annotLookup)
annotLookup1$ensembl_gene_id=NULL
annotLookup2=annotLookup1[,c("external_gene_name","logFC","logCPM","F","PValue","FDR")]
colnames(annotLookup2)=c("Gene Name","logFC","logCPM","F","PValue","FDR")

write.csv(annotLookup2, file = "Early_Day-Mid_Day-final_EdgeR.csv",row.names=TRUE)



annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id", "external_gene_name"),
  filter="ensembl_gene_id",
  values=rownames(full3),
  uniqueRows=TRUE)
annotLookup1 <- data.frame(full3[match(annotLookup$ensembl_gene_id, row.names(full3)),],
                           annotLookup)
annotLookup1$ensembl_gene_id=NULL
annotLookup2=annotLookup1[,c("external_gene_name","logFC","logCPM","F","PValue","FDR")]
colnames(annotLookup2)=c("Gene Name","logFC","logCPM","F","PValue","FDR")

write.csv(annotLookup2, file = "Early_Day-Mid_Night-final_EdgeR.csv",row.names=TRUE)



annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id", "external_gene_name"),
  filter="ensembl_gene_id",
  values=rownames(full4),
  uniqueRows=TRUE)
annotLookup1 <- data.frame(full4[match(annotLookup$ensembl_gene_id, row.names(full4)),],
                           annotLookup)
annotLookup1$ensembl_gene_id=NULL
annotLookup2=annotLookup1[,c("external_gene_name","logFC","logCPM","F","PValue","FDR")]
colnames(annotLookup2)=c("Gene Name","logFC","logCPM","F","PValue","FDR")

write.csv(annotLookup2, file = "Late_Day-Mid_Day-final_EdgeR.csv",row.names=TRUE)



annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id", "external_gene_name"),
  filter="ensembl_gene_id",
  values=rownames(full5),
  uniqueRows=TRUE)
annotLookup1 <- data.frame(full5[match(annotLookup$ensembl_gene_id, row.names(full5)),],
                           annotLookup)
annotLookup1$ensembl_gene_id=NULL
annotLookup2=annotLookup1[,c("external_gene_name","logFC","logCPM","F","PValue","FDR")]
colnames(annotLookup2)=c("Gene Name","logFC","logCPM","F","PValue","FDR")

write.csv(annotLookup2, file = "Late_Day-Mid_Night-final_EdgeR.csv",row.names=TRUE)



annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id", "external_gene_name"),
  filter="ensembl_gene_id",
  values=rownames(full6),
  uniqueRows=TRUE)
annotLookup1 <- data.frame(full6[match(annotLookup$ensembl_gene_id, row.names(full6)),],
                           annotLookup)
annotLookup1$ensembl_gene_id=NULL
annotLookup2=annotLookup1[,c("external_gene_name","logFC","logCPM","F","PValue","FDR")]
colnames(annotLookup2)=c("Gene Name","logFC","logCPM","F","PValue","FDR")

write.csv(annotLookup2, file = "Mid_Day-Mid_Night-final_EdgeR.csv",row.names=TRUE)