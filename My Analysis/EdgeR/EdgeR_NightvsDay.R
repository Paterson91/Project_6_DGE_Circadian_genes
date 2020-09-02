setwd("/Users/ap14958/OneDrive - University of Bristol/PhD/Low Salt Data/EdgeR Analysis/")
#Create DGEList data class
library("edgeR")
library("ggplot2")
require("biomaRt")
library(org.Rn.eg.db)

control = "Day"
test = "Night"

files <- read.table("QIAseqUltraplexRNA_2693_v2.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE, row.names = 1)
head(files)
rowClean <- function(x){ rownames(x) <- gsub("\\..*", "", rownames(x)); x } 
counts <- rowClean(files)

metadata = read.csv("SampleMetaData_2.csv", header = TRUE)
head(metadata)

order1_all = as.character(metadata$Sample)

iso_table_all = counts[, match(order1_all, colnames(counts))]


x = DGEList(iso_table_all, group = metadata$Group)
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

Comparison1 = makeContrasts(paste0(control,"-",test), levels = design)

qlf1 <- glmQLFTest(fit, contrast=Comparison1)
topTags(qlf1)
full1 = topTags(qlf1, n = Inf, adjust.method = "BH", sort.by = "PValue")

#columns(org.Rn.eg.db)
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

write.csv(annotLookup2, file = paste0(control,"_vs_",test, "-final_EdgeR.csv"),row.names=TRUE)













