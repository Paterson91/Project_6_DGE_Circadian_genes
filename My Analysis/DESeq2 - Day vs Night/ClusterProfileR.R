if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("biomaRt")


setwd("/Users/ap14958/OneDrive - University of Bristol/Genomics Facility Bioinformatics/Project #6 Hugh Piggins/My Analysis/DESeq2 - Day vs Night/")
library(clusterProfiler)
library(enrichplot)
require("biomaRt")
library(org.Mm.eg.db)

#LOAD GENES AND ANNOTATE

n=11

mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("mmusculus_gene_ensembl", mart)
ens = read.csv("Full_table.csv",header = T, row.names = 1)

if (exists("annotLookup_table")) {
print("Lookup Table already created")
} else {
annotLookup_table <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id", "external_gene_name","entrezgene_id"),
  filter="ensembl_gene_id",
  values=ens$gene,
  uniqueRows=TRUE)
annotLookup_match = merge(x=ens, y=annotLookup_table, by.x="gene", by.y = "ensembl_gene_id", all.x = TRUE)
}

head(ens)
head(annotLookup_table)
head(annotLookup_match)

subset_pvalue = annotLookup_match[annotLookup_match$pvalue<0.05,]
nrow(subset_pvalue)

GeneList_reads = subset_pvalue[,c(1,8:ncol(subset_pvalue))]
GeneList_reads$external_gene_name=NULL
GeneList_reads$gene_biotype=NULL
GeneList_reads$gene=NULL

#ordering <- order(rowMeans(as.numeric(GeneList_reads)),decreasing=T)[1:n]
#GeneList_reads_ordered = GeneList_reads[ordering,]
#head(GeneList_reads_ordered)


GeneList_fc = subset_pvalue[,c(1,3)]
GeneList_fc$log2FoldChange = 2^GeneList_fc$log2FoldChange
colnames(GeneList_fc)=c("Gene","FoldChange")
GeneList_fc_order = GeneList_fc[order(-GeneList_fc$FoldChange),]
head(GeneList_fc_order)
nrow(GeneList_fc_order)

### GO Classification

#Set levels
#Set choice of; CC - Cellular Component, MF - Molecular Function, BP - Biological Process

ont_selection = "CC"
level_selection = 7

ggo1 <- groupGO(gene     = as.character(GeneList_reads$entrezgene),
               keyType = "ENTREZID",
               OrgDb    = org.Mm.eg.db,
               ont      = ont_selection,
               level    = level_selection,
               readable = TRUE)
barplot(ggo1, drop=TRUE, showCategory=12, border = c(10,10))
dev.copy(png, width = 1000,height = 1000, paste0(ont_selection, "_lvl_", level_selection, "_Go.png"))
dev.off()

#Set levels
#Set choice of; CC - Cellular Component, MF - Molecular Function, BP - Biological Process

ont_selection = "CC"
level_selection = 7

de = GeneList_fc_order$Gene
edo <- enrichDGN(de)
  
  
  
barplot(ggo1, drop=TRUE, showCategory=12, border = c(10,10))
dev.copy(png, width = 1000,height = 1000, paste0(ont_selection, "_lvl_", level_selection, "_Go.png"))
dev.off()


### GO over-representation test


ego <- enrichGO(gene          = as.character(annotLookup_match$ensembl_gene_id),
                OrgDb         = org.Mm.eg.db,
                keyType = "ENSEMBL",
                ont           = ont_selection,
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)
barplot(ego, showCategory=8)
dev.copy(png, width = 1000,height = 1000, paste0(ont_selection, "_lvl_", level_selection, "_OverRep.png"))
dev.off()

### GSEA test

geneList = as.numeric(GeneList_fc_order[,2])
names(geneList) = as.character(GeneList_fc_order[,1])

ego3 <- gseGO(geneList     = geneList,
              OrgDb        = org.Mm.eg.db,
              ont          = "CC",
              keyType ="ENSEMBL",
              nPerm        = 1000,
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)

ridgeplot(ego3)

gseaplot(ego3, geneSetID = 2, title = ego3$Description[2])
gseaplot2(ego3, geneSetID = 1:6)


### Pubmed Trend of enriched terms

terms <- ego3$Description[1:3]
p <- pmcplot(terms, 2010:2017)
p2 <- pmcplot(terms, 2010:2017, proportion=FALSE)
plot_grid(p, p2, ncol=2)

pmcplot("Ifit1", 2000:2017)

# Visualisation

dotplot(ggo)
emapplot(ego)
goplot(ego)
gseaplot(kk2, geneSetID = "hsa04145")
## categorySize can be scaled by 'pvalue' or 'geneNum'
cnetplot(ego, categorySize="pvalue", foldChange=geneList)




write.table(sessionInfo(), sep = "\t", file = "Session Info")
