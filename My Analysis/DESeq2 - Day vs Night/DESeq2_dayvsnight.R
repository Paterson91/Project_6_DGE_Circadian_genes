setwd("/Users/ap14958/OneDrive - University of Bristol/Genomics Facility Bioinformatics/Project #6 Hugh Piggins/My Analysis/DESeq2 - Day vs Night/")

library("DESeq2")
library(viridis)


# Set the prefix for each output file name

sampleFiles<- c("QIAseqUltraplexRNA_2693_v2.csv")

countData <- as.matrix(read.csv("QIAseqUltraplexRNA_2693_v2.csv", header=TRUE,
                                  row.names = 1,
                                  as.is=TRUE))
head(countData)

sampleNames <- as.data.frame(colnames(countData))
colnames(sampleNames)="SampleNames"

metadata = read.csv("SampleMetaData_2.csv", header = TRUE)

order1_all = as.character(metadata$Sample)

iso_table_all = countData[,match(order1_all,colnames(countData))]
colnames(iso_table_all)

rownames(iso_table_all) = gsub("\\..*","",rownames(iso_table_all))

colData = merge(x=sampleNames, y=metadata, by.x="SampleNames", by.y="Sample", all.y=TRUE)

ddsHTSeq = DESeqDataSetFromMatrix(iso_table_all, colData, design = ~Group, ignoreRank = FALSE)

# Guts of DESeq2
dds <- DESeq(ddsHTSeq)
keep <- rowSums(counts(dds)) >= 1
ddskeep = dds[keep,]

res1 = results(ddskeep,contrast=c("Group","Night","Day"), alpha = 0.05)

# save data results and normalized reads to csv

resdata1 <- merge(as.data.frame(res1), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata1)[1] <- 'gene'
write.csv(resdata1, file = "Full_table.csv")

#order results by padj (p-adjusted) value (most significant to least)
#res1= subset(res1, padj<0.05)  #No need to remove all under 0.05

res1.1 <- res1[order(resdata1$padj),]

#source("http://bioconductor.org/biocLite.R")
#
library(biobroom)

res_tidy1 <- tidy.DESeqResults(res1.1)
head(res_tidy1)
library(dplyr)
library(annotables)

final_tidy1 = res_tidy1 %>%
  arrange(p.adjusted) %>% 
  inner_join(grcm38, by=c("gene"="ensgene")) %>%
  dplyr::select(gene, estimate, p.value, p.adjusted, symbol, description)

write.csv(final_tidy1, file = "final_tidied_DESeq2.csv",row.names=FALSE)

# send normalized counts to tab delimited file for GSEA, etc.
write.table(as.data.frame(counts(dds),normalized=T), file = "normalized_counts_DESeq2.txt", sep = '\t')

# Produce DataFrame of results of statistical tests
mcols(res1, use.names = T)
write.csv(as.data.frame(mcols(res1, use.name = T)),file = "Test-conditions_DESeq2.csv")

# transform raw counts into normalized values. variance stabilization - variance stabilization is very good for heatmaps, etc.
vsd <- varianceStabilizingTransformation(dds, blind=T)

#PCA Plots
library(ggplot2)
PCA=plotPCA(vsd, intgroup = "Group", ntop = 100000, returnData = FALSE)
ggsave("PCA_DESeq2.jpg", plot = PCA)

ColourStrategy_all = metadata$Group

library("RColorBrewer")
library("gplots")
library(rafalib)

# 100 top expressed genes with heatmap.2
select <- order(rowMeans(counts(dds,normalized=T)),decreasing=T)[1:1000]
cols <- palette(magma(256))[as.fumeric(as.character(ColourStrategy_all))]
heatmap.2(assay(vsd)[select,], col=viridis(256), Colv=TRUE,
          scale="row", key=T, keysize=1, symkey=T, labCol = ColourStrategy_all, ColSideColors=cols,lhei = c(1.5,5),
          density.info="none", trace="none",margins=c(6,5),
          cexCol=1.1, labRow=F,
          main="1000 Top Expressed Genes Heatmap")
dev.copy(png, width = 1000,height = 1000, "1000-Heatmap_DESeq2.png")
dev.off()

.rs.restartR()
