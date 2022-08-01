# Set directories & load libraries 
MAIN <- "/Users/tkay/Desktop/Work/IGE/data/"
FIGS <- "/Users/tkay/Desktop/Work/IGE/figures/"
library('jpeg'); library('DESeq2'); library('ggplot2'); library('EnhancedVolcano')

# Import data
setwd(MAIN)
counts   <- read.csv("counts.csv")[,-(2:6)]
metadata <- read.csv("metadata.csv", row.names = 1)
rownames(counts)   <- counts[,1]
counts   <- counts[,-1]

# Are transcripts mapped or counted at different rates between the two lineages?
summary(lm(metadata$Mapping ~ metadata$Focal)) # No: mapping p-value = 0.1
summary(lm(metadata$Counting ~ metadata$Focal)) # No: mapping p-value = 0.8

# Reorder counts data and metadata
metadata <- metadata[order(rownames(metadata)),]
counts <- counts[,order(colnames(counts))]

# DEA by time-point
ddsTP <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ TimePoint)
# Low count filtering
keep <- rowMeans(counts(ddsTP))>=10
ddsTP <- ddsTP[keep,]

# Differential expression analysis
ddsTP <- DESeq(ddsTP)
summary(results(ddsTP))

# vst transformation
vsdTP <- vst(ddsTP, blind = FALSE)

# PCA
pcaDataTP <- plotPCA(vsdTP, intgroup=c("TimePoint"), returnData=TRUE)
percentVarTP <- round(100*attr(pcaDataTP, "percentVar"))

# plot PCA
setwd(FIGS)
jpeg('PCA_Stage.jpg', width=6000, height=3000, unit='px')
par(mfrow=c(1,2), mar = c(1,1,1,1), mai = c(5,5,10,1), mgp = c(30, 5, 0))
PCA_TP <- ggplot(pcaDataTP, aes(PC1, PC2, color=TimePoint)) +
  geom_point(size=50) +
  xlab(paste0("PC1: ", percentVarTP[1], "% variance")) +
  ylab(paste0("PC2: ", percentVarTP[2], "% variance")) +
  coord_fixed() +
  theme(axis.text=element_text(size=0), axis.title=element_text(size=200,face="bold"),
        legend.text=element_text(size=150), legend.title=element_text(size=0))
PCA_TP
dev.off()

# Separate the samples by age
AdultMetaData  <- metadata[metadata$TimePoint == "Day_8",]
LarvaeMetaData <- metadata[metadata$TimePoint == "Day_1",]
AdultCounts    <- counts[,colnames(counts) %in% rownames(AdultMetaData)]
LarvaeCounts   <- counts[,colnames(counts) %in% rownames(LarvaeMetaData)]

# Day_8 samples
# DE by IGE
ddsAI <- DESeqDataSetFromMatrix(countData = AdultCounts,
                                colData = AdultMetaData,
                                design = ~ Env + Focal)
ddsAI <- DESeq(ddsAI[keep,])
summary(results(ddsAI))

# DE by DGE
ddsAD <- DESeqDataSetFromMatrix(countData = AdultCounts,
                                colData = AdultMetaData,
                                design = ~ Env + Focal)
ddsAD <- DESeq(ddsAD[keep,])
summary(results(ddsAD))

# Transform and plot
vsdAD <- vst(ddsAD, blind = FALSE)
pcaDataAD <- plotPCAD(vsdAD, intgroup=c("Focal"), returnData=TRUE)
percentVarAD <- round(100*attr(pcaDataAD, "percentVar"))
jpeg('PCA_day8.jpg', width=6000, height=3000, unit='px')
par(mfrow=c(1,1), mar = c(1,1,1,1), mai = c(5,5,10,1), mgp = c(30, 5, 0))
PCA_AD <- ggplot(pcaDataAD, aes(PC1, PC2, color=Focal, shape = as.character(AdultMetaData[,2]))) +
  geom_point(size=50) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  theme(axis.text=element_text(size=0), axis.title=element_text(size=200,face="bold"),
        legend.text=element_text(size=150), legend.title=element_text(size=0))
PCA_AD
dev.off()
jpeg('VolcanoPlot_day8.jpg', width=3000, height=3000, unit='px')
par(mfrow=c(1,1), mar = c(1,1,1,1), mai = c(5,5,10,1), mgp = c(30, 5, 0))
EnhancedVolcano(results(ddsAD),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1.5,
                lab = "",
                cutoffLineWidth = 2,
                pointSize = 50,
                title = 'Day_8') +
  theme(axis.text.x = element_text(face="bold", size=80),
        axis.text.y = element_text(face="bold", size=80),
        axis.line = element_line(size = 3),
        axis.title = element_text(size = 100),
        plot.title =  element_text(size = 150),
        legend.position = "none")
dev.off()

# Day_1 samples
# DE by IGE
ddsLI <- DESeqDataSetFromMatrix(countData = LarvaeCounts,
                                colData = LarvaeMetaData,
                                design = ~ Focal + Env)
ddsLI <- DESeq(ddsLI[keep,])
summary(results(ddsLI))

# DE by DGE
ddsLD <- DESeqDataSetFromMatrix(countData = LarvaeCounts,
                                colData = LarvaeMetaData,
                                design = ~ Env + Focal)
ddsLD <- DESeq(ddsLD[keep,])
summary(results(ddsLD))

# Transform and plot
vsdLD <- vst(ddsLD, blind = FALSE)
pcaDataLD <- plotPCA(vsdLD, intgroup=c("Focal"), returnData=TRUE)
percentVarLD <- round(100*attr(pcaDataLD, "percentVar"))
jpeg('PCA_day1.jpg', width=6000, height=3000, unit='px')
par(mfrow=c(1,2), mar = c(1,1,1,1), mai = c(5,5,10,1), mgp = c(30, 5, 0))
PCA_LD <- ggplot(pcaDataLD, aes(PC1, PC2, color=Focal, shape = as.character(LarvaeMetaData[,2]))) +
  geom_point(size=50) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  theme(axis.text=element_text(size=0), axis.title=element_text(size=200,face="bold"),
        legend.text=element_text(size=150), legend.title=element_text(size=0))
PCA_LD
dev.off()
jpeg('VolcanoPlot_day1.jpg', width=3000, height=3000, unit='px')
par(mfrow=c(1,1), mar = c(1,1,1,1), mai = c(5,5,10,1), mgp = c(30, 5, 0))
EnhancedVolcano(results(ddsL),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1.5,
                lab = "",
                cutoffLineWidth = 2,
                pointSize = 50,
                title = 'Day_1') +
  theme(axis.text.x = element_text(face="bold", size=80),
        axis.text.y = element_text(face="bold", size=80),
        axis.line = element_line(size = 3),
        axis.title = element_text(size = 100),
        plot.title =  element_text(size = 150),
        legend.position = "none")
dev.off()

# Calculate transcriptomic variance within day_8 and day_1 samples
countsddsAD <- counts(ddsAD, normalized = TRUE)
countsddsLD <- counts(ddsLD, normalized = TRUE)
day_8SD <- c()
day_1SD <- c()
for (row in 1:nrow(countsddsAD)){
  day_8SD <- c(day_8SD, sd(countsddsAD[row,]))
}
for (row in 1:nrow(countsddsLD)){
  day_1SD <- c(day_1SD, sd(countsddsLD[row,]))
}
# Compare variance between the two time-points
wilcox.test(day_1SD, day_8SD, paired = TRUE)

# DE using data and controlling for age
# DE by DGE
ddsAllD <- DESeqDataSetFromMatrix(countData = counts,
                                  colData = metadata,
                                  design = ~ TimePoint + Env + Focal)
ddsAllD <- DESeq(ddsAllD[keep,])
summary(results(ddsAllD))

# DE by IGE
ddsAllI <- DESeqDataSetFromMatrix(countData = counts,
                                  colData = metadata,
                                  design = ~ TimePoint + Env + Focal)
ddsAllI <- DESeq(ddsAllI[keep,])
summary(results(ddsAllI))

# Using all data and controlling for age and batch effects (extraction group and sequencing lane)
# DE by DGE
ddsAllCD <- DESeqDataSetFromMatrix(countData = counts,
                                   colData = metadata,
                                   design = ~ TimePoint + Lane + Extraction + Env + Focal)
ddsAllCD <- DESeq(ddsAllCD[keep,])
summary(results(ddsAllCD))

# DE by IGE
ddsAllCI <- DESeqDataSetFromMatrix(countData = counts,
                                   colData = metadata,
                                   design = ~ TimePoint + Lane + Extraction + Focal + Env)
ddsAllCI <- DESeq(ddsAllCI[keep,])
summary(results(ddsAllCI))

# Check for an interaction between DGE and IGE
ddsAllI <- DESeqDataSetFromMatrix(countData = counts,
                                  colData = metadata,
                                  design = ~ TimePoint + Lane + Extraction + Focal + Env + Focal:Env)
ddsAllI <- DESeq(ddsAllI[keep,])
summary(results(ddsAllI))

# Quantify batch effects
# Extraction group
ddsAllE <- DESeqDataSetFromMatrix(countData = counts,
                                  colData = metadata,
                                  design = ~ TimePoint + Lane + Focal + Env + Extraction)
ddsAllE <- DESeq(ddsAllE[keep,])
summary(results(ddsAllE))

# Sequencing lane
ddsAllL <- DESeqDataSetFromMatrix(countData = counts,
                                  colData = metadata,
                                  design = ~ TimePoint + Focal + Env + Extraction + Lane)
ddsAllL <- DESeq(ddsAllL[keep,])
summary(results(ddsAllL))