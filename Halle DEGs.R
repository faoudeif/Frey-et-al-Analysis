library( "DESeq2" )
library(ggplot2)
library(plyr); library(dplyr)
library(radiant)
library(pheatmap)
library(usethis)
library(devtools)
library(ashr)

log2fHalleLiver <- read.csv('DEGs_padj_0.05_fc_1.csv', header = TRUE, sep = ',')
countsWeizmanLiver <- read.csv('countsWeizmanLiver.csv', header = TRUE, sep = ",")
log2fRobinetteLiver <- read.csv('log2FoldRobinetteLiver.csv', header = TRUE, sep = ",")
countsWeizmanSpleen <- read.csv('countsWeizmanSpleen.csv', header = TRUE, sep = ",")
log2fRobinetteSpleen <- read.csv('log2FoldRobinetteSpleen.csv', header = TRUE, sep = ",")
log2fRobinetteSiLP <- read.csv('log2FoldRobinettesiLP.csv', header = TRUE, sep = ",")


log2fRobinetteSpleen <- log2fRobinetteSpleen[,1:2]
log2fRobinetteLiver <- log2fRobinetteLiver[,1:2]
log2fRobinetteSiLP <- log2fRobinetteSiLP[,1:2]
log2fHalleLiver <- log2fHalleLiver[,c("X","log2FoldChange")]

metadataWeizmanSpleen <- read.csv('metadataWeizmanSpleen.csv', header = TRUE, sep = ",")
metadataWeizmanLiver <- read.csv('metadataWeizmanLiver.csv', header = TRUE, sep = ",")


allUniqueWeizmanLiver = distinct(countsWeizmanLiver, X, .keep_all = TRUE)
allUniqueWeizmanSpleen = distinct(countsWeizmanSpleen, X, .keep_all = TRUE)


matrixWeizmanLiver <- allUniqueWeizmanLiver[,-1]
rownames(matrixWeizmanLiver) <- allUniqueWeizmanLiver[,1]

matrixWeizmanSpleen <- allUniqueWeizmanSpleen[,-1]
rownames(matrixWeizmanSpleen) <- allUniqueWeizmanSpleen[,1]


ddsWeizmanLiver <- DESeqDataSetFromMatrix(countData = matrixWeizmanLiver, colData = metadataWeizmanLiver, design = ~ cellType)
ddsWeizmanSpleen <- DESeqDataSetFromMatrix(countData = matrixWeizmanSpleen, colData = metadataWeizmanSpleen, design = ~ cellType)

keepWeizmanLiver <- rowSums(counts(ddsWeizmanLiver)) >= 10
ddsWeizmanLiver <- ddsWeizmanLiver[keepWeizmanLiver,]

keepWeizmanSpleen <- rowSums(counts(ddsWeizmanSpleen)) >= 10
ddsWeizmanSpleen <- ddsWeizmanSpleen[keepWeizmanSpleen,]

ddsWeizmanLiver$cellType <- relevel(ddsWeizmanLiver$cellType, ref = "NK")
ddsWeizmanSpleen$cellType <- relevel(ddsWeizmanSpleen$cellType, ref = "NK")


ddsWeizmanLiver <- DESeq(ddsWeizmanLiver)
resWeizmanLiver <- results(ddsWeizmanLiver)

ddsWeizmanSpleen <- DESeq(ddsWeizmanSpleen)
resWeizmanSpleen <- results(ddsWeizmanSpleen)


#resWeizmanLiver0.05 <- results(ddsWeizmanLiver, alpha = 0.05)
#resWeizmanSpleen0.05 <- results(ddsWeizmanSpleen, alpha = 0.05)


#Weizman Liver

resWeizmanLiverSigTable <- subset(resWeizmanLiver, padj < 0.05)

#Weizman Spleen

resWeizmanSpleenSigTable <- subset(resWeizmanSpleen, padj < 0.05)


#write.csv(as.data.frame(resWeizmanLiverSigTable), file="WeizmanLiverSigLog2F.csv")
#write.csv(as.data.frame(resWeizmanSpleenSigTable), file="WeizmanSpleenSigLog2F.csv")

WeizmanLiverSigL2f = data.frame(rownames(resWeizmanLiverSigTable), resWeizmanLiverSigTable$log2FoldChange)
WeizmanSpleenSigL2f = data.frame(rownames(resWeizmanSpleenSigTable), resWeizmanSpleenSigTable$log2FoldChange)


colnames(log2fHalleLiver) <- c("gene", "HalleLiver_log2FoldChange")
colnames(WeizmanLiverSigL2f) <- c("gene", "WeizmanLiver_log2FoldChange")
colnames(WeizmanSpleenSigL2f) <- c("gene", "WeizmanSpleen_log2FoldChange")
colnames(log2fRobinetteLiver) <- c("gene", "RobinetteLiver_log2FoldChange")
colnames(log2fRobinetteSpleen) <- c("gene", "RobinetteSpleen_log2FoldChange")
colnames(log2fRobinetteSiLP) <- c("gene", "RobinetteSiLP_log2FoldChange")


merged_df_Liver <- join_all(list(log2fHalleLiver,WeizmanLiverSigL2f,log2fRobinetteLiver), by = 'gene', type = 'inner')
merged_df_Spleen <- join_all(list(log2fHalleLiver,WeizmanSpleenSigL2f,log2fRobinetteSpleen), by = 'gene', type = 'inner')
merged_df_SiLP <- join_all(list(log2fHalleLiver,log2fRobinetteSiLP), by = 'gene', type = 'inner')

write.csv(as.data.frame(merged_df_Liver), file="merged_df_liver_updated.csv")
write.csv(as.data.frame(merged_df_SiLP), file="merged_df_SiLP_updated.csv")
write.csv(as.data.frame(merged_df_Spleen), file="merged_df_spleen_updated.csv")

