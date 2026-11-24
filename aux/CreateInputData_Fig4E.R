library(Seurat)
library(ComplexHeatmap)
library(matrixStats)
library(pheatmap)
library(reshape2)

y1 <- readRDS("./scRNA-seq/Y1_final.rds")

aggExp <- AverageExpression(y1,assays = "SCT",group.by = "CTs")
aggExp$SCTscaled <- t(apply(aggExp$SCT,1,scale))
colnames(aggExp$SCTscaled) <- colnames(aggExp$SCT)
selectedGenes <- c()
for(ct in colnames(aggExp$SCT)){
  rmax <- rowMaxs(aggExp$SCTscaled[,which(colnames(aggExp$SCT) != ct)])
  diffs <- aggExp$SCTscaled[,ct] - rmax
  diffs <- diffs[which(aggExp$SCT[,ct] > 1 & diffs > 0)]
  diffs <- sort(diffs,decreasing = T)
  selectedGenes <- c(selectedGenes,names(diffs)[1:10])
}

mat <- as.matrix(aggExp$SCTscaled)
p <- pheatmap::pheatmap(mat[selectedGenes,],
                        clustering_method = "ward.D2",
                        silent = T)

order <- c("Neuron",
           "Fibroblast-like",
           "OPC",
           "OLIG1+GFAP+ Cells",
           "Radial Glia",
           "Astrocyte",
           "Proliferating Microglia",
           "Microglia",
           "Activated Microglia")
geneordering <- p$tree_row$labels[p$tree_row$order]

inp <- melt(mat[selectedGenes,])
inp$Var1 <- factor(inp$Var1,levels = geneordering)
inp$Var2 <- factor(inp$Var2,levels = order)

saveRDS(inp,"./InputData/InputData_Fig4E.rds")
