###############################################################################
# This script processes the CellRanger output of the single-nuclei RNA-seq
# data and generates Figures S4D-F. Note: Label size has been adjusted in
# Biorender to increase readibility.
###############################################################################

library(Seurat)
library(DoubletFinder)
library(RColorBrewer)
library(ggplot2)
library(ggExtra)
library(cowplot)
library(dplyr)
library(cluster, quietly = TRUE)
library(ggstatsplot)
library(SCINA)
library(plyr)
library(biomaRt)
library(densvis)
library(monocle3)
library(singleseqgset)
library(pheatmap)
library(SeuratWrappers)
library(GSA)

## Functions ##
#seur: A Seurat object
#pass: An identifier of the number of times this function is running. Only for plotting.
#org: The organism the data is from. Only relevant if not all metadata info is there.
filterCells <- function(seur,pass = -1,org = "HUMAN",plot.path = "./"){
  seur@meta.data$Barcodes <- rownames(seur@meta.data)
  #Calculate percent.mito and percent.ribo metadata columns if they are not there
  if(!any(colnames(seur@meta.data) == "percent.mito")){
    if(org == "HUMAN"){
      seur[["percent.mito"]] <- PercentageFeatureSet(seur, pattern = "^MT-")
    } else if(org == "MOUSE"){
      seur[["percent.mito"]] <- PercentageFeatureSet(seur, pattern = "^mt-")
    } else {
      stop("The specified organism is not supported")
    }
  }
  if(!any(colnames(seur@meta.data) == "percent.ribo")){
    if(org == "HUMAN"){
      seur[["percent.ribo"]] <- PercentageFeatureSet(seur, pattern = "(^RPL|^RPS|^MRP)")
    } else if(org == "MOUSE"){
      seur[["percent.ribo"]] <- PercentageFeatureSet(seur, pattern = "(^Rpl|^Rps|^Mrp)")
    } else {
      stop("The specified organism is not supported")
    }
  }
  # Filtering cells based on percentage of mitochondrial transcripts
  cell.QC.stat <- seur@meta.data
  
  max.mito.thr <- median(cell.QC.stat$percent.mito) + 3*mad(cell.QC.stat$percent.mito)
  min.mito.thr <- median(cell.QC.stat$percent.mito) - 3*mad(cell.QC.stat$percent.mito)
  
  p1 <- ggplot(cell.QC.stat, aes(x=nFeature_RNA, y=percent.mito)) +
    geom_point() +
    geom_hline(aes(yintercept = max.mito.thr), colour = "red", linetype = 2) +
    geom_hline(aes(yintercept = min.mito.thr), colour = "red", linetype = 2) +
    annotate(geom = "text", label = paste0(as.numeric(table(cell.QC.stat$percent.mito > max.mito.thr | cell.QC.stat$percent.mito < min.mito.thr)[2])," cells removed\n",
                                           as.numeric(table(cell.QC.stat$percent.mito > max.mito.thr | cell.QC.stat$percent.mito < min.mito.thr)[1])," cells remain"), x = 6000, y = -10)
  
  p <- ggMarginal(p1, type = "histogram", fill="lightgrey", bins=100) 
  ggsave(paste0(plot.path,"Mitofilter_Marginal_Pass",pass,".png"),plot = p)
  
  cell.QC.stat <- cell.QC.stat %>% filter(percent.mito < max.mito.thr) %>% filter(percent.mito > min.mito.thr)
  
  # Filtering cells based on number of genes and transcripts detected
  # Set low and hight thresholds on the number of detected genes
  min.features.thr <- median(log10(cell.QC.stat$nFeature_RNA)) - 3*mad(log10(cell.QC.stat$nFeature_RNA))
  max.features.thr <- median(log10(cell.QC.stat$nFeature_RNA)) + 3*mad(log10(cell.QC.stat$nFeature_RNA))
  
  # Set hight threshold on the number of transcripts
  max.count.thr <- median(log10(cell.QC.stat$nCount_RNA)) + 3*mad(log10(cell.QC.stat$nCount_RNA))
  
  p1 <- ggplot(cell.QC.stat, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) +
    geom_point() +
    geom_smooth(method="lm") +
    geom_hline(aes(yintercept = min.features.thr), colour = "green", linetype = 2) +
    geom_hline(aes(yintercept = max.features.thr), colour = "green", linetype = 2) +
    geom_vline(aes(xintercept = max.count.thr), colour = "red", linetype = 2)
  
  p <- ggMarginal(p1, type = "histogram", fill="lightgrey")
  ggsave(paste0(plot.path,"FeatureAndCountFilter_1_Pass",pass,".png"),plot = p)
  
  # Filter cells base on both metrics
  cell.QC.stat <- cell.QC.stat %>% 
    filter(log10(nFeature_RNA) > min.features.thr) %>%
    filter(log10(nFeature_RNA) < max.features.thr) %>%
    filter(log10(nCount_RNA) < max.count.thr)
  
  lm.model <- lm(data = cell.QC.stat, formula = log10(nFeature_RNA) ~ log10(nCount_RNA))
  
  p2 <- ggplot(cell.QC.stat, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) +
    geom_point() +
    geom_smooth(method="lm") +
    geom_hline(aes(yintercept = min.features.thr), colour = "green", linetype = 2) +
    geom_hline(aes(yintercept = max.features.thr), colour = "green", linetype = 2) +
    geom_vline(aes(xintercept = max.count.thr), colour = "red", linetype = 2) +
    geom_abline(intercept = lm.model$coefficients[1] - 0.09 , slope = lm.model$coefficients[2], color="orange") +
    annotate(geom = "text", label = paste0(dim(cell.QC.stat)[1], " QC passed cells"), x = 4, y = 3.8)
  
  p <- ggMarginal(p2, type = "histogram", fill="lightgrey")
  ggsave(paste0(plot.path,"FeatureAndCountOutlier_2_Pass",pass,".png"),plot = p)
  
  # Cells to exclude lie below an intercept offset of -0.09
  cell.QC.stat$validCells <- log10(cell.QC.stat$nFeature_RNA) > (log10(cell.QC.stat$nCount_RNA) * lm.model$coefficients[2] + (lm.model$coefficients[1] - 0.09))
  
  p3 <- ggplot(cell.QC.stat, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) +
    geom_point(aes(colour = validCells), fill = "black",pch=21) +
    scale_color_manual(breaks = c("TRUE", "FALSE"), 
                       values = c("black","firebrick1")) +
    geom_smooth(method="lm") +
    geom_abline(intercept = lm.model$coefficients[1] - 0.09 , slope = lm.model$coefficients[2], color="orange") + 
    theme(legend.position="none") +
    annotate(geom = "text", label = paste0(as.numeric(table(cell.QC.stat$validCells)[2]), " QC passed cells\n",
                                           as.numeric(table(cell.QC.stat$validCells)[1]), " QC filtered"), x = 4, y = 3.8)
  
  p <- ggMarginal(p3, type = "histogram", fill="lightgrey")
  ggsave(paste0(plot.path,"FeatureAndCountOutlier_3_Pass",pass,".png"),plot = p)
  
  # Remove invalid cells
  cell.QC.stat <- cell.QC.stat %>% filter(validCells)
  
  seur <- subset(seur, subset = Barcodes %in% cell.QC.stat$Barcodes)
  return(seur)
}

findOptimalResolution <- function(seur,pcs,plot.path = "./"){
  emb <- Embeddings(object = seur@reductions$pca)[, pcs]
  resolutions <- seq(0.05,2,0.05)
  idx_list <- list()
  for(i in 1:length(resolutions)){
    print(paste0("Resolution: ",resolutions[i]))
    seur <- FindClusters(seur,
                         resolution = resolutions[i],
                         random.seed = 100,
                         method = "matrix",
                         algorithm = 2)
    part <- as.integer(levels(seur$seurat_clusters))[seur$seurat_clusters]
    idx_list[[i]] <- fpc::calinhara(emb,part)
    #intCriteria(traj = t(emb),
    #            part = part,
    #            crit = "Calinski_Harabasz")[[1]]
  }
  
  df <- data.frame(Resolution = resolutions,
                   Value = unlist(idx_list),
                   stringsAsFactors = F)
  
  p <- ggplot(data=df, aes(x=Resolution, y=Value)) +
    geom_line()+
    geom_point()+
    ylab("Calinski Harabasz Index")+
    xlab("Resolution")
  ggsave(paste0(plot.path,"Calinski_Harabasz_Index.png"),plot = p)
  return(seur)
}

annotateCells <- function(seur,plot.path = "./"){
  require(preprocessCore)
  inp_df <- as.matrix(seur@assays$RNA[,])
  inp_df <- log(inp_df+1)
  inp_df[] = normalize.quantiles(as.matrix(inp_df))
  
  cm_human <- read.table("./InputData/Cell_marker_Human.txt",header = T, sep = "\t",stringsAsFactors = F, quote = "",comment.char = "")
  cm_human <- cm_human[cm_human$tissue_class == "Brain",]
  cm_human <- cm_human[cm_human$cancer_type == "Normal",]
  cm_human <- cm_human[cm_human$cell_type != "Cancer cell",]
  cm_human <- cm_human[,c("cell_name","Symbol")]
  cm_human <- ddply(cm_human,.(cell_name,Symbol),nrow)
  cm_human <- cm_human[cm_human$Symbol != "",]
  cm_human <- cm_human[cm_human$V1 > 1,]
  cm_human <- cm_human[,c("cell_name","Symbol")]
  
  sigs <- lapply(unique(cm_human$cell_name),function(x){cm_human$Symbol[cm_human$cell_name == x]})
  names(sigs) <- unique(cm_human$cell_name)
  sigs <- lapply(sigs,function(x){intersect(x,rownames(inp_df))})
  sigs <- sigs[sapply(sigs,length) > 3]
  
  scina.results <- SCINA(inp_df,
                         sigs, 
                         max_iter = 100, 
                         convergence_n = 10, 
                         convergence_rate = 0.99, 
                         sensitivity_cutoff = 0.9, 
                         rm_overlap=F, 
                         allow_unknown=T
  )
  
  seur$SCINA_annot <- scina.results$cell_labels
  ggsave(paste0(plot.path,"UMAP_SCINA.png"),DimPlot(seur, reduction = "umap", group.by = "SCINA_annot"))
  return(seur)
}

runDensUMAP <- function(seur,pcs){
  require(densvis)
  emb <- Embeddings(object = seur@reductions$pca)[, pcs]
  fit <- densmap(emb,
                 dens_frac = 0.5,
                 dens_lambda = 0.5,
                 metric = "cosine",
                 n_epochs = 3000L,
                 random_state = 100L
  )
  colnames(fit) <- paste0("DensUMAP_", 1:2)
  rownames(fit) <- rownames(emb)
  seur[["densumap"]] <- CreateDimReducObject(embeddings = fit,
                                             key = "DensUMAP_",
                                             assay = DefaultAssay(seur)
  )
  return(seur)
}

findNumPCs <- function(seur){
  pct <- seur[["pca"]]@stdev / sum(seur[["pca"]]@stdev) * 100
  cumu <- cumsum(pct)
  
  co1 <- which(cumu > 90 & pct < 5)[1]
  
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.05), decreasing = T)[1] + 1
  
  return(min(co1,co2))
}

snuc <- Read10X(data.dir = "./sNuc-seq/CellRanger-Count-10X/")
snuc <- CreateSeuratObject(counts = snuc, project = "snuc", min.cells = 3, min.features = 200)
numCellsRecovered <- ncol(snuc)
snuc

plot.path = "./sNuc-seq/"

# Percent of mitochondrial counts
snuc[["percent.mito"]] <- PercentageFeatureSet(snuc, pattern = "^MT-")

# Percent of mitochondrial ribosomal
snuc[["percent.ribo"]] <- PercentageFeatureSet(snuc, pattern = "(^RPL|^RPS|^MRP)")

p1 <- ggplot(snuc@meta.data, aes(x=nCount_RNA, y=nFeature_RNA)) + geom_point() + geom_smooth(method="lm")
p1 <- ggMarginal(p1, type = "histogram", fill="lightgrey")

p2 <- ggplot(snuc@meta.data, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) + geom_point() + geom_smooth(method="lm")
p2 <- ggMarginal(p2, type = "histogram", fill="lightgrey")

p <- plot_grid(plotlist = list(p1,p2), ncol=2, align='h', rel_widths = c(1, 1))
ggsave(paste0(plot.path,"Count_Feature_Marginal_1stPass.png"),plot = p)

#Filter cells based on percent.mito, percent.ribo and count/feature relationships.
snuc_filtered <- filterCells(snuc,pass = 1, org = "HUMAN",plot.path = plot.path)

#Log-normalize data
snuc_filtered <- NormalizeData(snuc_filtered)

#Find variable features
snuc_filtered <- FindVariableFeatures(snuc_filtered, selection.method = "vst", nfeatures = 3000)

#Scale data
snuc_filtered <- ScaleData(snuc_filtered, features = rownames(snuc_filtered))

#Run PCA
snuc_filtered <- RunPCA(snuc_filtered, 
                        features = VariableFeatures(object = snuc_filtered),
                        npcs = 100)

cellcyclegenes <- read.csv("./InputData/CellCycleGenes_Human.csv")
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
ensToSymb <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=cellcyclegenes$geneID,mart= mart)
cellcyclegenes <- merge(cellcyclegenes,ensToSymb,by.x = 2, by.y = 1, all.x = T)
g2m_genes <- cellcyclegenes$hgnc_symbol[cellcyclegenes$phase == "G2/M"]
s_genes <- cellcyclegenes$hgnc_symbol[cellcyclegenes$phase == "S"]

snuc_filtered <- CellCycleScoring(snuc_filtered,
                                  s.features = s_genes,
                                  g2m.features = g2m_genes,
                                  set.ident = F)

ggsave(paste0(plot.path,"PCA_CellCycleGenes.png"),DimPlot(snuc_filtered,group.by = "Phase"))

snuc_filtered$CC.Difference <- snuc_filtered$S.Score - snuc_filtered$G2M.Score

snuc_filtered <- SCTransform(snuc_filtered, vars.to.regress = "CC.Difference",vst.flavor = "v2", verbose = T)

snuc_filtered <- RunPCA(snuc_filtered,
                        npcs = 100,
                        verbose = T)

ggsave(paste0(plot.path,"ElbowPlot.png"),ElbowPlot(snuc_filtered,ndims = 100)) #Take 50

numPCs <- findNumPCs(snuc_filtered)

#Remove outliers
sweep.list <- paramSweep_v3(snuc_filtered,
                            PCs = 1:numPCs,
                            sct = T)
sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

##Estimate expected percentage of doublets from 10X Genomics estimates from 3' Gene Expression v3.1 assay##
estDoublets <- c(0.4,0.8,1.6,2.4,3.2,4,4.8,5.6,6.4,7.2,8)
numCellsRec <- c(500,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000)
scatter.smooth(numCellsRec,estDoublets) #Looks linear
lm_doublets <- lm(estDoublets ~ numCellsRec)
summary(lm_doublets) #Perfect linear relationship r2 = 1

nExp <- round(ncol(snuc_filtered) * (unname(predict(lm_doublets,data.frame(numCellsRec = ncol(snuc_filtered))))/100))
pK <- as.numeric(levels(bcmvn$pK)[bcmvn$BCmetric == max(bcmvn$BCmetric)])
snuc_filtered <- doubletFinder_v3(seu = snuc_filtered,
                                  PCs = 1:numPCs,
                                  pN = 0.25, #default
                                  pK = pK,
                                  nExp = nExp,
                                  reuse.pANN = FALSE,
                                  sct = TRUE,
                                  annotations = NULL)

snuc_filtered$Doublets <- "Singlet"
snuc_filtered$Doublets[snuc_filtered[[paste0("pANN_0.25_",pK,"_",nExp)]] >= 0.5] <- "Doublet"
snuc_filtered$Doublets <- factor(snuc_filtered$Doublets,levels = c("Doublet","Singlet"))

p <- ggplot(snuc_filtered@meta.data, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) +
  geom_point(aes(colour = Doublets), fill = "black",pch=21) +
  scale_color_manual(breaks = c("Singlet", "Doublet"), 
                     values = c("black","firebrick1")) +
  geom_smooth(method="lm") +
  theme(legend.position="none") +
  annotate(geom = "text", label = paste0(as.numeric(table(snuc_filtered@meta.data$Doublets)[2]), " Singlets\n",
                                         as.numeric(table(snuc_filtered@meta.data$Doublets)[1]), " Doublets"), x = 4, y = 3.8)

ggsave(paste0(plot.path,"Doublets.png"),plot = p)

snuc_filtered_nodups <- subset(snuc_filtered, subset = Doublets == "Singlet") #There are no duplicates

# Run UMAP and tSNE
snuc_filtered_nodups <- RunUMAP(snuc_filtered_nodups,
                                dims = 1:numPCs,
                                n.neighbors = 20)

snuc_filtered_nodups <- RunTSNE(snuc_filtered_nodups,
                                dims = 1:numPCs)

#Let's see how it looks
ggsave(paste0(plot.path,"UMAP.png"),DimPlot(snuc_filtered_nodups, reduction = "umap"))
ggsave(paste0(plot.path,"tSNE.png"),DimPlot(snuc_filtered_nodups, reduction = "tsne"))

#Annotate
snuc_filtered_nodups <- annotateCells(snuc_filtered_nodups,plot.path = plot.path)
ggsave(paste0(plot.path,"UMAP_SCINA.png"),DimPlot(snuc_filtered_nodups, reduction = "umap",group.by = "SCINA_annot"))

#Densmap
snuc_filtered_nodups <- runDensUMAP(snuc_filtered_nodups,pcs = 1:numPCs)
ggsave(paste0(plot.path,"densMap.png"),DimPlot(snuc_filtered_nodups, reduction = "densumap"))

rm(snuc_filtered)
rm(snuc)

saveRDS(snuc_filtered_nodups,"./sNuc-seq/sNuc_SeuratObj.rds")

#Read scNuc data
scNuc <- readRDS("./sNuc-seq/sNuc_SeuratObj.rds")

cds <- as.cell_data_set(scNuc)
cds <- cluster_cells(cds,k = 20,
                     reduction_method = "UMAP",
                     cluster_method = "leiden",
                     resolution = 0.001,
                     num_iter = 5)

scNuc$seurat_clusters <- clusters(cds)
Idents(scNuc) <- clusters(cds)

subct_map <- c("1" = "Astrocyte",
            "2" = "Astrocyte",
            "3" = "OLIG1+GFAP+ Cells",
            "4" = "Astrocyte",
            "5" = "Astrocyte",
            "6" = "Astrocyte",
            "7" = "Astrocyte",
            "8" = "Astrocyte",
            "9" = "Radial Glia",
            "10" = "Astrocyte",
            "11" = "Astrocyte",
            "12" = "Astrocyte",
            "13" = "Radial Glia",
            "14" = "Astrocyte",
            "15" = "Microglia",
            "16" = "OPC",
            "17" = "Radial Glia",
            "18" = "Neuron",
            "19" = "Fibroblast-like")

ct_map <- c("1" = "Astrocyte",
            "2" = "Astrocyte",
            "3" = "OLIG1+GFAP+ Cells",
            "4" = "Astrocyte",
            "5" = "Astrocyte",
            "6" = "Astrocyte",
            "7" = "Astrocyte",
            "8" = "Astrocyte",
            "9" = "Radial Glia",
            "10" = "Astrocyte",
            "11" = "Astrocyte",
            "12" = "Astrocyte",
            "13" = "Radial Glia",
            "14" = "Astrocyte",
            "15" = "Microglia",
            "16" = "OPC",
            "17" = "Radial Glia",
            "18" = "Neuron",
            "19" = "Fibroblast-like")
scNuc$CTs <- unname(ct_map[scNuc$seurat_clusters])
scNuc$SubCTs <- unname(subct_map[scNuc$seurat_clusters])

saveRDS(scNuc,file = "./sNuc-seq/sNuc_SeuratObj.rds")


