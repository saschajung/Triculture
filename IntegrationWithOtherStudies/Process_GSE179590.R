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
library(stringr)

## Functions ##
#seur: A Seurat object
#pass: An identifier of the number of times this function is running. Only for plotting.
#org: The organism the data is from. Only relevant if not all metadata info is there.
filterCells <- function(seur,pass = -1,z = 3,org = "HUMAN",plot.path = "./"){
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

  max.mito.thr <- median(cell.QC.stat$percent.mito) + z*mad(cell.QC.stat$percent.mito)
  min.mito.thr <- median(cell.QC.stat$percent.mito) - z*mad(cell.QC.stat$percent.mito)

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
  min.features.thr <- median(log10(cell.QC.stat$nFeature_RNA)) - z*mad(log10(cell.QC.stat$nFeature_RNA))
  max.features.thr <- median(log10(cell.QC.stat$nFeature_RNA)) + z*mad(log10(cell.QC.stat$nFeature_RNA))

  # Set hight threshold on the number of transcripts
  max.count.thr <- median(log10(cell.QC.stat$nCount_RNA)) + z*mad(log10(cell.QC.stat$nCount_RNA))

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

  cm_human <- read.table("./Cell_marker_Human.txt",header = T, sep = "\t",stringsAsFactors = F, quote = "",comment.char = "")
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

paths <- list.dirs("./IntegrationWithOtherStudies/GSE179590/")
paths <- paths[-1]

for(path in paths){
  print(path)
  plot.path <- paste0(path,"/")
  require(stringr)
  samp <- str_extract(path,"S[0-9]+")
  if(file.exists(paste0(plot.path,samp,"_SeuratObj.rds"))){
    next
  }

  y1 <- Read10X(data.dir = path)
  y1 <- CreateSeuratObject(counts = y1, project = samp, min.cells = 3, min.features = 200)
  numCellsRecovered <- ncol(y1)
  y1

  # Percent of mitochondrial counts
  y1[["percent.mito"]] <- PercentageFeatureSet(y1, pattern = "^MT-")

  # Percent of mitochondrial ribosomal
  y1[["percent.ribo"]] <- PercentageFeatureSet(y1, pattern = "(^RPL|^RPS|^MRP)")

  p1 <- ggplot(y1@meta.data, aes(x=nCount_RNA, y=nFeature_RNA)) + geom_point() + geom_smooth(method="lm")
  p1 <- ggMarginal(p1, type = "histogram", fill="lightgrey")

  p2 <- ggplot(y1@meta.data, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) + geom_point() + geom_smooth(method="lm")
  p2 <- ggMarginal(p2, type = "histogram", fill="lightgrey")

  p <- plot_grid(plotlist = list(p1,p2), ncol=2, align='h', rel_widths = c(1, 1))
  ggsave(paste0(plot.path,"Count_Feature_Marginal_1stPass.png"),plot = p)

  #Filter cells based on percent.mito, percent.ribo and count/feature relationships.
  y1_filtered <- filterCells(y1,pass = 1, z= 10, org = "HUMAN",plot.path = plot.path)

  #Log-normalize data
  y1_filtered <- NormalizeData(y1_filtered)

  #Find variable features
  y1_filtered <- FindVariableFeatures(y1_filtered, selection.method = "vst", nfeatures = 3000)

  #Scale data
  y1_filtered <- ScaleData(y1_filtered, features = rownames(y1_filtered))

  #Run PCA
  y1_filtered <- RunPCA(y1_filtered,
                        features = VariableFeatures(object = y1_filtered),
                        npcs = 100)

  cellcyclegenes <- read.csv("./InputData/CellCycleGenes_Human.csv")
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  ensToSymb <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=cellcyclegenes$geneID,mart= mart)
  cellcyclegenes <- merge(cellcyclegenes,ensToSymb,by.x = 2, by.y = 1, all.x = T)
  g2m_genes <- cellcyclegenes$hgnc_symbol[cellcyclegenes$phase == "G2/M"]
  s_genes <- cellcyclegenes$hgnc_symbol[cellcyclegenes$phase == "S"]

  y1_filtered <- CellCycleScoring(y1_filtered,
                                  s.features = s_genes,
                                  g2m.features = g2m_genes,
                                  set.ident = F)

  #y1_filtered <- RunPCA(y1_filtered, features = c(s_genes, g2m_genes))
  ggsave(paste0(plot.path,"PCA_CellCycleGenes.png"),DimPlot(y1_filtered,group.by = "Phase"))
  #No dramatic confounding. Nevertheless, remove!

  y1_filtered$CC.Difference <- y1_filtered$S.Score - y1_filtered$G2M.Score

  y1_filtered <- SCTransform(y1_filtered, vars.to.regress = "CC.Difference",vst.flavor = "v2", verbose = T)

  y1_filtered <- RunPCA(y1_filtered,
                        npcs = 100,
                        verbose = T)

  ggsave(paste0(plot.path,"ElbowPlot.png"),ElbowPlot(y1_filtered,ndims = 100)) #Take 50

  numPCs <- 30

  #Remove outliers
  sweep.list <- paramSweep_v3(y1_filtered,
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

  nExp <- round(ncol(y1_filtered) * (unname(predict(lm_doublets,data.frame(numCellsRec = ncol(y1_filtered))))/100))
  pK <- as.numeric(levels(bcmvn$pK)[bcmvn$BCmetric == max(bcmvn$BCmetric)])
  y1_filtered <- doubletFinder_v3(seu = y1_filtered,
                                  PCs = 1:numPCs,
                                  pN = 0.25, #default
                                  pK = pK,
                                  nExp = nExp,
                                  reuse.pANN = FALSE,
                                  sct = TRUE,
                                  annotations = NULL)

  y1_filtered$Doublets <- "Singlet"
  y1_filtered$Doublets[y1_filtered[[paste0("pANN_0.25_",pK,"_",nExp)]] >= 0.5] <- "Doublet"
  y1_filtered$Doublets <- factor(y1_filtered$Doublets,levels = c("Doublet","Singlet"))

  p <- ggplot(y1_filtered@meta.data, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) +
    geom_point(aes(colour = Doublets), fill = "black",pch=21) +
    scale_color_manual(breaks = c("Singlet", "Doublet"),
                       values = c("black","firebrick1")) +
    geom_smooth(method="lm") +
    theme(legend.position="none") +
    annotate(geom = "text", label = paste0(as.numeric(table(y1_filtered@meta.data$Doublets)[2]), " Singlets\n",
                                           as.numeric(table(y1_filtered@meta.data$Doublets)[1]), " Doublets"), x = 4, y = 3.8)

  ggsave(paste0(plot.path,"Doublets.png"),plot = p)

  y1_filtered_nodups <- subset(y1_filtered, subset = Doublets == "Singlet") #There are no duplicates

  # Run UMAP and tSNE
  y1_filtered_nodups <- RunUMAP(y1_filtered_nodups,
                                dims = 1:numPCs,
                                n.neighbors = 20)

  y1_filtered_nodups <- RunTSNE(y1_filtered_nodups,
                                dims = 1:numPCs)

  cds <- as.cell_data_set(y1_filtered_nodups)
  cds <- cluster_cells(cds,k = 20,
                       reduction_method = "UMAP",
                       cluster_method = "leiden",
                       resolution = 0.001,
                       num_iter = 5)

  y1_filtered_nodups$seurat_clusters <- clusters(cds,reduction_method = "UMAP")
  Idents(y1_filtered_nodups) <- clusters(cds,reduction_method = "UMAP")

  #Let's see how it looks
  ggsave(paste0(plot.path,"UMAP.png"),DimPlot(y1_filtered_nodups, reduction = "umap"))
  ggsave(paste0(plot.path,"tSNE.png"),DimPlot(y1_filtered_nodups, reduction = "tsne"))

  #Annotate
  y1_filtered_nodups <- annotateCells(y1_filtered_nodups,plot.path = plot.path)

  #Densmap
  y1_filtered_nodups <- runDensUMAP(y1_filtered_nodups,pcs = 1:numPCs)
  ggsave(paste0(plot.path,"densMap.png"),DimPlot(y1_filtered_nodups, reduction = "densumap"))

  saveRDS(y1_filtered_nodups,paste0(plot.path,samp,"_SeuratObj.rds"))
}

#Integration
files <- list.files("./OtherStudies/GSE179590",pattern = ".rds",recursive = T,full.names = T)
obj.list <- list()
for(f in files){
  samp <- str_extract(f,"S[0-9]+")
  print(samp)
  obj.list[[samp]] <- readRDS(f)
}

plot.path <- "./OtherStudies/GSE179590/"

common_genes <- Reduce(intersect, lapply(obj.list, rownames))

features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 3000, assay = rep("SCT",length(obj.list)))
obj.list <- PrepSCTIntegration(object.list = obj.list,
                                  anchor.features = features,
                                  verbose = T)

anchors <- FindIntegrationAnchors(object.list = obj.list,
                                  normalization.method = "SCT",
                                  anchor.features = features,
                                  dims = 1:30,
                                  verbose = T,
                                  k.filter = NA,
                                  reduction = "rpca")

combined.sct <- IntegrateData(anchorset = anchors,
                              normalization.method = "SCT",
                              dims = 1:30)

rm(anchors)
rm(obj.list)

combined.sct <- RunPCA(combined.sct, npcs = 100, verbose = TRUE)

ggsave(paste0(plot.path,"ElbowPlot.png"),ElbowPlot(combined.sct,ndims = 100)) #Take 30

numPCs <- 30

DefaultAssay(combined.sct) <- "integrated"
combined.sct <- RunUMAP(combined.sct,
                        reduction = "pca",
                        dims = 1:numPCs,
                        n.neighbors = 20,
                        verbose = T)
combined.sct <- RunTSNE(combined.sct, reduction = "pca", dims = 1:numPCs, verbose = TRUE)

options(Matrix.warnDeprecatedCoerce = 1)
cds <- as.cell_data_set(combined.sct)
cds <- cluster_cells(cds,k = 10,
                     reduction_method = "UMAP",
                     cluster_method = "leiden",
                     resolution = 0.00005,
                     num_iter = 5,
                     verbose = F)

combined.sct$seurat_clusters <- clusters(cds,reduction_method = "UMAP")
Idents(combined.sct) <- clusters(cds,reduction_method = "UMAP")

ggsave(paste0(plot.path,"UMAP.png"),DimPlot(combined.sct, reduction = "umap"))
ggsave(paste0(plot.path,"tSNE.png"),DimPlot(combined.sct, reduction = "tsne"))

ggsave(paste0(plot.path,"UMAP_Split.png"),DimPlot(combined.sct, reduction = "umap", split.by = "orig.ident"))
ggsave(paste0(plot.path,"tSNE_Split.png"),DimPlot(combined.sct, reduction = "tsne", split.by = "orig.ident"))

ggsave(paste0(plot.path,"UMAP_Experiment.png"),DimPlot(combined.sct, reduction = "umap", group.by = "orig.ident"))
ggsave(paste0(plot.path,"tSNE_Experiment.png"),DimPlot(combined.sct, reduction = "tsne", group.by = "orig.ident"))

combined.sct <- runDensUMAP(combined.sct,pcs = 1:numPCs)
ggsave(paste0(plot.path,"densMap.png"),DimPlot(combined.sct, reduction = "densumap"))

#Annotate
combined.sct <- annotateCells(combined.sct,plot.path = plot.path)

#Read metadata
meta <-read.table("./OtherStudies/GSE179590/GSE179590_MS_brains_8samples_cell_lvl_metadata.csv", sep = ",", header = T, stringsAsFactors = F, quote = "",comment.char = "")
table(meta$disease.status,meta$sample.ID)
combined.sct$DiseaseStatus <- "SPMS"
combined.sct$DiseaseStatus[combined.sct$orig.ident %in% c("S1","S10","S5")] <- "RRMS"

meta <- meta[,c("sex","age","RIN","broad.label","cell_barcode")]

combined.sct$CellLabel <- paste0(combined.sct$orig.ident,"_",unname(gsub("-[0-9]+","",combined.sct$Barcodes)))

m <- merge(combined.sct@meta.data,meta,by.x = "CellLabel",by.y = "cell_barcode",all.x = T)
combined.sct@meta.data <- m
rownames(combined.sct@meta.data) <- names(combined.sct$CellLabel)

table(combined.sct$seurat_clusters,combined.sct$broad.label)

combined.sct$FinalLabel <- ""
combined.sct$FinalLabel[combined.sct$seurat_clusters == "1"] <- "Oligodendrocyte"
combined.sct$FinalLabel[combined.sct$seurat_clusters == "2"] <- "Oligodendrocyte"
combined.sct$FinalLabel[combined.sct$seurat_clusters == "3"] <- "Oligodendrocyte"
combined.sct$FinalLabel[combined.sct$seurat_clusters == "4"] <- "Excitatory Neuron"
combined.sct$FinalLabel[combined.sct$seurat_clusters == "5"] <- "Oligodendrocyte"
combined.sct$FinalLabel[combined.sct$seurat_clusters == "6"] <- "Oligodendrocyte"
combined.sct$FinalLabel[combined.sct$seurat_clusters == "7"] <- "Oligodendrocyte"
combined.sct$FinalLabel[combined.sct$seurat_clusters == "8"] <- "Oligodendrocyte"
combined.sct$FinalLabel[combined.sct$seurat_clusters == "9"] <- "Astrocyte"
combined.sct$FinalLabel[combined.sct$seurat_clusters == "10"] <- "Microglia"
combined.sct$FinalLabel[combined.sct$seurat_clusters == "11"] <- "Oligodendrocyte"
combined.sct$FinalLabel[combined.sct$seurat_clusters == "12"] <- "Oligodendrocyte"
combined.sct$FinalLabel[combined.sct$seurat_clusters == "13"] <- "Oligodendrocyte Progenitor"
combined.sct$FinalLabel[combined.sct$seurat_clusters == "14"] <- "Inhibitory Neuron"
combined.sct$FinalLabel[combined.sct$seurat_clusters == "15"] <- "Astrocyte"
combined.sct$FinalLabel[combined.sct$seurat_clusters == "16"] <- "Excitatory Neuron"
combined.sct$FinalLabel[combined.sct$seurat_clusters == "17"] <- "Excitatory Neuron"
combined.sct$FinalLabel[combined.sct$seurat_clusters == "18"] <- "Inhibitory Neuron"
combined.sct$FinalLabel[combined.sct$seurat_clusters == "19"] <- "Inhibitory Neuron"
combined.sct$FinalLabel[combined.sct$seurat_clusters == "20"] <- "Oligodendrocyte"
combined.sct$FinalLabel[combined.sct$seurat_clusters == "21"] <- "Endothelial cell"
combined.sct$FinalLabel[combined.sct$seurat_clusters == "22"] <- "Excitatory Neuron"
combined.sct$FinalLabel[combined.sct$seurat_clusters == "23"] <- "Inhibitory Neuron"
combined.sct$FinalLabel[combined.sct$seurat_clusters == "24"] <- "Excitatory Neuron"
combined.sct$FinalLabel[combined.sct$seurat_clusters == "25"] <- "Excitatory Neuron"
combined.sct$FinalLabel[combined.sct$seurat_clusters == "26"] <- "Inhibitory Neuron"
combined.sct$FinalLabel[combined.sct$seurat_clusters == "27"] <- "Excitatory Neuron"
combined.sct$FinalLabel[combined.sct$seurat_clusters == "28"] <- "Inhibitory Neuron"
combined.sct$FinalLabel[combined.sct$seurat_clusters == "29"] <- "Excitatory Neuron"
combined.sct$FinalLabel[combined.sct$seurat_clusters == "30"] <- "Oligodendrocyte Progenitor"
combined.sct$FinalLabel[combined.sct$seurat_clusters == "31"] <- "Excitatory Neuron"

ggsave(paste0(plot.path,"UMAP_FinalLabels.png"),DimPlot(combined.sct, reduction = "umap",group.by = "FinalLabel"))
ggsave(paste0(plot.path,"TSNE_FinalLabels.png"),DimPlot(combined.sct, reduction = "tsne",group.by = "FinalLabel"))

saveRDS(combined.sct,file = "./IntegrationWithOtherStudies/GSE179590/Integrated.rds")

# Integration with Triculture
scNuc <- readRDS("./sNuc-seq/sNuc_SeuratObj.rds")
obj <- readRDS("./IntegrationWithOtherStudies/GSE179590/Integrated.rds")
obj$CTs <- obj$FinalLabel
obj$CTs[obj$CTs == "Oligodendrocyte Progenitor"] <- "OPC"
DefaultAssay(obj) <- "SCT"
obj <- NormalizeData(obj)

obj.list <- SplitObject(obj, split.by = "orig.ident")
rm(obj)
gc()
gc()
gc()

obj.list[["snuc"]] <- scNuc
rm(scNuc)
gc()
gc()
gc()

obj.list <- lapply(obj.list,function(x){
  DefaultAssay(x) <- "RNA"
  x = DietSeurat(x,assays = "RNA")
  x = SCTransform(x, vars.to.regress = "CC.Difference",vst.flavor = "v2", verbose = T)
  DefaultAssay(x) <- "SCT"
  x = DietSeurat(x,assays = c("RNA","SCT"),scale.data = T)
  x = RunPCA(x, npcs = 30,verbose = FALSE)
  return(x)
})

gc()
gc()
gc()
gc()

features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 3000, assay = rep("SCT",length(obj.list)))
obj.list <- PrepSCTIntegration(object.list = obj.list,
                               anchor.features = features,
                               verbose = T)

anchors <- FindIntegrationAnchors(object.list = obj.list,
                                  normalization.method = "SCT",
                                  anchor.features = features,
                                  dims = 1:30,
                                  reduction = "rpca",
                                  verbose = T)
rm(obj.list)
gc()
gc()
gc()

combined.sct <- IntegrateData(anchorset = anchors,
                              normalization.method = "SCT",
                              dims = 1:30,
                              verbose = T)

rm(anchors)

combined.sct <- RunPCA(combined.sct, npcs = 100, verbose = TRUE)
plot.path <- "./IntegrationWithOtherStudies/GSE179590/"
dir.create(plot.path)
ggsave(paste0(plot.path,"ElbowPlot.png"),ElbowPlot(combined.sct,ndims = 100)) #Take 50

numPCs <- 50

DefaultAssay(combined.sct) <- "integrated"
combined.sct <- RunUMAP(combined.sct,
                        reduction = "pca",
                        dims = 1:numPCs,
                        verbose = T)
combined.sct <- RunTSNE(combined.sct, reduction = "pca", dims = 1:numPCs, verbose = TRUE)

combined.sct$CTs[combined.sct$orig.ident == "snuc"] <- unname(ct_map[combined.sct$seurat_clusters[combined.sct$orig.ident == "snuc"]])

combined.sct$Origin <- "Brain"
combined.sct$Origin[combined.sct$orig.ident == "snuc"] <- "This study"

combined.sct$CTOrigin <- paste0(combined.sct$CTs," (",combined.sct$Origin,")")

combined.sct$Clusts <- combined.sct$seurat_clusters
combined.sct$Clusts[combined.sct$orig.ident == "snuc"] <- paste0(combined.sct$Clusts[combined.sct$orig.ident == "snuc"]," (This study)")
combined.sct$Clusts[combined.sct$orig.ident != "snuc"] <- paste0(combined.sct$Clusts[combined.sct$orig.ident != "snuc"]," (Brain)")

save(combined.sct,file = "./IntegrationWithOtherStudies/GSE179590/sNucIntegration/Combined_SeuratObj.RData",compression_level = 9)


