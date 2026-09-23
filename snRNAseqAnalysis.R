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
library(scCustomize)
library(viridis)
library(ggpubr)
library(ggsci)

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


setwd("./")

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

#DimPlot(scNuc,reduction = "umap",label =T)

#scNuc <- PrepSCTFindMarkers(scNuc)
#markers <- FindAllMarkers(scNuc,assay = "SCT")

#scNuc <- annotateCells(scNuc)

#1: Astrocyte
#2: Astrocyte
#3: Oligodendrocytes
#4: Astrocyte
#5: Astrocyte
#6: Astrocyte
#7: Astrocyte
#8: Astrocyte
#9: Radial Glia
#10: Astrocyte
#11: Astrocyte
#12: Astrocyte
#13: Astrocyte_1
#14: Astrocyte_2
#15: Microglia
#16: OPC
#17: Astrocyte_3
#18: Neuron
#19: Fibroblast

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
scNuc <- readRDS(file = "./sNuc-seq/sNuc_SeuratObj.rds")

#scNuc <- subset(scNuc, subset = CTs %in% c("Astrocyte","Microglia","Neuron","Radial Glia","OLIG1+GFAP+ Cells"))

plot_cor_full <- round(x = cor(x = scNuc@meta.data[,"nCount_RNA"], y = scNuc@meta.data[, "nFeature_RNA"]),digits = 2)
df <- scNuc@meta.data[,c("nCount_RNA","nFeature_RNA")]
df$nCount_RNA <- (df$nCount_RNA)
df$nFeature_RNA <- (df$nFeature_RNA)
p <- ggscatter(df, x = "nCount_RNA", y = "nFeature_RNA",
               color = "black", shape = 19, size = 1, # Points color, shape and size
               add = "reg.line",  # Add regressin line
               add.params = list(color = "blue"), # Customize reg. line
               conf.int = TRUE, # Add confidence interval
               cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
               cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n"),
               yscale = "log10",
               xscale = "log10",
               format.scale = T,
               xlab = "UMIs per Nucleus",
               ylab = "Genes per Nucleus",
               cor.coef.size = 6
) +
  font("xlab", size = 20, color = "black")+
  font("ylab", size = 20, color = "black")+
  font("xy.text", size = 16, color = "black")

ggsave("./FinalPlots_sNuc/GenesvsUMIs_rev.png",p)

cols <- c("Astrocyte" = "#CD534CFF",
          "Fibroblast-like" = "#868686FF",
          "Microglia" = "#0073C2FF",
          "Activated Microglia" = "#7AA6DCFF",
          "Proliferating Microglia" = "#003C67FF",
          "Neuron" = "#EFC000FF",
          "OLIG1+GFAP+ Cells" = "#8F7700FF",
          "OPC" = "#3B3B3BFF",
          "Radial Glia" = "#A73030FF")#pal_jco()(n = length(unique(scNuc$SubCTs)))
#names(cols) <- unique(scNuc$SubCTs)
p <- DimPlot_scCustom(seurat_object = scNuc,
                      repel = T,
                      group.by = "SubCTs",
                      figure_plot = TRUE,
                      label = F,
                      colors_use = cols,
                      prop_plot_percent = T,
                      add_prop_plot = T,
                      prop_plot_label = T,
                      label.size = 16)
p[[1]][[1]] <- p[[1]][[1]] + ggtitle("") + scale_color_manual(values = cols)
p[[2]][[1]] <- p[[2]][[1]] + scale_fill_manual(values = cols)

ggsave("./FinalPlots_sNuc/UMAP_rev.png",p,width = 12, height = 5)

opc_markers <- c("PDGFRA","OLIG1")
radialGlia_markers <- c("FABP7","SOX2","HES1","HES5")
neuron_markers <- c("RBFOX3","SNAP25","MAP2","STMN2","TUBB3")
fibroblast <- c("COL1A1","COL1A2","DCN","LUM","KRT18")
microglia <- c("P2RY12","CX3CR1","CSF1R","SALL1","CCL2","CCL3","CCL4","MKI67")
astro <- c("GFAP","AQP4","SLC1A3","ALDH1L1")
featuresToPlot <- c(opc_markers,radialGlia_markers,
                    neuron_markers,fibroblast,microglia,astro)

png(filename = "./FinalPlots_sNuc/Dotplot_Markers_rev.png",res = 300,width = 3000, height = 2400)
p <- Clustered_DotPlot(seurat_object = scNuc,
                 features = featuresToPlot,
                 group.by = "SubCTs",
                 assay = "SCT",
                 k = 6,
                 colors_use_idents = cols[c(1,7,9,3,8,6,2)],
                 row_label_size = 16,
                 column_label_size = 0,
                 legend_label_size = 16,
                 legend_title_size = 20)
dev.off()

Cytokines_All1<-c("IL7","IL1B","IL18")
Cytokines_All2<-c("EDA","LTB","TNF","TNFSF10","TNFSF12","TNFSF13","TNFSF13B","TNFSF15","FGF1", "FGF2", "FGF7", "FGF9", "FGF11","FGF12","FGF13","FGF14", "FGFBP3")
Cytokines_All3<-c("PDGFA","PDGFB","PDGFC","PDGFD","ANGPT1","ANGPT2","ANGPTL1","ANGPTL2","VEGFB","IGF1","SCG2", "CCL2", "CCL3","CCL4")
Cytokines_All4<-c("CXCL12","CXCL16", "CXCL17")
Cytokines_All5<-c("CSF1", "TGFB1","TGFB2","BMP7", "EGF","LRIG1","NRG1" ,"NRG3" ,"NRG4" ,"TMEFF1")
Cytokines_All <- unique(c(Cytokines_All1,Cytokines_All2,Cytokines_All3,Cytokines_All4,Cytokines_All5))

Cytokines_All <- setdiff(Cytokines_All,c("CXCL12","CXCL16")) #Those are not expressed in any cell in the SCT assay

png(filename = "./FinalPlots_sNuc/Dotplot_Cytokines_rev.png",res = 300,width = 2500, height = 3500)
Clustered_DotPlot(seurat_object = scNuc,
                  features = Cytokines_All,
                  group.by = "SubCTs",
                  assay = "SCT",
                  k = 3,
                  colors_use_idents = cols[c(1,7,9,3,8,6,2)],
                  row_label_size = 16,
                  column_label_size = 16,
                  legend_label_size = 16,
                  legend_title_size = 20)
dev.off()


#Plotting
# DotPlot(scNuc, features = c("PAX6","EOMES","NEUROD6","GAD1","RELN","SLC17A6","AQP4","OLIG1","CX3CR1","FN1","HBB"))
# DotPlot(scNuc, features = c("LINC00943","ADAMTS3","DOCK9","FZD8","NPR3"))
# DotPlot(scNuc, features = c("MAP2","OLIG1","DCX","STMN1","GFAP","VIM","PAX6","SLC1A3","FABP7","TNC","CDH2","NES","SOX2"))+theme(axis.text.x = element_text(angle = 90))
 DotPlot(scNuc, features = c("MAP2","OLIG1","DCX","STMN1","GFAP","AQP4","PDGFRA","NES","SOX2","NOTCH1","HES1","OCLN","CDH1","SOX10","PROM1","CD24","SOX1","SOX9","MKI67","SPARC","SNAP25"))+theme(axis.text.x = element_text(angle = 90))

####Plotting

#Integrate with GSE179590
obj <- readRDS("./OtherStudies/GSE179590/Integrated.rds")
obj$CTs <- obj$FinalLabel
obj$CTs[obj$CTs == "Oligodendrocyte Progenitor"] <- "OPC"
DefaultAssay(obj) <- "SCT"
obj <- NormalizeData(obj)
DefaultAssay(obj) <- "integrated"

Sys.setenv("R_MAX_MEM_SIZE"="50G")
anchors <- FindTransferAnchors(reference = obj, query = scNuc, normalization.method = "SCT", dims = 1:50)

predictions <- TransferData(anchorset = anchors, refdata = obj$CTs, dims = 1:50)
scNuc <- AddMetaData(scNuc, metadata = predictions)

astro_scores <- unname(scNuc$prediction.score.Astrocyte[scNuc$CTs == "Astrocyte"])
neuron_scores <- unname(rowMaxs(as.matrix(scNuc@meta.data[,c("prediction.score.Inhibitory.Neuron","prediction.score.Excitatory.Neuron")]))[scNuc$CTs == "Neuron"])
opc_scores <- unname(scNuc$prediction.score.OPC[scNuc$CTs == "OPC"])
microglia_scores <- unname(scNuc$prediction.score.Microglia[scNuc$CTs == "Microglia"])

preds_df <- data.frame(CT = c(rep("Astrocyte",length(astro_scores)),
                              rep("Neuron",length(neuron_scores)),
                              rep("OPC",length(opc_scores)),
                              rep("Microglia",length(microglia_scores))
                              ),
                       Scores = c(astro_scores,
                                  neuron_scores,
                                  opc_scores,
                                  microglia_scores),
                       stringsAsFactors = F)

saveRDS(preds_df,"./PlotsRevision/InputData_RidgeBars_GSE179590_vs_ours.rds")

predIds <- scNuc$predicted.id
predIds[which(predIds %in% c("Inhibitory Neuron","Excitatory Neuron"))] <- "Neuron"
inp_bar <- as.data.frame(prop.table(table(scNuc$CTs,predIds),margin = 1))
inp_bar <- inp_bar[which(inp_bar$Var1 %in% c("Astrocyte","Neuron","OPC","Microglia")),]
inp_bar$Var1 <- factor(inp_bar$Var1, levels = rev(levels(factor(inp_bar$Var1))))

saveRDS(inp_bar,"./PlotsRevision/InputData_inp_bar_RidgeBars_GSE179590_vs_ours.rds")
saveRDS(predIds,"./PlotsRevision/InputData_predIds_RidgeBars_GSE179590_vs_ours.rds")

commonGenes <- intersect(rownames(obj@assays$SCT),rownames(scNuc@assays$SCT))

pseudo_ref <- PseudobulkExpression(object = obj,
                                   assay = "SCT",
                                   features = commonGenes,
                                   group.by = "CTs",
                                   method = "aggregate",
                                   layer = "counts")
pseudo_ref$SCT <- apply(pseudo_ref$SCT,2,function(x){(x*1e6)/sum(x)})
cnames <- colnames(pseudo_ref$SCT)
pseudo_ref$SCT <- t(apply(pseudo_ref$SCT,1,function(x){as.numeric(scale(x))}))
colnames(pseudo_ref$SCT) <- cnames

cor(pseudo_ref$SCT,use = "complete.obs")

pseudo_query <- PseudobulkExpression(object = scNuc,
                                   assay = "SCT",
                                   features = commonGenes,
                                   group.by = "CTs",
                                   method = "aggregate",
                                   layer = "counts")
pseudo_query$SCT <- apply(pseudo_query$SCT,2,function(x){(x*1e6)/sum(x)})
cnames <- colnames(pseudo_query$SCT)
pseudo_query$SCT <- t(apply(pseudo_query$SCT,1,function(x){as.numeric(scale(x))}))
colnames(pseudo_query$SCT) <- cnames

cormat <- cor(as.matrix(pseudo_query[["SCT"]]),as.matrix(pseudo_ref[["SCT"]]),method = "pearson",use = "complete.obs")

saveRDS(cormat, file = "./PlotsRevision/InputData_Heatmap_GSE179590_vs_ours.rds")

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
plot.path <- "./OtherStudies/GSE179590/sNucIntegration/"
dir.create(plot.path)
ggsave(paste0(plot.path,"ElbowPlot.png"),ElbowPlot(combined.sct,ndims = 100)) #Take 50

numPCs <- 50

DefaultAssay(combined.sct) <- "integrated"
combined.sct <- RunUMAP(combined.sct,
                        reduction = "pca",
                        dims = 1:numPCs,
                        #n.neighbors = 20,
                        verbose = T)
combined.sct <- RunTSNE(combined.sct, reduction = "pca", dims = 1:numPCs, verbose = TRUE)

combined.sct$CTs[combined.sct$orig.ident == "snuc"] <- unname(ct_map[combined.sct$seurat_clusters[combined.sct$orig.ident == "snuc"]])

combined.sct$Origin <- "Brain"
combined.sct$Origin[combined.sct$orig.ident == "snuc"] <- "This study"

combined.sct$CTOrigin <- paste0(combined.sct$CTs," (",combined.sct$Origin,")")

combined.sct$Clusts <- combined.sct$seurat_clusters
combined.sct$Clusts[combined.sct$orig.ident == "snuc"] <- paste0(combined.sct$Clusts[combined.sct$orig.ident == "snuc"]," (This study)")
combined.sct$Clusts[combined.sct$orig.ident != "snuc"] <- paste0(combined.sct$Clusts[combined.sct$orig.ident != "snuc"]," (Brain)")

save(combined.sct,file = "./OtherStudies/GSE179590/sNucIntegration/Combined_SeuratObj.RData",compression_level = 9)
load("./OtherStudies/GSE179590/sNucIntegration/Combined_SeuratObj.RData")

p <- DimPlot_scCustom(seurat_object = combined.sct,
                      repel = T,
                      group.by = "CTs",
                      figure_plot = TRUE,
                      label = T,
                      label.color = "white",
                      label.box	= T,
                      colors_use = c(pal_jco()(10),"orange"))
p[[1]] <- p[[1]] + ggtitle("") + theme(legend.position = "none")

ggsave(paste0(plot.path,"UMAP.png"),p)

p <- DimPlot_scCustom(seurat_object = combined.sct,
                      repel = T,
                      group.by = "Origin",
                      figure_plot = TRUE,
                      label = F,
                      label.color = "white",
                      label.box	= T,
                      colors_use = c(pal_jco()(10),"orange"))
p[[1]] <- p[[1]] + ggtitle("")

ggsave(paste0(plot.path,"UMAP_origin.png"),p)

cols <- c("Astrocyte (This study)" = "#CD534CFF",
          "Fibroblast-like (This study)" = "darkgrey",
          "Microglia (This study)" = "#0073C2FF",
          "Activated Microglia (This study)" = "darkgrey",
          "Proliferating Microglia (This study)" = "darkgrey",
          "Neuron (This study)" = "#EFC000FF",
          "OLIG1+GFAP+ Cells (This study)" = "darkgrey",
          "OPC (This study)" = "#3B3B3BFF",
          "Radial Glia (This study)" = "darkgrey",
          "Astrocyte (Brain)" = lighten("#CD534CFF",amount = 0.5),
          "Endothelial cell (Brain)" = "darkgrey",
          "Microglia (Brain)" = lighten("#0073C2FF",amount = 0.5),
          "Activated Microglia (Brain)" = "darkgrey",
          "Proliferating Microglia (Brain)" = "darkgrey",
          "Inhibitory Neuron (Brain)" = lighten("#EFC000FF",amount = 0.5),
          "Excitatory Neuron (Brain)" = darken("#EFC000FF",amount = 0.2),
          "OLIG1+GFAP+ Cells (Brain)" = "darkgrey",
          "OPC (Brain)" = lighten("#3B3B3BFF",amount = 0.2),
          "Radial Glia (Brain)" = "darkgrey",
          "Oligodendrocyte (Brain)" = "darkgrey")
p <- DimPlot_scCustom(seurat_object = combined.sct,
                      repel = T,
                      group.by = "CTOrigin",
                      figure_plot = TRUE,
                      label = T,
                      label.color = "white",
                      label.box	= T,
                      colors_use = cols,
                      label.size = 16/.pt)
p[[1]] <- p[[1]] + ggtitle("") + theme(legend.position = "none")
p
ggsave("./PlotsRevision/Fig4H_UMAP_CTsOrigin_rev.png",p,width = 10, height = 8)

#ggsave(paste0(plot.path,"UMAP_clusters.png"),DimPlot(combined.sct, reduction = "umap", group.by = "Clusts"))

sub <- subset(combined.sct, subset = CTOrigin %in% c("Microglia (This study)","Microglia (Brain)"))
DefaultAssay(sub) <- "SCT"
sub$replicate <- sub$orig.ident
set.seed(1234)
sub$replicate[which(sub$orig.ident == "snuc")] <- sample(paste0("sNuc",1:2),length(which(sub$orig.ident == "snuc")),replace = T)
microglia.diff <- Libra::run_de(input = sub,
                                meta = sub@meta.data,
                                cell_type_col = "CTs",
                                replicate_col = "replicate",
                                label_col = "Origin",
                                de_family = "pseudobulk",
                                de_method = "edgeR",
                                de_type = "LRT") #Negative FC means higher in this study

sub <- subset(combined.sct, subset = CTOrigin %in% c("Astrocyte (This study)","Astrocyte (Brain)"))
DefaultAssay(sub) <- "SCT"
sub$replicate <- sub$orig.ident
set.seed(1234)
sub$replicate[which(sub$orig.ident == "snuc")] <- sample(paste0("sNuc",1:2),length(which(sub$orig.ident == "snuc")),replace = T)
astrocyte.diff <- Libra::run_de(input = sub,
                                meta = sub@meta.data,
                                cell_type_col = "CTs",
                                replicate_col = "replicate",
                                label_col = "Origin",
                                de_family = "pseudobulk",
                                de_method = "edgeR",
                                de_type = "LRT") #Negative FC means higher in this study

sub <- subset(combined.sct, subset = CTOrigin %in% c("Neuron (This study)","Excitatory Neuron (Brain)"))
sub$CTs <- "Neuron"
DefaultAssay(sub) <- "SCT"
sub$replicate <- sub$orig.ident
set.seed(1234)
sub$replicate[which(sub$orig.ident == "snuc")] <- sample(paste0("sNuc",1:2),length(which(sub$orig.ident == "snuc")),replace = T)
excitNeuron.diff <- Libra::run_de(input = sub,
                                meta = sub@meta.data,
                                cell_type_col = "CTs",
                                replicate_col = "replicate",
                                label_col = "Origin",
                                de_family = "pseudobulk",
                                de_method = "edgeR",
                                de_type = "LRT") #Negative FC means higher in this study

sub <- subset(combined.sct, subset = CTOrigin %in% c("Neuron (This study)","Inhibitory Neuron (Brain)"))
sub$CTs <- "Neuron"
DefaultAssay(sub) <- "SCT"
sub$replicate <- sub$orig.ident
set.seed(1234)
sub$replicate[which(sub$orig.ident == "snuc")] <- sample(paste0("sNuc",1:2),length(which(sub$orig.ident == "snuc")),replace = T)
inhibNeuron.diff <- Libra::run_de(input = sub,
                                  meta = sub@meta.data,
                                  cell_type_col = "CTs",
                                  replicate_col = "replicate",
                                  label_col = "Origin",
                                  de_family = "pseudobulk",
                                  de_method = "edgeR",
                                  de_type = "LRT") #Negative FC means higher in this study

sub <- subset(combined.sct, subset = CTOrigin %in% c("OPC (This study)","OPC (Brain)"))
sub$CTs <- "OPC"
DefaultAssay(sub) <- "SCT"
sub$replicate <- sub$orig.ident
set.seed(1234)
sub$replicate[which(sub$orig.ident == "snuc")] <- sample(paste0("sNuc",1:2),length(which(sub$orig.ident == "snuc")),replace = T)
opc.diff <- Libra::run_de(input = sub,
                                  meta = sub@meta.data,
                                  cell_type_col = "CTs",
                                  replicate_col = "replicate",
                                  label_col = "Origin",
                                  de_family = "pseudobulk",
                                  de_method = "edgeR",
                                  de_type = "LRT") #Negative FC means higher in this study

save(astrocyte.diff,microglia.diff,inhibNeuron.diff,excitNeuron.diff,opc.diff,file = "./OtherStudies/GSE179590/sNucIntegration/DifferentialExpressionResults.RData",compression_level = 6)
load("./OtherStudies/GSE179590/sNucIntegration/DifferentialExpressionResults.RData")
library(WebGestaltR)
inp_microglia_list <- list("This study" = microglia.diff$gene[microglia.diff$p_val_adj <= 0.01 & microglia.diff$avg_logFC < -0.5],
                           "Brain" = microglia.diff$gene[microglia.diff$p_val_adj <= 0.01 & microglia.diff$avg_logFC > 0.5])

micBrain <- WebGestaltR(enrichMethod = "ORA",
            organism = "hsapiens",
            enrichDatabase = "geneontology_Biological_Process",
            interestGene = inp_microglia_list$Brain,
            interestGeneType = "genesymbol",
            referenceGene = rownames(combined.sct@assays$RNA),
            referenceGeneType = "genesymbol",
            sigMethod = "fdr",
            fdrMethod = "BH",
            fdrThr = 0.01,
            isOutput = F
            )

micThis <- WebGestaltR(enrichMethod = "ORA",
                        organism = "hsapiens",
                        enrichDatabase = "geneontology_Biological_Process",
                        interestGene = inp_microglia_list$`This study`,
                        interestGeneType = "genesymbol",
                        referenceGene = rownames(combined.sct@assays$RNA),
                        referenceGeneType = "genesymbol",
                        sigMethod = "fdr",
                        fdrMethod = "BH",
                        fdrThr = 0.01,
                        isOutput = F
)

#Excitatory Neuron
inp_excitneuron_list <- list("This study" = excitNeuron.diff$gene[excitNeuron.diff$p_val_adj <= 0.01 & excitNeuron.diff$avg_logFC < -0.5],
                         "Brain" = excitNeuron.diff$gene[excitNeuron.diff$p_val_adj <= 0.01 & excitNeuron.diff$avg_logFC > 0.5])

excitNeuBrain <- WebGestaltR(enrichMethod = "ORA",
                        organism = "hsapiens",
                        enrichDatabase = "geneontology_Biological_Process",
                        interestGene = inp_excitneuron_list$Brain,
                        interestGeneType = "genesymbol",
                        referenceGene = rownames(combined.sct@assays$RNA),
                        referenceGeneType = "genesymbol",
                        sigMethod = "fdr",
                        fdrMethod = "BH",
                        fdrThr = 0.01,
                        isOutput = F
)

excitNeuThis <- WebGestaltR(enrichMethod = "ORA",
                       organism = "hsapiens",
                       enrichDatabase = "geneontology_Biological_Process",
                       interestGene = inp_excitneuron_list$`This study`,
                       interestGeneType = "genesymbol",
                       referenceGene = rownames(combined.sct@assays$RNA),
                       referenceGeneType = "genesymbol",
                       sigMethod = "fdr",
                       fdrMethod = "BH",
                       fdrThr = 0.01,
                       isOutput = F
)

#Inhibitory Neuron
inp_inhibneuron_list <- list("This study" = inhibNeuron.diff$gene[inhibNeuron.diff$p_val_adj <= 0.01 & inhibNeuron.diff$avg_logFC < -0.5],
                         "Brain" = inhibNeuron.diff$gene[inhibNeuron.diff$p_val_adj <= 0.01 & inhibNeuron.diff$avg_logFC > 0.5])

inhibNeuBrain <- WebGestaltR(enrichMethod = "ORA",
                         organism = "hsapiens",
                         enrichDatabase = "geneontology_Biological_Process",
                         interestGene = inp_inhibneuron_list$Brain,
                         interestGeneType = "genesymbol",
                         referenceGene = rownames(combined.sct@assays$RNA),
                         referenceGeneType = "genesymbol",
                         sigMethod = "fdr",
                         fdrMethod = "BH",
                         fdrThr = 0.01,
                         isOutput = F
)

inhibNeuThis <- WebGestaltR(enrichMethod = "ORA",
                        organism = "hsapiens",
                        enrichDatabase = "geneontology_Biological_Process",
                        interestGene = inp_inhibneuron_list$`This study`,
                        interestGeneType = "genesymbol",
                        referenceGene = rownames(combined.sct@assays$RNA),
                        referenceGeneType = "genesymbol",
                        sigMethod = "fdr",
                        fdrMethod = "BH",
                        fdrThr = 0.01,
                        isOutput = F
)

#Astrocyte
inp_astro_list <- list("This study" = astrocyte.diff$gene[astrocyte.diff$p_val_adj <= 0.01 & astrocyte.diff$avg_logFC < -0.5],
                       "Brain" = astrocyte.diff$gene[astrocyte.diff$p_val_adj <= 0.01 & astrocyte.diff$avg_logFC > 0.5])

astroBrain <- WebGestaltR(enrichMethod = "ORA",
                         organism = "hsapiens",
                         enrichDatabase = "geneontology_Biological_Process",
                         interestGene = inp_astro_list$Brain,
                         interestGeneType = "genesymbol",
                         referenceGene = rownames(combined.sct@assays$RNA),
                         referenceGeneType = "genesymbol",
                         sigMethod = "fdr",
                         fdrMethod = "BH",
                         fdrThr = 0.01,
                         isOutput = F
)

astroThis <- WebGestaltR(enrichMethod = "ORA",
                        organism = "hsapiens",
                        enrichDatabase = "geneontology_Biological_Process",
                        interestGene = inp_astro_list$`This study`,
                        interestGeneType = "genesymbol",
                        referenceGene = rownames(combined.sct@assays$RNA),
                        referenceGeneType = "genesymbol",
                        sigMethod = "fdr",
                        fdrMethod = "BH",
                        fdrThr = 0.01,
                        isOutput = F
)

#OPC
inp_opc_list <- list("This study" = opc.diff$gene[opc.diff$p_val_adj <= 0.01 & opc.diff$avg_logFC < -0.5],
                       "Brain" = opc.diff$gene[opc.diff$p_val_adj <= 0.01 & opc.diff$avg_logFC > 0.5])

opcBrain <- WebGestaltR(enrichMethod = "ORA",
                          organism = "hsapiens",
                          enrichDatabase = "geneontology_Biological_Process",
                          interestGene = inp_opc_list$Brain,
                          interestGeneType = "genesymbol",
                          referenceGene = rownames(combined.sct@assays$RNA),
                          referenceGeneType = "genesymbol",
                          sigMethod = "fdr",
                          fdrMethod = "BH",
                          fdrThr = 0.01,
                          isOutput = F
)

opcThis <- WebGestaltR(enrichMethod = "ORA",
                         organism = "hsapiens",
                         enrichDatabase = "geneontology_Biological_Process",
                         interestGene = inp_opc_list$`This study`,
                         interestGeneType = "genesymbol",
                         referenceGene = rownames(combined.sct@assays$RNA),
                         referenceGeneType = "genesymbol",
                         sigMethod = "fdr",
                         fdrMethod = "BH",
                         fdrThr = 0.01,
                         isOutput = F
)

save(astroBrain,astroThis,micBrain,micThis,excitNeuBrain,excitNeuThis,inhibNeuBrain,inhibNeuThis,opcThis,opcBrain,file = "./OtherStudies/GSE179590/sNucIntegration/GODifferences.RData")

load("./OtherStudies/GSE179590/sNucIntegration/DifferentialExpressionResults.RData",verbose = T)
write.table(microglia.diff[which(abs(microglia.diff$avg_logFC) > 0.5 & microglia.diff$p_val_adj <= 0.01),],"./OtherStudies/GSE179590/sNucIntegration/DiffExp_Microglia.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(astrocyte.diff[which(abs(astrocyte.diff$avg_logFC) > 0.5 & astrocyte.diff$p_val_adj <= 0.01),],"./OtherStudies/GSE179590/sNucIntegration/DiffExp_Astrocyte.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(inhibNeuron.diff[which(abs(inhibNeuron.diff$avg_logFC) > 0.5 & inhibNeuron.diff$p_val_adj <= 0.01),],"./OtherStudies/GSE179590/sNucIntegration/DiffExp_InhibNeuron.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(excitNeuron.diff[which(abs(excitNeuron.diff$avg_logFC) > 0.5 & excitNeuron.diff$p_val_adj <= 0.01),],"./OtherStudies/GSE179590/sNucIntegration/DiffExp_ExcitNeuron.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(opc.diff[which(abs(opc.diff$avg_logFC) > 0.5 & opc.diff$p_val_adj <= 0.01),],"./OtherStudies/GSE179590/sNucIntegration/DiffExp_OPC.txt",sep = "\t",row.names = F, col.names = T, quote = F)

load("./OtherStudies/GSE179590/sNucIntegration/GODifferences.RData",verbose = T)
write.table(astroBrain,file = "./OtherStudies/GSE179590/sNucIntegration/GO_AstrocyteBrain.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(astroThis,file = "./OtherStudies/GSE179590/sNucIntegration/GO_AstrocyteThis.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(micBrain,file = "./OtherStudies/GSE179590/sNucIntegration/GO_MicrogliaBrain.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(micThis,file = "./OtherStudies/GSE179590/sNucIntegration/GO_MicrogliaThis.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(excitNeuBrain,file = "./OtherStudies/GSE179590/sNucIntegration/GO_ExcitNeuBrain.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(excitNeuThis,file = "./OtherStudies/GSE179590/sNucIntegration/GO_ExcitNeuThis.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(inhibNeuBrain,file = "./OtherStudies/GSE179590/sNucIntegration/GO_InhibNeuBrain.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(inhibNeuThis,file = "./OtherStudies/GSE179590/sNucIntegration/GO_InhibNeuThis.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(opcBrain,file = "./OtherStudies/GSE179590/sNucIntegration/GO_OPCBrain.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(opcThis,file = "./OtherStudies/GSE179590/sNucIntegration/GO_OPCThis.txt",sep = "\t",row.names = F, col.names = T, quote = F)


inp.list <- list("astroThis" = astroThis,
                 "astroBrain" = astroBrain,
                 "excitNeuBrain" = excitNeuBrain,
                 "excitNeuThis" = excitNeuThis,
                 "inhibNeuBrain" = inhibNeuBrain,
                 "inhibNeuThis" = inhibNeuThis,
                 "microgliaThis" = micThis,
                 "microgliaBrain" = micBrain,
                 "opcThis" = opcThis,
                 "opcBrain" = opcBrain)

title_map <- c("astroThis" = "Astrocyte (This study)",
               "astroBrain" = "Astrocyte (Brain)",
               "excitNeuBrain" = "Excitatory Neuron (Brain)",
               "excitNeuThis" = "Neuron (This study)",
               "inhibNeuBrain" = "Inhibitory Neuron (Brain)",
               "inhibNeuThis" = "Neuron (This study)",
               "microgliaThis" = "Microglia (This study)",
               "microgliaBrain" = "Microglia (Brain)",
               "opcThis" = "OPC (This study)",
               "opcBrain" = "OPC (Brain)")

plot.path <- "./OtherStudies/GSE179590/sNucIntegration/"
#Plot top 20

enrichment_df <- do.call("rbind.data.frame",lapply(names(inp.list),function(n){
  if(is.null(inp.list)){
    return(data.frame("GO" = character(0),
                      "FDR" = numeric(0),
                      "EnrRatio" = numeric(0),
                      "CellType" = character(0),
                      stringsAsFactors = F))
  }
  goterms <- inp.list[[n]]$description
  fdr <- inp.list[[n]]$FDR
  enrRatio <- inp.list[[n]]$enrichmentRatio
  df <- data.frame("GO" = goterms,
             "FDR" = fdr,
             "EnrRatio" = enrRatio,
             "CellType" = rep(unname(title_map[n]),length(goterms)),
             stringsAsFactors = F)
  if(nrow(df) > 10){
    df <- df[1:10,]
  }
  return(df)
}))

enrichment_df$GO <- factor(enrichment_df$GO,levels = rev(enrichment_df$GO[!duplicated(enrichment_df$GO)]))
enrichment_df$CellType <- factor(enrichment_df$CellType,levels = c("Astrocyte (Brain)",
                                                                   "Astrocyte (This study)",
                                                                   "Excitatory Neuron (Brain)",
                                                                   "Inhibitory Neuron (Brain)",
                                                                   "Neuron (This study)",
                                                                   "Microglia (Brain)",
                                                                   "Microglia (This study)",
                                                                   "OPC (Brain)",
                                                                   "OPC (This study)"))   
library(viridis)
p <- ggballoonplot(enrichment_df,
              size = "EnrRatio",
              fill = "FDR",
              x = "CellType",
              y = "GO") +
  scale_fill_viridis(option = "C", direction = -1, begin = 0.15, end = 0.85) +
  ylab("") +
  xlab("") + 
  font("legend.title",size = 20) +
  font("legend.text",size = 16) +
  font("xy",size = 20) +
  font("xy.text",size = 16)
p$labels$size = "Enrichment Ratio"

ggsave(paste0(plot.path,"GOEnrichment_CTsOrigin.png"),p,width = 15, height = 15)

avgExp <- AverageExpression(combined.sct,assay = c("RNA","SCT"),group.by = "CTOrigin")
avgExp$SCT <- avgExp$SCT[rowSds(avgExp$SCT) > 0,]
avgExp$RNA <- avgExp$RNA[rowSds(avgExp$RNA) > 0,]

avgExp$SCT <- avgExp$SCT[,c("Astrocyte (Brain)", "Astrocyte (This study)",
                            "Excitatory Neuron (Brain)","Inhibitory Neuron (Brain)",
                            "Microglia (Brain)","Microglia (This study)","Neuron (This study)")]
avgExp$SCT <- avgExp$SCT[rowSds(avgExp$SCT) > 0,]

library(pheatmap)
pheatmap(avgExp$SCT,
         scale = "row",
         show_rownames = F,
         clustering_method = "ward.D2",
         clustering_distance_cols = "correlation")
         filename = "./Documents/Microglia_Stanford/OtherStudies/GSE179590/sNucIntegration/AvgExpression.png")

##########
library(InterCom)
rm(combined.sct)
scNuc <- readRDS("./sNuc-seq/sNuc_SeuratObj.rds")
inpmat <- as.matrix(scNuc@assays$RNA@counts)
anno.tbl <- data.frame(cell.name = rownames(scNuc@meta.data),
                       cell.type = unname(scNuc$CTs),
                       stringsAsFactors = F)
anno.tbl$cell.type <- gsub(" ","\\.",anno.tbl$cell.type)
anno.tbl$cell.type <- gsub("-","\\.",anno.tbl$cell.type)
anno.tbl$cell.type <- gsub("\\+","\\.",anno.tbl$cell.type)
anno.tbl$cell.name <- gsub("-","\\.",anno.tbl$cell.name)
colnames(inpmat) <- gsub("-","\\.",colnames(inpmat))
inpmat[1:5,1:5]
InterCom(data=inpmat,
         anno.tbl = anno.tbl,
         species = "HUMAN",
         ncores = 4,
         sig.cutoff = 0,
         z.score.cutoff = 0,
         tissue.name = "forPaper_RNA",
         out.path = "./sNuc-seq/Intercom_forPaper_RNA/")


load("./sNuc-seq/Intercom_forPaper_RNA/output_forPaper_RNA.RData",verbose = T)
ints <- output$final
ints <- ints[ints$score >= 0.9,]
pop_map <- c("OLIG1.GFAP..Cells" = "OLIG1+GFAP+ Cells",
             "Intermediate.Astrocyte.Progenitor" = "Intermediate Astrocyte Progenitor",
             "Radial.Glia" = "Radial Glia",
             "Fibroblast.like" = "Fibroblast-like",
             "Astrocyte" = "Astrocyte",
             "Neuron" = "Neuron",
             "OPC" = "OPC",
             "Microglia" = "Microglia")
ints$Lig.pop <- pop_map[ints$Lig.pop]
ints$Rec.pop <- pop_map[ints$Rec.pop]

write.table(ints,file = "./sNuc-seq/Intercom_forPaper_RNA/FinalInteractions.txt",sep = "\t",row.names = F,col.names = T,quote = F)

library(CCPlotR)
inp <- data.frame(source = ints$Lig.pop,
                  target = ints$Rec.pop,
                  ligand = ints$Ligand,
                  receptor = ints$Receptor,
                  score = log(ints$score),
                  stringsAsFactors = F)
inp <- inp[which(inp$source %in% c("Neuron","Astrocyte","Microglia") & inp$target %in% c("Neuron","Astrocyte","Microglia")),]
CCPlotR::cc_sigmoid(inp,n_top_ints = nrow(ints),colours = cols[c("Neuron","Astrocyte","Microglia")])
p <- CCPlotR::cc_dotplot(inp,n_top_ints = nrow(ints),option = "CellPhoneDB") +
  scale_fill_viridis(option = "C") +
  guides(fill = "none",size = guide_legend(override.aes = list(fill = viridis(option = "C",n = 4,direction = -1))))
p$scales$scales[[1]]$name <- "Log-Score"

ggsave("./FinalPlots_sNuc/CCI.png",p,dpi = 300,height = 8)

ints_sub <- ints[,c(1,5)]
df <- aggregate(Num ~., data=transform(ints_sub,Num=1), length)
write.table(df,file = "./sNuc-seq/Intercom_forPaper_RNA/NumPopInteractions.txt",sep = "\t",row.names = F,col.names = T,quote = F)

library(circlize)
library(reshape2)
library(ggalluvial)
inp <- df
inp$Rec.pop <- paste0(inp$Rec.pop,"BLBLBL")
inp <- data.frame(Type = c("Sender","Receiver","Sender","Receiver","Sender","Receiver","Sender","Receiver","Sender","Receiver","Sender","Receiver","Sender","Receiver","Sender","Receiver","Sender","Receiver"),
                  CT = c("Astrocyte","Astrocyte","Astrocyte","Microglia","Astrocyte","Neuron","Microglia","Astrocyte","Microglia","Microglia","Microglia","Neuron","Neuron","Astrocyte","Neuron","Microglia","Neuron","Neuron"),
                  ID = c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9),
                  Num = c(10,10,16,16,12,12,7,7,14,14,8,8,11,11,15,15,13,13),
                  stringsAsFactors = F)
inp$Type <- factor(inp$Type,levels = c("Sender","Receiver"))
p <- ggplot(inp,
       aes(x = Type,stratum = CT, alluvium = ID,
           y = Num,
           fill = CT, label = CT)) +
  geom_stratum(alpha = 0.5) +
  scale_fill_viridis(discrete = T,option = "C") +
  geom_flow() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            reverse = FALSE) +
  theme_void() +
  theme(legend.position = "none",
        axis.text.x = element_text())

ggsave("./FinalPlots_sNuc/InteractionFlow.png",p,dpi = 300)


ints_sub <- ints[,c(2,4)]
df <- aggregate(Num ~., data=transform(ints_sub,Num=1), length)
write.table(df,file = "./sNuc-seq/Intercom_forPaper_RNA/NumInteractions.txt",sep = "\t",row.names = F,col.names = T,quote = F)
library(reshape2)
inp <- dcast(df, Ligand~Receptor, value.var="Num")
inp[is.na(inp)] <- 0
rownames(inp) <- inp$Ligand
inp$Ligand <- NULL

library(pheatmap)
pheatmap(inp,
         scale = "none",
         height = 12,
         filename = "./sNuc-seq/Intercom_forPaper_RNA/NumInteractions.png")



