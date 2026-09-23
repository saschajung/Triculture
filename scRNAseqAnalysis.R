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

y1 <- readRDS("./scRNA-seq/CellRanger-Count-10X/Analysis_Y1/scRNA_Y1_SeuratObj.rds")

cds <- as.cell_data_set(y1)
cds <- cluster_cells(cds,k = 20,
                     reduction_method = "UMAP",
                     cluster_method = "leiden",
                     resolution = 0.02,
                     num_iter = 5)

y1$seurat_clusters <- clusters(cds)
Idents(y1) <- clusters(cds)

DimPlot(y1,reduction = "umap", group.by = "seurat_clusters", label = T, repel = T)+theme(legend.position="none")

opc_markers <- c("PDGFRA","OLIG1")
radialGlia_markers <- c("FABP7","SOX2","HES1","HES5")
neuron_markers <- c("RBFOX3","SNAP25","MAP2","STMN2","TUBB3")
fibroblast <- c("COL1A1","COL1A2","DCN","LUM","KRT18")
microglia <- c("P2RY12","CX3CR1","CSF1R","SALL1","CCL2","CCL3","CCL4","MKI67")
astro <- c("GFAP","AQP4","SLC1A3","ALDH1L1")
featuresToPlot <- c(opc_markers,radialGlia_markers,
                    neuron_markers,fibroblast,microglia,astro)

DotPlot(y1,features = featuresToPlot,assay = "RNA")

ct_map <- c("1" = "Astrocyte",
            "2" = "Microglia",
            "3" = "Microglia",
            "4" = "Proliferating Microglia",
            "5" = "Microglia",
            "6" = "Microglia",
            "7" = "Microglia",
            "8" = "Radial Glia",
            "9" = "Microglia",
            "10" = "Microglia",
            "11" = "Astrocyte",
            "12" = "Activated Microglia",
            "13" = "Microglia",
            "14" = "Astrocyte",
            "15" = "Microglia",
            "16" = "Microglia",
            "17" = "Microglia",
            "18" = "Microglia",
            "19" = "Microglia",
            "20" = "Microglia",
            "21" = "Microglia",
            "22" = "Microglia",
            "23" = "Microglia",
            "24" = "Microglia",
            "25" = "Fibroblast-like",
            "26" = "Microglia",
            "27" = "Astrocyte",
            "28" = "Astrocyte",
            "29" = "Radial Glia",
            "30" = "Microglia",
            "31" = "Microglia",
            "32" = "Radial Glia",
            "33" = "OLIG1+GFAP+ Cells",
            "34" = "Fibroblast-like",
            "35" = "OPC",
            "36" = "Microglia")

y1$CTs <- unname(ct_map[y1$seurat_clusters])
y1@meta.data[names(which(y1@reductions$umap@cell.embeddings[,1] < -10 & y1@reductions$umap@cell.embeddings[,2] < -1)),"CTs"] <- "Neuron"

saveRDS(y1,"./scRNA-seq/CellRanger-Count-10X/Analysis_Y1/scRNA_Y1_Clustered_SeuratObj.rds")

DimPlot(y1,reduction = "umap", group.by = "CTs", label = T, repel = T)+theme(legend.position="none")

plot_cor_full <- round(x = cor(x = y1@meta.data[,"nCount_RNA"], y = y1@meta.data[, "nFeature_RNA"]),digits = 2)
df <- y1@meta.data[,c("nCount_RNA","nFeature_RNA")]
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
p
ggsave("./FinalPlots_scRNA/GenesvsUMIs_rev.png",p)

cols <- c("Astrocyte" = "#CD534CFF",
          "Fibroblast-like" = "#868686FF",
          "Microglia" = "#0073C2FF",
          "Activated Microglia" = "#7AA6DCFF",
          "Proliferating Microglia" = "#003C67FF",
          "Neuron" = "#EFC000FF",
          "OLIG1+GFAP+ Cells" = "#8F7700FF",
          "OPC" = "#3B3B3BFF",
          "Radial Glia" = "#A73030FF")#cols <- pal_jco()(n = length(unique(y1$CTs)))
#names(cols) <- unique(y1$CTs)
p <- DimPlot_scCustom(seurat_object = y1,
                      repel = T,
                      group.by = "CTs",
                      figure_plot = TRUE,
                      label = F,
                      colors_use = cols,
                      prop_plot_percent = T,
                      add_prop_plot = T,
                      prop_plot_label = T,
                      label.size = 16)
p[[1]][[1]] <- p[[1]][[1]] + ggtitle("") + scale_color_manual(values = cols) + theme(legend.text = element_text(size = 16))
p[[2]][[1]] <- p[[2]][[1]] + scale_fill_manual(values = cols) + theme(axis.text = element_text(size = 16,color = "black"), axis.title = element_text(size = 20))
ggsave("./FinalPlots_scRNA/UMAP_rev.png",p, width = 12, height = 5)


png(filename = "./FinalPlots_scRNA/Dotplot_Markers.png",res = 300,width = 3000, height = 2400)
Clustered_DotPlot(seurat_object = y1,
                  features = featuresToPlot,
                  group.by = "CTs",
                  assay = "SCT",
                  k = 6,
                  colors_use_idents = cols[c(2,1,3,5,4,7,8,9,6)],
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

png(filename = "./FinalPlots_scRNA/Dotplot_Cytokines.png",res = 300,width = 2500, height = 3500)
p <- Clustered_DotPlot(seurat_object = y1,
                       features = Cytokines_All,
                       group.by = "CTs",
                       assay = "SCT",
                       k = 6,
                       colors_use_idents = cols[c(2,1,3,5,4,7,8,9,6)],
                       row_label_size = 16,
                       column_label_size = 16,
                       legend_label_size = 16,
                       legend_title_size = 20)
dev.off()

library(InterCom)
inpmat <- as.matrix(y1@assays$RNA@counts)
anno.tbl <- data.frame(cell.name = rownames(y1@meta.data),
                       cell.type = unname(y1$CTs),
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
         tissue.name = "forPaperFinal",
         out.path = "./scRNA-seq/Y1_forPaperFinal/")

load("./scRNA-seq/Y1_forPaperFinal/output_forPaperFinal.RData",verbose = T)
ints <- output$final
ints <- ints[ints$score >= 0.9,]
pop_map <- c("OLIG1.GFAP..Cells" = "OLIG1+GFAP+ Cells",
             "Intermediate.Astrocyte.Progenitor" = "Intermediate Astrocyte Progenitor",
             "Radial.Glia" = "Radial Glia",
             "Fibroblast.like" = "Fibroblast-like",
             "Astrocyte" = "Astrocyte",
             "Proliferating.Microglia" = "Proliferating Microglia",
             "Activated.Microglia" = "Activated Microglia",
             "Neuron" = "Neuron",
             "OPC" = "OPC",
             "Microglia" = "Microglia")
ints$Lig.pop <- pop_map[ints$Lig.pop]
ints$Rec.pop <- pop_map[ints$Rec.pop]

write.table(ints,file = "./scRNA-seq/Y1_forPaperFinal/FinalInteractions.txt",sep = "\t",row.names = F,col.names = T,quote = F)

library(Libra)
for(ct in unique(y1$CTs)){
  y1$replicate <- ""
  y1$replicate[which(y1$CTs == ct)] <- sample(paste0(ct,1:3),length(which(y1$CTs == ct)),replace = T)
  y1$replicate[which(y1$CTs != ct)] <- sample(paste0("None",1:3),length(which(y1$CTs != ct)),replace = T)
  y1$label <- ct
  y1$label[which(y1$CTs != ct)] <- "Background"
  y1$ctcol <- ct
  diff <- Libra::run_de(input = y1,
                        cell_type_col = "ctcol",
                        replicate_col = "replicate",
                        label_col = "label",
                        de_family = "pseudobulk",
                        de_method = "edgeR",
                        de_type = "LRT")
}

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

library(ComplexHeatmap)
mat <- as.matrix(aggExp$SCTscaled)
p <- pheatmap::pheatmap(mat[selectedGenes,],
         clustering_method = "ward.D2")

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

saveRDS(inp,"./PlotsRevision/InputData_Fig4E.rds")
inp <- readRDS("./PlotsRevision/InputData_Fig4E.rds")
p <- ggplot(inp, aes(x = Var2, y = Var1, fill= value)) + 
  geom_tile() +
  scale_fill_gradient2(low = "#F0F921FF", high = '#0D0887FF', mid = "#CC4678FF") +
  xlab("") +
  ylab("") +
  labs(fill = "Z-Score")+
  theme_pubclean() +
  scale_y_discrete(position = "right")+
  theme(axis.text = element_text(size = 16),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 20))

ggsave("./FinalPlots_scRNA/Heatmap_Top10Markers.png",p,width = 5, height = 20)



ints_sub <- ints[,c(1,5)]
df <- aggregate(Num ~., data=transform(ints_sub,Num=1), length)
write.table(df,file = "./scRNA-seq/Y1_Intercom/NumPopInteractions.txt",sep = "\t",row.names = F,col.names = T,quote = F)

library(circlize)
library(reshape2)
inp <- dcast(df, Lig.pop~Rec.pop, value.var="Num")
inp[is.na(inp)] <- 0
rownames(inp) <- inp$Lig.pop
inp$Lig.pop <- NULL

png(filename = "./scRNA-seq/Y1_Intercom/NumPopInteractions.png",width = 7, height = 7, units = "in", res = 600)
chordDiagram(inp, transparency = 0.5)
dev.off()

ints_sub <- ints[,c(2,4)]
df <- aggregate(Num ~., data=transform(ints_sub,Num=1), length)
write.table(df,file = "./scRNA-seq/Y1_Intercom/NumInteractions.txt",sep = "\t",row.names = F,col.names = T,quote = F)
library(reshape2)
inp <- dcast(df, Ligand~Receptor, value.var="Num")
inp[is.na(inp)] <- 0
rownames(inp) <- inp$Ligand
inp$Ligand <- NULL

library(pheatmap)
pheatmap(inp,
         scale = "none",
         height = 12,
         filename = "./scRNA-seq/Y1_Intercom/NumInteractions.png")


#### Aux analyses ####
library(viridis)
clustProps <- as.data.frame(prop.table(table(y1$seurat_clusters,y1$orig.ident),margin = 2))
colnames(clustProps) <- c("Cluster","Sample","Percent")
p <- ggplot(clustProps, aes(fill=Sample, y=Percent, x=Cluster)) + 
  geom_bar(position="dodge", stat="identity") +
  scale_fill_viridis(option = "B",discrete = T,begin = 0.2, end=0.75) +
  theme(axis.title = element_text(size = rel(2)),
        axis.text = element_text(size = rel(2)),
        legend.text = element_text(size = rel(2)),
        legend.title = element_text(size = rel(2)))

ggsave("./scRNA-seq/Y1_ClusterProportions.png",p,width = 15)

####Organoid integration
rm(combined.sct)
files <- list.files("./OtherStudies/scOrganoids",pattern = ".rds",recursive = T, full.names = T)
obj.list <- list()
for(f in files){
  n <- rev(strsplit(f,"/")[[1]])[2]
  obj.list[[n]] <- readRDS(f)
  obj.list[[n]]$CTs <- obj.list[[n]]$CellType
}

obj.list[["y1"]] <- y1
gc()
gc()
gc()

rm(y1)

obj.list <- lapply(obj.list,function(x){
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
plot.path <- "./OtherStudies/scOrganoids/Y1Integration/"
dir.create(plot.path)
ggsave(paste0(plot.path,"ElbowPlot.png"),ElbowPlot(combined.sct,ndims = 100)) #Take 50

numPCs <- 30

DefaultAssay(combined.sct) <- "integrated"
combined.sct <- RunUMAP(combined.sct,
                        reduction = "pca",
                        dims = 1:numPCs,
                        verbose = T)
combined.sct <- RunTSNE(combined.sct, reduction = "pca", dims = 1:numPCs, verbose = TRUE)
save(combined.sct,file = "./OtherStudies/scOrganoids/Y1Integration/Combined_SeuratObj.RData",compression_level = 6)

DimPlot(combined.sct,group.by = "CTs")

##### Plot
load("./OtherStudies/scOrganoids/Y1Integration/Combined_SeuratObj.RData",verbose = T)
combined.sct$Origin <- "This study"
combined.sct$Origin[combined.sct$orig.ident != "Y1"] <- "Organoid"

combined.sct$CTOrigin <- paste0(combined.sct$CTs," (",combined.sct$Origin,")")

combined.sct$Clusts <- paste0(combined.sct$seurat_clusters," (",combined.sct$Origin,")")

sub <- subset(combined.sct,subset = orig.ident != "Y1")
DefaultAssay(sub) <- "RNA"
sub[["SCT"]] <- NULL
sub[["integrated"]] <- NULL
options(future.globals.maxSize = 4000 * 1024^2)
sub <- SCTransform(sub, vars.to.regress = "CC.Difference",vst.flavor = "v2", verbose = T)

commonGenes <- intersect(rownames(sub@assays$SCT),rownames(y1@assays$SCT))

pseudo_ref <- PseudobulkExpression(object = sub,
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

pseudo_query <- PseudobulkExpression(object = y1,
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

saveRDS(cormat, file = "./PlotsRevision/InputData_Heatmap_scOrganoid_vs_ours.rds")


load("./OtherStudies/scOrganoids/Y1Integration/Combined_SeuratObj.RData",verbose = T)
combined.sct$Origin <- "This study"
combined.sct$Origin[combined.sct$orig.ident != "Y1"] <- "Organoid"

combined.sct$CTOrigin <- paste0(combined.sct$CTs," (",combined.sct$Origin,")")

combined.sct$Clusts <- paste0(combined.sct$seurat_clusters," (",combined.sct$Origin,")")

cols <- c("Astrocyte (This study)" = "#CD534CFF",
          "Fibroblast-like (This study)" = "darkgrey",
          "Microglia (This study)" = "darkgrey",
          "Activated Microglia (This study)" = "darkgrey",
          "Proliferating Microglia (This study)" = "darkgrey",
          "Neuron (This study)" = "#EFC000FF",
          "OLIG1+GFAP+ Cells (This study)" = "darkgrey",
          "OPC (This study)" = "darkgrey",
          "Radial Glia (This study)" = "darkgrey",
          "Astroglia (Organoid)" = lighten("#CD534CFF",amount = 0.5),
          "oRG/Astroglia (Organoid)" = "darkgrey",
          "oRG (Organoid)" = "darkgrey",
          "Ventral Precursors (Organoid)" = "darkgrey",
          "Immature CPNs (Organoid)" = lighten("#EFC000FF",amount = 0.5),
          "CPNs (Organoid)" = darken("#EFC000FF",amount = 0.2),
          "Immature Interneurons (Organoid)" = lighten("#EFC000FF",amount = 0.3),
          "Cycling (Organoid)" = "darkgrey",
          "IPCs (Organoid)" = "darkgrey",
          "Immature PNs (Organoid)" = lighten("#EFC000FF",amount = 0.4),
          "RG (Organoid)" = "darkgrey",
          "Unknown (Organoid)" = "darkgrey")
p <- DimPlot_scCustom(seurat_object = combined.sct,
                      repel = T,
                      group.by = "CTOrigin",
                      figure_plot = TRUE,
                      label = F,
                      label.color = "white",
                      label.box	= T,
                      colors_use = cols,
                      label.size = 16/.pt)
#p[[1]] <- p[[1]] + ggtitle("") + theme(legend.position = "none")

# Define which cell types to label
highlight_types <- c("Astrocyte (This study)", 
                     "Neuron (This study)",
                     "Astroglia (Organoid)",
                     "Immature PNs (Organoid)",
                     "Immature CPNs (Organoid)",
                     "Immature Interneurons (Organoid)",
                     "CPNs (Organoid)")

# Compute centroids from plot data
centroids <- p[[1]]$data |>
  dplyr::group_by(CTOrigin) |>
  dplyr::summarise(x = median(.data[[names(p[[1]]$data)[1]]]),
                   y = median(.data[[names(p[[1]]$data)[2]]]))

# Filter to only highlighted types
centroids_highlight <- centroids |> 
  dplyr::filter(CTOrigin %in% highlight_types)

# Add labels only for highlighted types
p[[1]] <- p[[1]] + 
  geom_label_repel(
    data        = centroids_highlight,
    aes(x = x, y = y, label = CTOrigin, fill = CTOrigin),
    color       = "white",
    force       = 15,
    max.overlaps = Inf,
    box.padding = 0.5,
    size        = 16/.pt,
    show.legend = FALSE,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = cols) +  # match your existing color scheme
  ggtitle("") + 
  theme(legend.position = "none")

p
ggsave("./PlotsRevision/Fig4I_UMAP_CTsOrigin_rev.png",p,width = 8, height = 8)


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

library(Libra)
load("./OtherStudies/scOrganoids/Y1Integration/Combined_SeuratObj.RData",verbose = T)
combined.sct$Origin <- "This study"
combined.sct$Origin[combined.sct$orig.ident != "Y1"] <- "Organoid"

combined.sct$CTOrigin <- paste0(combined.sct$CTs," (",combined.sct$Origin,")")

combined.sct$Clusts <- paste0(combined.sct$seurat_clusters," (",combined.sct$Origin,")")

sub <- subset(combined.sct, subset = CTOrigin %in% c("Astrocyte (This study)","Astroglia (Organoid)"))
DefaultAssay(sub) <- "RNA"
sub$CTs <- "Astrocyte"
sub$replicate <- sub$orig.ident
set.seed(1234)
sub$replicate[which(sub$orig.ident == "11a.6mon")] <- sample(paste0("11a.6mon",1:2),length(which(sub$orig.ident == "11a.6mon")),replace = T)
sub$replicate[which(sub$orig.ident == "PGP1.6mon")] <- sample(paste0("PGP1.6mon",1:2),length(which(sub$orig.ident == "PGP1.6mon")),replace = T)
sub$replicate[which(sub$orig.ident == "Y1")] <- sample(paste0("Y1",1:4),length(which(sub$orig.ident == "Y1")),replace = T)
for(assay_name in Assays(sub)) {
  sub[[assay_name]] <- as(sub[[assay_name]], "Assay5")
}
astrocyte.diff <- Libra::run_de(input = sub,
                                cell_type_col = "CTs",
                                replicate_col = "replicate",
                                label_col = "Origin",
                                de_family = "pseudobulk",
                                de_method = "edgeR",
                                de_type = "QLF") #Negative FC means higher in organoid

sub <- subset(combined.sct, subset = CTOrigin %in% c("Neuron (This study)","CPNs (Organoid)"))
DefaultAssay(sub) <- "RNA"
sub$CTs <- "Neuron"
sub$replicate <- sub$orig.ident
set.seed(1234)
sub$replicate[which(sub$orig.ident == "11a.6mon")] <- sample(paste0("11a.6mon",1:2),length(which(sub$orig.ident == "11a.6mon")),replace = T)
sub$replicate[which(sub$orig.ident == "PGP1.6mon")] <- sample(paste0("PGP1.6mon",1:2),length(which(sub$orig.ident == "PGP1.6mon")),replace = T)
sub$replicate[which(sub$orig.ident == "Y1")] <- sample(paste0("Y1",1:4),length(which(sub$orig.ident == "Y1")),replace = T)
for(assay_name in Assays(sub)) {
  sub[[assay_name]] <- as(sub[[assay_name]], "Assay5")
}
neuron.diff <- Libra::run_de(input = sub,
                                cell_type_col = "CTs",
                                replicate_col = "replicate",
                                label_col = "Origin",
                                de_family = "pseudobulk",
                                de_method = "edgeR",
                                de_type = "QLF") #Negative FC means higher in organoid

sub <- subset(combined.sct, subset = CTOrigin %in% c("Radial Glia (This study)","RG (Organoid)"))
DefaultAssay(sub) <- "RNA"
sub$CTs <- "Radial Glia"
sub$replicate <- sub$orig.ident
set.seed(1234)
sub$replicate[which(sub$orig.ident == "11a.6mon")] <- sample(paste0("11a.6mon",1:2),length(which(sub$orig.ident == "11a.6mon")),replace = T)
sub$replicate[which(sub$orig.ident == "PGP1.6mon")] <- sample(paste0("PGP1.6mon",1:2),length(which(sub$orig.ident == "PGP1.6mon")),replace = T)
sub$replicate[which(sub$orig.ident == "Y1")] <- sample(paste0("Y1",1:4),length(which(sub$orig.ident == "Y1")),replace = T)
for(assay_name in Assays(sub)) {
  sub[[assay_name]] <- as(sub[[assay_name]], "Assay5")
}
radialglia.diff <- Libra::run_de(input = sub,
                             cell_type_col = "CTs",
                             replicate_col = "replicate",
                             label_col = "Origin",
                             de_family = "pseudobulk",
                             de_method = "edgeR",
                             de_type = "QLF") #Negative FC means higher in organoid

sub <- subset(combined.sct, subset = CTOrigin %in% c("Radial Glia (This study)","RG (Organoid)"))


save(astrocyte.diff,radialglia.diff,neuron.diff,file = "./OtherStudies/scOrganoids/Y1Integration/DiffExp_rev.RData",compression_level = 6)
load("./OtherStudies/scOrganoids/Y1Integration/DiffExp_rev.RData", verbose = T)

write.table(astrocyte.diff[which(abs(astrocyte.diff$avg_logFC) > 0.5 & astrocyte.diff$p_val_adj <= 0.01),],"./OtherStudies/scOrganoids/Y1Integration/DiffExp_Astrocyte.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(neuron.diff[which(abs(neuron.diff$avg_logFC) > 0.5 & neuron.diff$p_val_adj <= 0.01),],"./OtherStudies/scOrganoids/Y1Integration/DiffExp_Neuron.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(radialglia.diff[which(abs(radialglia.diff$avg_logFC) > 0.5 & radialglia.diff$p_val_adj <= 0.01),],"./OtherStudies/scOrganoids/Y1Integration/DiffExp_RadialGlia.txt",sep = "\t",row.names = F, col.names = T, quote = F)


library(WebGestaltR)
# Neuron
inp_neuron_list <- list("This study" = neuron.diff$gene[neuron.diff$p_val_adj <= 0.01 & neuron.diff$avg_logFC > -0.5],
                             "Brain" = neuron.diff$gene[neuron.diff$p_val_adj <= 0.01 & neuron.diff$avg_logFC < 0.5])

neuronBrain <- WebGestaltR(enrichMethod = "ORA",
                             organism = "hsapiens",
                             enrichDatabase = "geneontology_Biological_Process",
                             interestGene = inp_neuron_list$Brain,
                             interestGeneType = "genesymbol",
                             referenceGene = neuron.diff$gene,
                             referenceGeneType = "genesymbol",
                             sigMethod = "fdr",
                             fdrMethod = "BH",
                             fdrThr = 0.01,
                             isOutput = F
)

neuronThis <- WebGestaltR(enrichMethod = "ORA",
                            organism = "hsapiens",
                            enrichDatabase = "geneontology_Biological_Process",
                            interestGene = inp_neuron_list$`This study`,
                            interestGeneType = "genesymbol",
                            referenceGene = neuron.diff$gene,
                            referenceGeneType = "genesymbol",
                            sigMethod = "fdr",
                            fdrMethod = "BH",
                            fdrThr = 0.01,
                            isOutput = F
)

# Astrocyte
inp_astro_list <- list("This study" = astrocyte.diff$gene[astrocyte.diff$p_val_adj <= 0.01 & astrocyte.diff$avg_logFC > -0.5],
                        "Brain" = astrocyte.diff$gene[astrocyte.diff$p_val_adj <= 0.01 & astrocyte.diff$avg_logFC < 0.5])

astroBrain <- WebGestaltR(enrichMethod = "ORA",
                           organism = "hsapiens",
                           enrichDatabase = "geneontology_Biological_Process",
                           interestGene = inp_astro_list$Brain,
                           interestGeneType = "genesymbol",
                           referenceGene = astrocyte.diff$gene,
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
                          referenceGene = astrocyte.diff$gene,
                          referenceGeneType = "genesymbol",
                          sigMethod = "fdr",
                          fdrMethod = "BH",
                          fdrThr = 0.01,
                          isOutput = F
)

# RG
inp_rg_list <- list("This study" = radialglia.diff$gene[radialglia.diff$p_val_adj <= 0.01 & radialglia.diff$avg_logFC > -0.5],
                       "Brain" = radialglia.diff$gene[radialglia.diff$p_val_adj <= 0.01 & radialglia.diff$avg_logFC < 0.5])

rgBrain <- WebGestaltR(enrichMethod = "ORA",
                          organism = "hsapiens",
                          enrichDatabase = "geneontology_Biological_Process",
                          interestGene = inp_rg_list$Brain,
                          interestGeneType = "genesymbol",
                          referenceGene = radialglia.diff$gene,
                          referenceGeneType = "genesymbol",
                          sigMethod = "fdr",
                          fdrMethod = "BH",
                          fdrThr = 0.01,
                          isOutput = F
)

rgThis <- WebGestaltR(enrichMethod = "ORA",
                         organism = "hsapiens",
                         enrichDatabase = "geneontology_Biological_Process",
                         interestGene = inp_rg_list$`This study`,
                         interestGeneType = "genesymbol",
                         referenceGene = radialglia.diff$gene,
                         referenceGeneType = "genesymbol",
                         sigMethod = "fdr",
                         fdrMethod = "BH",
                         fdrThr = 0.01,
                         isOutput = F
)

write.table(astroBrain,file = "./OtherStudies/scOrganoids/Y1Integration/GO_AstrocyteOrganoid_rev.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(astroThis,file = "./OtherStudies/scOrganoids/Y1Integration/GO_AstrocyteThis_rev.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(rgBrain,file = "./OtherStudies/scOrganoids/Y1Integration/GO_RadialGliaOrganoid_rev.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(rgThis,file = "./OtherStudies/scOrganoids/Y1Integration/GO_RadialGliaThis_rev.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(neuronBrain,file = "./OtherStudies/scOrganoids/Y1Integration/GO_NeuronOrganoid_rev.txt",sep = "\t",row.names = F, col.names = T, quote = F)
write.table(neuronThis,file = "./OtherStudies/scOrganoids/Y1Integration/GO_NeuronThis_rev.txt",sep = "\t",row.names = F, col.names = T, quote = F)

syngo <- read.table("./syngo_ontologies_BP_GO.txt",header = F)
neuronBrain_SynGO <- neuronBrain[which(neuronBrain$geneSet %in% syngo$V1),]
neuronThis_SynGO <- neuronThis[which(neuronThis$geneSet %in% syngo$V1),]

saveRDS(neuronThis_SynGO,"./PlotsRevision/InputData_Fig4J_rev.rds")

neuronThis_SynGO <- 
p <- ggbarplot(neuronThis_SynGO,
          x = "description",
          y = "enrichmentRatio",
          fill = "FDR",
          sort.val = "asc") +
  guides(fill = guide_colorbar(barheight = 1, barwidth = 15)) +
  scale_fill_material("red",reverse = T)+
  coord_flip() +
  xlab("") + 
  font("legend.title",size = 20) +
  font("legend.text",size = 16) +
  font("xy",size = 20) +
  font("xy.text",size = 16)

ggsave(paste0(plot.path,"GOEnrichment_Neuron_ThisStudy.png"),p,width = 10, height = 10)

write.table(astroBrain,file = "./OtherStudies/scOrganoids/Y1Integration/AstroBrain.txt",row.names = F,col.names = T,quote =F, sep = "\t")
write.table(astroThis,file = "./OtherStudies/scOrganoids/Y1Integration/AstroThisn.txt",row.names = F,col.names = T,quote =F, sep = "\t")
write.table(neuronBrain,file = "./OtherStudies/scOrganoids/Y1Integration/NeuronBrain.txt",row.names = F,col.names = T,quote =F, sep = "\t")
write.table(neuronThis,file = "./OtherStudies/scOrganoids/Y1Integration/NeuronThis.txt",row.names = F,col.names = T,quote =F, sep = "\t")
write.table(rgBrain,file = "./OtherStudies/scOrganoids/Y1Integration/RGBrain.txt",row.names = F,col.names = T,quote =F, sep = "\t")
write.table(rgThis,file = "./OtherStudies/scOrganoids/Y1Integration/RGThis.txt",row.names = F,col.names = T,quote =F, sep = "\t")

inp.list <- list("astroThis" = astroThis,
                 "astroOrg" = astroOrg,
                 "rgOrg" = rgOrg,
                 "rgThis" = rgThis,
                 "cpnOrg" = cpnOrg,
                 "cpnThis" = cpnThis,
                 "interOrg" = interOrg,
                 "interThis" = interThis)

title_map <- c("astroThis" = "Astrocyte (This study)",
               "astroOrg" = "Astrocyte (Organoid)",
               "rgOrg" = "RG (Organoid)",
               "rgThis" = "Radial Glia (This study)",
               "cpnOrg" = "Immature CPNs (Organoid)",
               "cpnThis" = "Neuron (This study)",
               "interOrg" = "Immature Interneurons (Organoid)",
               "interThis" = "Neuron (This study)")



avgExp <- AverageExpression(combined.sct,assay = c("RNA","SCT"),group.by = "CTOrigin")

avgExp$SCT <- avgExp$SCT[rowSds(avgExp$SCT) > 0,]
avgExp$RNA <- avgExp$RNA[rowSds(avgExp$RNA) > 0,]

library(pheatmap)
pheatmap(avgExp$RNA,
         scale = "row",
         show_rownames = F,
         clustering_distance_cols = "correlation",
         filename = "./OtherStudies/scOrganoids/Y1Integration/AvgExpression.png")

