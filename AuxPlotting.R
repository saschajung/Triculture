setwd("./sNuc-seq")
library(Seurat)
library(DoubletFinder)
library(RColorBrewer)
library(ggplot2)
library(ggExtra)
library(cowplot)
library(dplyr)
library(cluster, quietly = TRUE)
library(ggstatsplot)
library(plyr)
library(biomaRt)
library(densvis)
library(monocle3)
library(singleseqgset)
library(pheatmap)
library(SeuratWrappers)
library(GSA)
library(ggh4x)
library(viridis)
library(corto)
library(Nebulosa)
library(ggheatmap)
library(ggtree)
library(ggridges)
library(patchwork)
library(scCustomize)
library(ggpubr)

pal_jco_cont <- colorRampPalette(c("#CD534CFF","#EFC000FF","#0073C2FF"))
pal_jco_heatmap <- colorRampPalette(c("#CD534CFF","#FFFFFFFF","#0073C2FF"))
pal_jco_size <- colorRampPalette(c("#CD534CFF","#0073C2FF"))

colorMap <- c("Astrocyte" = "#CD534CFF",
              "Fibroblast-like" = "#868686FF",
              "Microglia" = "#0073C2FF",
              "Activated Microglia" = "#7AA6DCFF",
              "Proliferating Microglia" = "#003C67FF",
              "Neuron" = "#EFC000FF",
              "OLIG1+GFAP+ Cells" = "#8F7700FF",
              "OPC" = "#3B3B3BFF",
              "Radial Glia" = "#A73030FF")

### Fig4D:
y1 <- readRDS(file = "./scRNA-seq/CellRanger-Count-10X/Analysis_Y1/scRNA_Y1_Clustered_SeuratObj.rds")

opc_markers <- c("PDGFRA","OLIG1")
radialGlia_markers <- c("FABP7","SOX2","HES1","HES5")
neuron_markers <- c("RBFOX3","SNAP25","MAP2","STMN2","TUBB3")
fibroblast <- c("COL1A1","COL1A2","DCN","LUM","KRT18")
microglia <- c("P2RY12","CX3CR1","CSF1R","SALL1","CCL2","CCL3","CCL4","MKI67")
astro <- c("GFAP","AQP4","SLC1A3","ALDH1L1")
featuresToPlot <- c(opc_markers,radialGlia_markers,
                    neuron_markers,fibroblast,microglia,astro)

png(filename = "./PlotsRevision/Fig4D_Dotplot_Markers_rev.png",res = 300,width = 8, height = 8,units = "in",pointsize = 16)
Clustered_DotPlot(seurat_object = y1,
                  features = featuresToPlot,
                  group.by = "CTs",
                  assay = "SCT",
                  colors_use_exp = pal_jco_heatmap(256),
                  k = 6,
                  colors_use_idents = colorMap[c(1,3,5,9,4,2,7,6,8)],
                  plot_km_elbow = F,
                  row_label_size = 16,
                  column_label_size = 0,
                  legend_label_size = 16,
                  legend_title_size = 20)
dev.off()

### Fig4C:
scNuc <- readRDS(file = "./sNuc-seq/sNuc_SeuratObj.rds")

opc_markers <- c("PDGFRA","OLIG1")
radialGlia_markers <- c("FABP7","SOX2","HES1","HES5")
neuron_markers <- c("RBFOX3","SNAP25","MAP2","STMN2","TUBB3")
fibroblast <- c("COL1A1","COL1A2","DCN","LUM","KRT18")
microglia <- c("P2RY12","CX3CR1","CSF1R","SALL1","CCL2","CCL3","CCL4","MKI67")
astro <- c("GFAP","AQP4","SLC1A3","ALDH1L1")
featuresToPlot <- c(opc_markers,radialGlia_markers,
                    neuron_markers,fibroblast,microglia,astro)

png(filename = "./PlotsRevision/Fig4C_Dotplot_Markers_rev.png",res = 300,width = 8, height = 8,units = "in",pointsize = 16)
Clustered_DotPlot(seurat_object = scNuc,
                 features = featuresToPlot,
                 group.by = "SubCTs",
                 assay = "SCT",
                 colors_use_exp = pal_jco_heatmap(256),
                 k = 6,
                 colors_use_idents = colorMap[c(1,7,9,3,8,6,2)],
                 plot_km_elbow = F,
                 row_label_size = 16,
                 column_label_size = 0,
                 legend_label_size = 16,
                 legend_title_size = 20)
dev.off()

### Fig4E
inp <- readRDS("./PlotsRevision/InputData_Fig4E.rds")

library(ComplexHeatmap)
library(circlize)
library(tibble)

# Reshape your data back to wide matrix if it's in long format
mat <- inp %>%
  pivot_wider(names_from = Var2, values_from = value) %>%
  column_to_rownames("Var1") %>%
  as.matrix()

# Recreate your diverging color scale
col_fun <- colorRamp2(
  c(min(mat), 0, max(mat)),
  pal_jco_heatmap(3)   # your existing palette
)

all_idx <- 1:nrow(mat)
groups <- ceiling(all_idx / 10)
right_idx <- all_idx[groups %% 2 == 1]  # groups 1, 3, 5... → right
left_idx  <- all_idx[groups %% 2 == 0]  # groups 2, 4, 6... → left

# Create annotations
right_anno <- rowAnnotation(
  genes = anno_mark(
    at = right_idx,
    labels = rownames(mat)[right_idx],
    labels_gp = gpar(fontsize = 16),
    side = "right",
    link_width = unit(10, "mm")
  )
)

left_anno <- rowAnnotation(
  genes = anno_mark(
    at = left_idx,
    labels = rownames(mat)[left_idx],
    labels_gp = gpar(fontsize = 16),
    side = "left",
    link_width = unit(10, "mm")
  )
)

ht <- Heatmap(
  mat,
  name = "Z-Score",
  col = col_fun,
  right_annotation = right_anno,
  left_annotation = left_anno,
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 20, fontface = "bold"),
    labels_gp = gpar(fontsize = 16)
  ),
  show_row_names = FALSE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_names_rot = 90,
  column_names_gp = gpar(fontsize = 16)
)

png("./PlotsRevision/Fig4E_Heatmap_Top10Markers_rev.png", width = 6.5, height = 12.5,units = "in", pointsize = 16, res = 300)
draw(ht)
dev.off()

###Fig4H
load("./OtherStudies/GSE179590/sNucIntegration/Combined_SeuratObj.RData")

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
                      label = F,
                      label.color = "white",
                      label.box	= T,
                      colors_use = cols,
                      label.size = 16/.pt)
#p[[1]] <- p[[1]] + ggtitle("") + theme(legend.position = "none")

# Define which cell types to label
highlight_types <- c("Astrocyte (This study)", 
                     "Neuron (This study)",
                     "OPC (This study)",
                     "Microglia (This study)",
                     "Astrocyte (Brain)", 
                     "Inhibitory Neuron (Brain)",
                     "Excitatory Neuron (Brain)",
                     "OPC (Brain)",
                     "Microglia (Brain)",
                     "Oligodendrocyte (Brain)")

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
ggsave("./PlotsRevision/Fig4H_UMAP_CTsOrigin_rev.png",p,width = 8, height = 8)


### Fig4I
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


### Fig4J: Organoid Comparison with enriched SynGO terms
library(colorspace)
neuronThis_SynGO <- readRDS("./PlotsRevision/InputData_Fig4J_rev.rds")
neuronThis_SynGO$description[which(neuronThis_SynGO$description == "modulation of chemical synaptic transmission")] <- "modulation of chemical\nsynaptic transmission" 
neuronThis_SynGO$description[which(neuronThis_SynGO$description == "regulation of synapse organization")] <- "regulation of\nsynapse organization" 

  p <- ggbarplot(neuronThis_SynGO,
                 x = "description",
                 y = "enrichmentRatio",
                 fill = "FDR",
                 sort.val = "asc") +
  guides(fill = guide_colorbar(barheight = 1, barwidth = 15)) +
 #scale_fill_material("red",reverse = T) +
    scale_fill_gradient(high = lighten("#EFC000FF", amount = 0.8),low="#EFC000FF")+
    coord_flip() +
  xlab("") + 
  ylab("Enrichment Ratio") +
    theme(legend.title = element_text(size = 24),
          legend.text = element_text(size = 22),
          axis.text = element_text(size = 22),
          axis.title = element_text(size = 24))

ggsave("./PlotsRevision/Fig4J_GOEnrichment_Neuron_ThisStudy.png",p,dpi = 300,width = 8, height = 8)

### scRNAseq: Plot QC 
y1 <- readRDS("./scRNA-seq/CellRanger-Count-10X/Analysis_Y1/scRNA_Y1_Clustered_SeuratObj.rds")
plot_cor_full <- round(x = cor(x = y1@meta.data[,"nCount_RNA"], y = y1@meta.data[, "nFeature_RNA"]),digits = 2)
df <- y1@meta.data[,c("nCount_RNA","nFeature_RNA")]
df$nCount_RNA <- (df$nCount_RNA)
df$nFeature_RNA <- (df$nFeature_RNA)
p <- ggscatter(df, x = "nCount_RNA", y = "nFeature_RNA",
               color = "black", shape = 19, size = 2, # Points color, shape and size
               add = "reg.line",  # Add regressin line
               add.params = list(color = "#0073C2FF", size = 2), # Customize reg. line
               conf.int = TRUE, # Add confidence interval
               cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
               cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n"),
               yscale = "log10",
               xscale = "log10",
               format.scale = T,
               xlab = "UMIs per Nucleus",
               ylab = "Genes per Nucleus",
               cor.coef.size = 16/.pt,
) +
  font("xlab", size = 20, color = "black")+
  font("ylab", size = 20, color = "black")+
  font("xy.text", size = 16, color = "black")
p
ggsave("./PlotsRevision/GenesvsUMIs_scRNAseq_rev.png",p,dpi = 300, width = 6, height = 6)

### snRNAseq: Plot QC 
sNuc <- readRDS("./sNuc-seq/sNuc_SeuratObj.rds")
plot_cor_full <- round(x = cor(x = sNuc@meta.data[,"nCount_RNA"], y = sNuc@meta.data[, "nFeature_RNA"]),digits = 2)
df <- sNuc@meta.data[,c("nCount_RNA","nFeature_RNA")]
df$nCount_RNA <- (df$nCount_RNA)
df$nFeature_RNA <- (df$nFeature_RNA)
p <- ggscatter(df, x = "nCount_RNA", y = "nFeature_RNA",
               color = "black", shape = 19, size = 2, # Points color, shape and size
               add = "reg.line",  # Add regressin line
               add.params = list(color = "#0073C2FF", size = 2), # Customize reg. line
               conf.int = TRUE, # Add confidence interval
               cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
               cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n"),
               yscale = "log10",
               xscale = "log10",
               format.scale = T,
               xlab = "UMIs per Nucleus",
               ylab = "Genes per Nucleus",
               cor.coef.size = 16/.pt,
) +
  font("xlab", size = 20, color = "black")+
  font("ylab", size = 20, color = "black")+
  font("xy.text", size = 16, color = "black")
p
ggsave("./PlotsRevision/GenesvsUMIs_snRNAseq_rev.png",p,dpi = 300, width = 6, height = 6)

### CSF1 expression in Astrocytes
sNuc <- readRDS("./sNuc-seq/sNuc_SeuratObj.rds")
sub <- subset(sNuc,subset = CTs %in% c("Astrocyte","Radial Glia","OPC","OLIG1+GFAP+ Cells"))
p <- DotPlot_scCustom(sub,
                 features = c("CSF1"),
                 group.by = "CTs",
                 colors_use = pal_jco_size(256))
sNuc_data <- p@data
sNuc_data$Dataset <- "snRNA-seq"

y1 <- readRDS("./scRNA-seq/CellRanger-Count-10X/Analysis_Y1/scRNA_Y1_Clustered_SeuratObj.rds")
sub <- subset(y1,subset = CTs %in% c("Astrocyte","Radial Glia","OPC","OLIG1+GFAP+ Cells"))
p <- DotPlot_scCustom(sub,
                      features = c("CSF1"),
                      group.by = "CTs",
                      colors_use = pal_jco_size(256))
sc_data <- p@data
sc_data$Dataset <- "scRNA-seq"

data <- rbind(sNuc_data,sc_data)

brks <- (scales::pretty_breaks(n = 4))(c(0, 12))
p <- ggplot(data, aes(x = Dataset, y = id, size = pct.exp, color = avg.exp.scaled)) +
  geom_point() +
  scale_color_gradientn(colors = pal_jco_heatmap(256),values = scales::rescale(c(min(data$avg.exp.scaled),max(data$avg.exp.scaled)),to = c(0,1) )) +
  scale_size(range = c(1, 12),breaks = brks) +
  theme_minimal() +
  theme(axis.text = element_text(size = 16, color = "black"),
        legend.text = element_text(size = 16, hjust = 1),
        legend.title = element_text(size = 20)
  ) +
  labs(color = "Scaled\nExpression", size = "Percent\nof Cells", x = "", y = "")

ggsave("./PlotsRevision/Fig_CSF1_Expression_rev.png",p,width = 8, height = 8, dpi = 300)


###Fig4 Kampmann vs ours
load("./OtherStudies/Kampmann/Combined_SeuratObj.RData",verbose = T)
combined.sct$Origin <- "This study"
combined.sct$Origin[combined.sct$orig.ident != "snuc"] <- "Assembloid"

combined.sct$CTOrigin <- paste0(combined.sct$CTs," (",combined.sct$Origin,")")

cols <- c("Astrocyte (This study)" = "#CD534CFF",
          "Fibroblast-like (This study)" = "darkgrey",
          "Microglia (This study)" = "#0073C2FF",
          "Neuron (This study)" = "#EFC000FF",
          "OLIG1+GFAP+ Cells (This study)" = "darkgrey",
          "OPC (This study)" = "darkgrey",
          "Radial Glia (This study)" = "darkgrey",
          "Astrocyte (Assembloid)" = lighten("#CD534CFF",amount = 0.5),
          "Microglia (Assembloid)" = lighten("#0073C2FF",amount = 0.5),
          "Neuron (Assembloid)" = lighten("#EFC000FF",amount = 0.5))
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
                     "Microglia (This study)",
                     "Astrocyte (Assembloid)", 
                     "Neuron (Assembloid)",
                     "Microglia (Assembloid)")

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
ggsave("./PlotsRevision/Fig_UMAP_CTsOrigin_Kampmann_vs_ours_rev.png",p,width = 8, height = 8)

#Integration Linnarsson data
combined.sct <- readRDS("./OtherStudies/Linnarsson/Combined_SeuratObj.rds")
commonGenes <- rownames(combined.sct)

obj <- subset(combined.sct,subset = donor_id %in% c("H18.30.002", "H19.30.001", "H19.30.002"))
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

scNuc <- subset(combined.sct,subset = donor_id %in% c("sNuc"))
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

saveRDS(cormat, file = "./PlotsRevision/InputData_Heatmap_Linnarsson_vs_ours.rds")

cormat <- readRDS("./PlotsRevision/InputData_Heatmap_Linnarsson_vs_ours.rds")
p <- ggheatmap(cormat,
               cluster_rows = T,
               cluster_cols = T,
               hclust_method = "ward.D2",
               dist_method = "maximum",
               border = "black")
p$plotlist[[1]] <- p$plotlist[[1]] + 
  geom_tile(color = "black",linewidth = 1) +
  scale_fill_gradientn(colors = pal_jco_heatmap(3),values = scales::rescale(c(min(cormat),0,max(cormat)),to = c(0,1) )) +
  theme(axis.text = element_text(size = 16, color = "black"),
        legend.text = element_text(size = 16, hjust = 1),
        legend.title = element_text(size = 20),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
  ) +
  labs(fill = "Correlation")
p$plotlist[[2]] <- p$plotlist[[2]] +
  geom_tree(linewidth = 1)
p$plotlist[[3]] <- p$plotlist[[3]] +
  geom_tree(linewidth = 1)
p

ggsave(filename = "./PlotsRevision/Heatmap_Correlation_Linnarsson_vs_ours.png",plot = p, dpi = 300, width = 12, height = 6)

sort(unique(combined.sct$CTOrigin),decreasing = F)
cols <- c("Astrocyte (This study)" = "#CD534CFF",
          "Astrocyte (Brain)" = lighten("#CD534CFF",amount = 0.5),
          "Fibroblast-like (This study)" = "#868686FF",
          "Fibroblast (Brain)" = lighten("#868686FF",amount = 0.5),
          "Microglia (This study)" = "#0073C2FF",
          "Microglia (Brain)" = lighten("#0073C2FF",amount = 0.5),
          "Macrophage (Brain)" = darken("#0073C2FF",amount = 0.5),
          "Neuron (This study)" = "#EFC000FF",
          "Glutamatergic Neuron (Brain)" = lighten("#EFC000FF",amount = 0.5),
          "Brush Cell (Brain)" = lighten("#EFC000FF",amount = 0.7),
          "OPC (This study)" = "#3B3B3BFF",
          "OPC (Brain)" = darken("#3B3B3BFF",amount = 0.7),
          "OLIG1+GFAP+ Cells (This study)" = "#8F7700FF",
          "B cell (Brain)" = "darkgrey",
          "Bergmann Glial Cell (Brain)" = "darkgrey",
          "Interneuron (Brain)" = "darkgrey",
          "Granule Cell (Brain)" = "darkgrey",
          "Cortical Neuron (Brain)" = "darkgrey",
          "Epithelial Cell (Brain)" = "darkgrey",
          "Medium Spiny Neuron (Brain)" = "darkgrey",
          "Endothelial Cell (Brain)" = "darkgrey",
          "Ependymal Cell (Brain)" = "darkgrey",
          "Hippocampal Neuron (Brain)" = "darkgrey",
          "Lamp5+Lhx6+ Neuron (Brain)" = "darkgrey",
          "Mammillary Body Neuron (Brain)" = "darkgrey",
          "GABAergic Neuron (Brain)" = "darkgrey",
          "Monocyte (Brain)" = "darkgrey",
          "NK Cell (Brain)" = "darkgrey",
          "Neuron (Brain)" = "darkgrey",
          "Oligodendrocyte (Brain)" = "darkgrey",
          "Pericyte (Brain)" = "darkgrey",
          "Purkinje Cell (Brain)" = "darkgrey",
          "T Cell (Brain)" = "darkgrey",
          "SMC (Brain)" = "darkgrey",
          "Radial Glia (This study)" = "darkgrey")

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
                     "Astrocyte (Brain)", 
                     "Neuron (This study)",
                     "Brush Cell (Brain)",
                     "Fibroblast (Brain)",
                     "Fibroblast-like (This study)",
                     "Glutamatergic Neuron (Brain)",
                     "Macrophage (Brain)",
                     "Microglia (Brain)",
                     "Microglia (This study)",
                     "OLIG1+GFAP+ Cells (This study)",
                     "OPC (This study)",
                     "OPC (Brain)"
                     )

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
ggsave("./PlotsRevision/Fig_UMAP_CTsOrigin_rev_Linnarsson.png",p,width = 8, height = 8)


### scRNAseq: Plot UMAP + fractions
y1 <- readRDS("./scRNA-seq/Y1_final.rds")
p <- DimPlot_scCustom(seurat_object = y1,
                      repel = T,
                      group.by = "CTs",
                      figure_plot = TRUE,
                      label = F,
                      colors_use = colorMap,
                      prop_plot_percent = T,
                      add_prop_plot = T,
                      prop_plot_label = T,
                      label.size = 16)
p[[1]][[1]] <- p[[1]][[1]] + ggtitle("") + scale_color_manual(values = colorMap) + theme(legend.text = element_text(size = 16))
p[[2]][[1]] <- p[[2]][[1]] + scale_fill_manual(values = colorMap) + theme(axis.text = element_text(size = 16,color = "black"), axis.title = element_text(size = 20))
ggsave("./PlotsRevision/Fig4B_UMAP_rev.png",p[[1]], width = 6, height = 6)
ggsave("./PlotsRevision/Fig4G_UMAP_rev.png",p[[2]], width = 6, height = 6)

### snRNAseq: Plot UMAP + fractions
sNuc <- readRDS("./sNuc-seq/sNuc_SeuratObj.rds")
p <- DimPlot_scCustom(seurat_object = sNuc,
                      repel = T,
                      group.by = "CTs",
                      figure_plot = TRUE,
                      label = F,
                      colors_use = colorMap,
                      prop_plot_percent = T,
                      add_prop_plot = T,
                      prop_plot_label = T,
                      label.size = 16)
p[[1]][[1]] <- p[[1]][[1]] + ggtitle("") + scale_color_manual(values = colorMap) + theme(legend.text = element_text(size = 16))
p[[2]][[1]] <- p[[2]][[1]] + scale_fill_manual(values = colorMap) + theme(axis.text = element_text(size = 16,color = "black"), axis.title = element_text(size = 20))
ggsave("./PlotsRevision/Fig4A_UMAP_rev.png",p[[1]], width = 6, height = 6)
ggsave("./PlotsRevision/Fig4F_UMAP_rev.png",p[[2]], width = 6, height = 6)

### Plot GSE179590 comparison (single-nuclei)
preds_df <- readRDS("../PlotsRevision/InputData_RidgeBars_GSE179590_vs_ours.rds")
inp_bar <- readRDS("../PlotsRevision/InputData_inp_bar_RidgeBars_GSE179590_vs_ours.rds")
predIds <- readRDS("../PlotsRevision/InputData_predIds_RidgeBars_GSE179590_vs_ours.rds")

p_ridge <- ggplot(preds_df, aes(x = Scores, y = CT, fill = CT)) +
  geom_density_ridges(scale = 1) +
  scale_fill_manual(values = colorMap) +
  theme_ridges() + 
  theme(legend.position = "none",
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16)) +
  ylab("")


p_bar <- ggplot(inp_bar, aes(y = "",x = Freq, fill = predIds)) +
  geom_col(position = "stack",orientation = "y",width = 0.5) +
  scale_x_continuous(labels = scales::percent) +
  facet_grid(Var1 ~ ., scales = "free_y") +
  guides(fill = guide_legend(position = "right")) +
  scale_fill_manual(values = c(colorMap,"Endothelial cell" = "grey50","Oligodendrocyte" = "grey90")) +
  theme_ridges() +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(color = "transparent"),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 16),
    axis.title.x = element_text(color = "transparent"),
    strip.text.y = element_blank(),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 20),
    panel.spacing = unit(0, "lines")
  ) +
  xlab(" ") +
  labs(fill = "Cell Types")

p <- p_ridge + p_bar + plot_layout(ncol = 2,axes = "collect") & theme(axis.title.x = element_text(hjust = 0.5))
ggsave("./OtherStudies/GSE179590/CTPrediction_Comparison.png",p,dpi = 300, width = 10, height = 4)

### Correlation plot scOrganoid vs ours
cormat <- readRDS("./PlotsRevision/InputData_Heatmap_scOrganoid_vs_ours.rds")
p <- ggheatmap(cormat,
               cluster_rows = T,
               cluster_cols = T,
               hclust_method = "ward.D2",
               dist_method = "maximum",
               border = "black")
p$plotlist[[1]] <- p$plotlist[[1]] + 
  geom_tile(color = "black",linewidth = 1) +
  scale_fill_gradientn(colors = pal_jco_heatmap(3),values = scales::rescale(c(min(cormat),0,max(cormat)),to = c(0,1) )) +
  theme(axis.text = element_text(size = 16, color = "black"),
        legend.text = element_text(size = 16, hjust = 1),
        legend.title = element_text(size = 20),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
  ) +
  labs(fill = "Correlation")
p$plotlist[[2]] <- p$plotlist[[2]] +
  geom_tree(linewidth = 1)
p$plotlist[[3]] <- p$plotlist[[3]] +
  geom_tree(linewidth = 1)
p

ggsave(filename = "./PlotsRevision/Heatmap_Correlation_scOrganoid_vs_ours.png",plot = p, dpi = 300, width = 8, height = 6)

### Correlation plot Kampmann vs ours
cormat_2w <- readRDS("./PlotsRevision/InputData_Heatmap_Kampmann_iAssembloid2w_vs_ours.rds")
colnames(cormat_2w) <- paste0(colnames(cormat_2w)," (2 weeks)")
cormat_4w <- readRDS("./PlotsRevision/InputData_Heatmap_Kampmann_iAssembloid4w_vs_ours.rds")
colnames(cormat_4w) <- paste0(colnames(cormat_4w)," (4 weeks)")

cormat <- cbind(cormat_2w,cormat_4w)

p <- ggheatmap(cormat,
               cluster_rows = T,
               cluster_cols = T,
               hclust_method = "ward.D2",
               dist_method = "maximum",
               border = "black")
p$plotlist[[1]] <- p$plotlist[[1]] + 
  geom_tile(color = "black",linewidth = 1) +
  scale_fill_gradientn(colors = pal_jco_heatmap(3),values = scales::rescale(c(min(cormat),0,max(cormat)),to = c(0,1) )) +
  theme(axis.text = element_text(size = 16, color = "black"),
        legend.text = element_text(size = 16, hjust = 1),
        legend.title = element_text(size = 20),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
  ) +
  labs(fill = "Correlation")
p$plotlist[[2]] <- p$plotlist[[2]] +
  geom_tree(linewidth = 1)
p$plotlist[[3]] <- p$plotlist[[3]] +
  geom_tree(linewidth = 1)
p

ggsave(filename = "./PlotsRevision/Heatmap_Correlation_Kampann_vs_ours.png",plot = p, dpi = 300, width = 8, height = 6)

### Correlation plot GSE179590 vs ours
cormat <- readRDS("./PlotsRevision/InputData_Heatmap_GSE179590_vs_ours.rds")
p <- ggheatmap(cormat,
               cluster_rows = T,
               cluster_cols = T,
               hclust_method = "ward.D2",
               dist_method = "maximum",
               border = "black")
p$plotlist[[1]] <- p$plotlist[[1]] + 
  geom_tile(color = "black",linewidth = 1) +
  scale_fill_gradientn(colors = pal_jco_heatmap(3),values = scales::rescale(c(min(cormat),0,max(cormat)),to = c(0,1) )) +
  theme(axis.text = element_text(size = 16, color = "black"),
        legend.text = element_text(size = 16, hjust = 1),
        legend.title = element_text(size = 20),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
  ) +
  labs(fill = "Correlation")
p$plotlist[[2]] <- p$plotlist[[2]] +
  geom_tree(linewidth = 1)
p$plotlist[[3]] <- p$plotlist[[3]] +
  geom_tree(linewidth = 1)
p

ggsave(filename = "./PlotsRevision/Heatmap_Correlation_GSE179590_vs_ours.png",plot = p, dpi = 300, width = 8, height = 6)

### Plot Intercom results
#Fig 5A & 5B
ints_sc <- read.table("./scRNA-seq/Y1_forPaperFinal/FinalInteractions.txt",header = T, sep = "\t")
ints_sn <- read.table("./sNuc-seq/Intercom_forPaper_RNA/FinalInteractions.txt",header = T, sep = "\t")

library(CCPlotR)
inp_sc <- data.frame(source = ints_sc$Lig.pop,
                  target = ints_sc$Rec.pop,
                  ligand = ints_sc$Ligand,
                  receptor = ints_sc$Receptor,
                  score = log(ints_sc$score),
                  stringsAsFactors = F)
inp_sc$source[which(inp_sc$source == "Activated Microglia")] <- "Microglia"
inp_sc$target[which(inp_sc$target == "Activated Microglia")] <- "Microglia"
inp_sc$source[which(inp_sc$source == "Proliferating Microglia")] <- "Microglia"
inp_sc$target[which(inp_sc$target == "Proliferating Microglia")] <- "Microglia"
inp_sc <- aggregate(score~.,inp_sc,max)
inp_sc$ID <- paste0(inp_sc$source,"_",inp_sc$target,"_",inp_sc$ligand,"_",inp_sc$receptor)
inp_sc <- inp_sc[which(inp_sc$source %in% c("Microglia","Neuron","Astrocyte") & inp_sc$target %in% c("Microglia","Neuron","Astrocyte")),]

inp_sn <- data.frame(source = ints_sn$Lig.pop,
                     target = ints_sn$Rec.pop,
                     ligand = ints_sn$Ligand,
                     receptor = ints_sn$Receptor,
                     score = log(ints_sn$score),
                     stringsAsFactors = F)
inp_sn$ID <- paste0(inp_sn$source,"_",inp_sn$target,"_",inp_sn$ligand,"_",inp_sn$receptor)
inp_sn <- inp_sn[which(inp_sn$source %in% c("Microglia","Neuron","Astrocyte") & inp_sn$target %in% c("Microglia","Neuron","Astrocyte")),]

input_df <- inp_sc %>% mutate(source = factor(source), 
                             target = factor(target)) %>% group_by(source, target, 
                                                                   .drop = FALSE) %>% tally()
brks <- (scales::pretty_breaks(n = 5))(c(round_any(min(input_df$n),10, f = floor), round_any(max(input_df$n), 10, f = ceiling)))
p_sc <- cc_dotplot(inp_sc, option = "A") +
  scale_fill_gradientn(colors = pal_jco_size(6)) +
  scale_size(range = c(10, 20), limits = c(round_any(min(input_df$n),10, f = floor), round_any(max(input_df$n), 10,f = ceiling)), breaks = rev(brks), name = "Number of\ninteractions") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(angle=45,hjust=1,vjust = 1),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 20))
p_sc$guides$guides$size$params$override.aes$fill = rev(pal_jco_size(6))#viridis(5,direction = 1, option = "C")


input_df <- inp_sn %>% mutate(source = factor(source), 
                              target = factor(target)) %>% group_by(source, target, 
                                                                    .drop = FALSE) %>% tally()
brks <- (scales::pretty_breaks(n = 5))(c(round_any(min(input_df$n),10, f = floor), round_any(max(input_df$n), 10, f = ceiling)))
p_sn <- cc_dotplot(inp_sn, option = "A") +
  scale_fill_gradientn(colors = pal_jco_size(5)) +
  scale_size(range = c(10, 20), limits = c(round_any(min(input_df$n),10, f = floor), round_any(max(input_df$n), 10,f = ceiling)), breaks = rev(brks), name = "Number of\ninteractions") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(angle=45,hjust=1,vjust = 1),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 20))
p_sn$guides$guides$size$params$override.aes$fill = rev(pal_jco_size(5))#viridis(5,direction = 1, option = "C")

ggsave("./PlotsRevision/Fig5A_rev.png",p_sn,dpi=300,width = 6, height = 6)
ggsave("./PlotsRevision/Fig5B_rev.png",p_sc,dpi=300,width = 6, height = 6)

commonInts <- intersect(inp_sn$ID,inp_sc$ID)
inp_common <- merge(inp_sc,inp_sn,by = "ID")
inp_common$score <- rowMaxs(as.matrix(inp_common[,c("score.x","score.y")]))
inp_common <- inp_common[,c("source.x","target.x","ligand.x","receptor.x","score")]
colnames(inp_common) <- c("source","target","ligand","receptor","score")

library(ggraph)
library(ggsci)
p <- cc_network(inp_common,option = "B",n_top_ints = nrow(inp_common),colours = colorMap,layout = "kk",label_size = 16/.pt, node_size = 2.5) +
  scale_fill_manual(values = colorMap,name = "Cell Type") +
  theme(legend.text=element_text(size=16),
        legend.title = element_text(size = 20))

ggsave("./PlotsRevision/Fig5C_rev.png",p,dpi=300,width = 8, height = 8)

library(ggpubr)

df <- as.data.frame(table(inp_sc$target,inp_sc$receptor))
df <- df[which(df$Freq > 0),]
df$Var1 <- factor(as.character(df$Var1),levels = c("Neuron","Astrocyte","Microglia"))
df$Var2 <- factor(as.character(df$Var2), levels = unique(as.character(df$Var2)))
df$Var2_unique <- paste0(df$Var2, "_", seq_len(nrow(df)))
rownames(df) <- NULL

oneCT_sc <- names(which(table(df$Var2)==1))
multInts_sc <- names(which(table(inp_sc$receptor) > 1))
oneInt_sc <- names(which(table(inp_sc$receptor) == 1))
oneCTmultInts_sc <- intersect(oneCT_sc,multInts_sc)
oneCToneInt_sc <- intersect(oneCT_sc,oneInt_sc)

p <- ggbarplot(df, x = "Var2_unique", y = "Freq",
          fill = "Var1",               # change fill color by cyl
          color = "white",            # Set bar border colors to white
          palette = "jco",            # jco journal color palett. see ?ggpar
          sort.val = "desc",          # Sort the value in dscending order
          sort.by.groups = T,     # Don't sort inside each group
          x.text.angle = 90,           # Rotate vertically x axis texts
          xlab = "Receptor"
) +
  font("legend.title",size = 20) +
  font("legend.text",size = 16) +
  font("xy.title", size = 20) +
  font("xy.text",size = 16) +
  scale_fill_manual(values = colorMap,name = "Cell Type")
pb <- ggplot_build(p)

# Access x-axis data (after sorting)
ordered_levels <- pb$data[[1]]$x
ordered_labels <- pb$layout$panel_params[[1]]$x$get_labels()
new_labels <- gsub("_[0-9]+","",ordered_labels)

p <- p + scale_x_discrete(labels = new_labels)
p

ggsave("./PlotsRevision/Fig5E_rev.png",p,dpi=300,width = 8, height = 6)

df <- as.data.frame(table(inp_sn$target,inp_sn$receptor))
df <- df[which(df$Freq > 0),]
df$Var1 <- factor(as.character(df$Var1),levels = c("Neuron","Astrocyte","Microglia"))
df$Var2 <- factor(as.character(df$Var2), levels = unique(as.character(df$Var2)))
df$Var2_unique <- paste0(df$Var2, "_", seq_len(nrow(df)))
rownames(df) <- NULL

oneCT_sn <- names(which(table(df$Var2)==1))
multInts_sn <- names(which(table(inp_sn$receptor) > 1))
oneInt_sn <- names(which(table(inp_sn$receptor) == 1))
oneCTmultInts_sn <- intersect(oneCT_sn,multInts_sn)
oneCToneInt_sn <- intersect(oneCT_sn,oneInt_sn)

p <- ggbarplot(df, x = "Var2_unique", y = "Freq",
               fill = "Var1",               # change fill color by cyl
               color = "white",            # Set bar border colors to white
               palette = "jco",            # jco journal color palett. see ?ggpar
               sort.val = "desc",          # Sort the value in dscending order
               sort.by.groups = T,     # Don't sort inside each group
               x.text.angle = 90,           # Rotate vertically x axis texts
               xlab = "Receptor"
) +
  font("legend.title",size = 20) +
  font("legend.text",size = 16) +
  font("xy.title", size = 20) +
  font("xy.text",size = 16) +
  scale_fill_manual(values = colorMap,name = "Cell Type")
pb <- ggplot_build(p)

# Access x-axis data (after sorting)
ordered_levels <- pb$data[[1]]$x
ordered_labels <- pb$layout$panel_params[[1]]$x$get_labels()
new_labels <- gsub("_[0-9]+","",ordered_labels)

p <- p + scale_x_discrete(labels = new_labels)
p

ggsave("./PlotsRevision/Fig5D_rev.png",p,dpi=300,width = 8, height = 6)

paste0(unique(sort(c(oneCTmultInts_sc,oneCTmultInts_sn))),collapse = ", ")
paste0(unique(sort(c(oneCToneInt_sc,oneCToneInt_sn))),collapse = ", ")

