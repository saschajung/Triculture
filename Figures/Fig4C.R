library(Seurat)
library(ggplot2)
library(scCustomize)

source("./aux/colorDefinitions.R")

scNuc <- readRDS(file = "./sNuc-seq/sNuc_SeuratObj.rds")

opc_markers <- c("PDGFRA","OLIG1")
radialGlia_markers <- c("FABP7","SOX2","HES1","HES5")
neuron_markers <- c("RBFOX3","SNAP25","MAP2","STMN2","TUBB3")
fibroblast <- c("COL1A1","COL1A2","DCN","LUM","KRT18")
microglia <- c("P2RY12","CX3CR1","CSF1R","SALL1","CCL2","CCL3","CCL4","MKI67")
astro <- c("GFAP","AQP4","SLC1A3","ALDH1L1")
featuresToPlot <- c(opc_markers,radialGlia_markers,
                    neuron_markers,fibroblast,microglia,astro)

if(!dir.exists("./PaperPlots/")){
  dir.create("./PaperPlots/")
}

png(filename = "./PaperPlots/Fig4C.png",res = 300,width = 8, height = 8,units = "in",pointsize = 16)
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
