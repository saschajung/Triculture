library(Seurat)
library(ggplot2)
library(scCustomize)

source("./aux/colorDefinitions.R")

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

if(!dir.exists("./PaperPlots/")){
  dir.create("./PaperPlots/")
}
ggsave("./PaperPlots/Fig4B.png",p[[1]], width = 6, height = 6)
