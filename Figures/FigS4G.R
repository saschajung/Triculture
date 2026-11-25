library(Seurat)
library(ggplot2)
library(scCustomize)
library(colorspace)
library(ggrepel)

combined.sct <- readRDS("./IntegrationWithOtherStudies/Linnarsson/Combined_SeuratObj.rds")

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
ggsave("./PaperPlots/FigS4G.png",p,width = 8, height = 8)
