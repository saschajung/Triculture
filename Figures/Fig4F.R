library(Seurat)
library(scCustomize)
library(colorspace)
library(ggrepel)

#### NOTE ######
# The label "Inhibitory Neuron (Brain)" of the lower-left, light yellow cluster
# has been manually added and is not reproduced with this script! 
################

source("./aux/colorDefinitions.R")

load("./IntegrationWithOtherStudies/GSE179590/Combined_SeuratObj.RData",verbose = T)

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
ggsave("./PaperPlots/Fig4F.png",p,width = 8, height = 8)
