library(Seurat)
library(scCustomize)
library(colorspace)
library(ggrepel)

source("./aux/colorDefinitions.R")

load("./IntegrationWithOtherStudies/scOrganoids/Combined_SeuratObj.RData",verbose = T)
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
ggsave("./PaperPlots/Fig4G.png",p,width = 8, height = 8)
