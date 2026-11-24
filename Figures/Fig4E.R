library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(tibble)
library(tidyr)

source("./aux/colorDefinitions.R")

if(!file.exists("./InputData/InputData_Fig4E.rds")){
  source("./aux/CreateInputData_Fig4E.R")
}

inp <- readRDS("./InputData/InputData_Fig4E.rds")

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

png("./PaperPlots/Fig4E.png", width = 6.5, height = 12.5,units = "in", pointsize = 16, res = 300)
draw(ht)
dev.off()
