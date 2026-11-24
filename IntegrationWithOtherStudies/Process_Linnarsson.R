library(anndataR)
library(Seurat)
library(scCustomize)
library(ggsci)

###### PYTHON #####
### Downsample original h5ad file to 1000 cells per cell type

# import os
# import glob
# import gc
# import numpy as np
# import scanpy as sc
# 
# files = glob.glob(os.path.expanduser("~/Downloads/*.h5ad"))
# output_dir = os.path.expanduser("./IntegrationWithOtherStudies/Linnarsson/")
# 
# rng = np.random.default_rng(42)
# 
# for f in files:
#   print(f)
# adata = sc.read_h5ad(f, backed="r")  # metadata only in memory; X stays on disk
# 
# cell_types = adata.obs["cell_type"].unique()
# 
# for ct in cell_types:
#   print(f"   {ct}")
# idx = np.where(adata.obs["cell_type"].values == ct)[0]
# n = min(len(idx), 1000)
# chosen_idx = np.sort(rng.choice(idx, size=n, replace=False))  # sorted: required for backed h5py indexing
# 
# adata_sub = adata[chosen_idx].to_memory()   # only now does data get read into RAM, and only for these cells
# adata_sub.layers["counts"] = adata_sub.X.copy()  # cheap now — subset is tiny
# 
# safe_ct = "".join(c if c.isalnum() else "_" for c in str(ct))
# out_path = os.path.join(output_dir, f"{safe_ct}.h5ad")
# adata_sub.write_h5ad(out_path)
# 
# del adata_sub
# gc.collect()
# 
# adata.file.close()  # release the backed file handle
# del adata
# gc.collect()

##### R
files <- list.files(path = "./IntegrationWithOtherStudies/Linnarsson",pattern = ".h5ad",full.names = T)

obj_list <- list()
for(f in files){
  print(f)
  seurat_obj <- read_h5ad(
    f,
    as = "Seurat"
  )
  seurat_obj[["RNA"]]$counts <- LayerData(seurat_obj[["RNA"]], layer = "X")
  seurat_obj[["RNA"]]$X <- NULL
  
  n <- tools::file_path_sans_ext(gsub("./IntegrationWithOtherStudies/Linnarsson/","",f,fixed = T))
  obj_list[[n]] <- seurat_obj
  
}

seurat_obj <- Reduce(merge,obj_list)
seurat_obj[["RNA"]] <- JoinLayers(seurat_obj[["RNA"]])

options(future.globals.maxSize = 10000 * 1024^2)
mem.maxVSize(100000)

donor_list <- SplitObject(seurat_obj, split.by = "donor_id")
donor_list[["scNuc"]] <- readRDS("./sNuc-seq/sNuc_SeuratObj.rds")
donor_list[["scNuc"]]$cell_type <- donor_list[["scNuc"]]$CTs
donor_list[["scNuc"]] <- UpdateSeuratObject(donor_list[["scNuc"]])
DefaultAssay(donor_list[["scNuc"]]) <- "RNA"

library(org.Hs.eg.db)
library(AnnotationDbi)

# identify which donor(s) use gene symbols vs Ensembl IDs
symbol_donor <- "scNuc"  # adjust to whichever one is the mismatch

symbols <- rownames(donor_list[[symbol_donor]][["RNA"]])  # rename BEFORE SCTransform if possible

ensembl_ids <- mapIds(
  org.Hs.eg.db,
  keys = symbols,
  column = "ENSEMBL",
  keytype = "SYMBOL",
  multiVals = "first"   # takes first match when a symbol maps to >1 Ensembl ID
)

# drop genes that didn't map, and any resulting duplicate Ensembl IDs
keep <- !is.na(ensembl_ids) & !duplicated(ensembl_ids)
new_names <- ensembl_ids[keep]

donor_list[[symbol_donor]][["RNA"]] <- as(donor_list[[symbol_donor]][["RNA"]], Class = "Assay5")
seurat_obj <- donor_list[[symbol_donor]][keep, ]
dimnames(seurat_obj[["RNA"]])[[1]] <- unname(new_names)
donor_list[[symbol_donor]] <- seurat_obj

donor_list <- lapply(donor_list,function(x){
  DefaultAssay(x) <- "RNA"
  x = DietSeurat(x,assays = "RNA")
  x = SCTransform(x,vst.flavor = "v2", verbose = T)
  DefaultAssay(x) <- "SCT"
  x = DietSeurat(x,assays = c("RNA","SCT"),scale.data = T)
  x = RunPCA(x, npcs = 30,verbose = FALSE)
  return(x)
})

gc()
gc()
gc()
gc()

features <- SelectIntegrationFeatures(object.list = donor_list, nfeatures = 3000, assay = rep("SCT",length(donor_list)))
donor_list <- PrepSCTIntegration(object.list = donor_list,
                               anchor.features = features,
                               verbose = T)

anchors <- FindIntegrationAnchors(object.list = donor_list,
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
plot.path <- "./IntegrationWithOtherStudies/Linnarsson/"
ggsave(paste0(plot.path,"ElbowPlot.png"),ElbowPlot(combined.sct,ndims = 100)) #Take 50

numPCs <- 50

DefaultAssay(combined.sct) <- "integrated"
combined.sct <- RunUMAP(combined.sct,
                        reduction = "pca",
                        dims = 1:numPCs,
                        #n.neighbors = 20,
                        verbose = T)

combined.sct <- RunTSNE(combined.sct, reduction = "pca", dims = 1:numPCs, verbose = TRUE)

combined.sct$donor_id[is.na(combined.sct$donor_id)] <- "sNuc"

combined.sct$Origin <- "Brain"
combined.sct$Origin[combined.sct$donor_id == "sNuc"] <- "This study"

ctMap <- c("amygdala excitatory neuron" = "Glutamatergic Neuron",
  "astrocyte" = "Astrocyte",
  "B cell" = "B cell",
  "Bergmann glial cell" = "Bergmann Glial Cell",
  "caudal ganglionic eminence derived interneuron" = "Interneuron",
  "central nervous system macrophage" = "Macrophage",
  "cerebellar granule cell" = "Granule Cell",
  "cerebellar inhibitory GABAergic interneuron" = "Interneuron",
  "cerebral cortex neuron" = "Cortical Neuron",
  "choroid plexus epithelial cell" = "Epithelial Cell",
  "committed oligodendrocyte precursor" = "OPC",
  "dentate gyrus granule cell" = "Granule Cell",
  "eccentric medium spiny neuron" = "Medium Spiny Neuron",
  "endothelial cell" = "Endothelial Cell",
  "ependymal cell" = "Ependymal Cell",
  "fibroblast" = "Fibroblast",
  "hippocampal CA1-3 neuron" = "Hippocampal Neuron",
  "hippocampal CA4 neuron" = "Hippocampal Neuron",
  "intratelencephalic-projecting glutamatergic cortical neuron" = "Glutamatergic Neuron",
  "Lamp5 Lhx6 neuron" = "Lamp5+Lhx6+ Neuron",
  "mammillary body neuron" = "Mammillary Body Neuron",
  "medial ganglionic eminence derived interneuron" = "Interneuron",
  "medium spiny neuron" = "Medium Spiny Neuron",
  "microglial cell" = "Microglia",
  "midbrain-derived inhibitory neuron" = "GABAergic Neuron",
  "monocyte" = "Monocyte",
  "natural killer cell" = "NK Cell",
  "near-projecting glutamatergic cortical neuron" = "Glutamatergic Neuron",
  "neuron" = "Neuron",
  "oligodendrocyte precursor cell" = "OPC",
  "oligodendrocyte" = "Oligodendrocyte",
  "pericyte" = "Pericyte",
  "Purkinje cell" = "Purkinje Cell",
  "pvalb chandelier GABAergic interneuron" = "Interneuron",
  "sst chodl GABAergic interneuron" = "Interneuron",
  "T cell" = "T Cell",
  "thalamic excitatory neuron" = "Glutamatergic Neuron",
  "unipolar brush cell" = "Brush Cell",
  "vascular associated smooth muscle cell" = "SMC",
  "Astrocyte" = "Astrocyte",
  "Radial Glia" = "Radial Glia",
  "OLIG1+GFAP+ Cells" = "OLIG1+GFAP+ Cells",
  "OPC" = "OPC",
  "Microglia" = "Microglia",
  "Fibroblast-like" = "Fibroblast-like",
  "Neuron" = "Neuron")

combined.sct$CTs <- unname(ctMap[combined.sct$cell_type])

combined.sct$CTOrigin <- paste0(combined.sct$CTs," (",combined.sct$Origin,")")

sub <- subset(combined.sct, subset = CTs == "OPC")
DimPlot(sub,reduction = "umap", group.by = "Origin")

saveRDS(combined.sct,file = "./IntegrationWithOtherStudies/Linnarsson/Combined_SeuratObj.rds")



