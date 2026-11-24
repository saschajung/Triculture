library(Seurat)
library(scCustomize)
library(colorspace)
library(ggrepel)

combined.sct <- readRDS("./IntegrationWithOtherStudies/Linnarsson/Combined_SeuratObj.rds")
commonGenes <- rownames(combined.sct)

#Remove donor H18.30.001 due to low cell number
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

saveRDS(cormat, file = "./InputData/InputData_FigS4H.rds")

rm(combined.sct,obj)

